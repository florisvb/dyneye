import numpy as np
import pickle
import matplotlib.pyplot as plt
import adskalman.adskalman as adskalman
from optparse import OptionParser

import data_fit
import copy
import fly_plot_lib.plot as fpl



import pykalman

def kalman_smoother(data, F, H, Q, R, initx, initv, plot=False):
    os = H.shape[0]
    ss = F.shape[0]
    
        
    xsmooth,Vsmooth = adskalman.kalman_smoother(data,F,H,Q,R,initx,initv)
    print xsmooth.shape
    if plot:
        fig = plt.figure()
        plt.plot(data[:], '.')
        plt.plot(xsmooth[:,0])

    return xsmooth,Vsmooth


def fix_array_size(data):
    l = []
    for key, item in data.items():
        l.append(len(item))
    l = np.max(l)
    for key, item in data.items():
        while len(item) < l:
            item.append(item[-1])
        item = np.array(item)
        data[key] = item
    return data
    
    
def dyneye_observer(control, rdes, fps):
    posest = control / rdes**2
    return posest
    
def kalman_smooth(estimates_data, stepper_data):
    
    estimates_data, stepper_data = interpolate_to_same_time(estimates_data, stepper_data)    
    measurements = dyneye_observer(estimates_data['control'], -0.1, 50.)
    
    A = np.array([[1,1/(2*np.pi*50.)],[0,1]])
    C = np.array([[1,0]])
    
    optic_flow_noise = get_optic_flow_noise_estimate(estimates_data, stepper_data)
    R = np.eye(1)
    R[0,0] = 1/optic_flow_noise**2
    Q = np.eye(2)*1e-4#np.array([[1, 0], [0, 1]])
    initial_covariance = 100*np.eye(2)
    
    #xsmooth, vsmooth = kalman_smoother(measurements[400:]*10, A, C, Q, R, np.array([0,0]), 100*np.eye(2), plot=True)
    
    #return xsmooth
    
    kf = pykalman.KalmanFilter(transition_matrices = A, observation_matrices = C, observation_covariance = R, transition_covariance = Q, initial_state_covariance=initial_covariance)
    #kf = pykalman.KalmanFilter(transition_matrices = A, observation_matrices = C, observation_covariance = [measurement_noise], transition_covariance = process_noise, initial_state_covariance=initial_covariance)
    (smoothed_state_means, smoothed_state_covariances) = kf.filter(measurements)
    
    return stepper_data['time'], smoothed_state_means
            
def interpolate_to_same_time(estimates_data, stepper_data):
    
    estimates_data = fix_array_size(estimates_data)
    estimates_data['time'] -= estimates_data['time'][0]
    stepper_data = fix_array_size(stepper_data)
    stepper_data['time'] -= stepper_data['time'][0]
    
    t_interp_start = np.max([stepper_data['time'][0], estimates_data['time'][0]])
    t_interp_stop = np.min([stepper_data['time'][-1], estimates_data['time'][-1]])
    t_interp = np.arange(t_interp_start, t_interp_stop, 0.005)
    
    def interp_data(data, t_interp):
        data_time = copy.copy(data['time'])
        data_interp = {}
        for key, item in data.items():
            if key == 'time':
                data_interp.setdefault('time', t_interp)
            else:
                interp_item = np.interp(t_interp, data_time, item)
                data_interp.setdefault(key, interp_item)
        return data_interp
        
    stepper_interp = interp_data(stepper_data, t_interp)
    estimates_interp = interp_data(estimates_data, t_interp)

    return estimates_interp, stepper_interp

def plot(estimates_data, stepper_data=None):
    
    fig = plt.figure()
    axpos = fig.add_subplot(141)
    axvel = fig.add_subplot(142)
    axof = fig.add_subplot(143)
    axctrl = fig.add_subplot(144)    
    
    fps = 50.
    rdes = -.1
    
    estimates_data, stepper_data = interpolate_to_same_time(estimates_data, stepper_data)    
    
    # plot raw data
    axpos.plot(stepper_data['time'], stepper_data['position'], 'blue')
    axvel.plot(stepper_data['time'], stepper_data['velocity'], 'blue')
    axof.plot(stepper_data['time'], stepper_data['velocity']/stepper_data['position'], 'blue')
    axctrl.plot(stepper_data['time'], stepper_data['control'], 'blue')
    
    # plot estimates
    posest = estimates_data['control']*(fps/100.) / rdes**2
    #axpos.plot(estimates_data['time'], posest, '.', color='green')
    measurements = dyneye_observer(estimates_data['control'], -0.1, 50.)
    axpos.plot(estimates_data['time'], measurements, 'green', zorder=-10)

    axof.plot(estimates_data['time'], estimates_data['opticflow'], 'red')

    optic_flow_noise = get_optic_flow_noise_estimate(estimates_data, stepper_data)
    print optic_flow_noise
    
    # 
    
    t, smoothed_state_means = kalman_smooth(estimates_data, stepper_data)
    axpos.plot(t,smoothed_state_means[:,0], 'red')
    axvel.plot(t,-1*smoothed_state_means[:,1], 'red')

    plt.show()
    
    
def get_optic_flow_noise_estimate(estimates_data, stepper_data):
    estimates_data, stepper_data = interpolate_to_same_time(estimates_data, stepper_data)   
    stepper_optic_flow = stepper_data['velocity']/stepper_data['position']
    
    difference = estimates_data['opticflow'] - stepper_optic_flow
    
    variance = np.mean(difference**2)
    
    return variance
    
if __name__=='__main__':

    parser = OptionParser()
    parser.add_option("--estimates", type="str", dest="estimates", default='',
                        help="filename to plot")
    parser.add_option("--stepper", type="str", dest="stepper", default='',
                        help="filename to plot") 
    (options, args) = parser.parse_args()
    
    f_estimates = open(options.estimates)
    estimates_data = pickle.load(f_estimates)
    f_estimates.close()
    
    if len(options.stepper) > 0:
        f_stepper = open(options.stepper)
        stepper_data = pickle.load(f_stepper)
        f_stepper.close()
    else:
        stepper_data = None
        
    plot(estimates_data, stepper_data)
