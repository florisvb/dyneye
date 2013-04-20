import numpy as np
import pickle
import matplotlib.pyplot as plt
from optparse import OptionParser

import data_fit
import copy
import fly_plot_lib.plot as fpl

from data_fit import least_squares

import pykalman

def fix_array_size(data):
    l = []
    for key, item in data.items():
        if type(item) is list:
            l.append(len(item))
    l = np.max(l)
    for key, item in data.items():
        if type(item) is list:
            while len(item) < l:
                item.append(item[-1])
            item = np.array(item)
            data[key] = item
    return data
    
'''
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
'''

def interpolate_to_same_time(estimates_data, stepper_data):
    
    estimates_data = fix_array_size(estimates_data)
    tstart = copy.copy(estimates_data['time'][0])
    estimates_data['time'] -= tstart
    stepper_data = fix_array_size(stepper_data)
    stepper_data['time'] -= tstart
    
    t_interp_start = np.max([stepper_data['time'][0], estimates_data['time'][0]])
    t_interp_stop = np.min([stepper_data['time'][-1], estimates_data['time'][-1]])
    t_interp = np.arange(t_interp_start, t_interp_stop, 0.02)
    
    def interp_data(data, t_interp):
        data_time = copy.copy(data['time'])
        data_interp = {}
        for key, item in data.items():
            print key, type(item)
            if key == 'time':
                data_interp.setdefault('time', t_interp)
            elif type(item) == np.ndarray:
                interp_item = np.interp(t_interp, data_time, item)
                data_interp.setdefault(key, interp_item)
            else:
                data_interp.setdefault(key, item)
        return data_interp
        
    stepper_interp = interp_data(stepper_data, t_interp)
    estimates_interp = interp_data(estimates_data, t_interp)

    return estimates_interp, stepper_interp
    

def sequential_least_squares_observer(estimates_data):
    optic_flow_noise = 0.006
    rdes = -.1
    fps = 50.
    
    LSE = least_squares.LeastSquaresEstimator()
    Hfunc = lambda t: np.exp(rdes*t)
    H = least_squares.H_LTI([Hfunc])
    
        
        
    if 1:
        LSE.initialize_with_guess(np.matrix([3]),np.matrix([1000]))
        
        measurements = estimates_data['control']*fps/100. / rdes**2
        err = np.abs(estimates_data['filteredopticflow'] - rdes)
        indices = np.arange(0,1000)
        best_guess = []
        best_d0_guess = 3
        for i, measurement in enumerate(measurements):
            if i in indices:
                t = estimates_data['time'][i]
                Ht = H(t)
                if i<100:
                    LSE.Pk = 100
                #LSE.Pk = 40*estimates_data['filteredopticflowdot'][i]**2
                LSE.update([measurement], Ht)
                best_guess.append(LSE.xhat[0,0]*Ht[0,0])
                best_d0_guess = LSE.xhat[0,0]
            else:
                try:
                    t = estimates_data['time'][i]
                    Ht = H(t)[0,0]
                    best_guess.append(best_d0_guess*Ht)
                except:
                    best_guess = [best_d0_guess]
        return np.array(best_guess)
        
        
    
    
def plot(estimates_data, stepper_data=None):
    
    fig = plt.figure()
    axpos = fig.add_subplot(151)
    axvel = fig.add_subplot(152)
    axof = fig.add_subplot(153)
    axofdot = fig.add_subplot(154)
    axctrl = fig.add_subplot(155)    
    
    fps = 50.
    rdes = -.1
    
    estimates_data, stepper_data = interpolate_to_same_time(estimates_data, stepper_data)    
    noise = get_optic_flow_noise_estimate(estimates_data, stepper_data)
    print noise
    
    # plot raw STEPPER data
    axpos.plot(stepper_data['time'], stepper_data['position'], 'blue')
    axvel.plot(stepper_data['time'], stepper_data['velocity'], 'blue')
    axof.plot(stepper_data['time'], stepper_data['velocity']/stepper_data['position'], 'blue')
    axctrl.plot(stepper_data['time'], stepper_data['control'], 'blue')
    
    # plot estimates data
    axof.plot(estimates_data['time'], estimates_data['filteredopticflow'], 'red')
    axofdot.plot(estimates_data['time'], estimates_data['filteredopticflowdot'], 'red')
    
    # plot raw guess
    guess = estimates_data['control']*fps/100. / rdes**2
    axpos.plot(estimates_data['time'], guess, '.', color='green')
    
    best_guess = sequential_least_squares_observer(estimates_data)

    axpos.plot(estimates_data['time'], best_guess, 'red')

    plt.show()
    
    
def get_optic_flow_noise_estimate(estimates_data, stepper_data):
    stepper_optic_flow = stepper_data['velocity']/stepper_data['position']
    
    difference = estimates_data['opticflow'] - stepper_optic_flow
    
    variance = np.mean(difference[200:-200]**2)
    
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
