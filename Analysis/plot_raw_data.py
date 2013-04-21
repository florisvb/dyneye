import os
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
    rdes = estimates_data['rdes']
    fps = estimates_data['fps']
    
    LSE = least_squares.LeastSquaresEstimator()
    Hfunc = lambda t: np.exp(rdes*t)
    H = least_squares.H_LTI([Hfunc])
    LSE.initialize_with_guess(np.matrix([3]),np.matrix([1000]))

    # observer    
    measurements = estimates_data['control']*fps/100. / rdes**2

    err = (estimates_data['filteredopticflow'] - rdes)**2
    position = []
    velocity = []
    smoothing = False
    for i, measurement in enumerate(measurements):
        t = estimates_data['time'][i]
        Ht = H(t)
        
        #########################################################
        # sketchy bit
        if np.abs(estimates_data['filteredopticflowdot'][i]) > 0.01 or i < 20:
            if smoothing is False:
                LSE.Pk = 100
                print LSE.Pk
        else:
            smoothing = True
        
        #LSE.Pk = 40*estimates_data['filteredopticflowdot'][i]**2
        #########################################################
        
        LSE.update([measurement], Ht)
        position.append(LSE.xhat[0,0]*Ht[0,0])
        velocity.append(LSE.xhat[0,0]*rdes*Ht[0,0])

    position = np.array(position)
    velocity = np.array(velocity)

    return position, velocity
        
        
    
    
def plot(estimates_data, stepper_data=None):
    
    fig = plt.figure()
    axpos = fig.add_subplot(151)
    axvel = fig.add_subplot(152)
    axof = fig.add_subplot(153)
    axofdot = fig.add_subplot(154)
    axctrl = fig.add_subplot(155)    
    axes = [axpos, axvel, axof, axofdot, axctrl]
    
    fps = estimates_data['fps']
    rdes = estimates_data['rdes']
    
    estimates_data, stepper_data = interpolate_to_same_time(estimates_data, stepper_data)    
    noise = get_optic_flow_noise_estimate(estimates_data, stepper_data)
    print 'optic flow noise estimate: ', noise
    
    print np.max(stepper_data['velocity'])
    index_stop = np.where(np.diff(stepper_data['position'])==0)[0][2]
    tstop = stepper_data['time'][index_stop]
    
    # plot raw STEPPER data
    axpos.plot(stepper_data['time'], stepper_data['position'], 'blue', linewidth=3)
    axvel.plot(stepper_data['time'], stepper_data['velocity'], 'blue', linewidth=3)
    axof.plot(stepper_data['time'], stepper_data['velocity']/stepper_data['position'], 'blue', linewidth=3)
    axctrl.plot(stepper_data['time'], stepper_data['control'], 'blue', linewidth=3)
    
    # plot estimates data
    axof.plot(estimates_data['time'], estimates_data['filteredopticflow'], 'red')
    axofdot.plot(estimates_data['time'], estimates_data['filteredopticflowdot'], 'red')
    
    # plot raw guess
    guess = estimates_data['control']*fps/100. / rdes**2
    axpos.plot(estimates_data['time'], guess, '.', color='green')
    
    lse_position, lse_velocity = sequential_least_squares_observer(estimates_data)

    axpos.plot(estimates_data['time'], lse_position, 'red')
    axvel.plot(estimates_data['time'], lse_velocity, 'red')

    # formatting
    for ax in axes:
        ax.set_xlim(0, tstop)
        
    axpos.set_ylim(-2,2.5)
    axof.set_ylim(rdes*(1+0.2), 0)

    plt.show()
    
    
def get_optic_flow_noise_estimate(estimates_data, stepper_data):
    stepper_optic_flow = stepper_data['velocity']/stepper_data['position']
    difference = estimates_data['opticflow'] - stepper_optic_flow
    variance = np.mean(difference[200:-200]**2)
    return variance**(1/2.)
    
if __name__=='__main__':

    parser = OptionParser()
    parser.add_option("--f", type="str", dest="directory", default='',
                        help="directory name that has the estimates and stepper file in it")
    parser.add_option("--estimates", type="str", dest="estimates", default='',
                        help="filename to plot")
    parser.add_option("--stepper", type="str", dest="stepper", default='',
                        help="filename to plot") 
    (options, args) = parser.parse_args()
    
    if len(options.directory) > 0:
        files = os.listdir(options.directory)
        for fname in files:
            if 'stepper' in fname:
                options.stepper = os.path.join(options.directory, fname)
            if 'estimates' in fname and 'image' not in fname:
                options.estimates = os.path.join(options.directory, fname)
    
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
