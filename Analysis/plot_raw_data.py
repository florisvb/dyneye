import os
import numpy as np
import scipy.stats
import scipy.io
import pickle
import matplotlib.pyplot as plt
from optparse import OptionParser

import data_fit
import copy
import fly_plot_lib.plot as fpl

from data_fit import least_squares

import pykalman
import pykalman.sqrt

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
    
    diffpos = np.diff(stepper_data['position'])
    middle_index = int(len(stepper_data['position'])/2.)
    index_stop = np.where(diffpos[middle_index:]==0)[0][0] + middle_index
    
    t_interp_stop = stepper_data['time'][index_stop] #np.min([stepper_data['time'][-1], estimates_data['time'][-1]])
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
    
    estimates_interp['control'] *= 2.

    return estimates_interp, stepper_interp
    
    

################################
    
    
class ukf_filter:
    def __init__(self, x0, P0):
        self.dt = 0.02
        
        self.xhat = x0
        self.P = P0
        self.xhistory = []
        self.rdes = -.1
        
        self.ukf = pykalman.AdditiveUnscentedKalmanFilter(self.f, self.g, n_dim_obs=1, n_dim_state=2, observation_covariance=0.006, transition_covariance=[[1e-1, 0], [0, 1e-1]])
        
        self.noise = scipy.stats.norm(0,1e-7)
    
    
    def f(self, x): # process
        x_new = [x[0] + self.dt*x[1], x[1] + self.dt*self.control]
        return x_new
    
    def g(self, x): # observation
        return x[1]/x[0] # + self.noise.rvs()
    
    def update(self, measurement, control):
        self.control = control
        self.xhat, self.P = self.ukf.filter_update(self.xhat, self.P, observation=measurement)
        self.P = self.P
        self.xhistory.append(self.xhat)
    
    def iterate_update(self, measurements, controls):
        for i, measurement in enumerate(measurements):
            self.update(measurement, controls[i])
        self.xhistory = np.array(self.xhistory)
    
    
    
#################################
    

def sequential_least_squares_observer(estimates_data):
    rdes = estimates_data['rdes']
    fps = estimates_data['fps']
    
    LSE = least_squares.LeastSquaresEstimator()
    Hfunc = lambda t: np.exp(rdes*t)
    H = least_squares.H_LTI([Hfunc])
    LSE.initialize_with_guess(np.matrix([5]),np.matrix([1000]))

    # observer    
    measurements = estimates_data['control']*fps/100./2. / rdes**2

    err = (estimates_data['filteredopticflow'] - rdes)
    position = []
    velocity = []
    smoothing = False
    for i, measurement in enumerate(measurements):
        t = estimates_data['time'][i]
        Ht = H(t)
        
        #########################################################
        # sketchy bit
        if np.abs(estimates_data['filteredopticflowdot'][i]) > 0.001 or i < 20:
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
    
    fig = plt.figure(figsize=(10,3))
    fig.subplots_adjust(wspace=1.3, bottom=0.2)
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
    
    # sequential least squares
    lse_position, lse_velocity = sequential_least_squares_observer(estimates_data)
    axpos.plot(estimates_data['time'], lse_position, 'red')
    axvel.plot(estimates_data['time'], lse_velocity, 'red')
    
    
    # UKF
    
    if 1:
        print 'starting ukf'
        ukf = ukf_filter([5,0], [[1,0], [0,1]])
        ukf.iterate_update(estimates_data['filteredopticflow'], estimates_data['control'])
        axpos.plot(estimates_data['time'], ukf.xhistory[:,0], color='purple')
        axvel.plot(estimates_data['time'], ukf.xhistory[:,1], color='purple')
        axof.plot(estimates_data['time'], ukf.xhistory[:,1]/ukf.xhistory[:,0], '*', color='purple', markeredgecolor='none')
        
    # formatting
    axpos.set_ylim(-15,25)
    axof.set_ylim(rdes*(1+0.2), 0)
    axofdot.set_ylim(-10,10)
    
    for ax in axes:
        #xticks = [0,5,10,15,20,25]
        xticks = ax.get_xticks()
        new_xticks = np.linspace(xticks[0], xticks[-1], 3)
        fpl.adjust_spines(ax, ['left', 'bottom'], xticks=new_xticks)
        xtickstrs = [str(tick)[0:4] for tick in new_xticks]
        ax.set_xticklabels(xtickstrs)
        
    axpos.set_xlabel('time, sec')
    
    axpos.set_ylabel('position, m')
    axvel.set_ylabel('velocity, m/s')
    axof.set_ylabel('optic flow, 1/s')
    axofdot.set_ylabel('d/dt optic flow, 1/s^2')
    axctrl.set_ylabel('control, m/s^2')    

    plt.show()
    
    
    
def plot_raw_stepper(stepper_data):
    fig = plt.figure()
    axpos = fig.add_subplot(121)
    axvel = fig.add_subplot(122)
    
    axpos.plot(stepper_data['time'], stepper_data['position'])
    axvel.plot(stepper_data['time'], stepper_data['velocity'])
    
    plt.show()
    
    
def get_optic_flow_noise_estimate(estimates_data, stepper_data):
    stepper_optic_flow = stepper_data['velocity']/stepper_data['position']
    difference = estimates_data['opticflow'] - stepper_optic_flow
    variance = np.mean(difference[200:-200]**2)
    return variance**(1/2.)
    
    
def export_to_matlab(estimates_data, stepper_data):
    estimates_data, stepper_data = interpolate_to_same_time(estimates_data, stepper_data)    
    
    time = estimates_data['time']
    control = estimates_data['control']
    optic_flow = estimates_data['filteredopticflow'] 
    position = stepper_data['position']
    velocity = stepper_data['velocity']
    
    mdict = {'time': time, 'control': control, 'optic_flow': optic_flow, 'position': position, 'velocity': velocity}
    
    scipy.io.savemat('data.mat', mdict=mdict)
    
    
    
if __name__=='__main__':

    parser = OptionParser()
    parser.add_option("--f", type="str", dest="directory", default='',
                        help="directory name that has the estimates and stepper file in it")
    parser.add_option("--estimates", type="str", dest="estimates", default='',
                        help="filename to plot")
    parser.add_option("--stepper", type="str", dest="stepper", default='',
                        help="filename to plot") 
    parser.add_option("--action", type="str", dest="action", default='analyze',
                        help='what kind of plot')
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
        
    if options.action == 'analyze':
        plot(estimates_data, stepper_data)
    elif options.action == 'raw':
        plot_raw_stepper(stepper_data)
    elif options.action == 'mat':
        export_to_matlab(estimates_data, stepper_data)
        
        
