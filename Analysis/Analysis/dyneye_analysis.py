from __future__ import division

import fly_plot_lib.set_params
fly_plot_lib.set_params.pdf()

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
import fly_plot_lib.text as fpltext

from data_fit import least_squares


DIVIDE = True

###########################################################################################################
# SLS estimator
###########################################################################################################  
    
def sequential_least_squares_observer(estimates_data):
    rdes = estimates_data['rdes']
    fps = estimates_data['fps']
    
    LSE = least_squares.LeastSquaresEstimator()
    Hfunc = lambda t: np.exp(rdes*t)
    H = least_squares.H_LTI([Hfunc])
    LSE.initialize_with_guess(np.matrix([5]),np.matrix([1000]))

    # observer    
    measurements = estimates_data['control'] / rdes**2

    err = (estimates_data['filteredopticflow'] - rdes)
    position = []
    velocity = []
    smoothing = False
    for i, measurement in enumerate(measurements):
        t = estimates_data['time'][i]
        Ht = H(t)
        
        #########################################################
        # artificially keep covariance high until optic flow reaches a steady state (ie. when r~=rd)  
        if np.abs(estimates_data['filteredopticflowdot'][i]) > 0.001 or i < 20:
            if smoothing is False:
                LSE.Pk = 20
                #print LSE.Pk
        else:
            smoothing = True
        #########################################################
        
        LSE.update([measurement], Ht)
        position.append(LSE.xhat[0,0]*Ht[0,0])
        velocity.append(LSE.xhat[0,0]*rdes*Ht[0,0])

    position = np.array(position)
    velocity = np.array(velocity)

    return position, velocity
    
    
###########################################################################################################
# Plotting functions
###########################################################################################################  
    
def plot(estimates_data, stepper_data=None, ukf_data=None):
    
    fig = plt.figure(figsize=(6,2.5))
    fig.subplots_adjust(wspace=1.3, bottom=0.2, right=0.95, left=0.1)
    axpos = fig.add_subplot(141)
    axvel = fig.add_subplot(142)
    axof = fig.add_subplot(143)
    #axofdot = fig.add_subplot(144)
    axctrl = fig.add_subplot(144)    
    axes = [axpos, axvel, axof, axctrl]
    
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
    #axofdot.plot(estimates_data['time'], estimates_data['filteredopticflowdot'], 'red')
    
    # plot raw guess
    guess = estimates_data['control'] / rdes**2
    
    indices = np.where(np.abs(guess)<5)[0]
    axpos.plot(estimates_data['time'][indices], guess[indices], '.', color='green', zorder=-10, markersize=1)
    
    # sequential least squares
    lse_position, lse_velocity = sequential_least_squares_observer(estimates_data)
    axpos.plot(estimates_data['time'], lse_position, 'red')
    axvel.plot(estimates_data['time'], lse_velocity, 'red')
    
    
    # UKF
    if ukf_data is not None:
        axpos.plot(ukf_data['time'], ukf_data['xf'][0,:], color='purple')
        axvel.plot(ukf_data['time'], ukf_data['xf'][1,:], color='purple')
        
    # formatting
    axpos.set_ylim(-2,4)
    axvel.set_ylim(-.5,.5)
    axof.set_ylim(-.15,0)
    axctrl.set_ylim(-.5,.5)
    pos_yticks = [-2,0,2,4]
    vel_yticks = [-.5,0,.5]
    of_yticks = [-.15,-.1,-.05,0]
    ctrl_yticks = [-.5,0,.5]
    #axof.set_ylim(rdes*(1+0.2), 0)
    
    for i, ax in enumerate(axes):
        #xticks = [0,5,10,15,20,25]
        xticks = ax.get_xticks()
        new_xticks = np.linspace(xticks[0], xticks[-1], 2)
        if i==0:
            fpl.adjust_spines(ax, ['left', 'bottom'], xticks=new_xticks, yticks=pos_yticks)
        elif i==1:
            fpl.adjust_spines(ax, ['left', 'bottom'], xticks=new_xticks, yticks=vel_yticks)
        elif i==2:
            fpl.adjust_spines(ax, ['left', 'bottom'], xticks=new_xticks, yticks=of_yticks)
        elif i==3:
            fpl.adjust_spines(ax, ['left', 'bottom'], xticks=new_xticks, yticks=ctrl_yticks)
        else:
            fpl.adjust_spines(ax, ['left', 'bottom'], xticks=new_xticks)
        xtickstrs = [str(tick)[0:4] for tick in new_xticks]
        ax.set_xticklabels(xtickstrs)
        
    axpos.set_xlabel('time, sec')
    
    axpos.set_ylabel('position, m')
    axvel.set_ylabel('velocity, m/s')
    axof.set_ylabel('optic flow, 1/s')
    #axofdot.set_ylabel('d/dt optic flow, 1/s^2')
    axctrl.set_ylabel('control, m/s2')    

    #plt.show()
    
    fpltext.set_fontsize(fig, 8)
    
    figname = 'dyneye_results_' + str(rdes) + '.pdf'
    fig.savefig(figname, format='pdf')
    
def plot_different_rvals(directory_list):
    if type(directory_list) is not list:
        directories = os.listdir(directory_list)
        for d, directory in enumerate(directories):
            directories[d] = os.path.join(directory_list, directory)
    else:
        directories = directory_list
        
    fig = plt.figure(figsize=(6,2.5))
    fig.subplots_adjust(wspace=1.3, bottom=0.2, right=0.95, left=0.1)
    n = len(directories)
    axes = [fig.add_subplot(1,n,i+1) for i in range(n)]
    
    # get sorting
    rdes_list = []
    for i, directory in enumerate(directories):
        stepper, estimates = get_stepper_and_estimates_filenames( directory )
        stepper_data, estimates_data = load_stepper_and_estimates_data(stepper, estimates)
        rdes = estimates_data['rdes']
        rdes_list.append(rdes)
    order = np.argsort(rdes_list)
    print np.array(rdes_list)[order]
    
    a = -1
    for i in order:
        directory = directories[i]
        a += 1
        ax = axes[a]
        
        stepper, estimates = get_stepper_and_estimates_filenames( directory )
        stepper_data, estimates_data = load_stepper_and_estimates_data(stepper, estimates)
        
        fps = estimates_data['fps']
        rdes = estimates_data['rdes']

        print a, rdes, rdes_list[i]
        
        estimates_data, stepper_data = interpolate_to_same_time(estimates_data, stepper_data)    
        noise = get_optic_flow_noise_estimate(estimates_data, stepper_data)
        print 'optic flow noise estimate: ', noise
    
        # plot raw STEPPER data
        ax.plot(stepper_data['time'], stepper_data['position'], 'blue', linewidth=3)
        
        # plot raw guess
        guess = estimates_data['control'] / rdes**2
        
        indices = np.where(np.abs(guess)<5)[0]
        ax.plot(estimates_data['time'][indices], guess[indices], '.', color='green', zorder=-10, markersize=1)
        
        # sequential least squares
        lse_position, lse_velocity = sequential_least_squares_observer(estimates_data)
        ax.plot(estimates_data['time'], lse_position, 'red')
        
        # UKF
        if 0:
            if ukf_data is not None:
                ax.plot(ukf_data['time'], ukf_data['xf'][0,:], color='purple')
                axvel.plot(ukf_data['time'], ukf_data['xf'][1,:], color='purple')
            
        # formatting
        ax.set_ylim(-2,4)
        yticks = [-2,0,2,4]
        ax.set_title(rdes)
    
        xticks = ax.get_xticks()
        new_xticks = np.linspace(xticks[0], xticks[-1], 2)
        
        if a==0:
            spines = ['left', 'bottom']
        else:
            spines = ['bottom']
        
        fpl.adjust_spines(ax, spines, xticks=new_xticks, yticks=yticks)
        
        xtickstrs = [str(tick)[0:4] for tick in new_xticks]
        ax.set_xticklabels(xtickstrs)
    
    
    ax = axes[0]
    ax.set_xlabel('time, sec')
    ax.set_ylabel('position, m')

    fpltext.set_fontsize(fig, 8)
    
    figname = 'dyneye_results_rvariation' + '.pdf'
    fig.savefig(figname, format='pdf')
    
def plot_raw_stepper(stepper_data):
    fig = plt.figure()
    axpos = fig.add_subplot(121)
    axvel = fig.add_subplot(122)
    
    axpos.plot(stepper_data['time'], stepper_data['position'])
    axvel.plot(stepper_data['time'], stepper_data['velocity'])
    
    plt.show()
    
    
###########################################################################################################
# Noise analysis
###########################################################################################################  
    
def get_optic_flow_noise_estimate(estimates_data, stepper_data, return_pdf=False):
    indices = [200,-200]
    stepper_optic_flow = stepper_data['velocity']/stepper_data['position']
    difference = estimates_data['opticflow'] - stepper_optic_flow
    variance = np.mean(difference[indices[0]:indices[-1]]**2)
    if not return_pdf:
        return variance
    else:
        return difference[indices[0]:indices[-1]]
        
###########################################################################################################
# UKF helper functions
###########################################################################################################   
    
def load_ukf_data(filename):
    matdata = scipy.io.loadmat(filename)
    return matdata   

def export_to_matlab(estimates_data, stepper_data):
    estimates_data, stepper_data = interpolate_to_same_time(estimates_data, stepper_data)    
    
    noise = get_optic_flow_noise_estimate(estimates_data, stepper_data)
    print 'optic flow variance: ', noise
    Q = np.matrix(np.eye(2))*1e-10
    
    time = estimates_data['time']
    control = estimates_data['control']
    optic_flow = estimates_data['filteredopticflow'] 
    position = stepper_data['position']
    velocity = stepper_data['velocity']
    
    dt = np.mean(np.diff(time))
    
    mdict = {'time': time, 'control': control, 'optic_flow': optic_flow, 'position': position, 'velocity': velocity, 'R': noise, 'Q': Q, 'dt': dt}
    
    scipy.io.savemat('../Matlab/data.mat', mdict=mdict)

###########################################################################################################
# Load data functions, interpolate to correct sizes
###########################################################################################################    

def get_stepper_and_estimates_filenames(directory):
    files = os.listdir(directory)
    for fname in files:
        if 'stepper' in fname:
            stepper = os.path.join(directory, fname)
        if 'estimates' in fname and 'image' not in fname:
            estimates = os.path.join(directory, fname)
    return stepper, estimates

def load_stepper_and_estimates_data(stepper, estimates):
    f_estimates = open(estimates)
    estimates_data = pickle.load(f_estimates)
    f_estimates.close()
    
    f_stepper = open(stepper)
    stepper_data = pickle.load(f_stepper)
    f_stepper.close()
    
    return stepper_data, estimates_data    

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
            #print key, type(item)
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
    
    if DIVIDE:
        estimates_interp['control'] /= 2.  # fixes bug in the original data collection runs 

    return estimates_interp, stepper_interp

    
if __name__=='__main__':

    parser = OptionParser()
    parser.add_option("--f", type="str", dest="directory", default='',
                        help="directory name that has the estimates and stepper file in it")
    parser.add_option("--d", type="str", dest="directory_list", default='',
                        help="directory list to make plot across different rvals")
    parser.add_option("--estimates", type="str", dest="estimates", default='',
                        help="filename to plot")
    parser.add_option("--stepper", type="str", dest="stepper", default='',
                        help="filename to plot") 
    parser.add_option("--action", type="str", dest="action", default='analyze',
                        help='what kind of plot')
    parser.add_option("--ukf", type="str", dest="ukf", default='',
                        help='destination of ukf mat data')  
                        
    (options, args) = parser.parse_args()
    
    # make plots for all the rvalues (Fig. 5C), no ukf plots
    if len(options.directory_list) > 0:
        plot_different_rvals(options.directory_list) 
        
    else:
        
        # load the data files from the directory (use --f to specify the directory)
        if len(options.directory) > 0:
            options.stepper, options.estimates = get_stepper_and_estimates_filenames(options.directory)
        stepper_data, estimates_data = load_stepper_and_estimates_data(options.stepper, options.estimates)
        
        # load the ukf results from the matlab code, if available
        if len(options.ukf) > 0:
            ukf_data = load_ukf_data(options.ukf)
        else:
            ukf_data = None
            
        # run the primary analysis - Fig. 5 A-D in the paper
        if options.action == 'analyze':
            plot(estimates_data, stepper_data, ukf_data)
            
        # plot the raw data for debugging
        elif options.action == 'raw':
            plot_raw_stepper(stepper_data)
            
        # save data to mat file, needed to run the matlab ukf
        elif options.action == 'mat':
            export_to_matlab(estimates_data, stepper_data)
            
            
            
            
            
            
            
            
            
            
            
