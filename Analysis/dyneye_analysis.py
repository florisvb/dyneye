import numpy as np
import pickle
import matplotlib.pyplot as plt

from optparse import OptionParser

import data_fit

import fly_plot_lib.plot as fpl

def plot(estimates_data, stepper_data=None):
    
    fig = plt.figure(figsize=(10,5))
    fig.subplots_adjust(wspace=1,hspace=1,bottom=0.2)
    ax_pos = fig.add_subplot(1,4,1) 
    ax_vel = fig.add_subplot(1,4,2) 
    ax_optic_flow = fig.add_subplot(1,4,3) 
    ax_control = fig.add_subplot(1,4,4) 
    
    # plot stepper data
    if stepper_data is not None:
        
        l = np.max([len(stepper_data['time']), len(stepper_data['position']), len(stepper_data['control']), len(stepper_data['velocity'])])
        def fix_arr_len(arr):
            while len(arr) < l:
                arr.append(arr[-1])
            return arr
        fix_arr_len(stepper_data['position'])
        fix_arr_len(stepper_data['velocity'])
        fix_arr_len(stepper_data['control'])
        fix_arr_len(stepper_data['time'])
        
        t = np.array(stepper_data['time'])
        t_start = t[0]
        t -= t_start
        
        stepper_pos = np.array(stepper_data['position'])
        stepper_vel = np.array(stepper_data['velocity'])
        stepper_control = np.array(stepper_data['control'])
        
        print t.shape, stepper_pos.shape, stepper_vel.shape, stepper_control.shape
        
        
        ax_pos.plot(t, stepper_pos, 'blue')
        ax_vel.plot(t, stepper_vel, 'blue')
        ax_control.plot(t, stepper_control, 'blue', zorder=10)
        
        ax_optic_flow.plot(t, stepper_vel/stepper_pos, 'blue', zorder=10)
    else:
        t_start = 0
    
    # plot estimates data
    t = np.array(estimates_data['time'])
    t -= t_start
    
    l = len(estimates_data['time'])
    def fix_arr_len(arr):
        while len(arr) < l:
            arr.append(arr[-1])
        return arr
    
    fix_arr_len(estimates_data['states'])
    fix_arr_len(estimates_data['control'])
    fix_arr_len(estimates_data['observations'])
    
    states = np.array(estimates_data['states'])
    control = np.array(estimates_data['control'])
    observations = np.array(estimates_data['observations'])
        
    
    estimates_pos = states[:,0]
    estimates_vel = states[:,1]
    estimates_optic_flow_est = states[:,3]
    estimates_optic_flow_real = observations
    estimates_control = control
    
    #ax_optic_flow.plot(t, estimates_optic_flow_real, 'blue')
    ax_optic_flow.plot(t, estimates_optic_flow_est, 'red')
    
    ax_control.plot(t, estimates_control, 'blue')
    ax_control.plot(t, states[:,2], 'red')
    
    ax_pos.plot(t, estimates_pos, 'red')
    ax_vel.plot(t, estimates_vel, 'red')
    
    
    # title
    ax_pos.set_title('position')
    ax_vel.set_title('velocity')
    ax_optic_flow.set_title('optic flow')
    ax_control.set_title('control signal')
    
    ax_pos.set_ylim(-.1,1.5)
    ax_vel.set_ylim(-1,1)
    ax_control.set_ylim(-1,.3)
    ax_optic_flow.set_ylim(-.3,0)
    
    xticks = [0,10,20,30]
    axes = [ax_pos, ax_vel, ax_optic_flow, ax_control]
    for ax in axes:
        fpl.adjust_spines(ax, ['left', 'bottom'], xticks=xticks)
    
    # label
    ax_pos.set_xlabel('time, seconds')
    ax_pos.set_ylabel('distance, m')
    ax_vel.set_ylabel('velocity, m/s')
    ax_optic_flow.set_ylabel('optic flow, 1/s')
    ax_control.set_ylabel('control: acceleration, m/s2')
    
    
    fig.savefig('dyneye_results.pdf', format='pdf')
    
    # v/d correlation
    
    if stepper_data is not None:
        t_stepper = np.array(stepper_data['time'])
        r_stepper = stepper_vel/(stepper_pos)
    
        t_estimates = np.array(estimates_data['time'])
        r_estimates = np.array(estimates_data['observations'])
        pos_estimates = np.array(estimates_data['states'])[:,0].reshape(len(r_estimates))
    
        t_interp_start = np.max([t_stepper[0], t_estimates[0]])
        t_interp_stop = np.min([t_stepper[-1], t_estimates[-1]])
        t_interp = np.arange(t_interp_start, t_interp_stop, 0.005)
        
        stepper_vel_interp = np.interp(t_interp, t_stepper, stepper_vel)
        stepper_pos_interp = np.interp(t_interp, t_stepper, stepper_pos)
        r_stepper_interp = np.interp(t_interp, t_stepper, r_stepper)
        r_estimates_interp = np.interp(t_interp, t_estimates, r_estimates)
        pos_estimates_interp = np.interp(t_interp, t_estimates, pos_estimates)
        control_interp = np.interp(t_interp, t_stepper, stepper_control)
        
        t_interp -= t_interp[0]

        indices = np.where( (t_interp<20) )[0]
       
        fig = plt.figure()
        
        ax = fig.add_subplot(131)
        ax.plot(t_interp, r_stepper_interp)
        ax.set_ylim(np.min(r_estimates_interp),np.max(r_estimates_interp))
        
        ax = fig.add_subplot(132)
        ax.plot(t_interp, r_estimates_interp)
        
        ax = fig.add_subplot(133)
        #ax.plot(r_estimates_interp, r_stepper_interp, '.', color='black')
        #ax.plot(r_estimates_interp[indices], r_stepper_interp[indices], '.', color='blue')
        
        fpl.scatter(ax, r_estimates_interp[indices], r_stepper_interp[indices], stepper_pos_interp, radius=0.001, use_ellipses=True, colornorm=(-.5,-.2))
        
        ax.set_xlabel('optic flow measurement')
        ax.set_ylabel('actual v/d')
        
        lm = data_fit.models.LinearModel()
        lm.fit(r_stepper_interp[indices], inputs=r_estimates_interp[indices])
        
        if 1:
            x = np.linspace(-1,0, 100)
            y = lm.get_val(x)
            ax.plot(x,y,'red')
            
        print lm.parameters
            
        ax.set_ylim(-1,.2)
    
    
        ##
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(states[:,0], states[:,1], '.')
    
    
    
    plt.show()
    
    
    
    
    
    
    
    
    
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
