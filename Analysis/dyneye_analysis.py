import numpy as np
import pickle
import matplotlib.pyplot as plt

from optparse import OptionParser

def plot(estimates_data, stepper_data):
    
    fig = plt.figure()
    ax_pos = fig.add_subplot(1,4,1) 
    ax_vel = fig.add_subplot(1,4,2) 
    ax_optic_flow = fig.add_subplot(1,4,3) 
    ax_control = fig.add_subplot(1,4,4)    
    
    # plot stepper data
    t = np.array(stepper_data['time'])
    stepper_pos = np.array(stepper_data['position'])
    stepper_vel = np.array(stepper_data['velocity'])
    ax_pos.plot(t, stepper_pos, 'blue')
    ax_vel.plot(t, stepper_vel, 'blue')
    
    # plot estimates data
    t = np.array(estimates_data['time'])
    states = np.array(estimates_data['states'])
    control = np.array(estimates_data['control'])
    observations = np.array(estimates_data['observations'])
    estimates_pos = states[:,0]
    estimates_vel = states[:,1]
    estimates_optic_flow_est = states[:,3]
    estimates_optic_flow_real = observations
    estimates_control = control
    
    ax_pos.plot(t, estimates_pos, 'red')
    ax_vel.plot(t, estimates_vel, 'red')
    ax_optic_flow.plot(t, estimates_optic_flow_est, 'blue')
    ax_optic_flow.plot(t, estimates_optic_flow_real, 'red')
    ax_control.plot(t, estimates_control, 'red')
    
    # title
    ax_pos.set_title('position')
    ax_vel.set_title('velocity')
    ax_optic_flow.set_title('optic flow')
    ax_control.set_title('control signal')
    
    # label
    ax_pos.set_xlabel('time, seconds')
    
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
    
    f_stepper = open(options.stepper)
    stepper_data = pickle.load(f_stepper)
    f_stepper.close()
    
    plot(estimates_data, stepper_data)
