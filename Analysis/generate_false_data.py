import scipy.stats
import numpy as np
import pickle

estimates_data = {}
estimates_data.setdefault('rdes', -0.1)
estimates_data.setdefault('fps', 50.)
estimates_data.setdefault('time', [])
estimates_data.setdefault('opticflow', [])
estimates_data.setdefault('filteredopticflow', [])
estimates_data.setdefault('filteredopticflowdot', [])
estimates_data.setdefault('control', [])
        
stepper_data = {}
stepper_data.setdefault('time', [0])
stepper_data.setdefault('position', [20])
stepper_data.setdefault('velocity', [-0.1])
stepper_data.setdefault('control', [-.1])

optic_flow_noise = scipy.stats.norm(0, 0.000006)

def step(stepper_data, estimates_data):
    pos = stepper_data['position'][-1]
    vel = stepper_data['velocity'][-1]
    dt = 1/float(estimates_data['fps'])
    rdes = estimates_data['rdes']
    optic_flow = vel/pos + optic_flow_noise.rvs()
    
    if len(estimates_data['opticflow']) > 0:
        optic_flow_dot = (optic_flow - estimates_data['opticflow'][-1]) / dt
    else:
        optic_flow_dot = optic_flow/dt
    
    control = -6*(optic_flow - rdes)
    new_pos = pos + dt*vel 
    new_vel = vel + dt*control
    new_t = stepper_data['time'][-1] + dt
    
    stepper_data['position'].append(new_pos)
    stepper_data['velocity'].append(new_vel)
    stepper_data['control'].append(control)
    stepper_data['time'].append(new_t)
    
    estimates_data['time'].append(new_t)
    estimates_data['opticflow'].append(optic_flow)
    estimates_data['filteredopticflow'].append(optic_flow)
    estimates_data['filteredopticflowdot'].append(optic_flow_dot)
    estimates_data['control'].append(control)
    
def empty_step(stepper_data, estimates_data):
    stepper_data['position'].append(stepper_data['position'][-1])
    stepper_data['velocity'].append(stepper_data['velocity'][-1])
    stepper_data['control'].append(stepper_data['control'][-1])
    stepper_data['time'].append(stepper_data['time'][-1])
    
    estimates_data['time'].append(estimates_data['time'][-1])
    estimates_data['opticflow'].append(estimates_data['opticflow'][-1])
    estimates_data['filteredopticflow'].append(estimates_data['filteredopticflow'][-1])
    estimates_data['filteredopticflowdot'].append(estimates_data['filteredopticflowdot'][-1])
    estimates_data['control'].append(estimates_data['control'][-1])
    
    
while stepper_data['position'][-1] > 0.01:
    step(stepper_data, estimates_data)
    print stepper_data['position'][-1], stepper_data['velocity'][-1]

for i in range(5):
    empty_step(stepper_data, estimates_data)
    
stepper_data_file = 'simulated_data/simulated_stepper_data.pickle'
f = open(stepper_data_file, 'w')
pickle.dump(stepper_data, f)    

estimates_data_file = 'simulated_data/simulated_estimates_data.pickle'  
f = open(estimates_data_file, 'w')
pickle.dump(estimates_data, f)    


  
