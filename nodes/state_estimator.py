#!/usr/bin/env python
import roslib
roslib.load_manifest('dyneye')
import sys, os
import rospy
import numpy as np
import pickle
import time

from dyneye.srv import *

from optparse import OptionParser

import matplotlib.pyplot as plt

class StateEstimator:
    def __init__(self):
    
        # states: distance, velocity, acceleration, optic flow, d/dt optic flow
        
        ## noise estimates
        # process
        self.Q_k = np.eye(5)*10
        self.Q_k[4,4] = 1e-10
        # observer
        self.R_k = np.eye(5)
        self.R_k[0,0] = 100000
        self.R_k[1,1] = 1000000
        self.R_k[2,2] = 1
        self.R_k[3,3] = 10000
        self.R_k[4,4] = 10000000
        
        ## dynamics models
        self.H_k = np.eye(5)
        self.H_k[4,4] = 0
        dt = 0.01
        self.F = np.matrix([[1,dt,0,0,0],
                            [0,1,dt,0,0],
                            [0,0,1,0,0],
                            [0,0,0,1,dt],
                            [0,0,0,0,1]])
                            
        ## initialization
        self.xhat = np.matrix([10,-1,0,-.1,0]).T
        self.cov = np.eye(5)*10
        self.r_desired = -0.1
        self.last_time = None
        
        self.use_optimal_kalman_gain = True
        
        #####################################################################
        # Controller
        self.gain_kp = 1
        self.gain_kd = 1e3
        
        #####################################################################
        # Save states
        fname = 'state_estimates_' + time.strftime("%Y%m%d_%s") + '.pickle'
        path = '~/'
        self.savepath = os.path.join(os.path.expanduser(path), fname)
        self.data = {}
        self.data.setdefault('time', [])
        self.data.setdefault('states', [])
        self.data.setdefault('observations', [])
        self.data.setdefault('control', [])
        rospy.on_shutdown(self.on_shutdown)
        
        #####################################################################
        # ROS stuff
        self.service_estimator = rospy.Service('state_estimator_service', StateEstimatorSrv, self.get_estimate)
        self.service_controller = rospy.Service('controller_service', ControllerSrv, self.get_control)
        rospy.spin()
        
    def on_shutdown(self):
        f = open(self.savepath, 'w')
        pickle.dump(self.data, f)
        f.close()
        
    def get_control(self, query=None):
        a_control = -1*self.gain_kp*(self.xhat[3,0] - self.r_desired) -1*self.gain_kd*(self.xhat[4,0] - 0) 
        return a_control
        
    def get_estimate(self, observations):
        a_control = self.get_control()
        t = rospy.get_time()
        if self.last_time is not None:
            dt = t - self.last_time
        else:
            dt = 0.01

        ## predict
        d_m = self.xhat[0,0] + self.xhat[1,0]*dt
        v_m = self.xhat[1,0] + self.xhat[2,0]*dt
        a_m = a_control
        r_m = self.xhat[3,0] + self.xhat[4,0]*dt #desired_r
        dr_m = 0
        xhat_k_k0 = np.matrix([d_m, v_m, a_m, r_m, dr_m]).T
        
        self.F[0,1] = dt
        self.F[1,2] = dt
        self.F[3,4] = dt
                        
        P_k_k0 = self.F*self.cov*self.F.T + self.Q_k

        ## observer
        a_o = a_control
        r_o = observations.r
        dr_o = (a_m/d_m - self.xhat[3,0]**2)
        d_o = self.xhat[2,0]/(self.xhat[3,0]**2 + self.xhat[4,0])
        v_o = self.xhat[3,0]*self.xhat[0,0]
        x_obs = np.matrix([d_o, v_o, a_o, r_o, dr_o]).T
        innovation = x_obs - xhat_k_k0
        
        ## update
        S_k = self.H_k*P_k_k0*self.H_k.T + self.R_k
        K_k_opt = P_k_k0*self.H_k.T*np.linalg.inv(S_k)
        if self.use_optimal_kalman_gain:
            K_k = K_k_opt
             
        xhat_k_k = xhat_k_k0 + K_k*innovation
        I = np.eye(5)
        if self.use_optimal_kalman_gain:
            P_k_k = (np.eye(5)-K_k*self.H_k)*P_k_k0
        else:
            P_k_k = (I-K_k*self.H_k)*P_k_k0*(I-K_k*self.H_k).T + K_k*R_k*K_k.T        

        self.xhat = xhat_k_k
        self.cov = P_k_k
        
        # save data
        self.data['time'].append(t)
        self.data['states'].append(self.xhat)
        self.data['observations'].append(observations.r)
        self.data['control'].append(a_control)
        self.last_time = t
        
        return [self.xhat[i,0] for i in range(len(self.xhat))]
        
def plot_estimates(filename):
    
    f = open(filename)
    data = pickle.load(f)
    f.close()
    
    fig = plt.figure()
    ax_pos = fig.add_subplot(1,4,1) 
    ax_vel = fig.add_subplot(1,4,2) 
    ax_optic_flow = fig.add_subplot(1,4,3)    
        
    t = np.array(data['time'])
    states = np.array(data['states'])

    ax_pos.plot(t, states[:,0])
    ax_vel.plot(t, states[:,1])
    ax_optic_flow.plot(t, states[:,3])
    
    plt.show()

if __name__=='__main__':

    parser = OptionParser()
    parser.add_option("--filename", type="str", dest="filename", default='',
                        help="filename to plot")
    (options, args) = parser.parse_args()

    
    if len(options.filename) == 0:
        rospy.init_node('state_estimator')
        state_estimator = StateEstimator()    
    else:
        print 'plotting: ', options.filename
        plot_estimates(options.filename)
    
    
    
    
    
    