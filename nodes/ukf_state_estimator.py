#!/usr/bin/env python
from __future__ import division
import roslib
roslib.load_manifest('dyneye')
import sys, os
import rospy
import numpy as np
import pickle
import time

from pykalman import AdditiveUnscentedKalmanFilter

from dyneye.srv import *
from std_msgs.msg import *

from optparse import OptionParser

import matplotlib.pyplot as plt

class StateEstimator:
    def __init__(self):
    
        # states: distance, velocity, acceleration, optic flow, d/dt optic flow
        
        ## noise estimates
        n_dim_state = 5
        n_dim_obs = 4

        transition_covariance = np.eye(n_dim_state)*1
        transition_covariance[0,0] = 1e-4
        transition_covariance[1,1] = 1e-4
        transition_covariance[2,2] = 1e-4
        transition_covariance[3,3] = 1e-4
        transition_covariance[4,4] = 1e-4

        observation_covariance = np.eye(n_dim_obs)*1
        observation_covariance[0,0] = 1e-2
        observation_covariance[1,1] = 1e-2
        observation_covariance[2,2] = 1
        observation_covariance[3,3] = 1e-2

        self.xhat = np.array([5, -.1, 0, 0, 0])
        self.covariance = np.eye(n_dim_state)*1e4
        self.covariance[0,0] = 1e4
        self.covariance[1,1] = 1e4

        self.ukf = AdditiveUnscentedKalmanFilter(self.f_ukf_dynamics, self.g_ukf_observer, transition_covariance=transition_covariance, observation_covariance=observation_covariance, initial_state_mean=self.xhat, initial_state_covariance=self.covariance)


        #####################################################################
        # Controller
        self.gain_kp = 6
        self.gain_kd = 20
        self.r_desired = -.1
        self.last_time = None
        self.fps = 50
        self.dt = 1/self.fps
                      
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
        self.data.setdefault('covariance', [])
        rospy.on_shutdown(self.on_shutdown)
        
        #####################################################################
        # ROS stuff
        self.publish_on_estimate = True
        self.use_optimal_kalman_gain = True
        self.service_estimator = rospy.Service('state_estimator_service', StateEstimatorSrv, self.get_estimate)
        self.service_controller = rospy.Service('controller_service', ControllerSrv, self.get_control)
        self.control_publisher = rospy.Publisher('dyneye_control', Float32)
        rospy.Subscriber("camnode/lk", Float32, self.optic_flow_callback)
        rospy.spin()
        
    def f_ukf_dynamics(self, state):
        x,xdot,xddot,r,of = state
        # u is control input, of is optic flow, for smoothing (decoupled from the rest)
        dt = self.dt
                
        new_x = x + xdot*dt
        new_xdot = xdot + xddot*dt #r*new_x
        new_r = new_xdot/new_x
        new_xddot = xddot
        new_of = of
        
        new_state = [new_x, new_xdot, new_xddot, new_r, of]
        return new_state
    
    def g_ukf_observer(self, state):
        return [state[3], state[4], state[0], state[2]]
    
    def optic_flow_callback(self, data):
        estimate = self.get_estimate(data.data)
        
    def on_shutdown(self):
        f = open(self.savepath, 'w')
        pickle.dump(self.data, f)
        f.close()
        
    def get_control(self, query=None):
        a_control = -1*self.gain_kp*(self.xhat[4] - self.r_desired)
        self.u = a_control
        return a_control
        
    def get_estimate(self, r):
        a_control = self.get_control()
        t = rospy.get_time()
        if self.last_time is not None:
            self.dt = t - self.last_time
        else:
            self.dt = 0.02
        
        if len(self.data['states']) > 0:
            rarr = np.array(self.data['states'])[:,4]
            f0 = np.max([0,len(rarr)-30])
            meanr = np.mean(rarr[f0:])
            if np.abs(meanr-self.r_desired) < np.abs(self.r_desired)*.25:
                xest = self.xhat[2]*(self.fps/100.) / self.r_desired**2
            else:
                xest = self.xhat[0] + self.xhat[1]*self.dt
        else:
            xest = self.xhat[0]
        
        observation = np.array([r, r, xest, self.u])
        new_xhat, new_cov = self.ukf.filter_update(self.xhat, self.covariance, observation)
        
        self.xhat = new_xhat
        self.covariance = new_cov
        
        # save data
        self.data['time'].append(t)
        self.data['states'].append(self.xhat)
        self.data['observations'].append(observation)
        self.data['control'].append(a_control)
        self.data['covariance'].append(self.covariance.diagonal().tolist()[0])
        self.last_time = t
        
        if self.publish_on_estimate:
            self.control_publisher.publish(self.get_control())
        
        return [self.xhat[i] for i in range(len(self.xhat))]
        
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
    
    
    
    
    
    
