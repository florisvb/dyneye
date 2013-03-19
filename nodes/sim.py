#!/usr/bin/env python
import roslib
roslib.load_manifest('dyneye')
import sys, os
import rospy
import numpy as np
import pickle
import time

import matplotlib.pyplot as plt

import scipy.stats

from dyneye.srv import *

class Simulation:
    def __init__(self):

        self.states = [[1,0,0,0]]
        self.control = [0]
        self.observations = [0]
        
        self.control_noise = 0.00001
        self.observation_noise = 0.1
        
        self.last_time = None
        
        # ROS stuff
        rospy.wait_for_service('state_estimator_service')
        rospy.wait_for_service('controller_service')
        rospy.on_shutdown(self.on_shutdown)
        
        self.estimator = rospy.ServiceProxy('state_estimator_service', StateEstimatorSrv)
        self.controller = rospy.ServiceProxy('controller_service', ControllerSrv)
        
    def step(self, control):
        t = rospy.get_time()
        if self.last_time is not None:
            dt = t - self.last_time
        else:
            dt = 0.01
        
        d,v,a,r = self.states[-1]
    
        d_ = d + v*dt
        v_ = v + a*dt
        a_ = control + scipy.stats.distributions.norm.rvs(0,self.control_noise)
        r_ = v_/d_
    
        new_states = [d_, v_, a_, r_]
        self.states.append(new_states)
        
        self.observations.append(r_)
        self.last_time = t
        
        delay = -1
        try:
            r_obs = self.observations[delay]
        except:
            r_obs = r_
        observations_noisy = r_obs + scipy.stats.distributions.norm.rvs(0,self.observation_noise)
        print d
        
        return observations_noisy
    
    def run_sim(self):
        
        r = rospy.Rate(60) # 10hz
        pos = 10
        while pos > 0.01:
            control = self.controller(1)   
            observations = self.step(control.val)
            state_estimates = self.estimator(observations)
            pos = self.states[-1][0]
            r.sleep()
        
        self.plot()
        
    def on_shutdown(self):
        self.plot()
        
    def plot(self):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        states = np.array(self.states)
        ax.plot(states[:,0])
        
        plt.show()
        
if __name__=='__main__':
    rospy.init_node('simulator')
    simulation = Simulation()
    simulation.run_sim()





    
