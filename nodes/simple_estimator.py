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
        

        #####################################################################
        # Controller
        self.gain_kp = 6
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
        self.data.setdefault('opticflow', [])
        self.data.setdefault('control', [])
        rospy.on_shutdown(self.on_shutdown)
        
        #####################################################################
        # ROS stuff
        self.publish_on_estimate = True
        self.use_optimal_kalman_gain = True
        self.control_publisher = rospy.Publisher('dyneye_control', Float32)
        rospy.Subscriber("camnode/lk", Float32, self.optic_flow_callback)
        rospy.spin()
        
    def optic_flow_callback(self, data):
        self.optic_flow = data.data
        self.get_estimate(self.optic_flow)
        
    def on_shutdown(self):
        f = open(self.savepath, 'w')
        pickle.dump(self.data, f)
        f.close()
        
    def get_control(self, query=None):
        a_control = -1*self.gain_kp*(self.optic_flow - self.r_desired)
        self.u = a_control
        return a_control
        
    def get_estimate(self, r):
        a_control = self.get_control()
        t = rospy.get_time()
        if self.last_time is not None:
            self.dt = t - self.last_time
        else:
            self.dt = 0.02
        
        # save data
        self.data['time'].append(t)
        self.data['opticflow'].append(self.optic_flow)
        self.data['control'].append(a_control)
        self.last_time = t
        
        if self.publish_on_estimate:
            self.control_publisher.publish(self.get_control())
        
        
        
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
    
    
    
    
    
    
