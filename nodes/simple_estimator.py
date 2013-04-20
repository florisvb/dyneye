#!/usr/bin/env python
from __future__ import division
import roslib
roslib.load_manifest('dyneye')
import sys, os
import rospy
import numpy as np
import pickle
import time

import pykalman

from dyneye.srv import *
from std_msgs.msg import *
from save_ros_movie.srv import *

from optparse import OptionParser

import matplotlib.pyplot as plt

class StateEstimator:
    def __init__(self):
    
        # states: distance, velocity, acceleration, optic flow, d/dt optic flow
        
        
        
        #####################################################################
        # Controller
        self.gain_kp = 6
        self.r_desired = -.4
        self.angfreq = 1 # for sin control
        self.last_time = None
        self.fps = 50
        self.dt = 1/self.fps
                      
        #####################################################################
        # Kalman
        
        optic_flow_noise = 3e-5
        R = np.eye(1)*optic_flow_noise
        Q = np.eye(2)*1
        Q[1,1] = 5e-1
        self.kf = pykalman.KalmanFilter(transition_matrices = [[1,self.dt],[0,1]], observation_matrices = [[1,0]], observation_covariance = R, transition_covariance = Q)
        self.filtered_optic_flow = 0
        self.filtered_optic_flow_dot = 0
        self.filtered_covariance = 100*np.eye(2)
                      
        #####################################################################
        # Save states
        basename = 'state_estimates_' + time.strftime("%Y%m%d_%s") + '_' + str(self.r_desired)[1:]
        fname = basename + '.pickle'
        videoname = basename + '_image_frames'
        path = '~/'
        self.savepath = os.path.join(os.path.expanduser(path), fname)
        self.data = {}
        self.data.setdefault('angfreq', self.angfreq)
        self.data.setdefault('rdes', self.r_desired)
        self.data.setdefault('fps', self.fps)
        self.data.setdefault('time', [])
        self.data.setdefault('opticflow', [])
        self.data.setdefault('filteredopticflow', [])
        self.data.setdefault('filteredopticflowdot', [])
        self.data.setdefault('control', [])
        rospy.on_shutdown(self.on_shutdown)
        
        #####################################################################
        # ROS stuff
        self.publish_on_estimate = True
        self.use_optimal_kalman_gain = True
        self.control_publisher = rospy.Publisher('dyneye_control', Float32)
        
        # start saving frames
        save_video_service_name = '/camnode/image_rect/save'
        rospy.wait_for_service(save_video_service_name)
        self.save_video_proxy = rospy.ServiceProxy(save_video_service_name, ToggleSavingVideo) 
        self.save_video_proxy(os.path.join(os.path.expanduser(path), videoname), 1)
        
        rospy.Subscriber("camnode/lk", Float32, self.optic_flow_callback)
        rospy.spin()
        
    def optic_flow_callback(self, data):
        self.optic_flow = data.data
        filtered_optic_flow, self.filtered_covariance = self.kf.filter_update([self.filtered_optic_flow,self.filtered_optic_flow_dot], self.filtered_covariance, observation=self.optic_flow)
        self.filtered_optic_flow = filtered_optic_flow[0]
        self.filtered_optic_flow_dot = filtered_optic_flow[1]
        
        self.get_estimate(self.filtered_optic_flow)
        
    def on_shutdown(self):
        f = open(self.savepath, 'w')
        pickle.dump(self.data, f)
        f.close()
        self.save_video_proxy('', 0)
        rospy.sleep(2)
        
        
    def get_control(self, query=None):
        a_control = -1*self.gain_kp*(self.filtered_optic_flow - self.r_desired)
        
        # sin control
        if 0:
            try:
                tstart = self.tstart
            except:
                self.tstart = time.time()
                tstart = self.tstart
            a_control = 1*np.cos(2*np.pi*self.angfreq*(time.time()-tstart)) -1*self.gain_kp*(self.filtered_optic_flow - -0.1)
        
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
        self.data['filteredopticflow'].append(self.filtered_optic_flow)
        self.data['filteredopticflowdot'].append(self.filtered_optic_flow_dot)
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
    
    
    
    
    
    
