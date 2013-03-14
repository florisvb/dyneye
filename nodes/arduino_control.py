#!/usr/bin/env python
import roslib
roslib.load_manifest('dyneye')
import sys, os
import rospy
import numpy as np
import pickle
import time

from dyneye.srv import *
from std_msgs.msg import *

import matplotlib.pyplot as plt
import time
import numpy as np
import arduino_stepper.arduino_stepper as arduino_stepper

class ArduinoControl:
    def __init__(self, port, baudrate):
        self.port = port
        self.baudrate = baudrate
        timeout = 1
        
        self.METERS_TO_STEPS = 20
        self.STEPS_TO_METERS = 1/float(self.METERS_TO_STEPS)
        
        # instantiate stepper class
        print 'Initiating arduino, allow a few seconds'
        self.astep = arduino_stepper.Arduino_Stepper(port=port,timeout=timeout, baudrate=baudrate)
        print
        
        # reset position to zero
        self.astep.reset_step_counter()
        self.last_time = None
        
        # Save states
        fname = 'arduino_stepper_data_' + time.strftime("%Y%m%d_%s") + '.pickle'
        path = '~/'
        self.savepath = os.path.join(os.path.expanduser(path), fname)
        self.data = []
        self.data.setdefault('time', [])
        self.data.setdefault('position_steps', [])
        self.data.setdefault('velocity_steps', [])
        self.data.setdefault('position', [])
        self.data.setdefault('velocity', [])
        rospy.on_shutdown(self.on_shutdown)
        
        ## ROS stuff
        rospy.Subscriber("dyneye_control", Float32, self.control)
        rospy.spin()

    def control(self, data):
        acceleration = data.data
        t = rospy.get_time()
        if self.last_time is not None:
            dt = t - self.last_time
        else:
            dt = 0.01
            
        new_vel = self.data['velocity'][-1] + acceleration*dt
        new_vel_in_steps_per_second = int(new_vel*self.METERS_TO_STEPS)
        self.astep.set_vel(new_vel_in_steps_per_second)
        
        pos = self.astep.get_pos()
        
        self.data['time'].append(t)
        self.data['position_steps'].append(pos)
        self.data['position'].append(pos*self.STEPS_TO_METERS)
        self.data['velocity'].append(new_vel)
        self.data['velocity_steps'].append(new_vel_in_steps_per_second)
        
        
    def on_shutdown(self):
        f = open(self.savepath, 'w')
        pickle.dump(self.data, f)
        f.close()
        
        
        
if __name__=='__main__':
    rospy.init_node('state_estimator')
    port = 'dev/tty/ACM0'
    baudrate = 19200
    arduino_control = ArduinoControl(port, baudrate)
        
        
        
        
        
        
        
        
