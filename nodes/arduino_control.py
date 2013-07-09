"""
-----------------------------------------------------------------------
dyneye
Copyright (C) Floris van Breugel, 2013.
  
florisvb@gmail.com

Released under the GNU GPL license, Version 3

This file is part of dyneye.

dyneye is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
    
dyneye is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details.

You should have received a copy of the GNU General Public
License along with dyneye.  If not, see <http://www.gnu.org/licenses/>.

------------------------------------------------------------------------
"""

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
        
        #INCHES_TO_STEPS = 1000./1.079
        self.METERS_TO_STEPS = 1000./0.02740
        self.STEPS_TO_METERS = 1/float(self.METERS_TO_STEPS)
        
        # instantiate stepper class
        print 'Initiating arduino, allow a few seconds'
        self.astep = arduino_stepper.Arduino_Stepper(port=port,timeout=timeout, baudrate=baudrate)
        print
        
        # reset position to zero
        #self.astep.reset_step_counter()
        self.find_home()
        pos = self.get_pos()
        #self.astep.disable_interrupts(0)
        self.last_time = None
        self.acceleration = 0
        self.t_start = None
        
        # Save states
        fname = 'arduino_stepper_data_' + time.strftime("%Y%m%d_%s") + '.pickle'
        path = '~/'
        self.savepath = os.path.join(os.path.expanduser(path), fname)
        self.data = {}
        self.data.setdefault('time', [rospy.get_time()])
        self.data.setdefault('position_steps', [pos])
        self.data.setdefault('velocity_steps', [0])
        self.data.setdefault('position', [pos*self.STEPS_TO_METERS])
        self.data.setdefault('velocity', [0])
        self.data.setdefault('control', [0])
        rospy.on_shutdown(self.on_shutdown)
        
        ## ROS stuff
        rospy.Subscriber("dyneye_control", Float32, self.save_control)
        rospy.spin()
        
        if 0:
            r = rospy.Rate(60) # 10hz
            while not rospy.is_shutdown():
                self.control()
                r.sleep()
                
    def get_pos(self):
        pos = self.astep.get_pos()+4171
        return pos
            
    def find_home(self):
        vel = 2000
        self.astep.set_interrupt_pins(0,1)
        self.astep.enable_interrupts()

        self.astep.set_vel(-1*np.abs(vel))
        interrupted = False
        
        while not interrupted:
            interrupt_0, interrupt_1 = self.astep.get_interrupt_states()
            if interrupt_0 == 1:
                interrupted = True
            time.sleep(.01)
        
        self.astep.reset_step_counter()
        self.astep.disable_interrupts(100)
        self.astep.go_to_pos(45000,5000)
        
    def save_control(self, data):
        self.acceleration = data.data
        self.control()
        
    def control(self):
    
        t = rospy.get_time()
        if self.last_time is not None:
            dt = t - self.last_time
        else:
            dt = 0.01
    
        if 1:
            acceleration = self.acceleration
            
            if 0:
                if self.t_start is None:
                    self.t_start = rospy.get_time()
                r = self.data['velocity'][-1]/self.data['position'][-1]
                r_des = -.1
                acceleration = -3*(r - r_des)# - 1*np.cos((rospy.get_time()-self.t_start)*np.pi)
            
                
            new_vel = self.data['velocity'][-1] + acceleration*dt
            new_vel_in_steps_per_second = int(new_vel*self.METERS_TO_STEPS)
            self.astep.set_vel(new_vel_in_steps_per_second)
        
        else:
            des_vel = -.3
            vel_err = (self.data['velocity'][-1] - des_vel)
            acceleration = -5*vel_err
            new_vel = self.data['velocity'][-1] + acceleration*dt
            new_vel_in_steps_per_second = int(new_vel*self.METERS_TO_STEPS)
            self.astep.set_vel(new_vel_in_steps_per_second)
        
        pos = self.get_pos()
        
        self.data['control'].append(acceleration)
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
    rospy.init_node('arduino_control')
    port = '/dev/ttyACM0'
    baudrate = 38400
    arduino_control = ArduinoControl(port, baudrate)
        
        
        
        
        
        
        
        
