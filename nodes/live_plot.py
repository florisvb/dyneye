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
import roslib; roslib.load_manifest('dyneye')
import rospy

import time

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy
np = numpy

import data_fit

import float_msg_to_numpy as fmtn

from optic_flow.srv import *
from optic_flow.msg import FloatArray

from std_msgs.msg import Float32

class Live_Plotter:
    
    def __init__(self, camera=""):
        self.image_source = 'camnode/image_raw'
        
        rospy.init_node('liveplotter', anonymous=True)
        self.image_sub = rospy.Subscriber(self.image_source,Image,self.image_callback)
        
        self.alphas = np.arange(0,260)
        
        # lines to animate
        self.lucas_kanade, = plt.plot(self.alphas, image_row, '.')
        self.line, = plt.plot(self.alphas, self.linearmodel.get_val(self.alphas), 'red')
                
        self.image_animation = animation.FuncAnimation(self.fig, self.update_line, self.alphas, init_func=self.init_plot, interval=20, blit=True)
        
        plt.xlim(self.alphas[0], self.alphas[-1])
        plt.ylim(-1, 1)
        plt.xlabel('pixels')
        plt.ylabel('pixel values')
        plt.title('test')
        plt.show()
        

    def update_line(self, i):
    
        data = self.return_optic_flow()
        image = fmtn.get_numpy_array(data.shape, data.data)
        image_row = np.mean(image, axis=0)
        
        indices_for_ransac = np.arange(200,400).tolist()
        self.linearmodel.ransac(image_row[indices_for_ransac], inputs=self.alphas[indices_for_ransac], min_data_vals=30, max_iterations=5, threshold=100.0, num_vals_req=30)
        self.ransac_pub.publish(self.linearmodel.parameters['slope'])
                
        self.lucas_kanade.set_data(self.alphas, image_row)
        self.line.set_data(self.alphas, self.linearmodel.get_val(self.alphas))

        return self.lucas_kanade, self.line,
    
        
    def init_plot(self): # required to start with clean slate
        self.lucas_kanade.set_data([],[])
        self.line.set_data([],[])
        return self.lucas_kanade, self.line,
        

if __name__ == '__main__':
    live_plotter = Live_Plotter(camera="camnode")
