#!/usr/bin/env python
import roslib
roslib.load_manifest('dyneye')
import sys
import rospy
import cv
from sensor_msgs.msg import Image
from cv_bridge import CvBridge, CvBridgeError
import numpy as np
import cv_numpy

import data_fit

from std_msgs.msg import *

class Optic_Flow_Calculator:

    def __init__(self, camera=None):
        if camera is None:
            camera = 'camera'
        self.image_source = camera + "/image_raw"
        
        # Initialize
        self.bridge = CvBridge()
        self.prev_image = None
        self.msg = None
        
        # Lucas Kanade Publisher
        pub_name = camera + "/optic_flow"
        self.optic_flow_pub = rospy.Publisher(pub_name, Float32)
        
        self.velx = None
        self.vely = None
        
        # line fit
        self.linear_model = data_fit.models.LinearModel()
        self.indices_for_fit = np.arange(0,260).tolist()
        self.alphas = np.arange(0,260)
        self.yindex=5
        self.unit_calibration = 4838.7
        
        # Raw Image Subscriber
        self.image_sub = rospy.Subscriber(self.image_source,Image,self.image_callback)
        
    def image_callback(self,image):
        try:
            curr_image = self.bridge.imgmsg_to_cv(image, "mono8")
            
            # For first loop
            if self.prev_image is not None:
                prev_image = self.prev_image
            else:
                prev_image = curr_image
                
            # get times
            secs = image.header.stamp.secs
            nsecs = image.header.stamp.nsecs
            curr_time = float(secs) + float(nsecs)*1e-9
            dt = curr_time - self.last_time
                
            if self.velx is None:
                self.velx = cv.CreateImage((curr_image.width, curr_image.height), cv.IPL_DEPTH_32F,1)
                self.vely = cv.CreateImage((curr_image.width, curr_image.height), cv.IPL_DEPTH_32F,1)
            
            winSize = (15,1)
            
            cv.CalcOpticalFlowLK(prev_image, curr_image, winSize, self.velx, self.vely)
            velx_np = cv_numpy.cv2array(self.velx)
            
            # fit line
            self.linear_model.fit(velx_np[self.indices_for_fit,self.yindex],inputs=self.alphas)
            
            self.optic_flow_pub.publish(self.linear_model.parameters['slope']*dt*self.unit_calibration)
            
            self.prev_image = curr_image
            self.last_time = curr_time
              
        except CvBridgeError, e:
            print e
  
def main(args):
  optic_flow_calculator = Optic_Flow_Calculator()
  try:
    rospy.spin()
  except KeyboardInterrupt:
    print "Shutting down"
  cv.DestroyAllWindows()

if __name__ == '__main__':
    rospy.init_node('optic_flow', anonymous=True)
    main(sys.argv)
