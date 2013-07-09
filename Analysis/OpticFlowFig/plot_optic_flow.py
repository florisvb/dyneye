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
import sys
sys.path.append('../nodes')
import cv
import numpy as np
import cv_numpy
from optparse import OptionParser
import data_fit
import pickle

import fly_plot_lib.plot as fpl
import fly_plot_lib.text as fpltext

import matplotlib.pyplot as plt


def draw_box_on_image(image, fname):   
    
    fig = plt.figure(figsize=(2.668,2), dpi=360)
    ax = fig.add_axes([0,0,1,1])
    
    im = cv_numpy.cv2array(image)
    
    ax.imshow(im[:,:,0], cmap='gray')
    
    ax.vlines([292,465], 230, 250, 'red')
    ax.hlines([230, 250], 292, 465, 'red')
    
    ax.set_frame_on(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    fig.savefig(fname, format='pdf')

def plot_optic_flow_sample(image1, image2, dt):
    fig = plt.figure(figsize=(2,2))
    fig.subplots_adjust(left=0.28, bottom=0.25)
    ax = fig.add_subplot(111)


    image1sub = cv.GetSubRect(image1, (0, 230, 659, 20))
    image2sub = cv.GetSubRect(image2, (0, 230, 659, 20))
        
    velx = cv.CreateImage((image1sub.width, image1sub.height), cv.IPL_DEPTH_32F,1)
    vely = cv.CreateImage((image1sub.width, image1sub.height), cv.IPL_DEPTH_32F,1)

    winSize = (15,3)

    cv.CalcOpticalFlowLK(image1sub, image2sub, winSize, velx, vely)
    velx_np = cv_numpy.cv2array(velx)

    velx_mean = np.mean(velx_np, axis=0)

    # fit line
    linear_model = data_fit.models.LinearModel()
    indices_for_fit = np.arange(292,465).tolist()
    alphas = np.arange(0,659)
    linear_model.ransac(velx_mean[indices_for_fit].T, inputs=alphas[indices_for_fit], min_data_vals=30, max_iterations=5, threshold=100.0, num_vals_req=30)
    unit_calibration = 4838.7*0.0069546813193345491*1.1*(1/50.)
    optic_flow = -1*linear_model.parameters['slope']*unit_calibration/dt


    # plot
    
    ax.plot(alphas, velx_mean, '.')
    y = linear_model.get_val(indices_for_fit)
    ax.plot(indices_for_fit, y, 'red', linewidth=2)
    
    ax.set_xlim(0,659)
    ax.set_ylim(-1,1)
    xticks = [0,659]
    yticks = [-1,0,1]
    
    fpl.adjust_spines(ax, ['left', 'bottom'], xticks=xticks, yticks=yticks)
    ax.set_xlabel('camera pixel')
    ax.set_ylabel('optic flow')
    
    fpltext.set_fontsize(fig, 8)
    
    figname = 'opticflowsample.pdf'
    fig.savefig(figname, format='pdf')
    
def load_image(image_filename):
    return cv.LoadImage(image_filename, False)
    
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--image1", type="str", dest="image1", default='',
                        help="filename of image1")
    parser.add_option("--image2", type="str", dest="image2", default='',
                        help="filename of image2")
    parser.add_option("--caminfo", type="str", dest="caminfo", default='',
                        help="pickle file of camera info")
    parser.add_option("--action", type="str", dest="action", default='plot',
                        help="action to perform")
                        
    (options, args) = parser.parse_args()

    if len(options.caminfo) > 0:
        f = open(options.caminfo)
        caminfo = pickle.load(f)
        
        fname1 = options.image1.split('/')[-1]
        frame1 = fname1.split('.')[0]
        
        fname2 = options.image2.split('/')[-1]
        frame2 = fname2.split('.')[0]
        
        time1 = float(caminfo[frame1]['secs']) + float(caminfo[frame1]['nsecs'])*1e-9
        time2 = float(caminfo[frame2]['secs']) + float(caminfo[frame2]['nsecs'])*1e-9
        dt = time2-time1
    
    else:
        dt = 0.02
        print "Guessing dt = 0.02!"
        
    image1 = load_image(options.image1)
    image2 = load_image(options.image2)
    
    if options.action == "plot":
        plot_optic_flow_sample(image1, image2, dt)
    elif options.action == "draw":
        draw_box_on_image(image1, 'image_box1.pdf')
        draw_box_on_image(image2, 'image_box2.pdf')
    
    
    
    
    
    
    
