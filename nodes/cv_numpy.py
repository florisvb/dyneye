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

import cv
import numpy as np

def cv2array(im):
    depth2dtype = {
    cv.IPL_DEPTH_8U: 'uint8',
    cv.IPL_DEPTH_8S: 'int8',
    cv.IPL_DEPTH_16U: 'uint16',
    cv.IPL_DEPTH_16S: 'int16',
    cv.IPL_DEPTH_32S: 'int32',
    cv.IPL_DEPTH_32F: 'float32',
    cv.IPL_DEPTH_64F: 'float64',
    }

    arrdtype=im.depth
    a = np.fromstring(
    im.tostring(),
    dtype=depth2dtype[im.depth],
    count=im.width*im.height*im.nChannels)
    a.shape = (im.height,im.width,im.nChannels)
    return a

def array2cv(a):
    dtype2depth = {
    'uint8': cv.IPL_DEPTH_8U,
    'int8': cv.IPL_DEPTH_8S,
    'uint16': cv.IPL_DEPTH_16U,
    'int16': cv.IPL_DEPTH_16S,
    'int32': cv.IPL_DEPTH_32S,
    'float32': cv.IPL_DEPTH_32F,
    'float64': cv.IPL_DEPTH_64F,
    }
    try:
        nChannels = a.shape[2]
    except:
        nChannels = 1
        cv_im = cv.CreateImageHeader((a.shape[1],a.shape[0]),
        dtype2depth[str(a.dtype)], nChannels)
        cv.SetData(cv_im, a.tostring(),a.dtype.itemsize*nChannels*a.shape[1])
    return cv_im
