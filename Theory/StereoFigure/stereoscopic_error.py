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

import fly_plot_lib.set_params
fly_plot_lib.set_params.pdf()
import fly_plot_lib.plot as fpl
import fly_plot_lib.text as fpltext

import numpy as np
import matplotlib.pyplot as plt

def stereoscopic_error(d,theta,l):
    return d**2*2*np.tan(theta/2.)/l

fig = plt.figure(figsize=(4,3))
fig.subplots_adjust(wspace=0.3, left=0.2, bottom=0.2, top=0.9)

ax1 = fig.add_subplot(1,1,1)

dhuman = np.linspace(1e-3,10,100)
theta = 0.003*np.pi/180.
l = 65e-3
human_error = stereoscopic_error(dhuman,theta,l)
ax1.plot(np.log(dhuman),human_error,'blue')

dfly = np.linspace(1e-3,.02,100)
theta = 5*np.pi/180.
l = 0.3e-3
fly_error = stereoscopic_error(dfly,theta,l)
ax1.plot(np.log(dfly),fly_error,'red')

human_xticks = np.log(np.array([1e-3,1e-2,1e-1,1,10]))
xticklabels = ['$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$1$', '$10$']

ax1.set_ylim(-.001,.12)
ax1.set_xlim(human_xticks[0], human_xticks[-1])

human_yticks = [0,0.03,0.06,0.09,0.12]

fpl.adjust_spines(ax1, ['left', 'bottom'], xticks=human_xticks, yticks=human_yticks)

ax1.set_xticklabels(xticklabels)

ax1.set_ylabel('Min. error in stereo depth estimate, m')
ax1.set_xlabel('Distance to object, m (log scale)')

fpltext.set_fontsize(fig, 8)

fig.savefig('stereo_error.pdf', format='pdf')
