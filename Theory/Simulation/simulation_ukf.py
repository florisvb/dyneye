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

from __future__ import division

# plotting functions
try:
    import fly_plot_lib
    fly_plot_lib.set_params.pdf()
    import fly_plot_lib.plot as fpl
    import fly_plot_lib.text as flytext
    nice_plots = True
except:
    print 'Cannot find FlyPlotLib to make nice plots!'
    nice_plots = False


# plotting library
import matplotlib.pyplot as plt

# matlab-like functions
import numpy as np
import scipy.stats
import scipy.io

# my kalman code
import kalman

dt = 0.02


def generate_synthetic_data(x0, nsamples, control='exp', accel=-0.1):
    
    x = [x0]
    y = []
    u = []
    
    drift = np.matrix([[0,1],[0,0]]) # A matrix
    gamma = np.matrix([[0],[1]])     # B matrix
    
    for k in range(nsamples):
    
        # calculate optic flow as v/d
        optic_flow = x[-1][1,0] / x[-1][0,0]
        
        # options for control
        if control == 'exp':
            optic_flow_des = -.1
            ctrl = -6*(optic_flow - optic_flow_des) # constant optic flow
        elif control == 'constantaccel':
            ctrl = accel
        else:
            ctrl = 0
            
        # propagate dynamics        
        x_new = x[-1] + (drift*x[-1] + gamma*ctrl)*dt
        y_new = optic_flow
        
        # stop simulation when getting very close to zero
        if x_new[0,0] <= 0.001:
            break
        
        # save new states and controls
        x.append(x_new)
        y.append(y_new)
        u.append(ctrl)
        
    return x[1:], y, u

def simple_ukf(x0, P0, Q, R, control='exp'):


    # generate synthetic data
    nsamples = 600
    xs,ys,us = generate_synthetic_data(x0, nsamples, control=control)

    # set initial condition for the filter to be some fraction of the real initial condition
    x0 = x0*0.8

    ###########################################################################################
    # UKF
    ###########################################################################################
    
    # dynamics function for UKF
    def f(x, u, t): 
        drift = np.matrix([[0,1],[0,0]])
        gamma = np.matrix([[0],[1]])
        x_new = x + (drift*x + gamma*u)*.02
        return x_new
    
    # measurement function for UKF
    def h(x, u, t): 
        y = x[1,0] / x[0,0]
        return np.matrix([[y]])
        
    # initialize UKF
    UKF = kalman.UnscentedKalmanFilter(x0=x0, P0=P0, f=f, h=h, Q=Q, R=R)

    t = 0       # time
    xhats = []  # estimates
    Psp = []    # covariance in the [0,0] entry (position)
    Psv = []    # covariance in the [1,1] entry (velocity)
    
    # run through all the observations
    for k in range(len(ys)):
        y = ys[k] # observation
        u = us[k] # control
        
        t += dt     
        xhat, Pk = UKF.update(y, u, t) # get new xhat and covariance from the UKF

        xhats.append(xhat)              # save new xhat
        Psp.append(Pk.diagonal()[0,0])  # save position covariance
        Psv.append(Pk.diagonal()[0,1])  # save veloicity covariance
            
            
    ###########################################################################################
    # Make Plots
    ###########################################################################################
    
    # reformat the data for plotting
    xhats = np.matrix([xhat.T.tolist()[0] for xhat in xhats])
    xs = np.matrix([x.T.tolist()[0] for x in xs])

    # calculate error between estiamtes and true values
    err0 = xhats[:,0] - xs[:,0]
    err1 = xhats[:,1] - xs[:,1]
    
    # plot error and covariance bounds
    fig = plt.figure(figsize=(7,8))
    axp = fig.add_subplot(411)
    axp.plot(xhats[:,0], 'red')
    axp.plot(xs[:,0], 'blue')
    axp.set_ylabel('position')
    
    axv = fig.add_subplot(412)
    axv.plot(xhats[:,1], 'red')
    axv.plot(xs[:,1], 'blue')
    axv.set_ylabel('velocity')
    
    axof = fig.add_subplot(413)
    axof.plot(xhats[:,1]/xhats[:,0], 'red')
    axof.plot(xs[:,1]/xs[:,0], 'blue')
    axof.set_ylim(-.2,.2)
    axof.set_ylabel('optic flow')

    axer = fig.add_subplot(414)
    axer.plot(xhats[:,0] - xs[:,0], 'black')
    
    print xhats[:,0] - xs[:,0]
    
    axer.plot(3*np.sqrt(np.array(Psp)), ':', color='black')
    axer.plot(-3*np.sqrt(np.array(Psp)), ':', color='black')
    axer.set_ylabel('position error')
    axer.set_ylim(-3,3)
    
    # save the figure to disk
    figname = 'test_ukf.png'
    fig.savefig(figname, format='png')
    
    return xhats
    
    
def export_to_matlab(x0, P0, Q, R, control='exp'):

    # generate synthetic data with these dynamics
    nsamples = 600
    xs,ys,us = generate_synthetic_data(x0, nsamples, control=control)
    
    def reformat_data(m):
        a = np.array([m[i].T.tolist()[0] for i in range(len(m))])
        return a
    xs,ys,us = generate_synthetic_data(x0, nsamples)
    xs = reformat_data(xs)
    
    time = np.arange(0,dt*nsamples,dt)
    control = us
    optic_flow = ys 
    position = xs[:,0].tolist()
    velocity = xs[:,1].tolist()
    
    mdict = {'time': time, 'control': control, 'optic_flow': optic_flow, 'position': position, 'velocity': velocity, 'Q': Q, 'R': R, 'dt': dt}
    
    scipy.io.savemat('data.mat', mdict=mdict)
    
    
def numerical_grammian(nsamples=200, control='exp', eps=1e-3, accel=-.1):

    x0 = np.matrix([[10+eps], [-1]])
    xs,ys,us = generate_synthetic_data(x0, nsamples, control, accel)
    ys_d_plus = np.array(ys)
    
    x0 = np.matrix([[10-eps], [-1]])
    xs,ys,us = generate_synthetic_data(x0, nsamples, control, accel)
    ys_d_minus = np.array(ys)
    
    x0 = np.matrix([[10], [-1+eps]])
    xs,ys,us = generate_synthetic_data(x0, nsamples, control, accel)
    ys_v_plus = np.array(ys)
    
    x0 = np.matrix([[10], [-1-eps]])
    xs,ys,us = generate_synthetic_data(x0, nsamples, control, accel)
    ys_v_minus = np.array(ys)
    
    cn = []
    
    for i in range(len(ys)):
        arr = np.matrix([[ys_d_plus[i]-ys_d_minus[i]], [ys_v_plus[i]-ys_v_minus[i]]])  # at tf
        arr0 = np.matrix([[ys_d_plus[0]-ys_d_minus[0]], [ys_v_plus[0]-ys_v_minus[0]]]) # at t0
        P = ((arr*arr.T)-(arr0*arr0.T))*i*0.02/(4*eps**2)
        u,s,v = np.linalg.svd(P)
        cn.append( np.max(s) / np.min(s) )
        
    def reformat_data(m):
        a = np.array([m[i].T.tolist()[0] for i in range(len(m))])
        return a
    x0 = np.matrix([[10], [-1]])
    xs = reformat_data(xs)
    
    return cn, xs[:,0]
    

def plot_numerical_grammian(eps=1e-5):
    fig = plt.figure(figsize=(3,3))
    ax = fig.add_axes([0.2,0.2,0.7,0.7])
    
    cn, distance = numerical_grammian(nsamples=2000, control='exp', eps=eps)
    ax.plot(distance, np.log(cn), 'red')
    
    cn, distance = numerical_grammian(nsamples=2000, control='constantaccel', eps=eps, accel=-0.1)
    print distance.shape, len(cn)
    ax.plot(distance, np.log(cn), 'blue')
    
    cn, distance = numerical_grammian(nsamples=2000, control='constantaccel', eps=eps, accel=-0.01)
    print distance.shape, len(cn)
    ax.plot(distance, np.log(cn), 'lightblue')
    
    cn, distance = numerical_grammian(nsamples=2000, control='none', eps=eps)
    ax.plot(distance, np.log(cn), 'black')
    
    yticks = np.log(np.array([1e0,1e4,1e8,1e12,1e16]))
    fpl.adjust_spines(ax, ['left', 'bottom'], yticks=yticks)
    ax.set_yticklabels(['$10^{0}$','$10^{4}$', '$10^{8}$', '$10^{12}$', '$10^{16}$'])
    
    ax.set_xlabel('distance to target, m')
    ax.set_ylabel('condition number')
    
    flytext.set_fontsize(fig, 8)
    
    fname = 'condition_numbers_' + str(np.log10(eps)) + '.pdf'
    fig.savefig(fname, format='pdf')
    
    
if __name__ == '__main__':
    
    # initial conditions
    P0 = np.matrix(np.eye(2))*10
    x0 = np.matrix([[10], [-1]])
    
    # define Q and R, note there is no actual noise in the simulated system
    Q = np.matrix(np.eye(2))*1e-5
    R = np.matrix(np.eye(1))*1e-5
    
    # choose type of control
    control = 'exp'
    
    simple_ukf(x0, P0, Q, R, control)
    export_to_matlab(x0, P0, Q, R, control)

    plot_numerical_grammian(eps=1e-0)
    plot_numerical_grammian(eps=1e-3)
    plot_numerical_grammian(eps=1e-6)








