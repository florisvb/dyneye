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
import numpy as np
import scipy.integrate.ode
import scipy.linalg
import copy

# class for running the discrete kalman filter algorithm from table 3.1
class DiscreteKalmanFilter(object):
    def __init__(self, x0=None, P0=None, phi=None, gamma=None, H=None, Q=None, R=None, gammaW=None):
    
        # initialize variables, to make class more user friendly
        assert phi is not None, "need to define phi"
        assert H is not None, "need to define H"
    
        self.nstates = phi.shape[0]
        self.nmeasurements = H.shape[0]
        
        if gamma is None:
            gamma = np.matrix(np.zeros([self.nstates, 1]))
        if Q is None:
            Q = np.matrix(np.eye(self.nstates))
        if R is None:
            R = np.matrix(np.eye(self.nmeasurements))
        if gammaW is None:
            gammaW = np.matrix(np.eye(self.nstates))
            
        if x0 is None:
            x0 = np.ones([self.nstates, 1])
        if P0 is None:
            P0 = np.matrix(np.eye(self.nstates))

        # save variables for later use
        self.xhat_apriori = x0
        self.xhat = x0
        self.P_apriori = P0
        self.P = P0
        self.phi = phi
        self.gamma = gamma
        self.Q = Q
        self.R = R
        self.H = H
        self.gammaW = gammaW
        
    # the meat of the algorithm
    def update(self, measurement, control=None):
    
        if control is None:
            control = 0
        
        # get a priori terms from memory
        xhat_apriori = self.xhat_apriori
        P_apriori = self.P_apriori
        phi = self.phi
        gamma = self.gamma
        Q = self.Q
        R = self.R
        H = self.H
        gammaW = self.gammaW
        
        # calculate kalman gain
        K = P_apriori*H.T*(H*P_apriori*H.T+R).I

        # update step
        xhat_aposteriori = xhat_apriori + K*(measurement - H*xhat_apriori)
        I = np.matrix(np.eye(self.nstates))
        P_aposteriori = (I-K*H)*P_apriori
        
        # propagate step
        xhat_new_apriori = phi*xhat_aposteriori + gamma*control
        P_new_apriori = phi*P_aposteriori*phi.T + gammaW*Q*gammaW.T
        
        # save new terms
        self.xhat_apriori = xhat_new_apriori
        self.P_apriori = P_new_apriori
        self.xhat = xhat_aposteriori
        self.P = P_aposteriori
        
        # return the current aposteriori estimates
        return self.xhat, self.P, K
        
        
        
        
        
        
################################################################################################


# class for running the continuous discrete kalman filter algorithm from table 3.7
class ContinuousDiscreteKalmanFilter(object):
    def __init__(self, x0=None, P0=None, f=None, gamma=None, H=None, Q=None, R=None, gammaW=None):
        '''
        f is a function: f(t, x, *fargs). See scipy.integrate.ode
        '''
    
        # initialize variables, to make class more user friendly
        assert f is not None, "need to define f"
        assert H is not None, "need to define H"
    
        self.nstates = x0.shape[0]
        self.nmeasurements = H.shape[0]
        
        if gamma is None:
            gamma = np.matrix(np.zeros([self.nstates, 1]))
        if Q is None:
            Q = np.matrix(np.eye(self.nstates))
        if R is None:
            R = np.matrix(np.eye(self.nmeasurements))
        if gammaW is None:
            gammaW = np.matrix(np.eye(self.nstates))
            
        if x0 is None:
            x0 = np.ones([self.nstates, 1])
        if P0 is None:
            P0 = np.matrix(np.eye(self.nstates))

        # save variables for later use
        self.xhat_apriori = x0
        self.xhat = x0
        self.P_apriori = P0
        self.P = P0
        self.f = f
        self.pf = pf
        self.gamma = gamma
        self.Q = Q
        self.R = R
        self.H = H
        self.gammaW = gammaW
        
        print self.xhat_apriori.shape
        
    # the meat of the algorithm
    def update(self, t, dt, measurement=None, control=None):
    
        if control is None:
            control = 0
        
        # get a priori terms from memory
        xhat_apriori = self.xhat_apriori
        P_apriori = self.P_apriori
        f = self.f
        pf = self.pf
        gamma = self.gamma
        Q = self.Q
        R = self.R
        H = self.H
        gammaW = self.gammaW
        
        # measurement update
        if measurement is not None:
            # calculate kalman gain
            K = P_apriori*H.T*(H*P_apriori*H.T+R).I 
            # update step
            xhat_aposteriori = xhat_apriori + K*(measurement - H*xhat_apriori)
            I = np.matrix(np.eye(self.nstates))
            P_aposteriori = (I-K*H)*P_apriori
        
        # propagate:
        
        # run integration on the dynamics and covariance
        integrator = scipy.integrate.ode(f).set_integrator('dopri5')
        
        initial_value = xhat_aposteriori.T.tolist()[0]
        initial_value.extend(P_aposteriori.reshape(1,self.nstates**2).tolist()[0])
        
        integrator.set_initial_value(initial_value, t)
        while integrator.successful() and integrator.t < t+dt:
            integrator.integrate(integrator.t+dt*1e-1)
        xhat_new_apriori = np.matrix(integrator.y[0:3]).T
        P_new_apriori = np.matrix(integrator.y[3:]).reshape(self.nstates,self.nstates)
        
        # save new terms
        self.xhat_apriori = xhat_new_apriori
        self.P_apriori = P_new_apriori
        self.xhat = xhat_aposteriori
        self.P = P_aposteriori
        
        # return the current aposteriori estimates
        return self.xhat, self.P
        
        
        
        
################################################################################################


# class for running the continuous discrete extended kalman filter algorithm from table 3.9
class ContinuousDiscreteExtendedKalmanFilter(object):
    def __init__(self, x0=None, P0=None, f=None, pf=None, h=None, gamma=None, H=None, Q=None, R=None, gammaW=None):
        '''
        h is a function: h(t, x, *fargs), returns np.matrix of observations
        H is a function: H(t, x), returns np.matrix that is the jacobian of h evaluated at xk
        f is a function: f(t, x, *fargs). See scipy.integrate.ode
        '''
    
        # initialize variables, to make class more user friendly
        assert f is not None, "need to define f"
        assert H is not None, "need to define H"
        assert h is not None, "need to define h"
    
        print x0
    
        self.nstates = x0.shape[0]
        self.nmeasurements = H(0, x0).shape[0]
        
        if gamma is None:
            gamma = np.matrix(np.zeros([self.nstates, 1]))
        if Q is None:
            Q = np.matrix(np.eye(self.nstates))
        if R is None:
            R = np.matrix(np.eye(self.nmeasurements))
        if gammaW is None:
            gammaW = np.matrix(np.eye(self.nstates))
            
        if x0 is None:
            x0 = np.ones([self.nstates, 1])
        if P0 is None:
            P0 = np.matrix(np.eye(self.nstates))

        # save variables for later use
        self.xhat_apriori = x0
        self.xhat = x0
        self.P_apriori = P0
        self.P = P0
        self.f = f
        self.h = h
        self.gamma = gamma
        self.Q = Q
        self.R = R
        self.H = H
        self.gammaW = gammaW
        
        
    # the meat of the algorithm
    def update(self, t, dt, measurement=None, control=None):
    
        if control is None:
            control = np.matrix(np.zeros([self.gamma.shape[1], self.nstates]))
        
        # get a priori terms from memory
        xhat_apriori = self.xhat_apriori
        P_apriori = self.P_apriori
        f = self.f
        h = self.h
        gamma = self.gamma
        Q = self.Q
        R = self.R
        gammaW = self.gammaW

        H = self.H(t, xhat_apriori)
        
        # measurement update
        if measurement is not None:
            # calculate kalman gain
            K = P_apriori*H.T*(H*P_apriori*H.T+R).I 
            # update step
            xhat_aposteriori = xhat_apriori + K*(measurement - h(t,xhat_apriori))
            I = np.matrix(np.eye(self.nstates))
            P_aposteriori = (I-K*H)*P_apriori
        
        # propagate:
        
        # integrate dynamics and covariance
        integrator = scipy.integrate.ode(f).set_integrator('dopri5')
        
        initial_value = xhat_aposteriori.T.tolist()[0]
        initial_value.extend(P_aposteriori.reshape(1,self.nstates**2).tolist()[0])
        
        integrator.set_initial_value(initial_value, t)
        while integrator.successful() and integrator.t < t+dt:
            integrator.integrate(integrator.t+dt*1e-1)
        xhat_new_apriori = np.matrix(integrator.y[0:self.nstates]).T
        P_new_apriori = np.matrix(integrator.y[self.nstates:]).reshape(self.nstates,self.nstates)
        
        print '***: ', xhat_new_apriori
        
        # save new terms
        self.xhat_apriori = xhat_new_apriori
        self.P_apriori = P_new_apriori
        self.xhat = xhat_aposteriori
        self.P = P_aposteriori
        
        # return the current aposteriori estimates
        return self.xhat, self.P        
        
        
###################################################################################


# class for running the discrete unscented kalman filter algorithm
class UnscentedKalmanFilter(object):
    def __init__(self, x0=None, P0=None, f=None, h=None, Q=None, R=None):
        '''
        h is a function: h(t, x, *fargs), returns np.matrix of observations
        f is a function: f(t, x, *fargs). See scipy.integrate.ode
        '''
    
        # initialize variables, to make class more user friendly
        assert f is not None, "need to define f"
        assert h is not None, "need to define h"
    
        self.nstates = x0.shape[0]
        self.nmeasurements = h(x0, 0, 0).shape[0]
        
        if x0 is None:
            x0 = np.ones([self.nstates, 1])
        if P0 is None:
            P0 = np.matrix(np.eye(self.nstates))

        # save variables for later use
        self.xhat_apriori = x0
        self.xhat = x0
        self.P_apriori = P0
        self.P = P0
        self.f = f
        self.h = h
        
        if Q is None:
            Q = np.matrix(np.eye(self.nstates))
        if R is None:
            R = np.matrix(np.eye(self.nmeasurements))
            
        self.Q = Q
        self.R = R
    
    def get_sigma_points(self, xhat, P):
        chol = np.matrix(scipy.linalg.cholesky(self.nstates*P))
        #chol = self.nstates*self.P_apriori
        x_squiggle_p = [chol[:,i] for i in range(chol.shape[0])]
        x_squiggle_n = [-1*chol[:,i] for i in range(chol.shape[0])]
        sigma_pts = [xhat + x_squiggle_p[i] for i in range(len(x_squiggle_p))]
        sigma_pts.extend( [xhat + x_squiggle_n[i] for i in range(len(x_squiggle_n))] )
        return sigma_pts
        
    def update(self, yk, uk, tk):
        
        sigma_pts = self.get_sigma_points(self.xhat, self.P)
        sigma_pts_propagated = [self.f(sigma_pt, uk, tk) for sigma_pt in sigma_pts]
        summed_sigma_pts = copy.copy(sigma_pts_propagated[0])
        for i in range(1,len(sigma_pts_propagated)):
            summed_sigma_pts += sigma_pts_propagated[i]
        xhat_apriori = summed_sigma_pts / float(len(sigma_pts_propagated))
        
        # apriori covariance
        P_apriori = (sigma_pts_propagated[0] - xhat_apriori)*(sigma_pts_propagated[0] - xhat_apriori).T
        for i in range(1,len(sigma_pts_propagated)):
            P_apriori += (sigma_pts_propagated[i] - xhat_apriori)*(sigma_pts_propagated[i] - xhat_apriori).T
        P_apriori /= float(len(sigma_pts_propagated))
        P_apriori += self.Q
        
        self.xhat_apriori = xhat_apriori
        self.P_apriori = P_apriori
        
        ### now measurement stuff ###
        sigma_pts = self.get_sigma_points(self.xhat_apriori, self.P_apriori)
        
        transformed_sigma_pts = [self.h(sigma_pt, uk, tk) for sigma_pt in sigma_pts]
        summed_transformed_sigma_pts = copy.copy(transformed_sigma_pts[0])
        for i in range(1,len(transformed_sigma_pts)):
            summed_transformed_sigma_pts += transformed_sigma_pts[i]
        yhat = summed_transformed_sigma_pts / float(len(transformed_sigma_pts))
        
        # measurement covariance
        Py = (transformed_sigma_pts[0] - yhat)*(transformed_sigma_pts[0] - yhat).T
        for i in range(1,len(transformed_sigma_pts)):
            Py += (transformed_sigma_pts[i] - yhat)*(transformed_sigma_pts[i] - yhat).T
        Py /= float(len(sigma_pts_propagated))
        Py += self.R
        
        # crosscovariance
        Pxy = (sigma_pts_propagated[0] - xhat_apriori)*(transformed_sigma_pts[0] - yhat).T
        for i in range(1,len(sigma_pts_propagated)):
            Pxy += (sigma_pts_propagated[i] - xhat_apriori)*(transformed_sigma_pts[i] - yhat).T
        Pxy /= float(len(sigma_pts_propagated))
        
        ### measurement update ###
        K = Pxy*Py.I
        
        xhat_aposteriori = xhat_apriori + K*(yk - yhat)
        P_aposteriori = P_apriori - K*Py*K.T
        
        self.xhat = xhat_aposteriori
        self.P = P_aposteriori
        
        return xhat_aposteriori, P_aposteriori
        
        
        
        
        

        
###################################################################################


# class for running the discrete unscented kalman filter algorithm
class FancyUnscentedKalmanFilter(object):
    def __init__(self, x0=None, P0=None, f=None, h=None, Q=None, R=None):
        '''
        h is a function: h(t, x, *fargs), returns np.matrix of observations
        f is a function: f(t, x, *fargs). See scipy.integrate.ode
        '''
    
        # initialize variables, to make class more user friendly
        assert f is not None, "need to define f"
        assert h is not None, "need to define h"
    
        self.nstates = x0.shape[0]
        self.nmeasurements = h(x0, 0, 0).shape[0]
        
        if x0 is None:
            x0 = np.ones([self.nstates, 1])
        if P0 is None:
            P0 = np.matrix(np.eye(self.nstates))

        # save variables for later use
        self.xhat_apriori = x0
        self.xhat = x0
        self.P_apriori = P0
        self.P = P0
        self.f = f
        self.h = h
        
        # fudge factos:
        self.L = self.nstates
        self.alpha = 1
        self.beta = 2
        self.kappa = 1
        self.lamb = self.alpha**2*(self.L+self.kappa)-self.L 
        self.fudgeN = np.sqrt(self.L+self.lamb)
        
        # weights
        self.W_mean_0 = self.lamb / (self.L+self.lamb)
        self.W_cov_0 = self.lamb / (self.L+self.lamb) + (1-self.alpha**2+self.beta)
        
        self.W_mean_i = 1 / (2*(self.L+self.lamb))
        self.W_cov_i = 1 / (2*(self.L+self.lamb))
        
        if Q is None:
            Q = np.matrix(np.eye(self.nstates))
        if R is None:
            R = np.matrix(np.eye(self.nmeasurements))
            
        self.Q = Q
        self.R = R
    
    def get_sigma_points(self, xhat, P):
        chol = np.matrix(self.fudgeN*scipy.linalg.cholesky(P))
        #chol = self.nstates*self.P_apriori
        x_squiggle_p = [chol[:,i] for i in range(chol.shape[0])]
        x_squiggle_n = [-1*chol[:,i] for i in range(chol.shape[0])]
        sigma_pts = [xhat]
        sigma_pts.extend( [xhat + x_squiggle_p[i] for i in range(len(x_squiggle_p))] )
        sigma_pts.extend( [xhat + x_squiggle_n[i] for i in range(len(x_squiggle_n))] )
        return sigma_pts
        
    def update(self, yk, uk, tk):
        
        sigma_pts = self.get_sigma_points(self.xhat, self.P)
        sigma_pts_propagated = [self.f(sigma_pt, uk, tk) for sigma_pt in sigma_pts]
        summed_sigma_pts = copy.copy(sigma_pts_propagated[0])*self.W_mean_0
        for i in range(1,len(sigma_pts_propagated)):
            summed_sigma_pts += sigma_pts_propagated[i]*self.W_mean_i
        xhat_apriori = summed_sigma_pts# / float(len(sigma_pts_propagated))
        
        # apriori covariance
        P_apriori = ((sigma_pts_propagated[0] - xhat_apriori)*(sigma_pts_propagated[0] - xhat_apriori).T)*self.W_cov_0
        for i in range(1,len(sigma_pts_propagated)):
            P_apriori += ((sigma_pts_propagated[i] - xhat_apriori)*(sigma_pts_propagated[i] - xhat_apriori).T)*self.W_cov_i
        #P_apriori /= float(len(sigma_pts_propagated))
        P_apriori += self.Q
        
        self.xhat_apriori = xhat_apriori
        self.P_apriori = P_apriori
        
        ### now measurement stuff ###
        sigma_pts = self.get_sigma_points(self.xhat_apriori, self.P_apriori)
        
        transformed_sigma_pts = [self.h(sigma_pt, uk, tk) for sigma_pt in sigma_pts]
        summed_transformed_sigma_pts = copy.copy(transformed_sigma_pts[0])*self.W_mean_0
        for i in range(1,len(transformed_sigma_pts)):
            summed_transformed_sigma_pts += transformed_sigma_pts[i]*self.W_mean_i
        yhat = summed_transformed_sigma_pts #/ float(len(transformed_sigma_pts))
        
        # measurement covariance
        Py = ((transformed_sigma_pts[0] - yhat)*(transformed_sigma_pts[0] - yhat).T)*self.W_cov_0
        for i in range(1,len(transformed_sigma_pts)):
            Py += ((transformed_sigma_pts[i] - yhat)*(transformed_sigma_pts[i] - yhat).T)*self.W_cov_i
        #Py /= float(len(sigma_pts_propagated))
        Py += self.R
        
        # crosscovariance
        Pxy = ((sigma_pts_propagated[0] - xhat_apriori)*(transformed_sigma_pts[0] - yhat).T)*self.W_cov_0
        for i in range(1,len(sigma_pts_propagated)):
            Pxy += ((sigma_pts_propagated[i] - xhat_apriori)*(transformed_sigma_pts[i] - yhat).T)*self.W_cov_i
        #Pxy /= float(len(sigma_pts_propagated))
        
        ### measurement update ###
        K = Pxy*Py.I
        
        xhat_aposteriori = xhat_apriori + K*(yk - yhat)
        P_aposteriori = P_apriori - K*Py*K.T
        
        self.xhat = xhat_aposteriori
        self.P = P_aposteriori
        
        return xhat_aposteriori, P_aposteriori

        
        
        
        
        
        
        
