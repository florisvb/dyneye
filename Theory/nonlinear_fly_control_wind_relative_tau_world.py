import Lie_Algebra as LA

import sympy as sp # this is a symbolic python package
sin = sp.sin
cos = sp.cos

import copy # this is for copying memory

import numpy as np

'''
v, theta, vw, thetadot, vdot = sp.symbols('v, theta, vw, thetadot, vdot')
eqn = 0.5*((v*cos(theta)-vw)**2+(v*sin(theta))**2)**(-.5)*(2*(v*cos(theta)-vw)*(vdot*cos(theta)-v*sin(theta)*thetadot)+2*(v*sin(theta))*(vdot*sin(theta)+v*cos(theta)*thetadot))
'''

'''
# states
vgx, vgy
vwx, vwy

# observations
airspeed: sqrt( (vgx-vwx)**2 + (vgy-vwy)**2 )
direction of flight: vgx/sqrt(vgx**2+vgy**2), vgy/sqrt(vgx**2+vgy**2)

'''

#x,y,z,vgx,vgy,vgz,vwx,vwy,vwz = sp.symbols('x,y,z,vgx,vgy,vgz,vwx,vwy,vwz')
#x,y,vgangle,vgmag,vwangle,vwmag = sp.symbols('x,y,vgangle,vgmag,vwangle,vwmag')
r, v, airspeed, theta, thetadot, vw = sp.symbols('r, v, airspeed, theta, thetadot, vw')


#dx = [x,y,z,vgx,vgy,vgz,vwx,vwy,vwz]
dx = [r, v, airspeed, theta, thetadot, vw]

if 1:
    f0 = sp.Matrix([-r**2, 0, 0, thetadot, 0, 0])
    f1 = sp.Matrix([0,0,v*vw*(v**2*sin(theta)**2 + v**2*cos(theta)**2 - 2*v*vw*cos(theta) + vw**2)**-0.5*sin(theta),0,1,0]) # thetadot control
    f2 = sp.Matrix([r/v**2,1,thetadot*v*vw*(v**2*sin(theta)**2 + v**2*cos(theta)**2 - 2*v*vw*cos(theta) + vw**2)**-0.5*sin(theta),0,0,0]) # acceleration
    
h = sp.Matrix([r, airspeed, thetadot])


loa = LA.little_observability_algebra(f0,h,dx,0)
dloa = loa.jacobian(dx)
rank = LA.column_rank(dloa)
print 'little obs. alg.'
print 'rank of d_loa: ', rank
print 'num states: ', len(dx)

# now add some terms from big observability matrix:
terms_to_add = []

if 0:
    f0f1 = LA.lie_bracket(f0,f1,dx)
    Lf0f1h = LA.directional_derivative(f0f1,h,dx)
    terms_to_add.extend(Lf0f1h)

if 1:
    f0f2 = LA.lie_bracket(f0,f2,dx)
    Lf0f2h = LA.directional_derivative(f0f2,h,dx)
    terms_to_add.extend(Lf0f2h)

f1f2 = LA.lie_bracket(f1,f2,dx)
#f1f2f0 = LA.lie_bracket(f1f2,f0,dx)

if 0:
    Lf0f2f1h = LA.directional_derivative(f0f2f1,h,dx)
    terms_to_add.extend(Lf0f2f1h)
    

for term in terms_to_add:
    if np.sum(np.abs(term)) != 0: # remove terms that are equal to zero
        loa = loa.col_join(sp.Matrix([term]))
dloa = loa.jacobian(dx)
rank = LA.column_rank(dloa)
print 'big obs. alg.'
print 'rank of d_loa: ', rank
print 'num states: ', len(dx)

        
LA_original, LA_complete, LA_rejected, LA_zero = LA.lie_algebra([f0,f1,f2,f1f2], dx, 0, show=False, ignore_linear_independence=True)
LA_mat = LA.get_matrix_from_LA(LA_complete)
rank = LA.column_rank(LA_mat)
print 'rank of lie algebra: ', rank
print 'num states: ', len(dx)

