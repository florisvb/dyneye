import Lie_Algebra as LA

import sympy as sp # this is a symbolic python package
import copy # this is for copying memory

import numpy as np

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
rf, rp, theta, vg_f, vg_p, wind_speed, thetadot = sp.symbols('rf, rp, theta, vg_f, vg_p, wind_speed, thetadot')


#dx = [x,y,z,vgx,vgy,vgz,vwx,vwy,vwz]
dx = [rf, rp, theta, vg_f, vg_p, wind_speed, thetadot]

# states: velocity in direction of travel, and perpendicular to direction of travel, angle relative to wind, wind speed

if 0:
    f0 = sp.Matrix([vgx,vgy,vgz,0,0,0,0,0,0])
    f1 = sp.Matrix([0,0,0,1,0,0,0,0,0])
    f2 = sp.Matrix([0,0,0,0,1,0,0,0,0])
    f3 = sp.Matrix([0,0,0,0,0,1,0,0,0])
if 1:
    f0 = sp.Matrix([-rf**2, -rp**2, thetadot, 0, 0, 0,0])
    f1 = sp.Matrix([0,0,0,0,0,0,1])
    f2 = sp.Matrix([rf/vg_f,0,0,1,0,0,0])
    f3 = sp.Matrix([0,rp/vg_p,0,0,1,0,0])
    
#h = sp.Matrix([vg_f*sp.cos(alpha)/sp.sqrt(x**2+y**2), vg_p*sp.sin(alpha)/sp.sqrt(x**2+y**2), sp.sqrt((vg_f*sp.cos(alpha)+vg_p*sp.cos(alpha)-wind_speed)**2 + (vg_f*sp.sin(alpha)+vg_p*sp.sin(alpha))**2)])
#h = sp.Matrix([vg_f*sp.cos(alpha)/x, vg_p*sp.sin(alpha)/y, vg_f*sp.cos(alpha)-wind_speed, vg_p*sp.sin(alpha)])
h = sp.Matrix([rf, rp, sp.sqrt((vg_f*sp.cos(theta)+vg_p*sp.cos(theta)-wind_speed)**2 + (vg_f*sp.sin(theta)+vg_p*sp.sin(theta))**2), thetadot])
#h = sp.Matrix([rf, rp, (vg_f*sp.cos(theta)+vg_p*sp.cos(theta)-wind_speed), (vg_f*sp.sin(theta)+vg_p*sp.sin(theta)), thetadot])


loa = LA.little_observability_algebra(f0,h,dx,0)
dloa = loa.jacobian(dx)
rank = LA.column_rank(dloa)
print 'little obs. alg.'
print 'rank of d_loa: ', rank
print 'num states: ', len(dx)

# now add some terms from big observability matrix:
terms_to_add = []

if 1:
    f0f1 = LA.lie_bracket(f0,f1,dx)
    Lf0f1h = LA.directional_derivative(f0f1,h,dx)
    terms_to_add.extend(Lf0f1h)

if 1:
    f0f2 = LA.lie_bracket(f0,f2,dx)
    Lf0f2h = LA.directional_derivative(f0f2,h,dx)
    terms_to_add.extend(Lf0f2h)

if 1:
    f0f3 = LA.lie_bracket(f0,f3,dx)
    Lf0f3h = LA.directional_derivative(f0f3,h,dx)
    terms_to_add.extend(Lf0f3h)
    
if 0:
    f0f2f1 = LA.lie_bracket(f0f2,f1,dx)
    Lf0f2f1h = LA.directional_derivative(f0f2f1,h,dx)
    terms_to_add.extend(Lf0f2f1h)
    
if 0:
    f0f3f1 = LA.lie_bracket(f0f3,f1,dx)
    Lf0f3f1h = LA.directional_derivative(f0f3f1,h,dx)
    terms_to_add.extend(Lf0f3f1h)
  
if 0:
    f0f2f3 = LA.lie_bracket(f0f2,f3,dx)
    Lf0f2f3h = LA.directional_derivative(f0f2f3,h,dx)
    terms_to_add.extend(Lf0f2f3h)


for term in terms_to_add:
    if np.sum(np.abs(term)) != 0: # remove terms that are equal to zero
        loa = loa.col_join(sp.Matrix([term]))
dloa = loa.jacobian(dx)
rank = LA.column_rank(dloa)
print 'big obs. alg.'
print 'rank of d_loa: ', rank
print 'num states: ', len(dx)

        
LA_original, LA_complete, LA_rejected, LA_zero = LA.lie_algebra([f0,f1,f2,f3,f0f2f1], dx, 0, show=False, ignore_linear_independence=True)
LA_mat = LA.get_matrix_from_LA(LA_complete)
rank = LA.column_rank(LA_mat)
print 'rank of lie algebra: ', rank
print 'num states: ', len(dx)

