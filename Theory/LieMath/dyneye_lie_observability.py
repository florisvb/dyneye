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

import Lie_Algebra as LA
import sympy as sp # this is a symbolic python package
import numpy as np

d,v = sp.symbols('d,v')

def r():
    return v/d

dx = [d,v]

f0 = sp.Matrix([v, 0])
f1 = sp.Matrix([0, 1])

h = sp.Matrix([v/d])

print 'h(x): '
print h
print 'rank: ', LA.column_rank(h)
print

f0h = LA.directional_derivative(f0,h,dx)
print 'f0h: '
print f0h
print

loa = LA.little_observability_algebra(f0,h,dx,0)

print 'little observability algebra (h, f0h):'
print loa
print 'rank: ', LA.column_rank(loa)
print 

# now add some terms from big observability matrix:
f0f1 = LA.lie_bracket(f0,f1,dx)
Lf0f1h = LA.directional_derivative(f0f1,h,dx)

for term in Lf0f1h:
    if np.sum(np.abs(term)) != 0: # remove terms that are equal to zero
        loa = loa.col_join(sp.Matrix([term]))

dloa = loa.jacobian(dx)
rank = LA.column_rank(dloa)

print 'big observability algebra: '
print loa
print 'rank: ', rank
print

print 'num states: ', len(dx)

'''
RESULTS:

h(x): 
[v/d]
rank:  1

f0h: 
[-v**2/d**2]

little observability algebra (h, f0h):
[v/d]
rank:  1

big observability algebra: 
[   v/d]
[v/d**2]
rank:  2

num states:  2
'''


