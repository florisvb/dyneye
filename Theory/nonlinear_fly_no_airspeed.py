import Lie_Algebra as LA
import sympy as sp # this is a symbolic python package


r,v = sp.symbols('r,v')
p = sp.symbols('p')

dx = [r,v]

f0 = sp.Matrix([-r**2, 0])
f1 = sp.Matrix([r/(v),1])

h = sp.Matrix([p*r])
loa = LA.little_observability_algebra(f0,h,dx,0)

# now add some terms from big observability matrix:
f0f1 = LA.lie_bracket(f0,f1,dx)
Lf0f1h = LA.directional_derivative(f0f1,h,dx)

for term in Lf0f1h:
    if np.sum(np.abs(term)) != 0: # remove terms that are equal to zero
        loa = loa.col_join(sp.Matrix([term]))

# Lf0h = d/dr(r)*rdot = rdot

dloa = loa.jacobian(dx)
rank = LA.column_rank(dloa)
print 'rank of d_loa: ', rank
print 'num states: ', len(dx)

LA_original, LA_complete, LA_rejected, LA_zero = LA.lie_algebra([f0,f1], dx, 3, show=False, ignore_linear_independence=True)
LA_mat = LA.get_matrix_from_LA(LA_complete)
rank = LA.column_rank(LA_mat)
print 'rank of lie algebra: ', rank
print 'num states: ', len(dx)









    





