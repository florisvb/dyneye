import sys
sys.path.insert(0,'/home/floris/src/sympy')
import sympy as sp
print sp.__version__
import numpy as np
import inspect
import copy

INF = np.inf

def lie_bracket(f1, f2, dx, simplify=True):
    
    n = len(dx)
    f1dx = sp.matrices.zeros([n,n])
    f2dx = sp.matrices.zeros([n,n])
    
    for i in range(n):
        for j in range(n):
            f1dx[i,j] = sp.diff(f1[i], dx[j])
            f2dx[i,j] = sp.diff(f2[i], dx[j])
            
    lb = (f2dx*f1 - f1dx*f2).expand()
    
    if simplify:
        for i, l in enumerate(lb):
            lb[i] = sp.simplify( l )
    
    return lb
    
def adfn(f0,f1,dx,n):
    i=0
    adi = f1
    while i<n:
       adi = lie_bracket(f0,adi,dx) 
       i += 1
    return adi
    
def ad_matrix(f0,f1,dx,n):
    admat = copy.copy(f1)
    for i in range(1, n+1):
        adi = adfn(f0,f1,dx,i)
        admat = admat.row_join(adi)
    return admat    
    
def delta_ad_matrix(f0, fc, dx, n, x0=None):
    if type(fc) is list:
        l = len(fc)
        fc_mat = fc[0]
        for i in range(1, l):
            fc_mat = fc_mat.row_join(fc[i])
        fc = fc_mat

    delta_ad_n = None
    for c in range(fc.shape[1]):
        admat_f = ad_matrix(f0, fc[:,c], dx, n)
        if delta_ad_n is None:
            delta_ad_n = admat_f
        else:
            delta_ad_n = delta_ad_n.row_join(admat_f)
    
    if x0 is not None:
        return delta_ad_n.subs(x0)
    else:
        return delta_ad_n

def column_rank(delta, method='nullspace', x0=None):

    if method == 'nullspace':
        if x0 is None:
            rank = delta.shape[1] - len(delta.nullspace())
            return rank
        else:
            rank = delta.shape[1] - len(delta.subs(x0).nullspace())
            return rank
        
    else:
        if x0 is not None:
            return len(delta.subs(x0).T.rref()[1])  
        else:
            return len(delta.T.rref()[1])  
        
def is_involutive(delta, dx, x0=None):
    rank_delta = column_rank(delta)
    
    result = True
    
    for i in range(delta.shape[1]):
        for j in range(i+1, delta.shape[1]):
            new_bracket = lie_bracket(delta[:,i], delta[:,j], dx)
            delta_tmp = delta.row_join(new_bracket)
            rank_delta_tmp = column_rank(delta_tmp, x0)
            if rank_delta_tmp > rank_delta:
                result = False
                
    return result

    
# calculate Lf0h n times
def directional_derivative(f0, h, dx, n=1):
    # f0 and h should both be sympy.matrices
    for iteration in range(n):
        hdx = sp.matrices.zeros([h.shape[0],len(dx)])
        for i in range(h.shape[0]):
            for j in range(len(dx)):
                hdx[i,j] = sp.diff(h[i], dx[j])
        f0h = (hdx*f0).expand()
        h = f0h
    return f0h
    
def little_observability_algebra(f0, h, dx, n=1):
    loa = copy.copy(h)
    for iteration in range(n):  
        new_terms = directional_derivative(f0,h,dx,iteration+1)
        for term in new_terms:
            if np.sum(np.abs(term)) != 0: # remove terms that are equal to zero
                loa = loa.col_join(sp.Matrix([term]))
    return loa
    
    
def lie_algebra(F, dx, n_iterations=INF, show=True, ignore_linear_independence=False):
    # F should either be a list of sp.Matrix elements corresponding to the vector fields, or a sp.Matrix where each column is a vector field
    # dx should be a list of the sp.symbols that correspond to the partial derivatives

    # Make a matrix of Vector Fields
    if type(F) is list:
        n = len(F)
        newF = F[0]
        for i in range(1, n):
            newF = newF.row_join(F[i])
        F = newF
    else:
        newF = F
        
    # number of original vector fields
    nF = F.shape[1]
        
    # original strings:
    F_strings = ['f'+str(i) for i in range(nF)]

    # initialize lie algebra dictionaries
    LA_zero = {} # lie brackets that equal zero
    LA_rejected = {} # lie brackets that are linearly dependent on other terms
    LA_complete = {} # the complete lie algebra
    LA_original = {} # the original vector fields
    tested_lie_brackets = [] # keep track of all the computed lie brackets in string form
        
    # match vector fields with their string names
    for i, s in enumerate(F_strings):
        LA_original.setdefault(F_strings[i], F[:,i])
    # initialize complete Lie Algebra with original
    LA_complete = copy.copy(LA_original)
    
    # temporary lie algebra, used for linear dependence calculations
    newF_tmp = newF
    
    # initialize    
    have_new_lie_brackets = True
    iterations = 0
    
    # run loop
    while have_new_lie_brackets is True and iterations < n_iterations:
        if show:
            print 'ITERATION: ', iterations, 'out of: ', n_iterations
        
        # re-initialize
        have_new_lie_brackets = False
        iterations += 1
        n = newF.shape[1]
        
        original_vector_field_strings = LA_original.keys()
        lie_algebra_strings = LA_complete.keys()
        
        for i in range(len(original_vector_field_strings)):
            for j in range(len(lie_algebra_strings)):
            
                k1 = original_vector_field_strings[i]
                k2 = lie_algebra_strings[j]
            
                if show:
                    print
                #print i,j
                
                lie_bracket_string_name = '[' + k1 + ', ' + k2 + ']'
                
                if show:
                    print 'calculating: '
                    print lie_bracket_string_name

                if lie_bracket_string_name in tested_lie_brackets:
                    if show:
                        print 'already computed'
                    continue
                else:
                    tested_lie_brackets.append(lie_bracket_string_name)
                
                # hack to make sure f0 term is on the left:
                if k2 == 'f0':
                    if show:
                        print 'skipping - will get covered another time'
                    continue
                
                lb = lie_bracket(LA_original[k1], LA_complete[k2], dx)
                if np.sum(np.abs(lb)) == 0:
                    if show:
                        print 'lie bracket equals 0'
                    LA_zero.setdefault(lie_bracket_string_name, lb)
                    continue
                    
                # temporarily add the new lie bracket to our lie algebra    
                newF_tmp = lb.row_join(newF)
                
                # check to see if there is any linear dependence of the first term on the other terms
                if ignore_linear_independence is False:
                    nullspace = newF_tmp.nullspace()
                    for ii in range(len(nullspace)):
                        if nullspace[ii][0] != 0:
                            nullspace[ii] *= nullspace[ii][0]**-1
                    nullspace_linear_components = [is_vec_linear(nullspace[ii], dx) for ii in range(len(nullspace))]
                    
                    
                    # if the nullspace is linear, then we do not have a linearly independent lie bracket, else we do
                    if np.sum(nullspace_linear_components) == 0:
                        newF = newF_tmp
                        have_new_lie_brackets = True
                        LA_complete.setdefault(lie_bracket_string_name, lb)
                        if show:    
                            print '* Added to Lie Algebra'
                    else:
                        LA_rejected.setdefault(lie_bracket_string_name, lb)
                        if show:
                            print '* Linearly Dependent on other terms'
                        
                else:
                    newF = newF_tmp
                    have_new_lie_brackets = True
                    LA_complete.setdefault(lie_bracket_string_name, lb)
                    
    if show:
        print '*'*25
        print 'showing results for ', n_iterations, ' iterations'
    
        print '*'*25
        print 'original vector fields: '
        print_LA_simple(LA_original)
                    
        print '*'*25
        print 'complete lie algebra: '
        print_LA_simple(LA_complete)
        
        print '*'*25
        print 'linearly dependent lie brackets: '
        print_LA_simple(LA_rejected)
        
        print '*'*25
        print 'lie brackets = zero: '
        print_LA_simple(LA_zero)
                    
    return LA_original, LA_complete, LA_rejected, LA_zero
    
def print_LA_simple(LA):
    # shows the values of the dictionary in a practical way

    keys = LA.keys()
    keys.sort(key = len) # sort shortest to longest
    
    for key in keys:
        print
        print key, ': '
        print LA[key]
        
def get_matrix_from_LA(LA):
    keys = LA.keys()
    keys.sort(key = len) # sort shortest to longest
    
    LA_matrix = None
    for key in keys:
        if LA_matrix is None:
            LA_matrix = copy.copy(LA[key])
        else:
            LA_matrix = LA_matrix.row_join(LA[key])
            
    return LA_matrix
    


    
def is_term_linear(m, var_list):
    # method: tries to typecast the term as a polynomial
    try:
        p = m.as_poly()
        
        if p is None:
            return False
        else:
            deg = p.degree
            if deg <= 1:
                return True
            else:
                return False 
    except (sp.SymbolsError, AttributeError):
        return True
    except:
        print 'help!'
        return False
            
def is_vec_linear(vec, var_list):
    typs = []
    for v in vec:
        typs.append(is_term_linear(v, var_list))
    if np.sum(typs) == len(typs):
        return True
    else:
        return False
    

    
if __name__ == "__main__":
    # example

    w1,w2,w3 = sp.symbols('w1,w2,w3')
    g1, g2, g3 = sp.symbols('g1,g2,g3')

    f0 = sp.Matrix([g1*w2*w3, g2*w1*w3, g3*w1*w2])
    f1 = sp.Matrix([1, 0, 0])
    f2 = sp.Matrix([0, 1, 0])
    dx = [w1, w2, w3]

    F = [f0, f1, f2]

    LA_original, LA_complete, LA_rejected, LA_zero = lbt.lie_algebra(F, dx, 2, show=True)



    
