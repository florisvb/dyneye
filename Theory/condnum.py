import fly_plot_lib.plot as fpl
import matplotlib.pyplot as plt
import sympy as sp
import numpy as np

r,t,d0,v0,w = sp.symbols('r,t,d0,v0,w')
tfs = np.linspace(0,10,100)
#rval = -.2
ds = np.linspace(10, 0.01, 100)
#v0val = -.1
#wval = 1


# const vel
def const_control():
    expA = sp.Matrix([[1,t],[0,1]])
    expAT = sp.Matrix([[1,0],[t,1]])
    C = sp.Matrix([[-v0/(v0*t+d0)**2, 1/(v0*t+d0)]])
    dP = expAT*C.T*C*expA
    P = sp.integrate(dP, t)

    vels = np.linspace(-.01,-100000,10)
    cn = np.zeros([len(ds), len(vels)])

    for j, vel in enumerate(vels):
        tfs = (ds-d0val)/vel
        for i, tf in enumerate(tfs):        
            Peval = (P.subs({d0:d0val,t:tf,v0:vel}) - P.subs({d0:d0val,t:0,v0:vel})).evalf()
            svals = Peval.singular_values()
            try:
                cond_number = sp.re(np.max(svals)) / sp.re(np.min(svals))
                cn[i,j] = cond_number.evalf()
            except:
                print vel, tf, sp.re(np.max(svals)), sp.re(np.min(svals))
                
    def replace_infs(cn):
        indices = np.where(cn==np.inf)
        for i in range(len(indices[0])):
            cn[indices[0][i], indices[1][i]] = 1e16
        return cn
        
    cn = replace_infs(cn)
    cn += 1e-16

    fig = plt.figure()
    ax = fig.add_subplot(111)
    norm = plt.Normalize(0, 10, clip=True)
    ax.imshow(np.log(cn), norm=norm, origin='lower', extent=[vels[0], vels[-1], ds[0], ds[-1]], interpolation='nearest')
    
    ax.set_aspect('auto')
    ax.set_xlabel('velocity setpoint, 1/s')
    ax.set_ylabel('distance to target, m')
    
    return cn

# control
def exp_control():


    expA = sp.Matrix([[1,t],[0,1]])
    expAT = sp.Matrix([[1,0],[t,1]])
    C = sp.Matrix([[-r*d0*sp.exp(r*t)/(d0**2*sp.exp(2*r*t)), 1/(d0*sp.exp(r*t))]])
    dP = expAT*C.T*C*expA
    P = sp.integrate(dP, t)

    rvals = np.linspace(-.01, -1, 20)
    cn_control = np.zeros([len(ds), len(rvals)])

    for j, rval in enumerate(rvals):
        tfs = np.log(ds/d0val)/rval
        #d = (d0val*sp.exp(rval*tfs))
        for i, tf in enumerate(tfs):        
            Peval = (P.subs({r:rval,d0:d0val,t:tf}) - P.subs({r:rval,d0:d0val,t:0})).evalf()
            svals = Peval.singular_values()
            try:
                cond_number = np.max(svals) / np.min(svals)
                cn_control[i,j] = cond_number.evalf()
            except:
                print rval, tf, svals
            
    fig = plt.figure()
    ax = fig.add_subplot(111)
    norm = plt.Normalize(0, 10, clip=True)
    ax.imshow(np.log(cn_control), norm=norm, origin='lower', extent=[rvals[0], rvals[-1], ds[0], ds[-1]])
    
    ax.set_aspect('auto')
    ax.set_xlabel('optic flow setpoint, 1/s')
    ax.set_ylabel('distance to target, m')

# sinusoidal
def sin_control():

    expA = sp.Matrix([[1,t],[0,1]])
    expAT = sp.Matrix([[1,0],[t,1]])
    
    def v():
        return 2*sp.pi*w*sp.cos(2*sp.pi*w*t)
        
    def d():
        return d0 + sp.sin(2*sp.pi*w*t)
        
        
    C = sp.Matrix([[-v()/d()**2, 1/d()]])
    dP = expAT*C.T*C*expA
    P = sp.integrate(dP, t)

    cn_sin = []
    d_sin = []

    for tf in tfs:

        Peval = (P.subs({r:rval,d0:d0val,t:tf,w:wval}) - P.subs({r:rval,d0:d0val,t:0,w:wval})).evalf()

        svals = Peval.singular_values()

        cond_number = np.max(svals) / np.min(svals)
        
        cn_sin.append(cond_number.evalf())
        d_sin.append( (d().subs({r:rval,d0:d0val,t:tf,w:wval})) )
        

    
    
