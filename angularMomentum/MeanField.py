import numpy as np
import scipy.optimize as opt

from numba import jit, prange, int32, uint32, uint64, int64, float64, complex128;

@jit(float64(int32,int32,float64,float64[:],float64[:]),
     nopython=True,nogil=True)
def efunc_summ(lmax,shift,g,cre,cim):
    # NUMBE USED TO SWEEP OVER SUMMATION INDEXES
    reM = 0
    imM = 0
    reN = 0
    imN = 0
    reP = 0
    imP = 0
    reQ = 0
    imQ = 0
    inter = 0
    kinect = 0
    for m in prange(-lmax,lmax+1,1):
        reM = cre[m+lmax]
        imM = cim[m+lmax]
        kinect += (reM*reM+imM*imM)*(m+shift)**2
        for n in prange(-lmax,lmax+1,1):
            for p in prange(-lmax,lmax+1,1):
                for q in prange(-lmax,lmax+1,1):
                    if (m+n != p+q): continue
                    reN = cre[n+lmax]
                    imN = cim[n+lmax]
                    reP = cre[p+lmax]
                    imP = cim[p+lmax]
                    reQ = cre[q+lmax]
                    imQ = cim[q+lmax]
                    inter += (reP*reQ-imP*imQ)*(reM*reN-imM*imN)+(reM*imN+reN*imM)*(reP*imQ+reQ*imP)
    return 0.5 * (kinect + g * inter)

@jit(float64(int32,int32,int32,float64,float64[:],float64[:]),
     nopython=True,nogil=True)
def efunc_grad_summ(k,lmax,shift,g,cre,cim):
    # NUMBE USED TO SWEEP OVER SUMMATION INDEXES
    ind = 0
    reM = 0.
    imM = 0.
    reN = 0.
    imN = 0.
    reIND = 0.
    imIND = 0.
    inter = 0.
    kinect = cre[k+lmax]*(k+shift)**2
    for m in prange(-lmax,lmax+1):
        reM = cre[m+lmax]
        imM = cim[m+lmax]
        for n in prange(-lmax,lmax+1):
            reN = cre[n+lmax]
            imN = cim[n+lmax]
            ind = m+n-k+lmax
            if (ind >= 0 and ind < 2*lmax+1):
                reIND = cre[ind]
                imIND = cim[ind]
                inter += reIND*(reM*reN-imM*imN)+imIND*(reM*imN+reN*imM)
    return kinect + 2 * g * inter

def efunc(c,lz,g):
    # COMPUTE ENERGY GIVEN THE REAL AND IMAGINARY PART OF THE COEFFICIENTS
    # OF THE CONDENSATE WAVE FUNCTION EXPANSION IN MOMENTUM STATES.  HERE,
    # THE DIMENSIONLESS UNITS IS ACHIEVED USING  THE  RADIUS  AS  UNIT  OF
    # DISTANCE WHEREAS THE SYSTEM HAS 2*PI*R PERIOD IN SPACE.
    # Input : c
    # 'c' has 2*lmax+1 numbers for real part concatenated to 2*lmax+1  for
    # imaginary part of coefficients
    # Input : lz
    # shift the corresponding momentum values for a better basis
    if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/2-1)/2)
    shift = round(lz)
    kinect = 0
    inter = 0
    cmat = c.reshape(2,2*lmax+1)
    cre = cmat[0]
    cim = cmat[1]
    return efunc_summ(lmax,shift,g,cre,cim)

def normCons(c):
    # CONSTRAINT CONDITION FOR A NORMALIZED STATE. MUST RETURN ZERO
    if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/2-1)/2)
    norm = 0
    cmat = c.reshape(2,2*lmax+1)
    cre = cmat[0]
    cim = cmat[1]
    for m in range(-lmax,lmax+1):
        re = cre[m+lmax]
        im = cim[m+lmax]
        norm = norm + re*re + im*im
    return norm - 1

def momentumCons(c,lz):
    # AVERAGE MOMENTUM PER PARTICLE SET TO lz
    if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/2-1)/2)
    shift = round(lz)
    mom = 0
    cmat = c.reshape(2,2*lmax+1)
    cre = cmat[0]
    cim = cmat[1]
    for m in range(-lmax,lmax+1):
        re = cre[m+lmax]
        im = cim[m+lmax]
        mom = mom + (re*re+im*im)*(m+shift)
    return mom - lz

def efunc_grad(c,lz,g):
    # MULTI-DIMENSIONAL GRADIENT WITH RESPECT TO INPUT COEFFICIENTS
    if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/2-1)/2)
    shift = round(lz)
    cmat = c.reshape(2,2*lmax+1)
    cre = cmat[0]
    cim = cmat[1]
    grad_re = np.zeros(cre.size)
    grad_im = np.zeros(cim.size)
    for k in range(-lmax,lmax+1):
        grad_re[k+lmax] = efunc_grad_summ(k,lmax,shift,g,cre,cim)
        grad_im[k+lmax] = efunc_grad_summ(k,lmax,shift,g,cim,cre)
    return np.concatenate([grad_re,grad_im])

def normCons_grad(c):
    # MULTI-DIMENSIONAL GRADIENT OF NORM CONSTRAINT
    if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/2-1)/2)
    cmat = c.reshape(2,2*lmax+1)
    cre = cmat[0]
    cim = cmat[1]
    grad_re = np.zeros(cre.size)
    grad_im = np.zeros(cim.size)
    for k in range(-lmax,lmax+1):
        grad_re[k+lmax] = 2*cre[k+lmax]
        grad_im[k+lmax] = 2*cim[k+lmax]
    return np.concatenate([grad_re,grad_im])

def momentumCons_grad(c,lz):
    # MULTI-DIMENSIONAL GRADIENT OF MOMENTUM CONSTRAINT
    if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/2-1)/2)
    shift = round(lz)
    cmat = c.reshape(2,2*lmax+1)
    cre = cmat[0]
    cim = cmat[1]
    grad_re = np.zeros(cre.size)
    grad_im = np.zeros(cim.size)
    for k in range(-lmax,lmax+1):
        grad_re[k+lmax] = 2*cre[k+lmax]*(k+shift)
        grad_im[k+lmax] = 2*cim[k+lmax]*(k+shift)
    return np.concatenate([grad_re,grad_im])

def yrast(g,lz_init=0,lz_final=2,dl=0.01):
    lmax = 7
    if (g > 1) : lmax = 10
    if (g > 5) : lmax = 15
    if (g > 10): lmax = 20
    if (g > 15): lmax = 25
    if (g > 20): lmax = 30
    c0 = np.random.random(2*(2*lmax+1)) - 0.5
    c0[int(c0.size/4)] = 2.
    c0 = c0 / np.sqrt(c0.dot(c0))
    lz_vals = np.arange(lz_init,lz_final+dl/2,dl)
    E = np.zeros(lz_vals.size)
    for i in range(lz_vals.size):
        lz = lz_vals[i]
        print("[{:4d}/{}] Working on lz = {:.3f}".format(i+1,lz_vals.size,lz))
        cons = [{'type':'eq',
         'fun' : normCons,
         'jac' : normCons_grad},
        {'type':'eq',
         'fun' : momentumCons,
         'jac' : momentumCons_grad,
         'args':(lz,)}]
        extra = {'ftol':1E-6,
        'disp':False,
        'maxiter':1000}
        res = opt.minimize(efunc,c0,method='SLSQP',jac=efunc_grad,
                           args=(lz,g),constraints=cons,options=extra)
        if (not res.success):
            raise RuntimeError("Fail to converge for lz = {.2f}".format(lz))
        E[i] = res.fun
    return lz_vals, E
