import numpy as np
import scipy.optimize as opt

from numba import jit, prange, int32, uint32, uint64, int64, float64, complex128;

@jit(float64(int32,int32,float64,float64,float64,float64,float64,float64[:],float64[:],float64[:],float64[:]),
     nopython=True,nogil=True)
def efunc_summ(lmax,shift,g1,g2,g12,Nratio,Mratio,cre_c1,cim_c1,cre_c2,cim_c2):
    # NUMBE USED TO SWEEP OVER SUMMATION INDEXES
    reM1 = 0
    imM1 = 0
    reN1 = 0
    imN1 = 0
    reP1 = 0
    imP1 = 0
    reQ1 = 0
    imQ1 = 0
    reM2 = 0
    imM2 = 0
    reN2 = 0
    imN2 = 0
    reP2 = 0
    imP2 = 0
    reQ2 = 0
    imQ2 = 0
    intra1 = 0
    intra2 = 0
    inter = 0
    kinect1 = 0
    kinect2 = 0
    kin_fac1 = 0.5*Nratio / (Nratio + 1)
    kin_fac2 = 0.5*Mratio / (Nratio + 1)
    # Compute energy of the system
    for m in prange(-lmax,lmax+1,1):
        reM1 = cre_c1[m+lmax]
        imM1 = cim_c1[m+lmax]
        reM2 = cre_c2[m+lmax]
        imM2 = cim_c2[m+lmax]
        kinect1 += kin_fac1*(reM1*reM1+imM1*imM1)*(m+shift)**2
        kinect2 += kin_fac2*(reM2*reM2+imM2*imM2)*(m+shift)**2
        for n in prange(-lmax,lmax+1,1):
            for p in prange(-lmax,lmax+1,1):
                for q in prange(-lmax,lmax+1,1):
                    if (m+n != p+q): continue
                    reN1 = cre_c1[n+lmax]
                    imN1 = cim_c1[n+lmax]
                    reP1 = cre_c1[p+lmax]
                    imP1 = cim_c1[p+lmax]
                    reQ1 = cre_c1[q+lmax]
                    imQ1 = cim_c1[q+lmax]
                    reN2 = cre_c2[n+lmax]
                    imN2 = cim_c2[n+lmax]
                    reP2 = cre_c2[p+lmax]
                    imP2 = cim_c2[p+lmax]
                    reQ2 = cre_c2[q+lmax]
                    imQ2 = cim_c2[q+lmax]
                    intra1 += (reP1*reQ1-imP1*imQ1)*(reM1*reN1-imM1*imN1)+\
                              (reM1*imN1+reN1*imM1)*(reP1*imQ1+reQ1*imP1)
                    intra2 += (reP2*reQ2-imP2*imQ2)*(reM2*reN2-imM2*imN2)+\
                              (reM2*imN2+reN2*imM2)*(reP2*imQ2+reQ2*imP2)
                    inter  += (reP1*reQ2-imP1*imQ2)*(reM1*reN2-imM1*imN2)+\
                              (reM1*imN2+reN2*imM1)*(reP1*imQ2+reQ2*imP1)
    return kinect1 + kinect2 + 0.5*(g1*intra1 + g2*intra2) + g12*inter

@jit(float64(int32,int32,int32,float64,float64,float64,float64[:],float64[:],float64[:],float64[:]),
     nopython=True,nogil=True)
def efunc_grad_summ(k,lmax,shift,kin_fac,g,g12,cre_c1,cim_c1,cre_c2,cim_c2):
    # NUMBE USED TO SWEEP OVER SUMMATION INDEXES
    ind = 0
    reM1 = 0.
    imM1 = 0.
    reN1 = 0.
    imN1 = 0.
    reM2 = 0.
    imM2 = 0.
    reN2 = 0.
    imN2 = 0.
    reIND1 = 0.
    imIND1 = 0.
    reIND2 = 0.
    imIND2 = 0.
    intra = 0.
    inter = 0.
    kinect = kin_fac*cre_c1[k+lmax]*(k+shift)**2
    for m in prange(-lmax,lmax+1):
        reM1 = cre_c1[m+lmax]
        imM1 = cim_c1[m+lmax]
        for n in prange(-lmax,lmax+1):
            reN1 = cre_c1[n+lmax]
            imN1 = cim_c1[n+lmax]
            ind = m+n-k+lmax
            if (ind >= 0 and ind < 2*lmax+1):
                reIND1 = cre_c1[ind]
                imIND1 = cim_c1[ind]
                reIND2 = cre_c2[ind]
                imIND2 = cim_c2[ind]
                intra += reIND1*(reM1*reN1-imM1*imN1)+imIND1*(reM1*imN1+reN1*imM1)
                inter += reIND2*(reM1*reN2-imM1*imN2)+imIND2*(reM1*imN2+reN2*imM1)
    return kinect + 2 * (g * intra + g12 * inter)

def efunc(c,lz,Nratio,Mratio,g1,g2,g12):
    # COMPUTE ENERGY GIVEN THE REAL AND IMAGINARY PART OF THE COEFFICIENTS
    # OF THE CONDENSATE WAVE FUNCTION EXPANSION IN MOMENTUM STATES.  HERE,
    # THE DIMENSIONLESS UNITS IS ACHIEVED USING  THE  RADIUS  AS  UNIT  OF
    # DISTANCE WHEREAS THE SYSTEM HAS 2*PI*R PERIOD IN SPACE.
    # Input : c
    # 'c' has 2*lmax+1 numbers for real part concatenated to 2*lmax+1  for
    # imaginary part of coefficients
    # Input : lz
    # shift the corresponding momentum values for a better basis
    # if ((c.size/4)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/4-1)/2)
    shift = round(lz)
    # shift = 0
    cmat = c.reshape(4,2*lmax+1)
    cre_c1 = cmat[0]
    cim_c1 = cmat[1]
    cre_c2 = cmat[2]
    cim_c2 = cmat[3]
    return efunc_summ(lmax,shift,g1,g2,g12,Nratio,Mratio,cre_c1,cim_c1,cre_c2,cim_c2)

def normCons1(c):
    # CONSTRAINT CONDITION FOR A NORMALIZED STATE. MUST RETURN ZERO
    # if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/4-1)/2)
    norm = 0
    cmat = c.reshape(4,2*lmax+1)
    cre = cmat[0]
    cim = cmat[1]
    for m in range(-lmax,lmax+1):
        re = cre[m+lmax]
        im = cim[m+lmax]
        norm = norm + re*re + im*im
    return norm - 1

def normCons2(c):
    # CONSTRAINT CONDITION FOR A NORMALIZED STATE. MUST RETURN ZERO
    # if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/4-1)/2)
    norm = 0
    cmat = c.reshape(4,2*lmax+1)
    cre = cmat[2]
    cim = cmat[3]
    for m in range(-lmax,lmax+1):
        re = cre[m+lmax]
        im = cim[m+lmax]
        norm = norm + re*re + im*im
    return norm - 1

def momentumCons(c,lz,Nratio):
    # AVERAGE MOMENTUM PER PARTICLE SET TO lz
    # if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/4-1)/2)
    shift = round(lz)
    # shift = 0
    mom1 = 0
    mom2 = 0
    cmat = c.reshape(4,2*lmax+1)
    cre_c1 = cmat[0]
    cim_c1 = cmat[1]
    cre_c2 = cmat[2]
    cim_c2 = cmat[3]
    for m in range(-lmax,lmax+1):
        re = cre_c1[m+lmax]
        im = cim_c1[m+lmax]
        mom1 = mom1 + (re*re+im*im)*(m+shift)
    for m in range(-lmax,lmax+1):
        re = cre_c2[m+lmax]
        im = cim_c2[m+lmax]
        mom2 = mom2 + (re*re+im*im)*(m+shift)
    return mom1*Nratio/(Nratio+1)+mom2/(Nratio+1) - lz

def efunc_grad(c,lz,Nratio,Mratio,g1,g2,g12):
    # MULTI-DIMENSIONAL GRADIENT WITH RESPECT TO INPUT COEFFICIENTS
    # if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/4-1)/2)
    shift = round(lz)
    # shift = 0
    cmat = c.reshape(4,2*lmax+1)
    cre_c1 = cmat[0]
    cim_c1 = cmat[1]
    cre_c2 = cmat[2]
    cim_c2 = cmat[3]
    grad1_re = np.zeros(2*lmax+1)
    grad1_im = np.zeros(2*lmax+1)
    grad2_re = np.zeros(2*lmax+1)
    grad2_im = np.zeros(2*lmax+1)
    kin_fac1 = Nratio/(1+Nratio)
    kin_fac2 = Mratio/(1+Nratio)
    for k in range(-lmax,lmax+1):
        grad1_re[k+lmax] = efunc_grad_summ(k,lmax,shift,kin_fac1,g1,g12,cre_c1,cim_c1,cre_c2,cim_c2)
        grad1_im[k+lmax] = efunc_grad_summ(k,lmax,shift,kin_fac1,g1,g12,cim_c1,cre_c1,cim_c2,cre_c2)
        grad2_re[k+lmax] = efunc_grad_summ(k,lmax,shift,kin_fac2,g2,g12,cre_c2,cim_c2,cre_c1,cim_c1)
        grad2_im[k+lmax] = efunc_grad_summ(k,lmax,shift,kin_fac2,g2,g12,cim_c2,cre_c2,cim_c1,cre_c1)
    return np.concatenate([grad1_re,grad1_im,grad2_re,grad2_im])

def normCons1_grad(c):
    # MULTI-DIMENSIONAL GRADIENT OF NORM CONSTRAINT
    # if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/4-1)/2)
    cmat = c.reshape(4,2*lmax+1)
    cre = cmat[0]
    cim = cmat[1]
    grad_re = 2*cre
    grad_im = 2*cim
    return np.concatenate([grad_re,grad_im,np.zeros(2*lmax+1),np.zeros(2*lmax+1)])

def normCons2_grad(c):
    # MULTI-DIMENSIONAL GRADIENT OF NORM CONSTRAINT
    # if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/4-1)/2)
    cmat = c.reshape(4,2*lmax+1)
    cre = cmat[2]
    cim = cmat[3]
    grad_re = 2*cre
    grad_im = 2*cim
    return np.concatenate([np.zeros(2*lmax+1),np.zeros(2*lmax+1),grad_re,grad_im])

def momentumCons_grad(c,lz,Nratio):
    # MULTI-DIMENSIONAL GRADIENT OF MOMENTUM CONSTRAINT
    # if ((c.size/2)%2 == 0): raise IOError("WRONG NUMBER OF MOMENTUM STATES")
    lmax = int((c.size/4-1)/2)
    shift = round(lz)
    # shift = 0
    cmat = c.reshape(4,2*lmax+1)
    cre_c1 = cmat[0]
    cim_c1 = cmat[1]
    cre_c2 = cmat[2]
    cim_c2 = cmat[3]
    grad1_re = np.zeros(cre_c1.size)
    grad1_im = np.zeros(cim_c1.size)
    grad2_re = np.zeros(cre_c2.size)
    grad2_im = np.zeros(cim_c2.size)
    for k in range(-lmax,lmax+1):
        grad1_re[k+lmax] = 2*cre_c1[k+lmax]*(k+shift)*Nratio/(Nratio+1)
        grad1_im[k+lmax] = 2*cim_c1[k+lmax]*(k+shift)*Nratio/(Nratio+1)
        grad2_re[k+lmax] = 2*cre_c2[k+lmax]*(k+shift)/(Nratio+1)
        grad2_im[k+lmax] = 2*cim_c2[k+lmax]*(k+shift)/(Nratio+1)
    return np.concatenate([grad1_re,grad1_im,grad2_re,grad2_im])

def yrast(Nratio,Mratio,g,lz_init=0,lz_final=3,dl=0.01):
    lmax = 8
    g1 = g[0]
    g2 = g[1]
    g12 = g[2]
    # random initial guess
    c0mat = np.random.random([2,2*(2*lmax+1)]) - 0.5
    # Renormalize each component coefficients
    c0mat[0] = c0mat[0] / np.sqrt(c0mat[0].dot(c0mat[0]))
    c0mat[1] = c0mat[1] / np.sqrt(c0mat[1].dot(c0mat[1]))
    c0 = np.concatenate([c0mat[0],c0mat[1]])
    c0_last = c0
    lz_vals = np.arange(lz_init,lz_final+dl/2,dl)
    E = np.zeros(lz_vals.size)
    extra = {'ftol':1E-6,'disp':False,'maxiter':3500}
    for i in range(lz_vals.size):
        lz = lz_vals[i]
        print("[{:4d}/{}] Working on lz = {:.3f}".format(i+1,lz_vals.size,lz))
        cons = [
        {'type':'eq',
         'fun' : normCons1,
         'jac' : normCons1_grad},
        {'type':'eq',
         'fun' : normCons2,
         'jac' : normCons2_grad},
        {'type':'eq',
         'fun' : momentumCons,
         'jac' : momentumCons_grad,
         'args':(lz,Nratio)}]
        res = opt.minimize(efunc,c0,method='SLSQP',jac=efunc_grad,
              args=(lz,Nratio,Mratio,g1,g2,g12),constraints=cons,options=extra)
        if (not res.success):
            res = opt.minimize(efunc,c0_last,method='SLSQP',jac=efunc_grad,
                  args=(lz,Nratio,Mratio,g1,g2,g12),constraints=cons,options=extra)
        if (not res.success):
            print(res)
            raise RuntimeError("Fail to converge for lz = {:.2f}".format(lz))
        E[i] = res.fun
        c0mat = (res.x*(1 + (np.random.random(c0.size)-0.5)/4)).reshape(2,2*(2*lmax+1))
        c0mat[0] = c0mat[0] / np.sqrt(c0mat[0].dot(c0mat[0]))
        c0mat[1] = c0mat[1] / np.sqrt(c0mat[1].dot(c0mat[1]))
        c0_last = np.concatenate([c0mat[0],c0mat[1]])
    return lz_vals, E
