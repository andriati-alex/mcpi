"""

    SCRIPT TO PLOT THE LIEB GAS GROUND STATE ENERGY
    ***********************************************

    The Lieb gas assume an anlytic formula for the energy spectrum through
    the solution of a nonlinear system of equations.  Here  this system of
    equation is solved to plot the ground state energy as function  of the
    interaction parameter strength.  The convention adopted here  for  the
    dimensionless Schrodinger equation describing the Lieb Gas is

    [ -(1/2) * (d/dx_j^2) + g * delta(x_i - x_j) ] \Psi(x_1,...,x_N)
                                   = E \Psi(x_1,...,x_N)

    where it is implied the sum over j, and for i > j.  This convention is
    the same convention adopted in equation (1) of ref [1] with  \hbar and
    mass being 1. It is different from  equation (3) of  ref [2] where the
    equivalence is achieved by multiplying this equation by 2 and defining
    2 * E = E_L(c) from their paper.

    The systems of nonlinear equations to be solved are given in  equation
    (39) of ref [1] and is implemented here.

    [1] 'Ground-state  energy  and excitation spectrum of the Lieb-Lininger
         model: accurate analytical results and conjectures about the exact
         solution', G. Lang, F. Hekking, Anna Minguzzi,
         doi : 10.21468/SciPostPhys.3.1.003
    [2] 'Exact ground state of finite Bose-Einstein condensates on a ring',
         K. Sakmann, A. Streltsov, O. Alon, L. Cederbaum,
         doi : 10.1103/PhysRevA.72.033613

"""

import sys;
import numpy as np;
import matplotlib as mpl;
import scipy.optimize as opt;
import matplotlib.pyplot as plt;
from mpl_toolkits.axes_grid1.inset_locator import inset_axes;

# figure parameters
mpl.rcParams['text.usetex'] = True;
mpl.rcParams['font.family'] = 'serif';
mpl.rcParams['font.serif'] = 'DejaVu Sans';
mpl.rcParams['font.size'] = 10;



def f(k,g,N):
    # Check ref [1] for more details
    comp = np.zeros(N);
    Ij = 0;
    for j in range(N):
        s = 0;
        Ij = (j + 1 - (N + 1) / 2);
        for l in range(N):
            s = s + np.arctan2(k[j] - k[l],g);
        comp[j] = 2 * np.pi * Ij - 2 * s - k[j];
    return comp;



N = int(5); # Number of particles
gmax = 1100.0; # Maximum interaction strength to sweep
gs = np.arange(0.0,gmax,gmax/1000); # values of interaction sweeped
x0 = 2 * np.random.random(N); # initial guess for root function
E = np.zeros(gs.size); # energy for each value of interaction
n = 0;
for g in gs:
    sol = opt.root(f,x0,args=(g,N));
    k = sol.x;
    # compute energy per particle
    E[n] = 0.5 * (k**2).sum() / N;
    n = n + 1;

# switch case. According to Tonks-Girardeau formula the fermi wave  function
# must be periodic/anti-periodic if the number of particles is odd/even. See
# 'Fermi-Bose mapping for one-dimensional Bose gases',  Yukalov & Girardeau,
# Laser Physics Letters, doi : 10.1002/lapl.200510011
if (N % 2 == 0) :
    k = np.pi * np.arange(-N+1,N-1+0.1,2);
    Efermi = 0.5 * (k**2).sum() / N;
else :
    k = 2 * np.pi * np.arange(-int(N/2),int(N/2)+0.1);
    Efermi = 0.5 * (k**2).sum() / N;

# Load numerical simulation data
Enum = np.loadtxt('liebGS.dat');
gnum = Enum[:,0];

fig = plt.figure(figsize=(7,6.2));
ax = plt.gca();

ax.plot(gs,E,lw=1,color='black',label='LL gas');
ax.plot(gnum,Enum[:,1],'bo',label='21 IPS',markerfacecolor='none');
ax.plot(gnum,Enum[:,2],'g+',label='11 IPS',);
ax.axhline(Efermi,0.7,1,lw=3,ls='--',color='red',label='TG gas');
ax = plt.gca();
ax.set_xlim(gs[0],gs[-1]);
ax.set_ylim(E[0],Enum[-1,2]*1.05);
ax.set_ylabel("$E_0 (m L^2 / \hbar^2)$",fontsize=12);
ax.set_xlabel("$g (m L^2 / \hbar^2)$",fontsize=12);
ax.legend(loc='lower right');
ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4));
ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4));

# plot zoom
axin = inset_axes(ax,width="50%",height="45%",loc=8,borderpad=2);
axin.plot(gs[:60],E[:60],'k-');
axin.plot(gnum,Enum[:,1],'bo',ms=6,markerfacecolor='none');
axin.plot(gnum,Enum[:,2],'g+',ms=6);
axin.tick_params(axis='both',which='major',labelsize=8);
axin.set_xlim(0,60);
axin.set_ylim(0,35);
axin.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2));
axin.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(4));

plt.savefig('LiebComparison.pdf',dpi=256,bbox_inches='tight');
