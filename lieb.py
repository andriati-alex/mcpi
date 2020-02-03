import sys;
import numpy as np;
import scipy.optimize as opt;
import matplotlib.pyplot as plt;

def f(k,g,N):
    comp = np.zeros(N);
    Ij = 0;
    for j in range(N):
        s = 0;
        Ij = (j + 1 - (N + 1) / 2);
        for l in range(N):
            s = s + np.arctan2(k[j] - k[l],g);
        comp[j] = 2 * np.pi * Ij - 2 * s - k[j];
    return comp;

N = int(5);
gmax = 1000.0;
x0 = 2 * np.random.random(N);
gs = np.arange(0.0,gmax,gmax/1000);
E = np.zeros(gs.size);
n = 0;
for g in gs:
    sol = opt.root(f,x0,args=(g,N));
    k = sol.x;
    # compute energy per particle
    for i in range(k.size): E[n] = E[n] + 0.5 * k[i] * k[i];
    E[n] = E[n]/N;
    n = n + 1;

Efermi = 0.5 * ((2 * np.pi * np.arange(-int(N/2),int((N+1)/2)))**2).sum() / N;
print(Efermi);

plt.plot(gs,E);
plt.show();
