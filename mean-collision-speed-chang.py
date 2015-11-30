from __future__ import print_function
import sys
from scipy.special import erf
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

def g(u):
    return np.exp(-u**2)/np.sqrt(np.pi) + (u + 1./(2*u))*erf(u)


figname = sys.argv[0].replace('.py', '.pdf')

fig, ax = plt.subplots(1, 1)
x = np.logspace(-1.0, 1.0, 200)
ax.plot(x, g(x))
ax.plot(x, x, lw=0.2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(0.1, 10)
ax.set_xlabel(r'$\xi = u / \hat{u}$')
ax.set_ylabel(r'Mean relative collision speed: $\langle g \rangle$')
fig.set_size_inches(4, 4)
fig.tight_layout()
fig.savefig(figname)
print(figname)
