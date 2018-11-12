from __future__ import print_function
import sys
from scipy.special import gamma, gammaln
import numpy as np
from numpy import exp, sqrt, pi
from matplotlib import pyplot as plt
import seaborn as sns

def A_kappa(kappa):
    if kappa < 100.0:
        gfactor = gamma(kappa+1)/gamma(kappa-0.5)
    else:
        # Avoid gamma function blowing up for large kappa
        gfactor = exp(gammaln(kappa+1) - gammaln(kappa-0.5))
    return gfactor/(kappa-1.5)**1.5


def f_M(x):
    '''Maxwellian phase-space density'''
    return exp(-x**2)/sqrt(pi**3)


def f_K(x, kappa):
    '''Kappa phase-space density'''
    return A_kappa(kappa) / sqrt(pi**3) / (1 + x**2/(kappa - 1.5))**(kappa + 1)


def x_K(x, kappa):
    return x/sqrt((1.0 + kappa)/(kappa - 1.5))


def f_Kcore(x, kappa):
    '''Kappa phase-space density, but with core fitted to Maxwellian'''
    return 1.0 / sqrt(pi**3) / (1 + x_K(x, kappa)**2/(kappa - 1.5))**(kappa + 1)


def f_MM(x, kappa):
    '''Maxwellian phase space density, but with same energy as f_Kcore'''
    return f_M(x_K(x, kappa))#*A_kappa(kappa)


def relative_difference(x, kappa):
    '''Relative to kappa difference between kappa and Maxwellian'''
    return (f_K(x, kappa) - f_M(x))/f_K(x, kappa)


def relative_difference_M(x, kappa):
    '''Relative to Maxwellian difference between kappa and Maxwellian'''
    return (f_K(x, kappa) - f_M(x))/f_M(x)


def fokker_planck_kappa(x, kappa):
    '''Multiply by f_K to get df/dt for e-e collisions'''
    top = 4*(kappa + 1)*(3.5 - x**2)
    bottom = x*(x**2 + kappa - 1.5)
    return top/bottom


def timescale_FP(x, kappa):
    '''Time for relaxation from kappa -> maxwell'''
    return -relative_difference(x, kappa)/fokker_planck_kappa(x, kappa)


fig, ax = plt.subplots(1, 1)
x = np.linspace(0.0, 100.0, 2000)
kappas = 2.0, 20.0, 200.0, 2000.0, 20000.0
for kappa in kappas:
    ax.plot(x, 4*pi*x**2*(f_Kcore(x, kappa) - f_M(x)),
            label=r'$\kappa = {:.0f}$'.format(kappa))
ax.set_color_cycle(None)
for kappa in kappas:
    ax.plot(x, 4*pi*x**2*(f_MM(x, kappa) - f_M(x)), '--', lw=0.5)
ax.plot(x, 4*pi*x**2*f_M(x), lw=4, color='k', alpha=0.3)
# ax.set_yscale('symlog', linthreshy=0.001, linscaley=0.1)
ax.set_xlim(0.1, 100.0)
ax.set_ylim(1e-12, 10.0)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('$u / \hat{u}$')
ax.set_ylabel('$4\pi (u / \hat{u})^{2} f$')
ax.legend()
figfile = sys.argv[0].replace('.py', '.pdf')
fig.savefig(figfile)
print(figfile)
