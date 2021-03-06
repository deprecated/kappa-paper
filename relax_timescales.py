from __future__ import print_function
import sys
from scipy.special import erf
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

logLambda = 20.0

def erfd(x):
    '''First derivative of error function'''
    return 2*np.exp(-x**2)/np.sqrt(np.pi)

def G(xi):
    '''Chandrasekhar big-G function'''
    return 0.5 * (erf(xi) - xi*erfd(xi)) / xi**2


def g(xi):
    '''Chandrasekhar little-g function for mu = 1 (equal masses)'''
    return 0.125 - 2.515*np.exp(-xi) + 4.369*np.exp(-2*xi)

def tau_F(xi):
    '''Velocity-space friction timescale'''
    return 0.5 / G(xi)


def tau_s(xi):
    '''Slowing down timescale'''
    return xi*tau_F(xi)


def tau_d(xi):
    '''Deflection timescale'''
    return 2*xi**5 / ((2*xi**2 - 1)*erf(xi) + xi*erfd(xi))


def tau_E(xi):
    '''Energy loss timescale'''
    return 0.25*xi**3 / (G(xi) + g(xi)/logLambda)


if __name__ == '__main__':
    fig, ax = plt.subplots(1, 1)
    xi = np.logspace(-1.5, 2.5, 500)
    styles = {'lw': 3, 'alpha': 0.8}
    ax.plot(xi, tau_F(xi), label=r'$\tau_F$', **styles)
    ax.plot(xi, tau_s(xi), label=r'$\tau_s$', **styles)
    ax.plot(xi, tau_d(xi), label=r'$\tau_d$', **styles)
    ax.plot(xi, tau_E(xi), label=r'$\tau_E$', **styles)
    ax.set_xlabel(r'$\xi = u\, /  \langle u \rangle$')
    ax.set_ylabel(r'$\tau \, / \, \tau_0 $')
    ax.set_xlim(xi[0], xi[-1])
    ax.set_ylim(1e-3, 1e8)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper left', fontsize='large')
    plotfile = sys.argv[0].replace('.py', '.pdf')
    fig.set_size_inches(3.5, 3.5)
    fig.tight_layout()
    fig.savefig(plotfile)
    print(plotfile)
