from __future__ import print_function
import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import astropy.constants as C
import astropy.units as u

# import k_B, Ryd, h, c
figname = sys.argv[0].replace('.py', '.pdf')
D = '../Storey-OIII/oiii/'
kT_ergs = 1e4*u.K*C.k_B.cgs 
Ryd_ergs = C.Ryd.cgs*C.h.cgs*C.c.cgs
ubar = np.sqrt(2*kT_ergs/C.m_e.cgs).value

def f_u(u):
    """Maxwellian distribution per unit velocity u"""
    xi = u / ubar
    return 4*xi**2 * np.exp(-xi**2) / (np.sqrt(np.pi)*ubar)

# ID, term, stat weight for ground level
lower = [
    ['1', '3P0', 1],
    ['2', '3P1', 3],
    ['3', '3P2', 5],
]

# ID, term, plotting scale factor for excited levels
upper = [
    ['4', '$^{1}$D', 5e-16],
    ['5', '$^{1}$S', 1e-17],
]

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False)
for id_upper, term_upper, scale in upper:
    label = '$^{3}$P - ' + term_upper
    E_stack = []
    sigma_stack = []
    E_common = []
    omegas = []
    for id_lower, term_lower, omega in lower:
        s = '{}_{}'.format(id_lower, id_upper)
        E_ryd, Omega = np.loadtxt(D + 'OMEGA_{}_OIII.dat'.format(s), unpack=True)
        E_ergs = E_ryd*Ryd_ergs
        sigma = 2*np.pi*(C.h.cgs/(2*np.pi))**2 * (Omega/omega) / (C.m_e.cgs*E_ergs)
        E_over_kT = E_ergs/kT_ergs
        omegas.append(omega)
        E_stack.append(E_over_kT.value)
        sigma_stack.append(sigma.value)
        E_common.extend(list(E_over_kT.value))
    E_common = np.array(sorted(list(set(E_common))))
    sigma_common = np.zeros_like(E_common)
    omega_sum = 0.0
    for E, sigma, omega in zip(E_stack, sigma_stack, omegas):
        sigma_common += omega*np.interp(E_common, E, sigma)
        omega_sum += omega
    sigma_common /= omega_sum
    ax1.plot(E_common*kT_ergs/Ryd_ergs, sigma_common/1e-15, lw=0.6, label=label)
    chi = E_common[0]
    E_final = E_common - chi
    u_final = np.sqrt(E_final)*ubar
    u_incident = np.sqrt(E_common)*ubar
    sigma_final = np.interp(E_final, E_common, sigma_common, left=0.0)
    usigf_final = u_final*sigma_final*f_u(u_final)
    usigf_incident = u_incident*sigma_common*f_u(u_incident)
    label2 = '({}) / {:.0e} cm$^{{-2}}$'.format(label, scale)
    ax2.plot(np.sqrt(E_final), (usigf_incident - usigf_final)/scale, label=label2)

ax1.set_xlabel('Electron energy, Rydbergs')
ax1.set_ylabel('Cross section, $10^{-15}$ cm$^2$')
ax1.set_xlim(0.0, 1.3)
ax1.set_ylim(0, 1.3)
ax1.set_ylim(3e-3, 3)
ax1.set_yscale('log')
ax1.legend()
xigrid = np.linspace(0.0, 4.0, 200)
fugrid = ubar*f_u(xigrid*ubar)
ax2.fill_between(xigrid, fugrid, -fugrid, color='r', alpha=0.2)
ax2.plot([], [], lw=10, solid_capstyle='butt', color='r', alpha=0.2, label="$\pm f^{*}_{u}$")
ax2.set_xlabel(r'Electron velocity: $\xi = u\, / \langle u \rangle$')
ax2.set_ylabel(r"Net $e^{-}$ gain-loss bracket: $\left[u'\sigma(u') f^{*}_{u'} - u\sigma(u) f^{*}_{u}\right]$")
ax2.set_xlim(0.0, 4.0)
ax2.set_ylim(-1.0, 1.0)
ax2.legend()
fig.set_size_inches(5, 8)
fig.tight_layout()
fig.savefig(figname)
print(figname)
