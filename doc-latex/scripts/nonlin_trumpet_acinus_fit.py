import numpy as np
from numpy import exp


import matplotlib.pyplot as plt

from scipy.optimize import fsolve


# parameters
TV = 0.0007
FRC = 0.0035
Ne = 300.

V0 = FRC/Ne

dP = 1500


# non-linear material law
def F(gamma, dP, dV, dP_tilde, dV_tilde):
    return dP_tilde - dP*(exp(gamma*dV_tilde) - 1.)/(exp(gamma*dV) - 1.)


def p_el(V, phi, bending):

    dV = phi*TV/Ne

    # middle fit point
    dP_tilde = 0.5*dP*(1. - bending)
    dV_tilde = 0.5*dV*(1. + bending)

    # find gamma
    gamma = fsolve(F, 1e03, args = (dP, dV, dP_tilde, dV_tilde))

    #print(gamma)

    # non-linear material law
    beta = dP/(exp(gamma*V0)*(exp(gamma*dV) - 1.))

    return beta*(exp(gamma*V) - exp(gamma*V0))

# domains
V = np.linspace(V0-TV/Ne,V0+2.*TV/Ne,1000)

N = 10
PHI = np.linspace(0.5, 1.5, N)
BEND = np.linspace(0.1, 0.9, N)
# calculate gamma


# plots
fig = plt.figure(figsize=(10, 4))

ax = fig.add_subplot(121, ylim=(0, 2000))

ax.axvspan(V0*1e06, (V0+TV/Ne)*1e06, color=[0.7,0.2,0.1], alpha=0.2)
ax.axhspan(0, dP, color=[0.1,0.2,0.8], alpha=0.2)
for phi in PHI:
    ax.plot(V*1e06,p_el(V, phi, 0.5), '-.k', lw=0.7)
ax.plot(V*1e06, p_el(V, 1.0, 0.5), 'k', lw=1.5)
ax.set_xlabel('lobule volume $V_\mathrm{lb}$ $[ml]$')
ax.set_ylabel('elastic pressure $p_\mathrm{el}(V)$ $[Pa]$')

plt.grid('on')



ax = fig.add_subplot(122, ylim=(0, 2000))

ax.axvspan(V0*1e06, (V0+TV/Ne)*1e06, color=[0.7,0.2,0.1], alpha=0.2)
ax.axhspan(0, dP, color=[0.1,0.2,0.8], alpha=0.2)
for bend in BEND:
    ax.plot(V*1e06,p_el(V, 1.0, bend), '-.k', lw=0.7)
ax.plot(V*1e06, p_el(V, 1.0, 0.5), 'k', lw=1.5)
ax.set_xlabel('lobule volume $V_\mathrm{lb}$ $[ml]$')

plt.grid('on')



plt.show()
plt.savefig('nonlinear_lobule.pdf')
