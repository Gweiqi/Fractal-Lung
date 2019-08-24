"""
    -------------------------------------------------------------------------
    This is a Python script
    -------------------------------------------------------------------------
    File      : pulsatile_transmissibility.py
    Institute : ARTORG
    Topic     : Airway modeling
    -------------------------------------------------------------------------
    Description:
    Comparison of the flow resistance (transmissibility) in a rigid tube for
    oscillatory and steady flow.

    References:
    Lecture notes biofluidmechanics
    [Kaczka DW 2011] Computational Anayses of Airway flwo and lung tissue
    dynamcis.

    -------------------------------------------------------------------------
    -------------------------------------------------------------------------
    """
# ----------
# imports
# ----------
import os

import numpy as np
from numpy import pi, sqrt

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import scipy.special

# ----------
# Paramter
# ----------

# fluid properties
mu  = 18e-06
rho = 1.22

# domain range
# !!! These values must correspond to the domain range variables hard-coded in the lung.cpp file of the 'FractalLung' source code !!!
N = 100

TBmax = 10.
TBmin = 0.1

dmax = 0.04
dmin = 0.0005




# ----------
# Law
# ----------

# transmissibility for oscillatory flow
def T(TB,a):
    omega = 2*pi/TB
    d = 2*a
    l = 10*d
    alpha = a*sqrt(omega*rho/mu)
    arg = (1j)**(3./2.)*alpha
    J0 = scipy.special.jv(0,arg)
    J1 = scipy.special.jv(1,arg)
    return np.real(-1j*(a**2*pi)/(omega*rho*l)*(1.-(2*J1)/((1j)**(3./2.)*alpha*J0)))

# transmissibility for steasy (Poiseuille) flow
def Tp(a):
    d = 2*a
    l = 10*d
    return (a**4*pi)/(8*mu*l)


# ----------
# ----------
N = 100

# frequency / period
TB = np.linspace(TBmin,TBmax,N)
TBm, none = np.meshgrid(TB,TB)

omegamin = 2*pi/TBmax
omega = np.linspace(omegamin,10*omegamin,N)

# radius
a = np.linspace(dmin/2.,dmax/2.,N)

womersley = a*sqrt(omega*rho/mu)

Omega, A = np.meshgrid(omega,a)
Diam = 2.0*A
Womersley = A*sqrt(Omega*rho/mu);

fig = plt.figure()

ax = fig.add_subplot(211, projection='3d', zlim=(0,1.1))

ax.plot_wireframe(TBm,Diam,T(TB,A)/Tp(A), rstride=5, cstride=5)
ax.set_xlabel('breath period $T_\mathrm{breath}$')
ax.set_ylabel('airway diameter $d$')
ax.set_zlabel('normalized transmissibility $\\frac{T}{T_p}$')

ax.view_init(40,-10)

ax = fig.add_subplot(212, projection='3d')

ax.plot_wireframe(Omega,Diam,Womersley, rstride=5, cstride=5)
ax.set_xlabel('frequency $\omega$')
ax.set_ylabel('airway diameter $d$')
ax.set_zlabel('Womersley number $\\alpha$')

ax.view_init(45,-10)

plt.show()

# ----------
# Save
# ----------
path = os.getcwd().replace('/scripts','/constant')
Output = T(TBm,A)/Tp(A)
np.savetxt(path+'/transFact',Output,'%f')


#plt.savefig('PulsatileTransmissibility.pdf', bbox_inches='tight', pad_inches=0)
