"""
-------------------------------------------------------------
This is a Python script
-------------------------------------------------------------
Author: David Hasler
Date:   25. 08. 2015
Topic:  Fractal lung - input generation
-------------------------------------------------------------
Description:
This script reads the flow channel of a B-file from a
multi-breath washout (MBW). The signal is low-
pass filtered and cut to a specific length corresponding to
'nbr' breaths. The processed signal, together with a table of
tidal volumes and breath periods serves as input for 'fractal
lung'
-------------------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob

from numpy import sin, pi

from numpy.matlib import repmat
from scipy.signal import cheby2, filtfilt
from scipy.interpolate import CubicSpline, BPoly
from scipy.integrate import trapz
from scipy.optimize import minimize

#--------------------
# PARAMETERS
#--------------------

# number of breaths to be analyzed
nbr = 3

# time-step => inverse of sampling frequency of 'InletFlow'
dt = 0.001

# breath period
tb = 3.2
TB = tb*np.ones(nbr)

# tidal volume
tv = 0.0005
TV = tv*np.ones(nbr)

# asymmetry in flow / pressure profile (inhalation fraction)
asym = 0.4;
zr_cr_slp = (tv*pi/tb)/(asym*tb/4.)

# slope of flow profile at zero-crossing

# pleural pressure parameters
ppl_os   = -1000.
ppl_ampl = 1500. # as absolute measure

# breath length
N = np.round(TB/dt)

# time domain
time = np.arange(0, nbr*tb, dt)
if time.size > int(np.round(nbr*tb/dt)):
    time = time[:int(np.round(nbr*tb/dt))]

#--------------------
# GENERATE FLOW INPUT
#--------------------

# control points inhalation (for one period, normalized)
def ctr_points(alpha, beta):
    t_ctr = np.array([0.00, 0.20*asym, 0.80*asym, 1.00*asym, 0.20*(1.00-asym) + asym, 0.80*(1.00-asym) + asym, 1.00])*tb
    f_ctr = np.array([0.00, 0.50,      0.50,      0.00,     -beta*0.50,              -beta*0.50,               0.00])*alpha
    return t_ctr, f_ctr

# construct cubic spline and evaluate on domain 't'
def curve(x, t):

    alpha, beta = x[0], x[1]

    x_ctr, y_ctr = ctr_points(alpha, beta)

    # B-spline
    cs = CubicSpline(x_ctr, y_ctr, bc_type='periodic')

    return cs(t)

# integral
def cost_fct(x):

    t_inh = np.arange(0., asym*tb+dt, dt)
    c = curve(x, t_inh)
    cost_inh = trapz(c, t_inh) - tv

    t_exh = np.arange(asym*tb, tb+dt, dt)
    c = curve(x, t_exh)
    cost_exh = trapz(c, t_exh) + tv

    return cost_inh**2 + cost_exh**2


# find amplitudes for inhalation / exhalation (one cycle)
x0 = np.array([tv*pi/tb, 0.5])
res = minimize(cost_fct, x0, method='Nelder-Mead')
alpha, beta = res.x

# flow profile
time_one = np.arange(0, tb, dt)
flow_one = curve(res.x, time_one)
flow = repmat(flow_one, 1, nbr)[0,:]

# check
print("total cumulative volume: "+str(np.round(trapz(flow, time),5)))


#--------------------
# GENERATE PLEURAL PRESSURE
#--------------------

# control points (for one period)
t_ctr = np.array([0.00, 0.70*asym, 1.00*asym, 0.30*(1.00-asym) + asym, 1.00])*tb
p_ctr = np.array([0.00, 0.50,      1.00,      0.60,                    0.00])*(-ppl_ampl) + ppl_os

# cubic spline
cs = CubicSpline(t_ctr, p_ctr, bc_type='periodic')
ppl_one = cs(time_one)

# pleural pressure profile
ppl = repmat(ppl_one, 1, nbr)[0,:]


#--------------------
# PLOT
#--------------------

fig = plt.figure()

ax = fig.add_subplot(211)

tf_ctr, f_ctr = ctr_points(alpha, beta)
ax.plot(time, flow, 'k')
ax.plot(tf_ctr, f_ctr, '.m')
ax.set_ylabel('flow rate at inlet $Q_\mathrm{in}(t)$ $[m^3/s]$')
ax.grid('on')

ax = fig.add_subplot(212)

ax.plot(time, ppl, 'k')
ax.plot(t_ctr, p_ctr, '.m')
ax.set_ylabel('pleural pressure $p_\mathrm{pl}(t)$ $[Pa]$')
ax.grid('on')


plt.show()


#--------------------
# SAVE
#--------------------
cpath = os.getcwd()
dpath = cpath.replace('scripts','data')

output = flow;
filename = 'inletFlow'
np.savetxt(dpath+'/'+filename, output, '%.8f')

output = ppl;
filename = 'pleuralPressure'
np.savetxt(dpath+'/'+filename, output, '%.8f')

output = np.array([TB, TV, N]).T;
filename = 'TBTVN'
np.savetxt(dpath+'/'+filename, output, ['%f', '%.8f', '%d'])

output = np.array([[nbr], [1./dt]]).T;
filename = 'nbfs'
np.savetxt(dpath+'/'+filename, output, ['%d', '%f'])
