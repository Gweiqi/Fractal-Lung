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
from scipy.interpolate import BPoly
from scipy.integrate import trapz
from scipy.optimize import minimize

#--------------------
# PARAMETERS
#--------------------

# number of breaths to be analyzed
nbr = 1

# time-step => inverse of sampling frequency of 'InletFlow'
dt = 0.001

# breath period
tb = 3.2
TB = tb*np.ones(nbr)

# tidal volume
tv = 0.0005
TV = tv*np.ones(nbr)

# slope of flow profile at zero-crossing

# breath length
N = np.round(TB/dt)

# time domain
time = np.arange(0, nbr*tb, dt)
if time.size > int(np.round(nbr*tb/dt)):
    time = time[:int(np.round(nbr*tb/dt))]

#--------------------
# GENERATE FLOW INPUT
#--------------------

# flow profile: simple sine function
flow = tv*pi/tb*sin(2.*pi/tb*time)


#--------------------
# PLOT
#--------------------

fig = plt.figure()

ax = fig.add_subplot(111)

ax.plot(time, flow, 'k')
ax.set_ylabel('flow rate at inlet $Q_\mathrm{in}(t)$ $[m^3/s]$')
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

output = np.array([TB, TV, N]).T;
filename = 'TBTVN'
np.savetxt(dpath+'/'+filename, output, ['%f', '%.8f', '%d'])

output = np.array([[nbr], [1./dt]]).T;
filename = 'nbfs'
np.savetxt(dpath+'/'+filename, output, ['%d', '%f'])
