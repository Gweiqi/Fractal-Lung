"""
-------------------------------------------------------------
This is a Python script
-------------------------------------------------------------
Author: David Hasler
Date:   24. 08. 2015
Topic:  Fractal lung - input generation
-------------------------------------------------------------
Description:
This script reads the flow channel of a A-file from a double
tracer gas single breath washout (DTG-SBW). The signal is low-
pass filtered and cut to a specific length corresponding to 
'nbr' breaths. The processed signal, together with a table of
tidal volumes and breath periods serves as input for 'fractal
lung'

It is assumed, that for single breath washouts (A-files) the 
washout always takes place in the last but two breaths.
-------------------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob

from scipy.signal import cheby2
from scipy.signal import filtfilt

#--------------------
# PARAMETERS
#--------------------

# number of breaths to be analyzed
nbr = 1

# sampling frequency of Afile
fs = 200.

# time-step => inverse of sampling frequency of 'InletFlow'
dt = 0.01

# file number
fnbr = 3;

#--------------------
# LOAD DATA
#--------------------

# directories
cpath = os.getcwd()

# get filelist
flist = glob.glob(cpath.replace('scripts','data/lufu/Afiles/*.txt'))
filename = flist[fnbr-1]

# load data
data = np.genfromtxt(filename, skip_header=1)

# extract
time = data[:,0]*1e-03
flow = data[:,1]*1e-03
MMss = data[:,4]
                     
#--------------------
# SIGNAL PROCESSING
#--------------------

# store raw
time_raw = np.copy(time)
flow_raw = np.copy(flow)

# low-pass filter flow signal
b, a = cheby2(4, 40, 0.05, btype='low')
flow_f = filtfilt(b, a, flow)

# re-sample
time = np.arange(time[0], time[-1], dt)

flow_f = np.interp(time, time_raw, flow_f)

# volume
volume = np.cumsum(flow_f)/fs

# zero crossings of flow signal
zr_cr = flow_f*np.roll(flow_f,-1) < 0
zr_cr_ind = np.arange(time.size)[zr_cr] + 1
zr_cr_ind = zr_cr_ind[:-1] # deleta spurious zero crossing from periodic 'roll'

# single breath washout indices
woa = zr_cr_ind[-5]; wob = zr_cr_ind[-3]; wom = zr_cr_ind[-4]

# cut
time_c = time[woa:wob+1]
time_c = time_c[1:] # start value is not input

flow_c = flow_f[woa:wob+1]
flow_c = flow_c - flow_c[0]
flow_c = flow_c[1:] # start value is not input

# tidal volume
TVisp = volume[wom] - volume[woa]
TVexp = volume[wob] - volume[wom]
TV = 0.5*(abs(TVisp) + abs(TVexp))

# breath periode
TB = time[wob] - time[woa]
TB = round(TB,3)

#--------------------
# PLOT
#--------------------

fig = plt.figure()

ax = fig.add_subplot(211)

ax.plot(time_raw, flow_raw, 'k')
ax.plot(time, flow_f, 'm')
ax.plot(time[zr_cr_ind], flow_f[zr_cr_ind], 'om')
ax.plot(time[woa:wob+1], flow_f[woa:wob+1], 'b')
ax.set_xlabel('time $t$ $[s]$')
ax.set_ylabel('flow rate at inlet $Q_\mathrm{in}(t)$ $[m^3/s]$')
ax.grid('on')

ax = fig.add_subplot(212)

ax.plot(time_c, flow_c, 'k')
ax.set_xlabel('time $t$ $[s]$')
ax.set_ylabel('flow rate at inlet $Q_\mathrm{in}(t)$ $[m^3/s]$')
ax.grid('on')

plt.show()


#--------------------
# SAVE
#--------------------
dpath = cpath.replace('scripts','data')

output = flow_c;
filename = 'InletFlow'
np.savetxt(dpath+'/'+filename, output, '%f')

output = np.array([[TB], [TV], [flow_c.size]]).T;
filename = 'TBTVN'
np.savetxt(dpath+'/'+filename, output, ['%f', '%f', '%d'])

output = np.array([[nbr], [1./dt]]).T;
filename = 'nbfs'
np.savetxt(dpath+'/'+filename, output, ['%d', '%f'])





