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
tidal volumes and breath periods serves as input for the 'flPROG'
application.
-------------------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob

from scipy.signal import cheby2
from scipy.signal import filtfilt
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d, UnivariateSpline

#--------------------
# PARAMETERS
#--------------------

# number of breaths to be analyzed
nbr = 27

# sampling frequency of Bfile
fs = 200

# time-step => inverse of sampling frequency of 'InletFlow'
dt = 0.001

# file number
fnbr = 1;

#--------------------
# LOAD DATA
#--------------------

# directories
cpath = os.getcwd()

# get filelist
flist = np.sort(glob.glob(cpath.replace('scripts','data/lufu/B-Files (for input)/*.txt')))
filename = flist[fnbr-1]

# load data
data = np.genfromtxt(filename, skip_header=1)

# extract
time = data[:,0]*1e-03
flow = data[:,1]*1e-03
N2 = data[:,4]

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

interpolant_flow = interp1d(time_raw, flow_f, kind='quadratic')
flow_f = interpolant_flow(time)
#flow_f = np.interp(time, time_raw, flow_f)

interpolant_N2 = interp1d(time_raw, N2, kind='quadratic')
N2 = interpolant_N2(time)
#N2 = np.interp(time, time_raw, N2)

# find start of washout
st_estim = min(time[N2 < 0.5*max(N2)]) - 0.5
st_estim_ind = np.argwhere(time==st_estim)[0,0]

# volume
volume = cumtrapz(flow_f, time)

# zero crossings of flow signal
zr_cr = flow_f*np.roll(flow_f,-1) < 0
zr_cr_ind = np.arange(time.size)[zr_cr] + 1
zr_cr_ind = zr_cr_ind[:-1] # delete spurious zero crossing from periodic 'roll'
zr_cr_ind = zr_cr_ind[zr_cr_ind > st_estim_ind]

# multi-breath breath washout sequence
woa = zr_cr_ind[0]; wob = zr_cr_ind[nbr*2];

# cut
time_c = time[woa:wob+1]
time_c = time_c[1:] # start value is not input

flow_c = flow_f[woa:wob+1]
flow_c = flow_c - flow_c[0]
flow_c = flow_c[1:] # start value is not input

N2_c = N2[woa:wob+1]
N2_c = N2_c[1:] # start value is not input

volume_c = volume[woa:wob+1]
volume_c = volume_c - volume_c[0]
volume_c = volume_c[1:] # start value is not input

# For signals longer than four breaths, correct flow for zero cumulative volume
spl = UnivariateSpline(time_c, volume_c)
if nbr > 4:
    volume_c_corr = volume_c - (spl(time_c) - spl(time_c[0]))

    N = np.size(time_c)
    flow_c_corr = np.zeros_like(time_c)
    for k in range(N):

        if k == 0:
            flow_c_corr[k] = (volume_c_corr[k+1] - volume_c_corr[k])/(dt)

        elif k == N-1:
            flow_c_corr[k] = (volume_c_corr[k] - volume_c_corr[k-1])/(dt)

        else:
            flow_c_corr[k] = (volume_c_corr[k+1] - volume_c_corr[k-1])/(2*dt)

    volume_c_corr2 = cumtrapz(flow_c_corr, time_c, initial=0)
else:
    volume_c_corr2 = volume_c_corr = volume_c
    flow_c_corr = flow_c


# tidal volume
TV = np.zeros(nbr)
for i in range(nbr):
    TV[i] = abs(volume_c_corr[zr_cr_ind[2*i + 1]-woa-1] - volume_c_corr[zr_cr_ind[2*i]-woa-1]) + abs(volume_c_corr[zr_cr_ind[2*i + 2]-woa-1] - volume_c_corr[zr_cr_ind[2*i + 1]-woa-1])
TV = 0.5*TV

# breath periode
TB = np.zeros(nbr)
Ts = np.zeros(nbr)
Te = np.zeros(nbr)
for i in range(nbr):
    Te[i] = time[zr_cr_ind[2*i + 2]]
    Ts[i] = time[zr_cr_ind[2*i]]
    TB[i] = Te[i] - Ts[i]

#TB[0] += dt
TB = np.round(TB,3)

# breath length
N = np.round(TB/dt)

#--------------------
# PLOT
#--------------------

fig = plt.figure()

ax = fig.add_subplot(411)

ax.plot(time_raw, flow_raw, 'k')
ax.plot(time, flow_f, 'm')
ax.plot(time[zr_cr_ind], flow_f[zr_cr_ind], 'om')
ax.plot(time[zr_cr_ind-1], flow_f[zr_cr_ind-1], 'or')
ax.plot(time[woa:wob+1], flow_f[woa:wob+1], 'b')
ax.set_ylabel('flow rate at inlet $Q_\mathrm{in}(t)$ $[m^3/s]$')
ax.grid('on')

ax = fig.add_subplot(412)

ax.plot(time, N2, 'b')
ax.set_xlabel('time $t$ $[s]$')
ax.set_ylabel('nitrogen $N_2$')
ax.grid('on')

ax = fig.add_subplot(413)

ax.plot(time_c, flow_c, 'k')
ax.plot(time_c, flow_c_corr, 'g')
ax.set_xlabel('time $t$ $[s]$')
ax.set_ylabel('flow rate at inlet $Q_\mathrm{in}(t)$ $[m^3/s]$')
ax.grid('on')

ax = fig.add_subplot(414)

ax.plot(time_c, volume_c, 'k')
ax.plot(time_c, spl(time_c), '--k')
ax.plot(time_c, volume_c_corr, 'g')
ax.plot(time_c, volume_c_corr2, '--g')
ax.set_xlabel('time $t$ $[s]$')
ax.set_ylabel('cum. vol. $V_\mathrm{c}(t)$ $[m^3]$')
ax.grid('on')


plt.show()

#--------------------
# SAVE
#--------------------
dpath = cpath.replace('scripts','data')

output = flow_c;
filename = 'inletFlow'
np.savetxt(dpath+'/'+filename, output, '%.8f')

output = np.array([TB, TV, N]).T;
filename = 'TBTVN'
np.savetxt(dpath+'/'+filename, output, ['%f', '%.8f', '%d'])

output = np.array([[nbr], [1./dt]]).T;
filename = 'nbfs'
np.savetxt(dpath+'/'+filename, output, ['%d', '%f'])

output = N2_c/N2_c.max();
filename = 'N2'
np.savetxt(dpath+'/'+filename, output, '%.8f')
