import numpy as np
import matplotlib.pyplot as plt
import os
import glob

from matplotlib import cm

write_new_BFile = False

# file number: (only important if an existing B-File is used as template)
fnbr = 1

# header
BFile_header = 'Time Flow_BTPS O2 CO2 N2 MMss SampleFlow MMms DelayO2 DelayMMss'

# nominal magnitude for results BF-File
Q_nominal     =  1e03
O2_nominal    = 20.9
CO2_nominal   =  4.5
N2_nominal    = 78.1
MMss_nominal  = 28.5
Qss_nominal   = -0.0035
MMms_nominal  = 28.9
dO2_nominal   =  0.56
dMMss_nominal =  0.64

# load data
data_path = os.getcwd().replace('scripts', 'data')
res_data  = np.genfromtxt(data_path+'/primary_results')
flow_data = np.genfromtxt(data_path+'/inletFlow')

# unpack
time = res_data[:,0]
sp1  = res_data[:,1]
sp2  = res_data[:,2]
ppl  = res_data[:,3]

Q    = flow_data[:]

# resample for B-File
fs = 200.
dt = 1./fs

time_BF = np.arange(time[0], time[-1], dt)
N2_BF   = np.interp(time_BF, time, sp1)*N2_nominal
Q_BF    = np.interp(time_BF, time, Q)*Q_nominal
time_BF *= 1e03

if write_new_BFile:
    # generate for BFile
    O2_BF    = np.ones_like(time_BF)*O2_nominal
    CO2_BF   = np.ones_like(time_BF)*CO2_nominal
    MMss_BF  = np.ones_like(time_BF)*MMss_nominal
    Qss_BF   = np.ones_like(time_BF)*Qss_nominal
    MMms_BF  = np.ones_like(time_BF)*MMms_nominal
    dO2_BF   = np.ones_like(time_BF)*dO2_nominal
    dMMss_BF = np.ones_like(time_BF)*dMMss_nominal

else:
    # read from existing B-File and shorten
    flist = np.sort(glob.glob(data_path+'/lufu/B-Files (for input)/*.txt'))
    filename = flist[fnbr-1]

    # load data
    data     = np.genfromtxt(filename, skip_header=1)

    O2_BF    = data[:,2][:time_BF.size]
    CO2_BF   = data[:,3][:time_BF.size]
    MMss_BF  = data[:,5][:time_BF.size]
    Qss_BF   = data[:,6][:time_BF.size]
    MMms_BF  = data[:,7][:time_BF.size]
    dO2_BF   = data[:,8][:time_BF.size]
    dMMss_BF = data[:,9][:time_BF.size]

# save B-File
output = np.array([time_BF, Q_BF, O2_BF, CO2_BF, N2_BF, MMss_BF, Qss_BF, MMms_BF, dO2_BF, dMMss_BF]).T;
filename = 'B-results.txt'
np.savetxt(data_path+'/lufu/B-Files (results)/'+filename, output, fmt=['%d', '%.8f', '%.5f', '%.5f', '%.5f', '%.5f', '%.8f', '%.5f', '%.5f', '%.5f'], header=BFile_header, comments='', delimiter=' ')

# plot
fig = plt.figure()
lw = 0.7

ax = fig.add_subplot(311, ylim=(0,1))

ax.plot(time, sp1, 'm', lw=lw, label='species 1')
ax.plot(time, sp2, 'r', lw=lw, label='species 2')
ax.set_ylabel('conc. $N_2$ (normalized)')

plt.legend()

ax = fig.add_subplot(312)

ax.plot(time, ppl, 'k', lw=lw)
ax.set_ylabel('pleural pressure [Pa]')
ax.set_xlabel('time $t$')

ax = fig.add_subplot(313)

ax.plot(time, Q*1e03, 'b', lw=lw)
ax.set_ylabel('flow [l/s]')
ax.set_xlabel('time $t$')

plt.show()

# save
save_path = os.getcwd().replace('scripts', 'images')
plt.savefig(save_path+'/results.pdf')
