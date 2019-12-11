import numpy as np
import matplotlib.pyplot as plt
import os
import glob

from matplotlib import cm
from scipy.interpolate import interp1d

# choose type:
# Use 'True' if a generic (e.g. sine function) flow profile was used for the simulation. 
# Use 'False' if the flow profile was read from an exisiting B-File
write_new_BFile = False

# file number: (only important if an existing B-File is used as template)
fnbr = 1

# header
BFile_header = 'Time Flow_BTPS O2 CO2 N2 MMss SampleFlow MMms DelayO2 DelayMMss (SPW_V3.2.1;Pre=0;Post=0)'

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
data_path   = os.getcwd().replace('scripts', 'data')
res_data    = np.genfromtxt(data_path+'/primary_results')
flow_data   = np.genfromtxt(data_path+'/inletFlow')
breath_data = np.genfromtxt(data_path+'/TBTVN')
nbfs_data   = np.genfromtxt(data_path+'/nbfs')

# unpack
time = res_data[:,0]
sp1  = res_data[:,1]
sp2  = res_data[:,2]
ppl  = res_data[:,3]

Q    = flow_data[:]

if len(breath_data.shape) < 2:
    TB1 = breath_data[0]
    TV1 = breath_data[1]
    N1  = int(breath_data[2])
else:
    TB1 = breath_data[0,0]
    TV1 = breath_data[0,1]
    N1  = int(breath_data[0,2])


# resample for B-File
fs = 200.
dt = 1./fs

# pre-phase
if write_new_BFile:

    # artificial pre-phase (one pulse)
    Q_pp  = Q[:N1]
    N2_pp = np.ones([1,N1])
    time_pp = time[:N1]

    time_wo_start = TB1

else:
    # read from existing B-File and shorten
    flist = np.sort(glob.glob(data_path+'/lufu/B-Files (for input)/*.txt'))
    filename = flist[fnbr]

    # load data
    data = np.genfromtxt(filename, skip_header=1)
    t  = data[:,0]*1e-03
    flow = data[:,1]*1e-03
    N2   = data[:,4]/N2_nominal

    # resample to simulation time resolution
    t_BF = np.copy(t)
    t = np.arange(t_BF[0], t_BF[-1], 1./nbfs_data[1])

    print(t_BF[-1])

    interpolant_flow = interp1d(t_BF, flow, kind='quadratic')
    flow = interpolant_flow(t)

    interpolant_N2 = interp1d(t_BF, N2, kind='quadratic')
    N2 = interpolant_N2(t)

    # estimate for start of washout
    time_wo_estim = min(t[N2 < 0.5*max(N2)]) - 0.5

    # zero crossings of flow signal
    zr_cr = flow*np.roll(flow,-1) < 0
    zr_cr_ind = np.arange(t.size)[zr_cr] + 1
    zr_cr_ind = zr_cr_ind[:-1] # delete spurious zero crossing from periodic 'roll'
    zr_cr_time = t[zr_cr_ind]
    time_wo_start = min(zr_cr_time[zr_cr_time > time_wo_estim])

    mask_pp = t < time_wo_start

    # data in prephase
    time_pp = t[mask_pp]
    Q_pp  = flow[mask_pp]
    flow = data[:,1]*1e-03
    #flow = Q;
    N2_pp = N2[mask_pp]

time_ppl = time+time_wo_start
time = np.append(time_pp, time+time_wo_start)
N2   = np.append(N2_pp, sp1)
Q    = np.append(Q_pp, Q)

time_BF = np.arange(time[0], time[-1], dt)
N2_BF   = np.interp(time_BF, time, N2)*N2_nominal
Q_BF    = np.interp(time_BF, time, Q)*Q_nominal
time_BF *= 1e03

if write_new_BFile:
    # generate for B-File
    O2_BF    = np.ones_like(time_BF)*O2_nominal
    CO2_BF   = np.ones_like(time_BF)*CO2_nominal
    MMss_BF  = np.ones_like(time_BF)*MMss_nominal
    Qss_BF   = np.ones_like(time_BF)*Qss_nominal
    MMms_BF  = np.ones_like(time_BF)*MMms_nominal
    dO2_BF   = np.ones_like(time_BF)*dO2_nominal
    dMMss_BF = np.ones_like(time_BF)*dMMss_nominal

else:
    # adapt B-File
    O2_BF    = data[:,2][:time_BF.size]
    CO2_BF   = data[:,3][:time_BF.size]
    MMss_BF  = data[:,5][:time_BF.size]
    Qss_BF   = data[:,6][:time_BF.size]
    MMms_BF  = data[:,7][:time_BF.size]
    dO2_BF   = data[:,8][:time_BF.size]
    dMMss_BF = data[:,9][:time_BF.size]

# save B-File
output = np.array([time_BF, Q_BF, O2_BF, CO2_BF, N2_BF, MMss_BF, Qss_BF, MMms_BF, dO2_BF, dMMss_BF]).T
filename = 'B-results.txt'
np.savetxt(data_path+'/lufu/B-Files (results)/'+filename, output, fmt=['%d', '%.8f', '%.5f', '%.5f', '%.5f', '%.5f', '%.8f', '%.5f', '%.5f', '%.5f'], header=BFile_header, comments='', delimiter=' ')

# plot
fig = plt.figure()
lw = 0.7

ax = fig.add_subplot(311, ylim=(-0.05,1.05), xlim=(time[0], time[-1]))

ax.plot(time, N2, 'm', lw=lw, label='N_2')
ax.set_ylabel('conc. $N_2$ (normalized)')

plt.legend()

ax = fig.add_subplot(312, xlim=(time[0], time[-1]))

ax.plot(time_ppl, ppl, 'k', lw=lw)
ax.set_ylabel('pleural pressure [Pa]')
ax.set_xlabel('time $t$')

ax = fig.add_subplot(313, xlim=(time[0], time[-1]))

ax.plot(time, Q*1e03, 'b', lw=lw)
ax.set_ylabel('flow [l/s]')
ax.set_xlabel('time $t$')

plt.show()

# save
save_path = os.getcwd().replace('scripts', 'images')
plt.savefig(save_path+'/results.pdf')
