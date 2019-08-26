import numpy as np
import matplotlib.pyplot as plt
import os
import glob

from matplotlib import cm

# load data
data_path = os.getcwd().replace('scripts', 'data')
res_data = np.genfromtxt(data_path+'/primary_results')

# unpack
time = res_data[:,0]
sp1  = res_data[:,1]
sp2  = res_data[:,2]
ppl  = res_data[:,3]

# plot
fig = plt.figure()
lw = 0.7

ax = fig.add_subplot(211, ylim=(0,1))

ax.plot(time, sp1, 'm', lw=lw, label='species 1')
# ax.plot(time, sp2, 'r', lw=lw, label='species 2')
ax.set_ylabel('concentration $N_2$')

plt.legend()

ax = fig.add_subplot(212)

ax.plot(time, ppl, 'k', lw=lw)
ax.set_ylabel('pleural pressure [Pa]')
ax.set_xlabel('time $t$')

plt.show()

# save
save_path = os.getcwd().replace('scripts', 'images')
plt.savefig(save_path+'/results.pdf')
