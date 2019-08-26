import numpy as np
import matplotlib.pyplot as plt

from numpy import log, tan, pi, sqrt, exp
from numpy.random import lognormal

import os

# length of mod-table. This corresponds with the total number of terminal duct / total number of trumpet lobules, and depends on the value set for d_Lim in the constant/systemProperties file. The total number of terminal ducts is prompted to the console/termial when the flPROG application is run.
N = 304
aw_ind = np.arange(N) # ID of airway unit

#---------------------
# initialize mod. parameter
#---------------------

# duct transmissibility mod. parameter
xi = np.ones_like(aw_ind, dtype='f')

# lobule residual volume mod. parameter
theta = np.ones_like(aw_ind, dtype='f')

# lobule compliance mod. parameter
phi = np.ones_like(aw_ind, dtype='f')

# lobule resistance mod. parameter
tau = np.ones_like(aw_ind, dtype='f')


#---------------------
# change mod. parameter: conctant increase
#---------------------
#stretch_width = 0.05
#phi = -2*stretch_width/N*aw_ind + 1 + stretch_width

#---------------------
# change mod. parameter: log-normal distribution
#---------------------
#mu, sigma = 0.0, 0.4
#phi = lognormal(mu, sigma, N)

#---------------------
# change mod. parameter: local, regional
#---------------------
mod_param = 1.0 # modification parameter value
f = 0.1        # fraction of modified airway units
for k in aw_ind:

    # either grouped (regional modification) ...
    
    if k < f*N:
        phi[k] *= mod_param

    if k > (1.-f)*N:
        phi[k] *= 1.+ (1. - mod_param)
    

    # ... or distributed (local modification)
    # if np.mod(k, np.int(1/f)) == 0:
    #    xi[k] *= mod_param


# generate table
mod_table = np.array([aw_ind, xi, theta, phi, tau]).T

# some plots
fig = plt.figure()

ax = fig.add_subplot(111, xlim=(0,2))

ax.hist(phi, round(N/5))
ax.set_xlabel('$\phi$')
ax.set_ylabel('prob')

plt.show()

# save
data_path = os.getcwd().replace('scripts', 'data')
np.savetxt(data_path+'/modifyLung', mod_table, fmt='%d %f %f %f %f')
