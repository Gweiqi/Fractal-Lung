import numpy as np
import matplotlib.pyplot as plt

from numpy import log, tan, pi, sqrt, exp
from numpy.random import lognormal

import os

# length of mod-table
N = 304
abs_ind = np.arange(N)

# duct transmissibility mod factor
duct_trans = np.ones_like(abs_ind)

# acinus volume mod factor
acin_vol = np.ones_like(abs_ind, dtype='f')

# acinus stretch value (phi)
#phi = np.ones_like(abs_ind, dtype='f')
stretch_width = 0.05
#phi = -2*stretch_width/N*abs_ind + 1 + stretch_width

mu, sigma = 0.0, 0.4
phi = lognormal(mu, sigma, N)

# acinus resistance mod factor
acin_res = np.ones_like(abs_ind, dtype='f')

fac = 1.0 # constant modification factor
f = 0.2 # fraction of modified values
for k in abs_ind:

    """
    # either grouped (lower index side) ...
    if k < f*N:
        acin_res[k] *= fac

    # ... or distributed
    if np.mod(k, np.int(1/f)) == 0:
        acin_res[k] *= fac
    """


# generate table
mod_table = np.array([abs_ind, duct_trans, acin_vol, phi, acin_res]).T

# some plots
fig = plt.figure()

ax = fig.add_subplot(111, xlim=(0,2))

ax.hist(phi, N/5)
ax.set_xlabel('$\phi$')
ax.set_ylabel('prob')

plt.show()

# save
data_path = os.getcwd().replace('scripts', 'data')
np.savetxt(data_path+'/modifyLung', mod_table, fmt='%d %f %f %f %f')
