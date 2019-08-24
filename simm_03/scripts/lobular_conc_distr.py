import numpy as np
import matplotlib.pyplot as plt
import os
import glob

from vtk import vtkUnstructuredGridReader
from vtk.util.numpy_support import vtk_to_numpy

from matplotlib import cm

# define name of lobule varibale
lob_var_name = 'concentrationI_acinus'

# directories
cpath = os.getcwd()

# get filelist
flist = glob.glob(cpath.replace('scripts','data/lobule/reg/*.vtk'))

# sampling frequeny
fs = 50.

# read lobular data
t_ind = -1
reader = vtkUnstructuredGridReader()
reader.SetFileName(flist[t_ind])
reader.Update()
output = reader.GetOutput()

# get number of lobules
nbrLobules = vtk_to_numpy(output.GetPointData().GetArray('absind_acinus')).max() + 1
nbrLobules = int(nbrLobules)

# get node coordinates
node_coords = output.GetPoints().GetData()
node_coords = vtk_to_numpy(node_coords)

# get absolute index
ind = output.GetPointData().GetArray('absind_acinus')
ind = vtk_to_numpy(ind)
ind, sel = np.unique(ind, return_index=True)
ind  = ind.astype(np.int32)

# get specific lobule variable
var = output.GetPointData().GetArray(lob_var_name)
var = vtk_to_numpy(var)
var_ = np.zeros_like(ind).astype(np.float32)
for s, i in zip(sel, range(ind.size)):
    var_[i] = var[s:s+51].mean()

varx = np.copy(var)
var = var_;


# plots
fig = plt.figure(1)
lw = 0.7
ms = 1.

ax = fig.add_subplot(111, ylim=(0,1))
ax.plot(ind, var, 'ok', ms=ms)
ax.set_xlabel('lobule index')
ax.set_ylabel('mean lob. conc.')

plt.show()

# save
filename = cpath.replace('scripts','images/lobular_conc_distr.pdf')
#plt.savefig(filename, bbox_inches='tight')
