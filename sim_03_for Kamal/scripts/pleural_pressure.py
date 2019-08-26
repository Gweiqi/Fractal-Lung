import numpy as np
import matplotlib.pyplot as plt
import os
import glob

from vtk import vtkUnstructuredGridReader
from vtk.util.numpy_support import vtk_to_numpy

from matplotlib import rc

# directories
cpath = os.getcwd()

# get filelist
flist = glob.glob(cpath.replace('scripts','data/acinus/*.vtk'))
flist = flist[:500]

t = np.zeros(np.size(flist))
v = np.zeros(np.size(flist))
i = 0
for filename in flist:

    # read dat
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    output = reader.GetOutput()

    # get nodes
    nodes = output.GetPoints().GetData()
    nodes = vtk_to_numpy(nodes)

    # get pleural pressure
    c1 = output.GetPointData().GetArray('pleural_pressure_acinus')
    c1 = vtk_to_numpy(c1)

    # pick value
    v[i] = c1[0]

    t[i] = (i+1)*0.001

    i = i+1


# N2 washout from lufu
filename = cpath.replace('scripts','data/N2')
N2 = np.genfromtxt(filename)

# use latex and select 'Helvetica'
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

# plots
fig = plt.figure()
lw = 0.5
fs = 18


ax = fig.add_subplot(111)

ax.plot(t,v,'k', linewidth=lw)
ax.set_xlabel('time $t$ $[s]$', fontsize=fs)
ax.set_ylabel('pleural pressure $p_\mathrm{pl}$ $[Pa]$', fontsize=fs)


plt.show()

# save
filename = cpath.replace('scripts','images/pleural_pressure.pdf')
plt.savefig(filename, bbox_inches='tight')
