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
flist = glob.glob(cpath.replace('scripts','data/duct/*.vtk'))
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

    # get concentration I
    c1 = output.GetPointData().GetArray('concentrationI_duct')
    c1 = vtk_to_numpy(c1)

    # get velocity
    v1 = output.GetCellData().GetArray('velocity_duct')
    v1 = vtk_to_numpy(v1)

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


ax = fig.add_subplot(111, ylim=(0,1))

ax.plot(t,v,'k', linewidth=lw)
ax.plot(np.linspace(t[0], t[-1], N2.size),N2,'m', linewidth=0.5)
ax.set_xlabel('time $t$ $[s]$', fontsize=fs)
ax.set_ylabel('inert gas washout curve $c_\mathrm{in}$ $[-]$', fontsize=fs)


plt.show()

# save
filename = cpath.replace('scripts','images/washout.pdf')
plt.savefig(filename, bbox_inches='tight')
