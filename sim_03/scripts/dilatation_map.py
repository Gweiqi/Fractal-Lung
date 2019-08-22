import numpy as np
import matplotlib.pyplot as plt
import os
import glob

from vtk import vtkUnstructuredGridReader
from vtk.util.numpy_support import vtk_to_numpy

from matplotlib import rc
from matplotlib import cm

from mpl_toolkits.mplot3d.axes3d import Axes3D

# directories
cpath = os.getcwd()

# get filelist
flist = glob.glob(cpath.replace('scripts','data/acinus/*.vtk'))

# get number of acini
reader = vtkUnstructuredGridReader()
reader.SetFileName(flist[0])
reader.Update()
output = reader.GetOutput()
nbrAcini = vtk_to_numpy(output.GetPointData().GetArray('absind_acinus')).max() + 1

# allocate
Ind = np.zeros([np.size(flist), nbrAcini])
Dil = np.zeros([np.size(flist), nbrAcini])
T   = np.zeros([np.size(flist), nbrAcini])
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

    # get absolute index
    ind = output.GetPointData().GetArray('absind_acinus')
    ind = vtk_to_numpy(ind)
    Ind[i,:], sel = np.unique(ind, return_index=True)

    # get dilatation
    dil = output.GetPointData().GetArray('dilatation_acinus')
    dil = vtk_to_numpy(dil)
    Dil[i,:] = dil[sel];
    
    # update time
    T[i,:] = np.ones(nbrAcini)*(i+1)*0.01
    
    i = i+1


# use latex and select 'Helvetica'
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

# plots
fig = plt.figure(1)
lw = 0.5
fs = 12

ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(T, Ind, Dil, rstride=10, cstride=10, alpha=0.1)
ax.set_xlabel('time $t$ $[s]$', fontsize=fs)
ax.set_ylabel('lobule', fontsize=fs)
ax.set_zlabel('lobule dilatation', fontsize=fs)

fig = plt.figure(2)

#v = np.linspace(-Dil.max(), Dil.max(), 21)
ax = fig.add_subplot(111)

cb = ax.contourf(T, Ind, Dil, 20, cmap = cm.cool)
ax.set_xlabel('time $t$ $[s]$', fontsize=fs)
ax.set_ylabel('lobule')
plt.colorbar(cb)


plt.show()

# save
filename = cpath.replace('scripts','images/dilatation_map.pdf')
plt.savefig(filename, bbox_inches='tight')



