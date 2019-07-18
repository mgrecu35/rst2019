fname='cm1_SquallLine.Fields2.nc'
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from numpy import *

f=Dataset(fname,'r')
kext=f['kext3d'][:,:,:,1]
z3d=f['z3d'][:,:,:,1]
qc=f['qc'][:,:,:].T
qs=f['qg'][:,:,:].T
plt.contourf(z3d[:,10,:].T,levels=arange(11)*4,vmin=0)
a=nonzero(z3d[:145,10,18:50]>0)
b=nonzero(qs[:145,10,18:50][a]>0.05)
plt.figure()
plt.scatter(z3d[:145,10,18:50][a][b],log10(kext[:145,10,18:50][a][b]))

a=nonzero(z3d[:145,10,:5]>0)
b=nonzero(qs[:145,10,:5][a]<0.05)
plt.figure()
plt.scatter(z3d[:145,10,:5][a][b],log10(kext[:145,10,:5][a][b]))
