from netCDF4 import Dataset
from numpy import *
import pyHB2 as pyHB2
pyHB2.initt()

fname='/home/grecu/cm1r19.8/run/cm1out.nc'

f=Dataset(fname,'r')
it=12
qv=f['qv'][it,:,:,:]    # water vapor
qr=f['qr'][it,:,:,:]     # rain mixing ratio
qs=f['qs'][it,:,:,:]     # snow mixing ratio
qc=f['qc'][it,:,:,:]    # cloud mixing ratio
qg=f['qg'][it,:,:,:]   # graupel mixing ratio
qh=f['qg'][it,:,:,:].copy()*0.   # graupel mixing ratio
ncs=f['ncs'][it,:,:,:]    # snow number concentration
ncg=f['ncg'][it,:,:,:]    # graupel number concentration
ncr=f['ncr'][it,:,:,:]    # rain number concentration
nch=ncg.copy()*0.
#z=f['z_coords'][:]/1000.             # height (km)
th=f['th'][it,:,:,:]   # potential temperature (K)
prs=f['prs'][it,:,:,:]  # pressure (Pa)
T=th*(prs/100000)**0.286  # Temperature
z=f['zf'][:]
w=f['w'][it,:,:,:]



R=287.058  #J*kg-1*K-1
rho=prs/(R*T)


wm=0.5*(w[1:,:,:]+w[:-1,:,:])
#wm=w[:,:,:]
nz,ny,nx=qv.shape
n0w=0.
nfreq=8
zKu_L=[]
zKa_L=[]
zW_L=[]
vKu_L=[]
vKa_L=[]
vW_L=[]
j=200



x1L=[]
from getNw import *
dr=0.25

fd=Dataset('nw_dm_dsd/dm_z.nc')
dm_z=fd['dm_z'][:,:]

fd=Dataset('nw_dm_dsd/doppler_2.nc')
vt_R=fd['rain'][:,:]
vt_S=fd['snow'][:,:]
ls=open('TablesN/BlueWhiteOrangeRed.rgb.txt','r').readlines()

#ls=open('BlueDarkRed18.rgb.txt','r').readlines()
rgbL=[]
for l in ls:
    v1=[float(v)/255. for v in l.split()[0:3]]
    rgbL.append(v1)
import matplotlib.pyplot as plt
import matplotlib.colors as col

#rgbL=array(rgbL)
from matplotlib.colors import LinearSegmentedColormap
cmp=LinearSegmentedColormap.from_list('ncar',rgbL[5:-5],N=125)


#z_clw = gcloud(freqy,t,clw) #
freqs=[13.8,35.5,94.]
#x1=f['x_coords'][nx1:nx2]/1000.
#y1=f['y_coords'][ny1:ny2]/1000.
import matplotlib.pyplot as plt
import matplotlib.colors as col




from fields3D_CM1 import *
nx1=120
nx2=320
ny1=5
ny2=25



z3d,kext3d,salb3d,asym3d,\
    v3d,dnr,dns,dng,rainDBL,snowDBL=radarFields_3d(nx1,nx2,ny1,ny2,qs,qg,qh,qr,qc,qv,T,prs,\
                                           ncs,ncg,nch,ncr,rho,wm,z,\
                                           dm_z,vt_R,vt_S,pyHB2,nz,freqs)

plt.figure()
rainDBL=array(rainDBL)
plt.semilogx(rainDBL[:,1],rainDBL[:,3],'*')

plt.figure()
snowDBL=array(snowDBL)
plt.semilogx(snowDBL[:,1],snowDBL[:,3],'*')

#stop
plt.figure()
#plt.contourf(arange(60),z[1:,60,60],z3d[:,34,:,0].T,vmin=0,levels=arange(13)*5,cmap='jet')

plt.contourf(arange(200),z[1:],z3d[:,0,:,0].T,levels=arange(30)*2+5,cmap='jet')
plt.colorbar()

plt.figure()
plt.contourf(arange(200),z[1:],dnr[:,0,:].T,cmap='RdBu')

import xarray as xr
z3dx=xr.DataArray(z3d,dims=['nx','ny','nz','nfreq'])
kext3dx=xr.DataArray(kext3d,dims=['nx','ny','nz','nfreq'])
salb3dx=xr.DataArray(salb3d,dims=['nx','ny','nz','nfreq'])
asym3dx=xr.DataArray(asym3d,dims=['nx','ny','nz','nfreq'])
v3dx=xr.DataArray(v3d,dims=['nx','ny','nz','nfreq'])

qrx=xr.DataArray(qr[:,ny1:ny2,nx1:nx2]*rho[:,ny1:ny2,nx1:nx2]*1e3,dims=['nz','ny','nx'])
qsx=xr.DataArray(qs[:,ny1:ny2,nx1:nx2]*rho[:,ny1:ny2,nx1:nx2]*1e3,dims=['nz','ny','nx'])
qgx=xr.DataArray(qg[:,ny1:ny2,nx1:nx2]*rho[:,ny1:ny2,nx1:nx2]*1e3,dims=['nz','ny','nx'])
qcx=xr.DataArray(qc[:,ny1:ny2,nx1:nx2]*rho[:,ny1:ny2,nx1:nx2]*1e3,dims=['nz','ny','nx'])

dnrx=xr.DataArray(dnr,dims=['nx','ny','nz'])
dnsx=xr.DataArray(dns,dims=['nx','ny','nz'])
dngx=xr.DataArray(dng,dims=['nx','ny','nz'])

zx=xr.DataArray(z,dims=['nz1'])
d=xr.Dataset({'z3d':z3dx,'kext3d':kext3dx,'salb3d':salb3dx,'asym3d':asym3dx,\
              'v3d':v3dx,'z':zx, 'qr':qrx, 'qs':qsx, 'qg':qgx, 'qc':qcx, 'dnr':dnrx, \
              'dns':dnsx, 'dng':dngx})
d.to_netcdf('cm1_SquallLine.Fields2.nc')
import pickle
pickle.dump([rainDBL,snowDBL],open('kz.pklz','wb'))
