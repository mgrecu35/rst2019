from netCDF4 import Dataset
from numpy import *
import pyHB2 as pyHB2
pyHB2.initt()

fname='/media/grecu/ExtraDrive1/IPHEX_0611/wrfout_d03_2014-06-11_19:00:00'

f=Dataset(fname,'r')
it=2
qv=f['QVAPOR'][it,:,:,:]    # water vapor
qr=f['QRAIN'][it,:,:,:]     # rain mixing ratio
qs=f['QSNOW'][it,:,:,:]     # snow mixing ratio
qc=f['QCLOUD'][it,:,:,:]    # cloud mixing ratio
qg=f['QGRAUP'][it,:,:,:]   # graupel mixing ratio
qh=f['QGRAUP'][it,:,:,:]*0.   # graupel mixing ratio
ncr=f['QNRAIN'][it,:,:,:]*5.     # rain mixing ratio
ncs=f['QNSNOW'][it,:,:,:]*5.     # snow mixing ratio
ncg=f['QNGRAUPEL'][it,:,:,:]*5.   # graupel mixing ratio
nch=f['QNGRAUPEL'][it,:,:,:]*0.   # graupel mixing ratio
#z=f['z_coords'][:]/1000.             # height (km)
th=f['T'][it,:,:,:]+300    # potential temperature (K)
prs=f['P'][it,:,:,:]+f['PB'][it,:,:,:]  # pressure (Pa)
T=th*(prs/100000)**0.286  # Temperature
#stop
z=(f['PHB'][it,:,:,:]+f['PH'][it,:,:,:])/9.81/1000.

R=287.058  #J*kg-1*K-1
rho=prs/(R*T)

w=f['W'][it,:,:,:]
#wm=0.5*(w[1:,:,:]+w[:-1,:,:])
wm=w[:,:,:]
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
nx2=220
nx1=90
ny1=150
ny2=180

#x1=f['x_coords'][nx1:nx2]/1000.
#y1=f['y_coords'][ny1:ny2]/1000.
import matplotlib.pyplot as plt
import matplotlib.colors as col
plt.pcolormesh(qr[0,ny1:ny2,nx1:nx2],norm=col.LogNorm())
plt.colorbar()

from fields3D_WRF2M import *

z3d,kext3d,salb3d,asym3d,\
    v3d,dnr,dns,dng,\
    rainDBL,snowDBL=radarFields_3d(nx1,nx2,ny1,ny2,qs,qg,\
                                   qh,qr,qc,qv,T,prs,ncs,ncg,nch,ncr,\
                                   rho,wm,z,dm_z,vt_R,vt_S,pyHB2,nz,freqs)


plt.figure()
#plt.contourf(arange(60),z[1:,60,60],z3d[:,34,:,0].T,vmin=0,levels=arange(13)*5,cmap='jet')
plt.contourf(arange(nx1,nx2),z[1:,60,60],z3d[:,10,:,0].T,vmin=10,levels=10+arange(21)*2,cmap='jet')
plt.ylim(0,12)
plt.colorbar()
#plt.figure()
#plt.plot(z3d[30,26,:,0],z[1:,60,60]/1000.)
#plt.xlim(0,60)
#plt.ylim(0,12)
#stop
import xarray as xr
z3dx=xr.DataArray(z3d,dims=['nx','ny','nz','nfreq'])
kext3dx=xr.DataArray(kext3d,dims=['nx','ny','nz','nfreq'])
salb3dx=xr.DataArray(salb3d,dims=['nx','ny','nz','nfreq'])
asym3dx=xr.DataArray(asym3d,dims=['nx','ny','nz','nfreq'])
v3dx=xr.DataArray(v3d,dims=['nx','ny','nz','nfreq'])
z1=z[:,ny1:ny2,nx1:nx2]
zx=xr.DataArray(z1,dims=['nz1','ny','nx'])
qrx=xr.DataArray(qr[:,ny1:ny2,nx1:nx2]*rho[:,ny1:ny2,nx1:nx2]*1e3,dims=['nz','ny','nx'])
qsx=xr.DataArray(qs[:,ny1:ny2,nx1:nx2]*rho[:,ny1:ny2,nx1:nx2]*1e3,dims=['nz','ny','nx'])
qgx=xr.DataArray(qg[:,ny1:ny2,nx1:nx2]*rho[:,ny1:ny2,nx1:nx2]*1e3,dims=['nz','ny','nx'])
qcx=xr.DataArray(qc[:,ny1:ny2,nx1:nx2]*rho[:,ny1:ny2,nx1:nx2]*1e3,dims=['nz','ny','nx'])

dnrx=xr.DataArray(dnr,dims=['nx','ny','nz'])
dnsx=xr.DataArray(dns,dims=['nx','ny','nz'])
dngx=xr.DataArray(dng,dims=['nx','ny','nz'])

d=xr.Dataset({'z3d':z3dx,'kext3d':kext3dx,'salb3d':salb3dx,'asym3d':asym3dx,\
              'v3d':v3dx,'z':zx, 'qr':qrx, 'qs':qsx, 'qg':qgx, 'qc':qcx, 'dnr':dnrx, \
              'dns':dnsx, 'dng':dngx})
d.to_netcdf('wrfJune11_2014.Fields_dn5.nc')
