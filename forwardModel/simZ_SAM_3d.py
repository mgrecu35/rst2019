from netCDF4 import Dataset
from numpy import *
import pyHB2 as pyHB2
pyHB2.initt()

f=Dataset('../../GATE_IDEAL_extracted_50m_0.5s_cu0.5_1296_0040-0600.nc',
          'r')

stop
qv=f['vapor'][-1,:,:,:]*1e-3    # water vapor
qr=f['rain'][-1,:,:,:]*1e-3     # rain mixing ratio
qs=f['snow'][-1,:,:,:]*1e-3     # snow mixing ratio
qag=f['aggregates'][-1,:,:,:]*1e-3     # aggregates mixing ratio
ncag=f['agg_concen_kg'][-1,:,:,:]      # aggregates number concentration
qc=f['cloud'][-1,:,:,:]*1e-3    # cloud mixing ratio
ncs=f['snow_concen_kg'][-1,:,:,:]          # snow number concentration
ncg=f['graup_concen_kg'][-1,:,:,:]          # graupel number concentration
nch=f['hail_concen_kg'][-1,:,:,:]          # graupel number concentration
ncr=f['rain_concen_kg'][-1,:,:,:]          # rain number concentration
qg=f['graupel'][-1,:,:,:]*1e-3   # graupel mixing ratio
qh=f['hail'][-1,:,:,:]*1e-3   # graupel mixing ratio
z=f['z_coords'][:]/1000.             # height (km)
th=f['theta'][-1,:,:,:]    # potential temperature (K)
prs=f['press'][-1,:,:,:]*100.  # pressure (Pa)
T=th*(prs/100000)**0.286  # Temperature

R=287.058  #J*kg-1*K-1
rho=prs/(R*T)

w=f['w'][-1,:,:,:]
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
nx2=330
nx1=130
ny1=190
ny2=210

x1=f['x_coords'][nx1:nx2]/1000.
y1=f['y_coords'][ny1:ny2]/1000.

from fields3D import *

z3d,kext3d,salb3d,asym3d,v3d=radarFields_3d(nx1,nx2,ny1,ny2,qs,qg,qh,qag,qr,qc,qv,T,prs,ncs,ncg,nch,ncag,ncr,rho,wm,z,\
                                 dm_z,vt_R,vt_S,pyHB2,nz,freqs)

import xarray as xr
z3dx=xr.DataArray(z3d,dims=['nx','ny','nz','nfreq'])
kext3dx=xr.DataArray(kext3d,dims=['nx','ny','nz','nfreq'])
salb3dx=xr.DataArray(salb3d,dims=['nx','ny','nz','nfreq'])
asym3dx=xr.DataArray(asym3d,dims=['nx','ny','nz','nfreq'])
v3dx=xr.DataArray(v3d,dims=['nx','ny','nz','nfreq'])
zx=xr.DataArray(z,dims=['nz'])
xx=xr.DataArray(f['x_coords'][nx1:nx2]/1000.,dims=['nx'])
yx=xr.DataArray(f['y_coords'][ny1:ny2]/1000.,dims=['ny'])
d=xr.Dataset({'z3d':z3dx,'kext3d':kext3dx,'salb3d':salb3dx,'asym3d':asym3dx,'v3d':v3dx,'z':zx,'x':xx,'y':yx})
d.to_netcdf('ramsZFields.nc')
