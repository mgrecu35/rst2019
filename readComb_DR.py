from netCDF4 import Dataset
import pyHB2 as pyHB2
pyHB2.initt()

fname='2B.GPM.DPRGMI.CORRA2018.20141013-S025456-E042728.003539.V06A.HDF5'
fname='2B.GPM.DPRGMI.CORRA2018.20181001-S002201-E015435.026079.V06A.HDF5'
fpath='/media/grecu/ExtraDrive1/DPR/'
#fpath='/gpmdata/2014/10/13/radar/'
fpath='/home/grecu/forwardSim/2BCMB/'
fh=Dataset(fpath+fname)
n1=12
n2=37
Lat=fh['NS']['Latitude'][:,n1:n2]
Lon=fh['NS']['Longitude'][:,n1:n2]
tb=fh['DiagGroup/dset1'][:,n1:n2,:]

import geopandas
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from numpy import *

#plt.plot(Lon[:,24],Lat[:,24])
a=nonzero((Lat[:,12]+2)*(Lat[:,12]-2)<0)
a=array([arange(1900,2010)])


tb=tb[a[0],:,:]
# = Basemap(width=1900000,height=1900000,
#           rsphere=(6378137.00,6356752.3142),\
#           resolution='l',area_thresh=1000.,projection='lcc',\
#           lat_1=-3,lat_2=3,lat_0=0,lon_0=95.)
m = Basemap(llcrnrlon=96.,llcrnrlat=-2.,urcrnrlon=100.,urcrnrlat=2.,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=2.,lon_0=98.,lat_ts=-2.)


#plt.scatter(x,y)

#stop
sfcRain=fh['MS/surfPrecipTotRate'][a[0],:]
#plt.contour(x,y,sfcRain,[0.01,0.1,1])
qv=fh['MS']['vaporDensity'][a[0],:,:]
qv_sfc=fh['MS/surfaceVaporDensity'][a[0],:]
simTb_CMB=fh['MS/simulatedBrightTemp'][a[0],:,:]
envNode=fh['MS/envParamNode'][a[0],:,:]
sfcTemp=fh['MS/skinTemperature'][a[0],:]
airTemp=fh['MS/airTemperature'][a[0],:,:]
press=fh['MS/airPressure'][a[0],:,:]
sfcEmiss=fh['MS/surfEmissivity'][a[0],:,:]
w10=fh['MS/tenMeterWindSpeed'][a[0],:]
pType=(fh['MS/Input/precipitationType'][a[0],:]/1e7).astype(int)
binNodes=fh['MS/phaseBinNodes'][a[0],:,:]
Nw=fh['MS/precipTotPSDparamLow'][a[0],:,:]-log10(8e6)
psdNodes=fh['MS/PSDparamLowNode'][a[0],:,:]
pRate=fh['MS/precipTotRate'][a[0],:,:]
emiss=fh['MS/surfEmissivity'][a[0],:,:]
z13=fh['MS/correctedReflectFactor'][a[0],:,:,0]

j=12
qvInt2D=zeros((110,88),float)-99
for i in range(110):
    bins=arange(envNode[i,j,0],envNode[i,j,-1]+1)
    qvInt=interp(bins,envNode[i,j,:],qv[i,j,:])
    qvInt2D[i,bins]=qvInt
#stop
m.drawcoastlines()
x,y=m(Lon[a[0],:],Lat[a[0],:])
plt.pcolormesh(x[00:,:],y[00:,:],tb[:,:,12],cmap='jet',vmin=136,vmax=270)

nx,ny=sfcRain.shape
# absair,abswv = gasabsr98(f,tk,rhowv,pa,ireturn)
#tb = radtran(umu,nlyr,btemp,lyrtemp,lyrhgt,kext,salb,asym,fisot,emis,ebar)
# emis,ebar = emit(f,npol,ts,w,umu)
# Marti,26Feb!2019
# Joi,21Mart!2019
f=190.3
umu=cos(53/180.*pi)
tbL=[]
tbL2=[]
tbL3=[]
freqs=[10.6,10.6,18.7,18.7,23.,37,37.,89,89.,166.,166.,186.3,190.3]
npol=[1,0,1,0,1,1,0,1,0,1,0,1,1]
iFreq=[0,0,1,1,2,3,3,4,4,5,5,6,6,7,7]
nfreq=8
tbSimM=zeros((110,n2-n1,13),float)
tbSimM0=zeros((110,n2-n1,13),float)
from rteModule import *
dfL=[]
s1=[]
s2=[]
kL=[]
for i in range(1,110):
    for j in range(ny):
        kextH=zeros((88,8),float)
        salbH=zeros((88,8),float)
        asymH=zeros((88,8),float)
        if pType[i,j]>=1:
            tbSim,s1L,s2L,kfL=rte(psdNodes,binNodes,Nw,pRate,z13,i,j,emiss,\
                      qv,airTemp,press,envNode,sfcTemp,umu,pyHB2,1.)

            tbSim_noad=rte_noad(psdNodes,binNodes,Nw,pRate,z13,i,j,emiss,\
                      qv,airTemp,press,envNode,sfcTemp,umu,pyHB2)
            dtb=array([array(tbSim)[-6:]-array(tbSim_noad)[-6:]])
            pinv=linalg.pinv(dot(dtb.T,dtb)+eye(6)*4)
            kgain=dot(dtb,pinv)
            e=tb[i,j,7:13]-array(tbSim_noad)[-6:]
            e[tb[i,j,7:13]<0]=0
            df=dot(kgain,e)[0]
            dfL.append(df)
            if df<-3:
                df=-3.
            if df>3:
                df=3.
            #df=2.
            #df=0.
            tbSim,s1L,s2L,kfL=rte(psdNodes,binNodes,Nw,pRate,z13,i,j,emiss,\
                      qv,airTemp,press,envNode,sfcTemp,umu,pyHB2,df)
            tbSimM[i,j,:]=tbSim
            tbSimM0[i,j,:]=tbSim_noad
            s1.extend(s1L)
            s2.extend(s2L)
            kL.extend(kfL)
            if binNodes[i,j,3]<binNodes[i,j,4] and z13[i,j,binNodes[i,j,3]]>10:
                s1.append(pRate[i,j,binNodes[i,j,3]])
                s2.append(pRate[i,j,binNodes[i,j,3]])
                kL.append(binNodes[i,j,3])
                
            if simTb_CMB[i,j,8]>0 and simTb_CMB[i,j,12]>0:
                tbL.append(tbSim)
                tbL2.append(tb[i,j,:13])
                tbL3.append(simTb_CMB[i,j,:13])
            
tbL=array(tbL)
tbL3=array(tbL3)
tbL2=array(tbL2)
s1=array(s1)
s2=array(s2)
kL=array(kL)
profL=[]
for k in range(kL.min(),kL.max()):
    a=nonzero(kL==k)
    print(s1[a].mean(), s2[a].mean(), k)
    profL.append([s1[a].mean(), s2[a].mean(), k])
    
print(corrcoef(tbL[:,12],tbL3[:,12]))

import xarray as xr
simTb=xr.DataArray(array(tbL),dims=['np','nf'])
simTb_x=xr.Dataset({'simTb_p':simTb})
simTb_x.to_netcdf('simTb.nc')
plt.figure()
tbSimMm=ma.array(tbSimM,mask=tbSimM<0.01)
m.drawcoastlines()
plt.pcolormesh(x[00:110,:],y[00:110,:],tbSimMm[00:,:,12],cmap='jet',vmin=136,vmax=270)
plt.colorbar()
plt.figure()
tbSimMm=ma.array(tbSimM0,mask=tbSimM0<0.01)
m.drawcoastlines()
plt.pcolormesh(x[00:110,:],y[00:110,:],tbSimMm[00:,:,12],cmap='jet',vmin=136,vmax=270)
plt.colorbar()
#a2=nonzero(array(tbL)[:,0]>0)

plt.figure()
tbSimMm=ma.array(simTb_CMB,mask=tbSimM0<0.01)
m.drawcoastlines()
plt.pcolormesh(x[00:110,:],y[00:110,:],tbSimMm[00:,:,12],cmap='jet',vmin=136,vmax=270)
plt.colorbar()
#a2=nonzero(array(tbL)[:,0]>0)
#print(corrcoef(array(tbL)[:,0][a2],array(tbL)[:,2][a2]))
#plt.figure()
#plt.scatter(array(tbL)[:,0][a2],array(tbL)[:,2][a2])
tbLm=array([190.54020808, 124.00596947, 231.75368474, 192.64629731,\
            264.70775687, 244.63579192, 221.3298183 , 236.30827453,\
            234.02888285, 211.69668184, 211.69595682, 230.49217778,\
            221.65313326])
#plt.figure()
#plt.scatter(tbL2[:,10],tbL3[:,10])
#plt.plot([100,250],[100,250])

plt.figure()
profL=array(profL)
plt.plot(profL[:,1],profL[:,2][::-1])
plt.plot(profL[:,0],profL[:,2][::-1])
