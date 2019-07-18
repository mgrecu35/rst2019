from netCDF4 import Dataset
import pyHB2 as pyHB2
pyHB2.initt()

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from numpy import *

fname='2B.GPM.DPRGMI.CORRA2018.20141013-S025456-E042728.003539.V06A.HDF5'
fpath='/media/grecu/ExtraDrive1/DPR/'
fpath='/gpmdata/2014/10/13/radar/'

def cmbFile(fpath,fname):
    fh=Dataset(fpath+fname)
    Lat=fh['NS']['Latitude'][:,:]
    Lon=fh['NS']['Longitude'][:,:]
    tb=fh['DiagGroup/dset1'][:,:,:]
    a=nonzero((Lon[:,24]+180)*(Lon[:,24]+0)<0)
    tb=tb[a[0],:,:]
    sfcRain=fh['NS/surfPrecipTotRate'][a[0],:]
    Lat=fh['NS']['Latitude'][a[0],:]
    Lon=fh['NS']['Longitude'][a[0],:]
    qv=fh['NS']['vaporDensity'][a[0],:,:]
    qv_sfc=fh['NS/surfaceVaporDensity'][a[0],:]
    simTb=fh['NS/simulatedBrightTemp'][a[0],:,:]
    envNode=fh['NS/envParamNode'][a[0],:,:]
    sfcTemp=fh['NS/skinTemperature'][a[0],:]
    airTemp=fh['NS/airTemperature'][a[0],:,:]
    press=fh['NS/airPressure'][a[0],:,:]
    sfcEmiss=fh['NS/surfEmissivity'][a[0],:,:]
    w10=fh['NS/tenMeterWindSpeed'][a[0],:]
    pType=(fh['NS/Input/precipitationType'][a[0],:]/1e7).astype(int)
    sfcType=fh['NS/Input/surfaceType'][a[0],:]
    nx,ny=sfcRain.shape
# Marti,26Feb!2019
# Joi,21Mart!2019
    umu=cos(53/180.*pi)
    tbL=[]
    freqs=[10.6,10.6,18.7,18.7,23.,37,37.,89,89.,166.,166.,186.3,190.3]
    npol=[1,0,1,0,1,1,0,1,0,1,0,1,1]
    
    map = Basemap(projection='merc', resolution = 'c', area_thresh=0.1)
    tbSim=[]
    for i in range(2,nx-3):
        for j in range(3,49-3):
            if pType[i,j]<1:
                if(sfcType[i-2:i+3,j-2:j+3].max()==0):
                    for (pol,f) in zip(npol,freqs):
                        kextL=[]
                        for q,tk,pa in zip(qv[i,j,:],airTemp[i,j,:],\
                                               press[i,j,:]):
                            absair,abswv=pyHB2.gasabsr98(f,tk,q*1e-3,pa*1e2,1)
                            kextL.append(absair+abswv)
                        bins=arange(envNode[i,j,0],envNode[i,j,-1]+1)
                        kextInt=interp(bins,envNode[i,j,:],kextL)
                        tLayer=list(interp(bins,envNode[i,j,:],airTemp[i,j,:]))
                        tLayer.append(sfcTemp[i,j])
                        salb=kextInt.copy()*0.
                        asym=kextInt.copy()*0.
                        emis,ebar=pyHB2.emit(f,pol,sfcTemp[i,j],w10[i,j],umu)
                        nL=kextInt.shape[0]
                        tb1=pyHB2.radtran(umu,nL,sfcTemp[i,j],tLayer[::-1],\
                                          arange(nL+1)*0.25,kextInt[::-1],\
                                          salb,asym,2.7,emis,ebar)
                        tbSim.append(tb1)
            
                    tbL.append(tbSim)
    return tbL

import xarray as xr
tbSim=cmbFile(fpath,fname)
simTb=xr.DataArray(array(tbL),dims=['np','nf'])
simTb_x=xr.Dataset({'simTb':simTb})
simTb_x.to_netcdf('simTb_180.nc')
