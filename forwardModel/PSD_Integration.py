from numpy import *
from bhmie import *
import pytmatrix.refractive
wl=[pytmatrix.refractive.wl_Ku,pytmatrix.refractive.wl_Ka,pytmatrix.refractive.wl_W]

#print(wl) # units are mm

print(dir(pytmatrix.refractive))

mw1= pytmatrix.refractive.m_w_0C[wl[0]]
mw2= pytmatrix.refractive.m_w_0C[wl[1]]
mw3= pytmatrix.refractive.m_w_0C[wl[2]]
ms1= pytmatrix.refractive.mi(wl[0],0.1)
ms1= pytmatrix.refractive.mi(wl[1],0.1)
ms2= pytmatrix.refractive.mi(wl[2],0.1)
sback13w=[]
sback13s=[]
sback37w=[]
sback95w=[]
nang=20

dbinsM=(0.025+arange(300)*0.05)
massW=(0.1*dbinsM)**3/6.*pi
dMax=10.*(massW/0.0061)**(1/2.05)

vM=(0.1*dbinsM)**3/6.*pi
vIntW=17.67*(0.1*dbinsM)**0.67  #Atlas and Ulbrich
D=0.1*dbinsM
vIntW=9.25-9.25*exp(-(6.8*D**2+4.88*D))
vInt=2.34*(0.1*dbinsM)**0.31
vInt=vIntW/4.
sfact=(0.92/0.1)**(0.333)
xW1=4*3.142/wl[0]*(0.25*dMax)
xf1=(1+0.159*xW1**2)/(1+(0.159+1./3)*xW1**2+0.164*xW1**4);
xW2=4*3.142/wl[1]*(0.25*dMax)
xf2=(1+0.159*xW2**2)/(1+(0.159+1./3)*xW2**2+0.164*xW2**4);
xW3=4*3.142/wl[2]*(0.25*dMax)
xf3=(1+0.159*xW3**2)/(1+(0.159+1./3)*xW3**2+0.164*xW3**4);

sback13s_rg=[]
sback37s_rg=[]
sback95s_rg=[]
mi1= pytmatrix.refractive.mi(wl[0],0.92)
mi2= pytmatrix.refractive.mi(wl[1],0.92)
mi3= pytmatrix.refractive.mi(wl[2],0.92)
m=array([mi1,mi2,mi3])
K2r=abs((m**2-1)/(m**2+2))**2
K2s=((m**2-1)/(m**2+2))**2
K2=0.93

for i,d in enumerate(dbinsM):
    x=d*pi/wl[0]
    s1,s2,qext13w,qsca13w,qback13w,gsca13w=bhmie(x,mw1,nang)
    sback13w.append(qback13w*d**2*pi/4.*wl[0]**4/pi**5) ## mm^6
    xs=d*pi/wl[0]*sfact
    s1s,s2s,qext13s,qsca13s,qback13s,gsca13s=bhmie(xs,ms1,nang)
    sback13s.append(qback13s*(d*sfact)**2*pi/4.*wl[0]**4/pi**5) ##m^6
    sback13s_rg.append(d**6*xf1[i]*K2r[0])
    sback37s_rg.append(d**6*xf2[i]*K2r[0])
    sback95s_rg.append(d**6*xf3[i]*K2r[0])
    x=d*pi/wl[1]
    s1,s2,qext37w,qsca37w,qback37w,gsca37w=bhmie(x,mw2,nang)
    sback37w.append(qback37w*(d)**2*pi/4.*wl[1]**4/pi**5)
    x=d*pi/wl[2]
    s1,s2,qext95w,qsca95w,qback95w,gsca95w=bhmie(x,mw3,nang)
    sback95w.append(qback95w*(d)**2*pi/4.*wl[2]**4/pi**5)

sback13w=array(sback13w)
sback37w=array(sback37w)
sback95w=array(sback95w)
K2=0.93
Dm=10.
dat1=[]
R=logspace(-3,3.5,500)


for r1 in R:
    r=r1+0.
    lambd=1.125*41.*r**(-0.21)
    nd=0.08*exp(-lambd*0.1*dbinsM)*0.005*(d/Dm)**0.5
    M=sum(nd*vM)*1e6  # sum(nd*vM) in g/cm^3 to g/m^3
    Z1w=log10(sum(nd*1e6*sback13w/K2))*10.
    Z1s=log10(sum(nd*1e6*sback13s/K2))*10.
    Z1s_rg=log10(sum(nd*1e6*sback13s_rg/K2))*10.
    Z2s_rg=log10(sum(nd*1e6*sback37s_rg/K2))*10.
    Z3s_rg=log10(sum(nd*1e6*sback95s_rg/K2))*10.
    Z2w=log10(sum(nd*1e6*sback37w/K2))*10.
    Z3w=log10(sum(nd*1e6*sback95w/K2))*10.
    dm=sum(nd*dbinsM**4)/sum(nd*dbinsM**3)
    v1W=sum(nd*1e6*sback13w*vIntW/K2)/sum(nd*1e6*sback13w/K2)
    v2W=sum(nd*1e6*sback37w*vIntW/K2)/sum(nd*1e6*sback37w/K2)
    v3W=sum(nd*1e6*sback95w*vIntW/K2)/sum(nd*1e6*sback95w/K2)
    v1S=sum(nd*1e6*sback13s_rg*vInt)/sum(nd*1e6*sback13s_rg)
    v2S=sum(nd*1e6*sback37s_rg*vInt)/sum(nd*1e6*sback37s_rg)
    v3S=sum(nd*1e6*sback95s_rg*vInt)/sum(nd*1e6*sback95s_rg)
    
    dat1.append([M,dm,Z1w,Z2w,Z3w,v1W,v2W,v3W,Z1s,v1S,v2S,v3S])
    
    print(Z1w,Z2w,Z3w,Z1s,Z1s_rg,Z2s_rg,Z3s_rg)

dat1=array(dat1)
d1=0.2+arange(191)*0.02
z1=interp(d1,dat1[:,1],dat1[:,2])
z2=interp(d1,dat1[:,1],dat1[:,3])
z3=interp(d1,dat1[:,1],dat1[:,4])
import xarray as xr
import numpy as np
dm_z=xr.DataArray(np.c_[d1,z1,z2,z3])
dm_zs=xr.Dataset({'dm_z':dm_z})
dm_zs.to_netcdf('dm_z.nc')
import matplotlib.pyplot as plt



zint=-12+arange(145)*0.5
v1=interp(zint,dat1[:,2],dat1[:,5])
v2=interp(zint,dat1[:,2],dat1[:,6])
v3=interp(zint,dat1[:,2],dat1[:,7])

zints=-12+arange(115)*0.5
v1s=interp(zints,dat1[:,8],dat1[:,9])
v2s=interp(zints,dat1[:,8],dat1[:,10])
v3s=interp(zints,dat1[:,8],dat1[:,11])

import xarray as xr

d1=xr.DataArray(np.c_[zint,v1,v2,v3],dims=['nbins','n4'])
d2=xr.DataArray(np.c_[zints,v1s,v2s,v3s],dims=['nbinsS','n4'])
doppler=xr.Dataset({'rain':d1,'snow':d2})
doppler.to_netcdf('doppler.nc')

