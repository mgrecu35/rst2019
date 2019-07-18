import matplotlib.pyplot as plt
from netCDF4 import Dataset
from numpy import *

f=Dataset("iphex.output10_freq1_0175_dn5.nc","r")
zMSKa=f['zMS'][:,:]
zSSKa=f['zSS'][:,:]
kext=f['kext'][:,:]
piaKa=f['piaKa'][:]
zMSKam=ma.array(zMSKa,mask=zMSKa<-10)
zSSKam=ma.array(zSSKa,mask=zSSKa<-10)

f=Dataset("iphex.output10_freq2_0175_dn5.nc","r")
zMSW=f['zMS'][:,:]
zSSW=f['zSS'][:,:]
kext=f['kext'][:,:]
piaW=f['piaKa'][:]
zMSWm=ma.array(zMSW,mask=zMSW<-10)
zSSWm=ma.array(zSSW,mask=zSSW<-10)

f=Dataset("iphex.output10_freq0_0175_dn5.nc","r")
zMSKu=f['zMS'][:,:]
zSSKu=f['zSS'][:,:]
kext=f['kext'][:,:]
piaKu=f['piaKa'][:]
zMSKum=ma.array(zMSKu,mask=zMSKu<-10)
zSSKum=ma.array(zSSKu,mask=zSSKu<-10)



import pyHB2 as pyHB2
pyHB2.initt()
hg=arange(80)*0.25+0.125
f=interp(hg,[0,4,5,20],[0.,0.,1,1.])

plt.plot(zMSKum[100,::-1],hg)
plt.plot(zMSKam[100,::-1],hg)
plt.plot(zMSWm[100,::-1],hg)
plt.legend(['Ku','Ka', 'W'])
plt.xlabel('dBZ')
plt.ylabel('Height (km)')
plt.savefig('Profile100_0175.png')
from gaussNewton import *
dn=-.5
piaKa=0
piaW=0.
zW_sim=zMSWm[100,:].copy()
for k in range(79,-1,-1):
    if zMSKam[100,79-k]>-10:
        wc,dpiaKa,zka_mod,dpiaW,zkaW=newton_raphsonS(f[k],dn,zMSKam[100,79-k]+piaKa,pyHB2)
        piaKa+=dpiaKa*2*0.25
        piaW+=dpiaW*2*0.25
stop
#f2=Dataset('cm1_SquallLine.Fields2.nc','r')
#qr=f2['qr'][:,9,:].T
#qs=f2['qs'][:,9,:].T
#qg=f2['qg'][:,9,:].T
#qc=f2['qc'][:,9,:].T1
#zTrue2=f2['z3d'][:,9,:,1]

nx=130
plt.figure(figsize=(8,11))
plt.subplot(311)
plt.pcolormesh(arange(nx)*1.3,arange(80)*0.25,zMSKum[::1,::-1].T,cmap='jet',vmax=50)
plt.ylabel('Height (km)')
plt.ylim(0,15)
plt.colorbar()
plt.subplot(312)
plt.suptitle('Ku/Ka/W-band reflectivity 0.175 beamwidth')
plt.pcolormesh(arange(nx)*1.3,arange(80)*0.25,zMSKam[::1,::-1].T,cmap='jet',vmax=50)
plt.ylabel('Height (km)')
plt.ylim(0,15)
plt.colorbar()
plt.subplot(313)
plt.pcolormesh(arange(nx)*1.3,arange(80)*0.25,zMSWm[:,::-1].T,cmap='jet',vmax=50)
plt.ylabel('Height (km)')
plt.xlabel('km')
plt.ylim(0,15)
plt.colorbar()
plt.savefig('KuKaWZ_0175.png')
plt.figure(figsize=(8,11))
plt.subplot(311)
plt.pcolormesh(arange(nx)*1.3,arange(80)*0.25,zSSKum[::1,:].T,cmap='jet',vmax=50)
plt.ylabel('Height (km)')
plt.ylim(0,15)
plt.colorbar()
plt.subplot(312)
plt.suptitle('Ku/Ka/W-band reflectivity no MSS')
plt.pcolormesh(arange(nx)*1.3,arange(80)*0.25,zSSKam[::1,:].T,cmap='jet',vmax=50)
plt.ylabel('Height (km)')
plt.ylim(0,15)
plt.colorbar()
plt.subplot(313)
plt.pcolormesh(arange(nx)*1.3,arange(80)*0.25,zSSWm[:,:].T,cmap='jet',vmax=50)
plt.ylabel('Height (km)')
plt.xlabel('km')
plt.ylim(0,15)
plt.colorbar()
plt.savefig('KuKaWZ.png')
plt.figure()
qr=f['qr'][:,:]
#qc=f['qr'][:,:]
plt.subplot(211)
plt.semilogy(arange(nx)*1.3,qr[:,0])
plt.xlim(0,nx*1.3)
plt.subplot(212)
plt.xlim(0,nx*1.3)
f2=Dataset('wrfJune11_2014.Fields.nc')
qr2=f2['qr'][:,9,:].T
qc2=f2['qc'][:,9,:].T
z1=f2['z'][:,9,:].T
zm=0.5*(z1[:,1:]+z1[:,:-1])
nz=80
hg=(arange(nz)+0.5)*0.25
qrL=[]
qcL=[]
for i in range(nx):
    qrL.append(interp(hg,zm[i,:],qr2[i,:]))
    qcL.append(interp(hg,zm[i,:],qc2[i,:]))
print(corrcoef(qr[:,0],array(qrL)[:,0]))

plt.figure()
plt.contour(arange(nx)*1.3,hg,array(qcL).T,levels=(arange(10)+0.25)*0.1)
plt.contour(arange(nx)*1.3,hg,array(qrL).T,levels=(arange(10)+0.25)*0.1)
stop
stop
vMS=f['vMS'][:,:]
vSS=f['vSS'][:,:]
zTrue=f['zTrue'][:,:]
vMSm=ma.array(vMS,mask=zMS<0)
vSSm=ma.array(vSS,mask=zSS<0)

plt.figure()
plt.subplot(211)
plt.suptitle('Doppler velocity')

plt.pcolormesh((arange(nx))*.5,arange(80)*.25,-vMSm[::,::-1].T,cmap='RdBu_r',vmin=-15,vmax=15)

plt.ylabel('Height (km)')
plt.colorbar()
plt.subplot(212)
plt.pcolormesh((arange(nx))*0.5,arange(80)*0.25,-vSSm[:,::].T,cmap='RdBu_r',vmin=-15,vmax=15)
plt.xlabel('km')
plt.ylabel('Height (km)')
plt.colorbar()
plt.savefig('KaDoppler.png')

piaKa_u=4.343*2*0.25*kext[:,:].sum(axis=1)

ix=150
plt.figure()
plt.plot(zSS[ix,:],arange(80)[:]*0.25)
plt.plot(zMS[ix,:],arange(80)[::-1]*0.25)
plt.plot(zTrue[ix,:],arange(80)[:]*0.25)
plt.plot(zTrue2[ix,:],arange(80)[:]*0.25)
plt.ylim(0,12)
plt.xlim(-10,50)
plt.figure()
plt.plot(qs[ix,:],arange(80)[:]*0.25,'*')
plt.plot(qg[ix,:],arange(80)[:]*0.25)
plt.plot(qr[ix,:],arange(80)[:]*0.25)
plt.plot(qc[ix,:],arange(80)[:]*0.25,'*')
plt.ylim(0,12)
