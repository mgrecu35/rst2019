from netCDF4 import Dataset
from numpy import *


f=Dataset('cm1out_2steps.nc','r')
qv=f['qv'][-1,:,:,:]    # water vapor
qr=f['qr'][-1,:,:,:]    # rain mixing ratio
qs=f['qs'][-1,:,:,:]*.95    # snow mixing ratio
qc=f['qc'][-1,:,:,:]    # cloud mixing ratio
ncs=f['ncs'][-1,:,:,:]    # snow number concentration
ncg=f['ncg'][-1,:,:,:]    # graupel number concentration
ncr=f['ncr'][-1,:,:,:]    # rain number concentration
qg=f['qg'][-1,:,:,:]*.95    # graupel mixing ratio
z=f['z'][:]             # height (km)
th=f['th'][-1,:,:,:]    # potential temperature (K)
prs=f['prs'][-1,:,:,:]  # pressure (Pa)
T=th*(prs/100000)**0.286  # Temperature
dbz=f['dbz'][-1,:,:,:]
R=287.058  #J*kg-1*K-1
rho=prs/(R*T)

w=f['w'][-1,:,:,:]
wm=0.5*(w[1:,:,:]+w[:-1,:,:])
nz,ny,nx=qv.shape
n0w=0.
nfreq=8
zKu_L=[]
zKa_L=[]
zW_L=[]
vKu_L=[]
vKa_L=[]
vW_L=[]
j=int(ny/2)+5
import pyHB2 as pyHB2
pyHB2.initt()
x1L=[]
from getNw import *
dr=0.25

fd=Dataset('dm_z.nc')
dm_z=fd['dm_z'][:,:]

fd=Dataset('doppler_2.nc')
vt_R=fd['rain'][:,:]
vt_S=fd['snow'][:,:]
ls=open('BlueWhiteOrangeRed.rgb.txt','r').readlines()

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
n2=350
for i in range(150,n2):
    zKu_1d=zeros((nz),float)-99
    zKa_1d=zeros((nz),float)-99
    zW_1d=zeros((nz),float)-99
    vW_1d=zeros((nz),float)-99
    vKa_1d=zeros((nz),float)-99
    vKu_1d=zeros((nz),float)-99
    piaKu=0.
    piaKa=0.
    piaW=0.
    for k in range(nz-1,-1,-1):
        iwc=(qs[k,j,i]+qg[k,j,i])*rho[k,j,i]*1e3
        pwc=qr[k,j,i]*rho[k,j,i]*1e3
        if iwc>1e-4:
            if qs[k,j,i]>0.001e-4:
                lams,nws=getn0s(qs[k,j,i],ncs[k,j,i])
                nws=nws*rho[k,j,i]*1e-8
            else:
                nws=0.08
            if qg[k,j,i]>0.00001e-4:
                lamg,nwg=getn0g(qg[k,j,i],ncg[k,j,i])
                nwg=nwg*rho[k,j,i]*1e-8
            else:
                nwg=0.08
            iwc1=(qs[k,j,i])*rho[k,j,i]*1e3
            iwc2=(qg[k,j,i])*rho[k,j,i]*1e3
            kextS1,salbS1,asymS1,dpiaKu_S1,dpiaKa_S2,zKu_S1,zKa_S1,zW_S1,dmS1 = pyHB2.getsnowp2(log10(nws/0.08),(iwc1),nfreq)
            kextS2,salbS2,asymS2,dpiaKu_S2,dpiaKa_S2,zKu_S2,zKa_S2,zW_S2,dmS2 = pyHB2.getsnowp2(log10(nwg/0.08),(iwc2),nfreq)
            zW_S11=getZ_W(dmS1,zKu_S1,nws,dm_z,zW_S1)
            #print(zW_S11,zW_S1)
            zW_S2=getZ_W(dmS2,zKu_S2,nwg,dm_z,zW_S2)
            zKu_S=log10(10.**(0.1*zKu_S1)+10.**(0.1*zKu_S2))*10
            zKa_S=log10(10.**(0.1*zKa_S1)+10.**(0.1*zKa_S2))*10
            zWS=log10(10.**(0.1*zW_S1)+10.**(0.1*zW_S2))*10
            i0s1=int((zKu_S1-10*log10(nws/0.08)+12.)/0.5)
            i0s1=max(i0s1,0)
            i0s1=min(114,i0s1)
            i0s2=int((zKu_S2-10*log10(nwg/0.08)+12.)/0.5)
            i0s2=max(i0s2,0)
            i0s2=min(114,i0s2)
            vs1=vt_S[i0s1,1:]
            vs2=vt_S[i0s2,1:]
            vS=(10**(0.1*zKu_S1)*vs1+10**(0.1*zKu_S2)*vs2)/\
                (10**(0.1*zKu_S1)+10**(0.1*zKu_S2))
            #vS=vS*0.
            dpiaW_S1=kextS1[4]*4.343
            dpiaW_S2=kextS2[4]*4.343
        else:
            zKu_S=-99
            zKa_S=-99
            zWS=-99
            vS=array([0.,0.,0.])
            dpiaKu_S1=0
            dpiaKu_S2=0
            dpiaKa_S1=0
            dpiaKa_S2=0
            dpiaW_S1=0
            dpiaW_S2=0
        n0w=0.
        if pwc>1e-9:
            lamr,nwr=getn0w(qr[k,j,i],ncr[k,j,i])
            nwr=nwr*rho[k,j,i]*1e-8
            kextR,salbR,asymR,dpiaKu_R,dpiaKa_R,zKu_R,zKa_R,zWR, dmR = pyHB2.getrainp2(log10(nwr/0.08),(pwc),nfreq)
            i0dm=int((dmR-0.2)/0.02)
            if i0dm>0 and i0dm<190:
                dZ=zKu_R-10*log10(nwr/0.08)-dm_z[i0dm,1]
                zWR=dm_z[i0dm,3]+dZ+10*log10(nwr/0.08)
            i0r=int((zKu_R-10*log10(nwr/0.08)+12)/0.5)
            i0r=max(i0r,0)
            i0r=min(124,i0r)
            vR=vt_R[i0r,1:]
            if zKu_R>300:
                stop
            dpiaW_R=kextR[4]*4.343
        else:
            zKu_R=-99
            zKa_R=-99
            zWR=-99
            dpiaKu_R=0
            dpiaKa_R=0.
            dpiaW_R=0.
            dmR=0
            vR=array([0,0.,0])
        a=2.65
        b=0.098
        dpiaL=[]
        ireturn=0
        for freq in freqs:
            cldw=qc[k,j,i]*1e3*rho[k,j,i]
            z_clw = pyHB2.gcloud(freq,T[k,j,i],cldw)
            absair,abswv = pyHB2.gasabsr98(freq,T[k,j,i],qv[k,j,i],prs[k,j,i],ireturn)
            dpiaL.append((z_clw+absair+abswv)*4.343)
        zKu_1d[k]=log10(10.**(0.1*zKu_S)+10.**(0.1*zKu_R))*10-piaKu-(dpiaKu_S1+dpiaKu_S2+dpiaKu_R+dpiaL[0])*dr
        zKa_1d[k]=log10(10.**(0.1*zKa_S)+10.**(0.1*zKa_R))*10-piaKa-(dpiaKa_S1+dpiaKa_S2+dpiaKa_R+dpiaL[1])*dr
        zW_1d[k]=log10(10.**(0.1*zWS)+10.**(0.1*zWR))*10.-piaW-(dpiaW_S1+dpiaW_S2+dpiaW_R+dpiaL[2])*dr
        piaKu+=2*(dpiaKu_S1+dpiaKu_S2+dpiaKu_R+dpiaL[0])*dr  # attenuation due to cloud and  water vapor
        piaKa+=2*(dpiaKa_S1+dpiaKa_S2+dpiaKa_R+dpiaL[1])*dr  # not included yet
        piaW+=2*(dpiaW_S1+dpiaW_S2+dpiaW_R+dpiaL[2])*dr
        vT=(10**(0.1*zKu_S)*vS+10**(0.1*zKu_R)*vR)/\
            (10**(0.1*zKu_S)+10**(0.1*zKu_R))
        vW_1d[k]=vT[2]-wm[k,j,i]
        vKu_1d[k]=vT[0]-wm[k,j,i]
        vKa_1d[k]=vT[1]-wm[k,j,i]
    zKu_L.append(zKu_1d)
    zKa_L.append(zKa_1d)
    zW_L.append(zW_1d)
    vW_L.append(vW_1d)
    vKu_L.append(vKu_1d)
    vKa_L.append(vKa_1d)
    print(piaKu,piaKa,piaW)
    #x1L.append(x
zKu_L=array(zKu_L)
zKa_L=array(zKa_L)
zW_L=array(zW_L)
vKu_L=array(vKu_L)
vKa_L=array(vKa_L)
vW_L=array(vW_L)
xh=f['xh'][:]
z=f['z'][:]
zKum=ma.array(zKu_L,mask=zKu_L<-20)
zKam=ma.array(zKa_L,mask=zKa_L<-20)
zWm=ma.array(zW_L,mask=zW_L<-20)

vKum=ma.array(vKu_L,mask=zKu_L<-20)
vKam=ma.array(vKa_L,mask=zKa_L<-20)
vWm=ma.array(vW_L,mask=zW_L<-20)

plt.figure(figsize=(8,11))
fig=plt.subplot(311)
plt.pcolormesh(xh[150:n2],z,zKum.T,cmap='jet',vmin=-20)
fig.axes.get_xaxis().set_visible(False)
plt.ylim(0,15)
plt.title('Simulated Ku (dBz)')
plt.ylabel('Height (km)')
plt.colorbar()
fig2=plt.subplot(312)
plt.pcolormesh(xh[150:n2],z,zKam.T,cmap='jet',vmin=-20)
plt.ylim(0,15)
plt.title('Simulated Ka (dBz)')
plt.ylabel('Height (km)')
fig2.axes.get_xaxis().set_visible(False)
plt.colorbar()
plt.subplot(313)
plt.pcolormesh(xh[150:n2],z,zWm.T,cmap='jet',vmin=-20)
plt.ylim(0,15)
plt.title('Simulated W (dBz)')
plt.xlabel('Distance (km)')
plt.ylabel('Height (km)')
plt.colorbar()
plt.savefig('squallLineZ.png')


plt.figure(figsize=(8,11))
fig=plt.subplot(311)
plt.pcolormesh(xh[150:n2],z,vKum.T,cmap=cmp,vmin=-12,vmax=12)
fig.axes.get_xaxis().set_visible(False)
plt.ylim(0,15)
plt.title('Simulated Ku-band Doppler (m/s)')
plt.ylabel('Height (km)')
plt.colorbar()
fig2=plt.subplot(312)
plt.pcolormesh(xh[150:n2],z,vKam.T,cmap=cmp,vmin=-12,vmax=12)
plt.ylim(0,15)
plt.title('Simulated Ka-band Doppler (m/s)')
plt.ylabel('Height (km)')
fig2.axes.get_xaxis().set_visible(False)
plt.colorbar()
plt.subplot(313)
plt.pcolormesh(xh[150:n2],z,vWm.T,cmap=cmp,vmin=-12,vmax=12)
plt.ylim(0,15)
plt.title('Simulated W-band Doppler (m/s)')
plt.xlabel('Distance (km)')
plt.ylabel('Height (km)')
plt.colorbar()
plt.savefig('squallLineDoppler.png')

#plt.figure()
#plt.pcolormesh(dbzm[:,j,150:n2],cmap=cmp,vmin=0)
#plt.colorbar()
