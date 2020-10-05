from math import *
import random
#from numpy import *


class photon:
    def __init__(self,xini,yini,zini):
        self.ilayer=-1
        self.muz=0.
        self.mux=0.
        self.muy=0.
        self.x=xini
        self.y=yini
        self.z=zini

class cloud:
    def __init__(self,n,dz):
        self.kext=list()
        self.salb=list()
        self.ebr=list()
        self.g=list()
        self.zTop=list()
        self.zBot=list()
        self.tauD=list()
        self.tauU=list()
        self.z=list()
        self.zobs=list()
        self.vterm=list()
        self.ztrue=list()
        self.n=n
        
from netCDF4 import Dataset
from numpy import *
def cloudprof(f,ix,iy,ifreq):
    c=cloud(80,dz)
    i=0
    nz=80
    #print(ix,iy)
    z1d=f['z3d'][ix,iy,:,ifreq]
    kext1d=f['kext3d'][ix,iy,:,ifreq]
    salb1d=f['salb3d'][ix,iy,:,ifreq]
    asym1d=f['asym3d'][ix,iy,:,ifreq]
    v1d=f['v3d'][ix,iy,:,ifreq]
    z=0.5*(f['z'][1:,iy,ix]+f['z'][:-1,iy,ix])
    hg=(arange(nz)+0.5)*0.25
    z1d_g=interp(hg,z,z1d)
    kext1d_g=interp(hg,z,kext1d)
    salb1d_g=interp(hg,z,salb1d)
    asym1d_g=interp(hg,z,asym1d)
    v1d_g=interp(hg,z,v1d)
    for i in range(nz):
        c.z.append(0.)
        c.zobs.append(-99.9)
        c.ztrue.append(z1d_g[i])
        c.kext.append(kext1d_g[i])
        c.salb.append(salb1d_g[i])
        c.g.append(asym1d_g[i])
        c.vterm.append(v1d_g[i])
        if(z1d_g[i]>-40):
            c.ebr.append(kext1d_g[i]/10.**(0.1*z1d_g[i]))
        else:
            c.ebr.append(0.)
        c.zTop.append((i+1)*dz)
        c.zBot.append((i)*dz)
        c.tauD.append(0.)
        c.tauU.append(0.)
        if(c.ebr[i]>0):
            c.z[i]=log10(c.kext[i]/c.ebr[i])*10.
        else:
            c.z[i]=-99.
        i=i+1

    c.tauD.append(0.)
    c.tauU.append(0.)
        
    dztop=380
    c.zTop[c.n-1]=c.zTop[c.n-1]+dztop
    
    for i in range(c.n-1,-1,-1):
        c.tauD[i]=c.tauD[i+1]+c.kext[i]*dz

    for i in range(c.n):
        c.tauU[i+1]=c.tauU[i]+c.kext[i]*dz
        
    for i in range(c.n-1,-1,-1):
        if(c.z[i]>-10):
            c.z[i]=c.z[i]+0*log10(exp(-(c.tauD[i+1]+c.tauD[i])))*10.
            c.zobs[i]=c.z[i]+log10(exp(-(c.tauD[i+1]+c.tauD[i])))*10.
                
    return c,log10(exp(-(c.tauD[0]+c.tauD[1])))*10.



def antenaF(ang,sigma):
    f=exp(-4*log(2)*(ang/sigma)*(ang/sigma))#/\
#        (2*pi)**.5/(sigma/2./(2*log(2))**.5)
    return f


dz=0.25





ifreq=1
#fname='cm1_SquallLine.Fields2.nc'


from numpy import *
c3d=[]

def readdata(fname,ifreq,nx,ny):
    ic=0
    f=Dataset(fname,'r')
    print(f)
    global c3d
    c3d=[]
    for ix in range(nx):
        c3d.append([])
        for jx in range(ny):
            c,pia=cloudprof(f,ix,jx,ifreq)
            c.pia=-pia
            c3d[ix].append(c)
            ic=ic+1
        #print(ix+1,jx+1,-pia)
    f.close()

def getkext(nx,ny):
    kext=zeros((nx,ny,80),float)
    for i in range(nx):
        for j in range(ny):
            kext[i,j,:]=array(c3d[i][j].kext)

    return kext

def getvterm(nx,ny):
    vterm=zeros((nx,ny,80),float)
    for i in range(nx):
        for j in range(ny):
            vterm[i,j,:]=array(c3d[i][j].vterm)

    return vterm

def getsalb(nx,ny):
    kext=zeros((nx,ny,80),float)
    for i in range(nx):
        for j in range(ny):
            kext[i,j,:]=array(c3d[i][j].salb)

    return kext

def getg(nx,ny):
    kext=zeros((nx,ny,80),float)
    for i in range(nx):
        for j in range(ny):
            kext[i,j,:]=array(c3d[i][j].g)

    return kext

def getz(nx,ny):
    kext=zeros((nx,ny,80),float)
    for i in range(nx):
        for j in range(ny):
            kext[i,j,:]=array(c3d[i][j].z)

    return kext

def getzobs(nx,ny):
    kext=zeros((nx,ny,80),float)
    for i in range(nx):
        for j in range(ny):
            kext[i,j,:]=array(c3d[i][j].zobs)

    return kext

def getebr(nx,ny):
    kext=zeros((nx,ny,80),float)
    for i in range(nx):
        for j in range(ny):
            kext[i,j,:]=array(c3d[i][j].ebr)

    return kext

def getztop(nx,ny):
    kext=zeros((nx,ny,80),float)
    for i in range(nx):
        for j in range(ny):
            kext[i,j,:]=array(c3d[i][j].zTop)

    return kext

def getzbot(nx,ny):
    kext=zeros((nx,ny,80),float)
    for i in range(nx):
        for j in range(ny):
            kext[i,j,:]=array(c3d[i][j].zBot)

    return kext

def gettauU(nx,ny):
    kext=zeros((nx,ny,81),float)
    for i in range(nx):
        for j in range(ny):
            kext[i,j,:]=array(c3d[i][j].tauU)

    return kext

def gettauD(nx,ny):
    kext=zeros((nx,ny,81),float)
    for i in range(nx):
        for j in range(ny):
            kext[i,j,:]=array(c3d[i][j].tauD)

    return kext

import xarray as xr
def saveMS(zMS,vMS,zObs,vTerm,pia,zTrue,iy,fname,fnameS):
    zMSx=xr.DataArray(zMS)
    vMSx=xr.DataArray(vMS)
    zObsx=xr.DataArray(zObs)
    zTruex=xr.DataArray(zTrue)
    vTermx=xr.DataArray(vTerm)
    f=Dataset(fnameS,'r')
    qgx=xr.DataArray(f['qg'][:,iy-1,:].T)
    qrx=xr.DataArray(f['qr'][:,iy-1,:].T)
    qsx=xr.DataArray(f['qs'][:,iy-1,:].T)
    qcx=xr.DataArray(f['qc'][:,iy-1,:].T)
    dngx=xr.DataArray(f['dng'][:,iy-1,:])
    dnrx=xr.DataArray(f['dnr'][:,iy-1,:])
    dnsx=xr.DataArray(f['dns'][:,iy-1,:])
    
    kextx=xr.DataArray(f['kext3d'][:,iy,:,1])
    print(zMSx.shape)
    print(vMSx.shape)
    print(zObsx.shape)
    print(zTruex.shape)
    print(vTermx.shape)
    print(qrx.shape)
    print(qsx.shape)
    nx=zMSx.shape[0]
    nz=80
    hg=(arange(nz)+0.5)*0.25
    qr1=zeros((nx,80),float)
    qg1=zeros((nx,80),float)
    qs1=zeros((nx,80),float)
    qc1=zeros((nx,80),float)
    dnr1=zeros((nx,80),float)
    dng1=zeros((nx,80),float)
    dns1=zeros((nx,80),float)
    kext1=zeros((nx,80),float)
    print(dnrx.shape)
    for ix in range(nx):
        z=0.5*(f['z'][1:,iy-1,ix]+f['z'][:-1,iy-1,ix])
        qr1[ix,:]=interp(hg,z,qrx[ix,:])
        qs1[ix,:]=interp(hg,z,qsx[ix,:])
        qg1[ix,:]=interp(hg,z,qgx[ix,:])
        qc1[ix,:]=interp(hg,z,qcx[ix,:])
        dnr1[ix,:]=interp(hg,z,dnrx[ix,:])
        dns1[ix,:]=interp(hg,z,dnsx[ix,:])
        dng1[ix,:]=interp(hg,z,dngx[ix,:])
        kext1[ix,:]=interp(hg,z,kextx[ix,:])
    qgx=xr.DataArray(qg1)
    qrx=xr.DataArray(qr1)
    qsx=xr.DataArray(qs1)
    qcx=xr.DataArray(qc1)
    dngx=xr.DataArray(dng1)
    dnrx=xr.DataArray(dnr1)
    dnsx=xr.DataArray(dns1)
    kextx=xr.DataArray(kext1)
    dnsx=xr.DataArray(dns1)
    d=xr.Dataset({'zMS':zMSx,'vMS':vMSx,'zSS':zObsx,'zTrue':zTruex,'vSS':vTermx, \
                  'qr':qrx,'qs':qsx,'qg':qgx, 'qc':qcx,\
                  'dnr':dnrx,'dns':dnsx,'dng':dngx,'kext':kextx,'piaKa':pia})
    d.to_netcdf(fname)
