from math import *
import random
#from numpy import


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
        self.n=n

def cloudprof(fname,il):
    c=cloud(80,dz)
    f=open(fname,'r')
    lines=f.readlines()
    i=0
    for line in lines[il*80:(il+1)*80]:
        s=line.split()
        s=s[1:]
        c.z.append(0.)
        c.kext.append(float(s[0]))
        c.salb.append(float(s[1]))
        c.g.append(float(s[2]))
        c.ebr.append(float(s[3]))
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
        if(c.z[i]>0):
            c.z[i]=c.z[i]+log10(exp(-(c.tauD[i+1]+c.tauD[i])))*10.
                
    return c,log10(exp(-(c.tauD[0]+c.tauD[1])))*10.


def antenaF(ang,sigma):
    f=exp(-4*log(2)*(ang/sigma)*(ang/sigma))#/\
#        (2*pi)**.5/(sigma/2./(2*log(2))**.5)
    return f


dz=0.25

fname='/home/grecu/FrontPacked/profScatP2.dat'
fname2='/home/grecu/FrontPacked/profScatP3.dat'

c,pia=cloudprof(fname,0)
c2,pia2=cloudprof(fname2,0)
fname2='/home/grecu/FrontPacked/scattProfs.dat'
cs=list()
for i in range(400):
    cs.append(cloudprof(fname2,i))

from procWTabs import *
f=open('res2.txt','w')
for i in range(57,288):
    if(abs(-cs[i][1]-20)<10):
        dZSS11,dZSS21,icr1=mainProc(c,cs[i][0])
        dZSS12,dZSS22,icr2=mainProc(cs[i][0],c)
        s1=''
        s2=''
        for k in range(4):
            s1=s1+'%6.2f '%dZSS21[k]
            s2=s2+'%6.2f '%dZSS12[k]
        s=s1+s2+'%6.2f '%-cs[i][1]
        f.write(s+'\n')
        f.flush()
        print s
        stop
