from numpy import *
#import pyHB2 as pyHB2
#pyHB2.initt()

def zka(wc,f,dn,pyHB2):
    nfreq=8
    swc=f*wc
    lwc=(1-f)*wc
    if swc>0.001:
        kextS,salbS,asymS,dpiaKu_S,dpiaKa_S,\
            zKu_S,zKa_S,zW_S,dmS = pyHB2.getsnowp2(dn,(swc),nfreq)
        #print(swc,dn,zKa_S)
        dpiaW_S=kextS[4]*4.343
    else:
        dpiaKa_S=0
        zKa_S=-99
        zW_S=-99
        dpiaW_S=0
    if lwc>0.0001:
        kextR,salbR,asymR,dpiaKu_R,dpiaKa_R,\
            zKu_R,zKa_R,zW_R,dmR = pyHB2.getrainp2(dn,(lwc),nfreq)
        dpiaW_R=kextR[4]*4.343
    else:
        dpiaKa_R=0
        zKa_R=-99
        zW_R=-99
        dpiaW_R=0
        
    zKa=log10(10.**(0.1*zKa_S)+10.**(0.1*zKa_R))*10.
    zW=log10(10.**(0.1*zW_S)+10.**(0.1*zW_R))*10.
    dpiaKa=dpiaKa_S+dpiaKa_R
    dpiaW=dpiaW_S+dpiaW_R
    return zKa, dpiaKa, zW, dpiaW    
def newton_raphsonS(f,dn,zka_obs,pyHB2):
    wc=0.1
    for it in range(25):
        zka_mod,dpiaKa,zW,dpiaW=zka(wc,f,dn,pyHB2)
        zka_mod1,dpiaKa1,zW1,dpiaW1=zka(wc+0.05,f,dn,pyHB2)
        if(it==14):
            print(it, zka_obs, zka_mod)
        dzdwc=(zka_mod1-zka_mod)/0.05
        wc+=(zka_obs-zka_mod)*(dzdwc)/(dzdwc**2+1e-3)
        if wc<0.005:
            wc=0.005
    return wc, dpiaKa, zka_mod, dpiaW, zW

#f=1.0
#dn=0
#wc,dpiaKa,zka_mod,dpiaW,zkaW=newton_raphsonS(f,dn,20.0,pyHB2)
