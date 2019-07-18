def radarFields_3d():
    for i in range(nx1,nx2):
        for j in range(ny1,ny2):
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
                iwc=(qs[k,j,i]+qg[k,j,i]+qh[k,j,i]+qag[k,j,i])*rho[k,j,i]*1e3
                pwc=qr[k,j,i]*rho[k,j,i]*1e3
                if k==nz-1:
                    dr=(z[nz-1]-z[nz-2])/1.
                else:
                    dr=(z[k+1]-z[k])/1.
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
                    if qh[k,j,i]>0.00001e-4:
                        lamh,nwh=getn0g(qh[k,j,i],nch[k,j,i])
                        nwh=nwh*rho[k,j,i]*1e-8
                    else:
                        nwh=0.08
                    if qag[k,j,i]>0.00001e-4:
                        lamag,nwag=getn0s(qag[k,j,i],ncag[k,j,i])
                        nwag=nwag*rho[k,j,i]*1e-8
                    else:
                        nwag=0.08
                    iwc1=(qs[k,j,i])*rho[k,j,i]*1e3
                    iwc2=(qg[k,j,i])*rho[k,j,i]*1e3
                    iwc3=(qh[k,j,i])*rho[k,j,i]*1e3
                    iwc4=(qag[k,j,i])*rho[k,j,i]*1e3
                    kextS1,salbS1,asymS1,dpiaKu_S1,dpiaKa_S2,zKu_S1,zKa_S1,zW_S1,dmS1 = pyHB2.getsnowp2(log10(nws/0.08),(iwc1),nfreq)
                    kextS2,salbS2,asymS2,dpiaKu_S2,dpiaKa_S2,zKu_S2,zKa_S2,zW_S2,dmS2 = pyHB2.getsnowp2(log10(nwg/0.08),(iwc2),nfreq)
                    kextS3,salbS3,asymS3,dpiaKu_S3,dpiaKa_S3,zKu_S3,zKa_S3,zW_S3,dmS3 = pyHB2.getsnowp2(log10(nwh/0.08),(iwc3),nfreq)
                    kextS4,salbS4,asymS4,dpiaKu_S4,dpiaKa_S4,zKu_S4,zKa_S4,zW_S4,dmS4 = pyHB2.getsnowp2(log10(nwag/0.08),(iwc4),nfreq)
                    zW_S1=getZ_W(dmS1,zKu_S1,nws,dm_z,zW_S1)
                    zW_S2=getZ_W(dmS2,zKu_S2,nwg,dm_z,zW_S2)
                    zW_S3=getZ_W(dmS3,zKu_S3,nwh,dm_z,zW_S3)
                    zW_S4=getZ_W(dmS4,zKu_S4,nwag,dm_z,zW_S4)
                    zKu_S=log10(10.**(0.1*zKu_S1)+10.**(0.1*zKu_S2)+10.**(0.1*zKu_S3)+10.**(0.1*zKu_S4))*10
                    zKa_S=log10(10.**(0.1*zKa_S1)+10.**(0.1*zKa_S2)+10.**(0.1*zKa_S3)+10.**(0.1*zKa_S4))*10
                    zWS=log10(10.**(0.1*zW_S1)+10.**(0.1*zW_S2)+10.**(0.1*zW_S3)+10.**(0.1*zW_S4))*10
    
                    i0s1=int((zKu_S1-10*log10(nws/0.08)+12.)/0.5)
                    i0s1=max(i0s1,0)
                    i0s1=min(114,i0s1)
                    i0s2=int((zKu_S2-10*log10(nwg/0.08)+12.)/0.5)
                    i0s2=max(i0s2,0)
                    i0s2=min(114,i0s2)
                    i0s3=int((zKu_S3-10*log10(nwh/0.08)+12.)/0.5)
                    i0s3=max(i0s2,0)
                    i0s3=min(114,i0s3)
                    vs1=vt_S[i0s1,1:]
                    vs2=vt_S[i0s2,1:]
                    vs3=vt_S[i0s3,1:]
                    z_S1=array([zKu_S1,zKa_S1,zW_S1])
                    z_S2=array([zKu_S2,zKa_S2,zW_S3])
                    z_S3=array([zKu_S3,zKa_S3,zW_S3])
                    vS=(10**(0.1*z_S1)*vs1+10**(0.1*z_S2)*vs2+10**(0.1*z_S3)*vs3)/\
                        (10**(0.1*z_S1)+10**(0.1*z_S2)+10**(0.1*z_S3))
                
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
                if piaKu<0:
                    stop
                vT=(10**(0.1*zKu_S)*vS+10**(0.1*zKu_R)*vR)/\
                    (10**(0.1*zKu_S)+10**(0.1*zKu_R))
                vW_1d[k]=vT[2]-wm[k,j,i]
                vKu_1d[k]=vT[0]-wm[k,j,i]
                vKa_1d[k]=vT[1]-wm[k,j,i]
                print(piaKu,piaKa,piaW)
