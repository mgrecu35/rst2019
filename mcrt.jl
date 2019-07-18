function mcrt(c,ix,iy,dx,nx,ny)
    nphot=3*50000
    nsscat=12
    SS=zeros(80)
    HOS=zeros(80)
    icrs=zeros(Int,80)
    wV=zeros(80,80)
    wZ=zeros(80,80)
    vdoppler=zeros(80)
    vdopplerR=zeros(80).-99.0
    z1d=zeros(80).-99.0
    bins=zeros(Int,nsscat+1)
    muzs=zeros(nsscat+1)
    surfReturn=0.0
    #dx=250.0
    #ix=100
    #iy=10
    println(" ",ix," ",iy," ",dx)
    dz=0.25
    dztop=380.0
    sigma=0.35/180*pi*0.5
    nz=80
    for i=1:nphot
        x0=(ix-0.5)*dx
        y0=(iy-0.5)*dx
        G=0.5 #  G  is the integral of gain^2 function
        #global nx, ny
        sx=sigma/2/(2*log(2))^0.5
        i00=Int(trunc(x0/dx))
        j00=Int(trunc(y0/dx))
        alpha1=randn()*sx
        beta1=randn()*2*pi
        r1=rand()
        r2=rand()
        eps=(sx*sqrt(-2*log(r1))*cos(2*pi*r2))
        eta=(sx*sqrt(-2*log(r1))*sin(2*pi*r2))
        x1=tan(eps)
        x2=tan(eta)
        mux1=x1/sqrt(1+x1*x1+x2*x2)
        muy1=x2/sqrt(1+x1*x1+x2*x2)
        muz1=-1.0/sqrt(1+x1*x1+x2*x2)
        tau1=-log(1-rand())
        i1,x1,y1,z1,l1,tauc,i01,j01=trace3D(c,muz1,mux1,muy1,tau1,x0,y0,dx,dz,nx,ny)
        att11, mux11,muy11,muz11,dtot1 = 
            attToRadar3D(c,x1,y1,z1,
                         x0,y0,i01,j01,i1,dx,dz,nx,ny)
        
        ndc=Int(0)
        if z1>0 && l1>dztop 
            nr=Int(trunc(l1/0.25))-1520
            if(nr<80 && nr>=0)
                icrs[nr+1]=icrs[nr+1]+1
                if(c.ebr[i01+1,j01+1,i1+1]>0) 
                    deltaP=1.0/c.ebr[i01+1,j01+1,i1+1]*
                    exp(-tau1)*antenna2F(acos(-muz1),sigma,0,sigma)/G
                    SS[nr+1]=SS[nr+1]+deltaP
                    vdoppler[nr+1]=vdoppler[nr+1]+c.vTerm[i01+1,j01+1,i1+1]*(muz11-muz1)/2.0*deltaP
                    wV[nr+1,nr+1]=wV[nr+1,nr+1]+(muz11-muz1)/2.0*deltaP
                end
                if nr+1==80
                    surfReturn=surfReturn+10.0^8.0/(c.kext[i01+1,j01+1,i1+1]+1e-5)*
                    exp(-tau1)*antenna2F(acos(-muz1),sigma,0,sigma)/G
                end
            end
            if(c.ebr[i01+1,j01+1,i1+1]>0)    
                weight=1.0/c.ebr[i01+1,j01+1,i1+1]/
                henyeym(c.g[i01+1,j01+1,i1+1],-1.0)
            else
                weight=0.
            end
            dshift=0.0
            for iscatt=1:nsscat
                i01b=Int(max(0,min(nx-1,Int(i01))))
                j01b=Int(max(0,min(ny-1,Int(j01))))
                i1b=Int(max(0,min(nz-1,Int(i1))))
                nmux,nmuy,nmuz,atheta=scattering2(mux1,muy1,muz1,
                                                  c.g[i01b+1,j01b+1,i1b+1])
                ndc=ndc+1
                bins[ndc]=Int(i1b)+1
                muzs[ndc]=(nmuz-muz1)/2.
                dshift=dshift+c.vTerm[i01b+1,j01b+1,i1b+1]*(nmuz-muz1)/2.
                tau=-log(rand())
                x21,y21,z21,i21,tauleft,i021,j021,l2=
                goFrom1To3D(x1,y1,z1,
                            i01,j01,i1,
                            nmux,nmuy,nmuz,c,tau,dx,dz,nx,ny)
                l2=sqrt((x21-x1)^2+(y21-y1)^2+(z21-z1)^2)
                tauc3,mux3,muy3,muz3,l3=
                attToRadar3D(c,x21,y21,z21,
                             x0,y0,i021,j021,i21,dx,dz,nx,ny)
                l3=sqrt((x21-x0)^2+(y21-y0)^2+(z21-400)^2)
                ctheta3=(nmux*mux3+nmuy*muy3+nmuz*muz3)
                i021b=Int(max(0,min(nx-1,i021)))
                j021b=Int(max(0,min(ny-1,j021)))
                i21b=Int(max(0,min(nz-1,i21)))
                t=weight*exp(-tauc3)*antenna2F(acos(muz3),sigma,0,sigma)/G*
                henyeym(c.g[i021b+1,j021b+1,i21b+1],ctheta3)
                weight=weight*c.salb[i021b+1,j021b+1,i21b+1]
                ltot=(l1+l2+l3)/2
                nr=Int(trunc(ltot/0.25))-1520
                
                if(nr<80 && nr>=0)
                    if(iscatt<nsscat)
                        HOS[nr+1]=HOS[nr+1]+t
                    end
                    wV[nr+1,80-i21b]=wV[nr+1,80-i21b]+(muz3-nmuz)/2.0*t
                    vdoppler[nr+1]=vdoppler[nr+1]+(dshift+(muz3-nmuz)/2*c.vTerm[i021b+1,j021b+1,i21b+1])*t
                    for indc=1:ndc
                        wV[nr+1,81-bins[indc]]=wV[nr+1,81-bins[indc]]+muzs[indc]*t
                        wZ[nr+1,81-bins[indc]]=wZ[nr+1,81-bins[indc]]+t
                    end
                end
                
                l1=l1+l2
                x1=x21
                y1=y21
                z1=z21
                i1=i21
                i01=i021
                j01=j021
                if(z1<0)
                    break
                end
                if(z1>20)
                    break
                end
                mux1=nmux
                muy1=nmuy
                muz1=nmuz
            end
        end
    end
    
    for i=1:80
        if SS[i] + HOS[i] >0
            z1d[i]=log10((SS[i]+HOS[i])/nphot/dz)*10.0
            vdopplerR[i]=vdoppler[i]/((SS[i]+HOS[i]))
            #println(i," ",log10(SS[i]/nphot/dz)*10.0," ",c.zObs[ix,iy,81-i]," ",log10((SS[i]+HOS[i])/nphot/dz)*10.0,
            #        " ",vdopplerR[i]," ",c.vTerm[ix,iy,81-i])
            #else
            
            #println(i," ",-99, " ", c.zObs[ix,iy,81-i])
        end
    end
    pia=80-10.0*log10(surfReturn/nphot/dz)
    println()
#exit()
    return z1d,vdopplerR,pia
end
