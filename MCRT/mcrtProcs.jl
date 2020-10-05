function trace3D(c,muz,mux,muy,tau,x0,y0,dx,dz,nx,ny)
    cmu=abs(muz)
    xin=x0
    yin=y0
    l=0.
    dstep=0.0125
    zin=400
    ilayin=79
    tauc=0.0
    
    dstep1=(380.)/abs(muz)
    xin=xin+dstep1*mux              
    yin=yin+dstep1*muy
    zin=zin+dstep1*muz
    i01=Int(trunc(xin/dx))
    j01=Int(trunc(yin/dx))
    #println("$tau $tauc")
    ilayn=79
    while tauc<tau
        #global ilayn
        xin=xin+dstep*mux              
        yin=yin+dstep*muy
        zin=zin+dstep*muz
        ilayn=Int(trunc(zin/dz))
        if(ilayn>79)
            ilayn=79
        end
        i0n=Int(trunc(xin/dx))
        j0n=Int(trunc(yin/dx))
        
        if(i0n>nx-1)
            i0n=nx-1
        end
        if(j0n>ny-1)
            j0n=ny-1
        end
        if(i0n<=0)
            i0n=0
        end
        if(j0n<=0)
            j0n=0
        end
        if(i01>nx-1)
            i01=nx-1
        end
        if(j01>ny-1)
            j01=ny-1
        end
        if(i01<=0)
            i01=0
        end
        if(j01<=0)
            j01=0
        end
        if(ilayn<0)
            break
        end
        
        tauc=tauc+dstep*0.5*(c.kext[i01+1,j01+1,ilayin+1]+
                             c.kext[i0n+1,j0n+1,ilayn+1])
        i01=i0n
        j01=j0n
        ilayin=ilayn
        
        l=sqrt((xin-x0)^2+(yin-y0)^2+(zin-400)^2)
    end
    return ilayn,xin,yin,zin,l,tauc,i01,j01
end 
    
function attToRadar3D(c,x0,y0,z0,xRad,yRad,i0,j0,ilay,dx,dz,nx,ny)
    xin=x0
    yin=y0
    zin=z0
    dtot=sqrt((400.0-z0)^2+(xRad-x0)^2+(yRad-y0)^2)
    muz=(400.0-z0)/dtot
    mux=(xRad-x0)/dtot
    muy=(yRad-y0)/dtot
    dstep=0.025

    tauc=0
    i01=Int(trunc(x0/dx))
    j01=Int(trunc(y0/dx))
    i01=Int(trunc(i0))
    j01=Int(trunc(j0))
    ilayin=Int(trunc(ilay))    
    if(ilayin<0)
        ilayin=0
    end
    while zin<20
        xin=xin+dstep*mux              
        yin=yin+dstep*muy
        zin=zin+dstep*muz
        ilayn=Int(trunc(zin/dz))

        if(ilayn>79)
            ilayn=79
        end
        if(ilayn<0)
            ilayn=0
        end
        i0n=Int(trunc(xin/dx))
        j0n=Int(trunc(yin/dx))
        
        if(i0n>nx-1)
            i0n=nx-1
        end
        if(j0n>ny-1)
            j0n=ny-1
        end
        if(i0n<=0)
            i0n=0
        end
        if(j0n<=0)
            j0n=0
        end
        if(ilayn<0)
            break
        end
        tauc=tauc+dstep*0.5*(c.kext[i01+1,j01+1,ilayin+1]+
                                 c.kext[i0n+1,j0n+1,ilayn+1])
        i01=i0n
        j01=j0n
        ilayin=ilayn
    end
    l=sqrt((xin-x0)^2+(yin-y0)^2+(zin-400)^2)
    return tauc, mux,muy,muz,dtot
end


function SIGN(x)
    if x>0
        return 1
    else
        return -1
    end
end

function scattering2(mux,muy,muz,g)
    PI=pi
    ONE_MINUS_COSZERO=1e-12
    rnd = rand()
    if (abs(g)< 0.1e-4)
        costheta = 2.0*rnd - 1.0
    else
        temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd);
        costheta = (1.0 + g*g - temp*temp)/(2.0*g)
    end
    sintheta = sqrt(1.0 - costheta*costheta)
    psi = 2.0*PI*rand()
    cospsi = cos(psi)
    if (psi < PI)
        sinpsi = sqrt(1.0 - cospsi*cospsi)
    else
        sinpsi = -sqrt(1.0 - cospsi*cospsi)
    end
    
    if (1 - abs(muz) <= ONE_MINUS_COSZERO)
        uxx = sintheta * cospsi;
        uyy = sintheta * sinpsi;
        uzz = costheta * SIGN(muz)
    else
        temp = sqrt(1.0 - muz * muz)
        uxx = sintheta * (mux * muz * cospsi - muy * sinpsi)/temp +
            mux * costheta
        uyy = sintheta * (muy * muz * cospsi + mux * sinpsi)/ temp +
            muy * costheta
        uzz = -sintheta * cospsi * temp + muz * costheta
    #print (uxx*mux+uyy*muy+uzz*muz), costheta
    end
    return uxx,uyy,uzz, (uxx*mux+uyy*muy+uzz*muz)
end
        
function goFrom1To3D(x,y,z,i0,j0,ilay,mux,muy,muz,c,tau,dx,dz,nx,ny)
    dstep=0.0125
    xin=x
    yin=y
    zin=z
    keepgoing=1
    ic=ilay
    tauc=0
    dz=0.25
    ilayin=Int(trunc(ilay))
    i01=Int(trunc(i0))
    j01=Int(trunc(j0))
   
    while(zin>0 && zin<20 && keepgoing>0)
        xin=xin+dstep*mux              
        yin=yin+dstep*muy
        zin=zin+dstep*muz
        ilayn=Int(trunc(zin/dz))
        i0n=Int(trunc(xin/dx))
        j0n=Int(trunc(yin/dx))

        if(i0n>nx-1)
            i0n=nx-1
        end
        if(j0n>ny-1)
            j0n=ny-1
        end
        if(i0n<=0)
            i0n=0
        end
        if(j0n<=0)
            j0n=0
        end

        if(ilayn>c.n-1)
            if(zin>40)
                break
            end
        else
            tauc=tauc+dstep*0.5*(c.kext[i01+1,j01+1,ilayin+1]+
                                 c.kext[i0n+1,j0n+1,ilayn+1])
        end
        ilayin=ilayn  
        if(tauc>tau)
            break
        end
        i01=i0n
        j01=j0n
    end
    l=sqrt((x-xin)^2+(y-yin)^2+(z-zin)^2)

    return xin, yin, zin, ilayin,tau-tauc, i01, j01, l
end
