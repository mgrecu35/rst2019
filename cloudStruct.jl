
struct cloud
    kext :: Array{Float64,3};
    salb :: Array{Float64,3};
    g    :: Array{Float64,3};
    ebr  :: Array{Float64,3};
    z    :: Array{Float64,3};
    zTop :: Array{Float64,3};
    zBot :: Array{Float64,3};
    tauU :: Array{Float64,3};
    tauD :: Array{Float64,3};
    vTerm ::  Array{Float64,3};
    zObs  :: Array{Float64,3};
    #zTrue :: Array{Float64,3};
    n :: Int64;
end

using PyCall
push!(PyVector(pyimport("sys")["path"]), "")
@pyimport pyMCRTPy6 as pymcrt

ifreq=1
function henyeym( g, mu)
    x=(1-g*g)/(1+g*g-2*g*mu)^1.5
    return x
end

function antenna2F(ang1,sigma1,ang2,sigma2)
    f=exp(-4*log(2.0)*((ang1/sigma1)*(ang1/sigma1)+
                     (ang2/sigma2)*(ang2/sigma2)))
    return f
end
include("mcrtProcs.jl")
include("mcrt.jl")
fname="wrfJune11_2014.Fields_dn5.nc"

for ifreq=0:0
    nx=130
    ny=30
    pymcrt.readdata(fname,ifreq,nx,ny)
    kext=pymcrt.getkext(nx,ny)
    g=pymcrt.getg(nx,ny)
    salb=pymcrt.getsalb(nx,ny)
    zTop=pymcrt.getztop(nx,ny)
    zBot=pymcrt.getzbot(nx,ny)
    ebr=pymcrt.getebr(nx,ny)
    z=pymcrt.getz(nx,ny)
    zobs=pymcrt.getzobs(nx,ny)
    tauD=pymcrt.gettauD(nx,ny)
    tauU=pymcrt.gettauU(nx,ny)
    vterm=pymcrt.getvterm(nx,ny)
    
    c=cloud(kext,salb,g,ebr,z,zTop,zBot,tauU,tauD,vterm,zobs,80)

    
    println("finished reading input data")
    rte_calc=1
    dx=1300.0
    nx,ny,nz=size(c.vTerm)
    nx1=nx
    zMS=zeros(nx1,nz).-99
    vMS=zeros(nx1,nz).-99
    piaKa=zeros(nx1).-99
    iy=10
    if (rte_calc==1)
        for ix=1:nx1
        #for ix=150:160  
            z1d,vdoppler,pia=mcrt(c,ix,iy,dx,nx,ny)
            zMS[ix,:]=z1d
            vMS[ix,:]=vdoppler
            piaKa[ix]=pia
            print(" $ix")
        end
    end
    
    z1=zobs[1:nx1,iy,:]
    z1true=z[1:nx1,iy,:]
    v1=vterm[1:nx1,iy,:]
 
    pymcrt.saveMS(copy(zMS),copy(vMS),copy(z1),copy(v1),copy(piaKa),
                  copy(z1true),iy,"iphex.output$(iy)_freq$(ifreq)_0175_dn5.nc",fname)
end



