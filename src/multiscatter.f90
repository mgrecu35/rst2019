subroutine multiscatterF(nrange,kext,salb,g,zTrue,zMS,dr,noMS,&
     alt,theta,freq,noNorm)
implicit none
integer :: nrange,noMS,i,noNorm
real :: extinct(nrange),salb(nrange),g(nrange),ext2bscatt(nrange), lambd
real , intent(out) :: zMS(nrange)
real :: bscatt(nrange)
real :: alt,dr,theta,freq

real :: kext(nrange),zTrue(nrange)!,noMS,freq,theta,dr,noNorm=inp
real :: pi4, pi, lamb4, Z(nrange)
!#  float lamb=0.00857;
!#  float kextFort[88],z35med[88], ext2bscattFort[88], salbFort[88], gFort[88]; 
!#  float pi4=97.409091034;
!#  float lamb4=5.39415e-9;
!#  float Z=pow(10.,0.1*z35med[i])*pi4/1e18/lamb4/4*0.93;
integer :: k
real :: dz
pi=3.1415
pi4=pi**4
lambd=300./freq/1e3
lamb4=lambd**4
Z=10.**(0.1*zTrue)*pi4/1e18/lamb4/4*0.93
dZ=log10(pi4/1e18/lamb4/4*0.93)*10
kext=kext*1e-3

ext2bscatt=1000.
if (noNorm==1) then
   salb=salb*1e-3
   salb=salb/kext
end if
extinct=kext
do k=1,nrange
   if (zTrue(k)>-10) then
      ext2bscatt(k)=kext(k)/Z(k)
   else
      ext2bscatt(k)=1000.
   endif
   !print*, ext2bscatt(k),k, salb(k), g(k), kext(k), zTrue(k)
end do

            
dr=dr*1000

call multiscatter2(nrange, extinct, ext2bscatt, salb, g, &
     bscatt, noMS)




do i=1,nrange
!   print*, bscatt(i)
enddo
zMS=log10(bscatt)*10-dZ

end subroutine multiscatterF

subroutine multiscatterf2(nrange,extinct,ext2bscatt,salb,g,bscatt,&
     lambd,noMS,&
     alt,dr,theta,freq)
  implicit none
  integer :: nrange,noMS,i
  real :: extinct(nrange),salb(nrange),g(nrange),ext2bscatt(nrange), lambd
  real , intent(out) :: bscatt(nrange)
  real :: alt,dr,theta,freq
  
!  call multiscatterF(nrange,extinct,ext2bscatt,salb,g,bscatt,lambd,noMS)
  
end subroutine multiscatterf2
