#include "scalelib.h"
!OCL SERIAL
module mod_scladv_slppm

  use scale_io
  use acc_amps
  use com_amps
  implicit none
  public

contains

! semi Lagrangian piecewise parabolic method
subroutine cal_ppmcoef(ns,n2,n2max,xx,CPPM)
  implicit none
  integer,intent(in) :: ns,n2,n2max
  real(MP_KIND), intent(in), dimension(n2) :: xx ! n2max -> n2
  real(MP_KIND), intent(out), dimension(11,n2max) :: CPPM ! how n2max -> n2

  real(MP_KIND),dimension(n2) :: dx
  integer :: i,j
  
  ! calculate the difference of grid interfaces
  do i=2,n2-1
     dx(i)=xx(i)-xx(i-1)
  end do
  dx(1)=dx(2)
  dx(n2)=dx(n2-1)

  ! coefficients for calculation of slope
  do i=2,n2-1
     !! for da1, these are unitless
     ! for a(i-1)
     CPPM(1,i)=-(dx(i)+2.0*dx(i+1))/(dx(i-1)+dx(i))
     ! for a(i)
     CPPM(2,i)=-(2.0*dx(i-1)+dx(i))/(dx(i+1)+dx(i))+(dx(i)+2.0*dx(i+1))/(dx(i-1)+dx(i))
     ! for a(i+1)
     CPPM(3,i)=(2.0*dx(i-1)+dx(i))/(dx(i+1)+dx(i))
     do j=1,3
        CPPM(j,i)=CPPM(j,i)*dx(i)/(dx(i-1)+dx(i)+dx(i+1))
     end do
     
     !! for da2, these have unit of L^(-2)
     ! for a(i-1)
     CPPM(4,i)=1.0/(dx(i)+dx(i-1))/(dx(i-1)+dx(i)+dx(i+1))
     ! for a(i)
     CPPM(5,i)=(-1.0/(dx(i+1)+dx(i))-1.0/(dx(i)+dx(i-1)))/(dx(i-1)+dx(i)+dx(i+1))
     ! for a(i+1)
     CPPM(6,i)=1.0/(dx(i+1)+dx(i))/(dx(i-1)+dx(i)+dx(i+1))
     
  end do
  
  
  ! coefficients for interpolation of interface values
  do i=ns,n2-2
     ! for a(i), dimensionless
     CPPM(7,i)=1.0-dx(i)/(dx(i)+dx(i+1))+&
          ((dx(i-1)+dx(i))/(2.0*dx(i)+dx(i+1))-(dx(i+2)+dx(i+1))/(2.0*dx(i+1)+dx(i)))*&
          (-1.0)*2.0*dx(i+1)*dx(i)/(dx(i)+dx(i+1))/&
          (dx(i-1)+dx(i)+dx(i+1)+dx(i+2))
     
     ! for a(i+1), dimensionless
     CPPM(8,i)=dx(i)/(dx(i)+dx(i+1))+&
          ( (dx(i-1)+dx(i))/(2.0*dx(i)+dx(i+1))-(dx(i+2)+dx(i+1))/(2.0*dx(i+1)+dx(i)) )*&
          (1.0)*2.0*dx(i+1)*dx(i)/(dx(i)+dx(i+1))/&
          (dx(i-1)+dx(i)+dx(i+1)+dx(i+2))
     
     
     ! for da(i), dimensionless
     CPPM(9,i)=dx(i+1)*(dx(i+1)+dx(i+2))/(dx(i)+2.0*dx(i+1))/&
            (dx(i-1)+dx(i)+dx(i+1)+dx(i+2))            
     ! for da(i+1), dimensionless
     CPPM(10,i)=(-1.0)*dx(i)*(dx(i-1)+dx(i))/(2.0*dx(i)+dx(i+1))/&
          (dx(i-1)+dx(i)+dx(i+1)+dx(i+2))
     
!!c       aLR(i+1,1)=a(i)+(a(i+1)-a(i))*dx(i)/(dx(i)+dx(i+1))+&
!!c            (&
!!c            ((dx(i-1)+dx(i))/(2.0*dx(i)+dx(i+1))-(dx(i+2)+dx(i+1))/(2.0*dx(i+1)+dx(i)))*&
!!c            (a(i+1)-a(i))*2.0*dx(i+1)*dx(i)/(dx(i)+dx(i+1))-&
!!c            da(i+1)*dx(i)*(dx(i-1)+dx(i))/(2.0*dx(i)+dx(i+1))+&
!!c            da(i)*dx(i+1)*(dx(i+1)+dx(i+2))/(dx(i)+2.0*dx(i+1))&
!!c            )/&
!!c            (dx(i-1)+dx(i)+dx(i+1)+dx(i+2))

     ! dimension is L^2
     CPPM(11,i)=-1.0/(0.5*dx(i+1)+dx(i)+0.5*dx(i-1))*&
                ( (0.5*dx(i)+0.5*dx(i-1))**3.0+(0.5*dx(i+1)+0.5*dx(i))**3.0)

  end do
  CPPM(7,n2-1)=CPPM(7,n2-2)
  CPPM(8,n2-1)=CPPM(8,n2-2)
  CPPM(9,n2-1)=CPPM(9,n2-2)
  CPPM(10,n2-1)=CPPM(10,n2-2)
  CPPM(11,n2-1)=CPPM(11,n2-2)
  
  
end subroutine cal_ppmcoef
!!c        call cal_ppmcoef_eta(nxp,nyp,nxpmax,nypmax,ngrid,maxgrds
!!c     *       ,zz,a(idz1),integera(ik1eta),ADVPPMZE)

subroutine cal_ppmcoef_eta(n1,zz,k1eta,CPPME)
  implicit none
  integer,intent(in) :: n1
  real(MP_KIND), intent(in), dimension(n1) :: zz
  integer, intent(in) :: k1eta
  real(MP_KIND),dimension(11) :: CPPME

  real(MP_KIND) :: dz1, dzkp1, dzkp2, dzkm1
  integer :: k
  

  dzkp1=zz(k1eta+1)-zz(k1eta)
  dzkp2=zz(k1eta+2)-zz(k1eta+1)
  dz1 = zz(2) - zz(1)
  dzkm1=dz1
  ! calculate the difference of grid interfaces
  ! coefficients for calculation of slope
  !! for da1
  ! for a(i-1)
  CPPME(1)=-(dz1+2.0_MP_KIND*dzkp1)/(dzkm1+dz1)
  ! for a(i)
  CPPME(2)=-(2.0_MP_KIND*dzkm1+dz1)/(dzkp1+dz1)+(dz1+2.0_MP_KIND*dzkp1)/(dzkm1+dz1)
  ! for a(i+1)
  CPPME(3)=(2.0_MP_KIND*dzkm1+dz1)/(dzkp1+dz1)
  do k=1,3
     CPPME(k)=CPPME(k)*dz1/(dzkm1+dz1+dzkp1)
  end do
        
  !! for da2
  ! for a(i-1)
  CPPME(4)=1.0_MP_KIND/(dz1+dzkm1)/(dzkm1+dz1+dzkp1)
  ! for a(i)
  CPPME(5)=(-1.0_MP_KIND/(dzkp1+dz1)-1.0_MP_KIND/(dz1+dzkm1))/(dzkm1+dz1+dzkp1)
  ! for a(i+1)
  CPPME(6)=1.0_MP_KIND/(dzkp1+dz1)/(dzkm1+dz1+dzkp1)
        
  ! for a(i)
  CPPME(7)=1.0_MP_KIND-dz1/(dz1+dzkp1)+&
       ((dzkm1+dz1)/(2.0_MP_KIND*dz1+dzkp1)-(dzkp2+dzkp1)/(2.0_MP_KIND*dzkp1+dz1))*&
       (-1.0_MP_KIND)*2.0*dzkp1*dz1/(dz1+dzkp1)/&
       (dzkm1+dz1+dzkp1+dzkp2)
  
  ! for a(i+1)
  CPPME(8)=dz1/(dz1+dzkp1)+&
       ((dzkm1+dz1)/(2.0_MP_KIND*dz1+dzkp1)-(dzkp2+dzkp1)/(2.0_MP_KIND*dzkp1+dz1))*&
       (1.0_MP_KIND)*2.0_MP_KIND*dzkp1*dz1/(dz1+dzkp1)/&
       (dzkm1+dz1+dzkp1+dzkp2)
        
        
  ! for da(i)
  CPPME(9)=dzkp1*(dzkp1+dzkp2)/(dz1+2.0_MP_KIND*dzkp1)/&
       (dzkm1+dz1+dzkp1+dzkp2)            
  ! for da(i+1)
  CPPME(10)=(-1.0_MP_KIND)*dz1*(dzkm1+dz1)/(2.0_MP_KIND*dz1+dzkp1)/&
       (dzkm1+dz1+dzkp1+dzkp2)
     
  CPPME(11)=-1.0/(0.5_MP_KIND*dzkp1+dz1+0.5_MP_KIND*dzkm1)*&
       ( (0.5_MP_KIND*dz1+0.5_MP_KIND*dzkm1)**3+(0.5_MP_KIND*dzkp1+0.5_MP_KIND*dz1)**3)

  return
end subroutine cal_ppmcoef_eta

subroutine cal_aLRa6_y(ns,ne,n3,n3max,lbin,istp,CPPM,a,aLR,a6)
  implicit none
  ! ns and ne is the starting and ending grid cell for aLR and a6
  integer,intent(in) :: ns,ne,n3,n3max,lbin,istp
  real(PS) :: a(lbin),a6(n3),aLR(n3,2)
  real(PS) :: CPPM(11,n3max) ! how about n3max -> n3
  real(PS),dimension(n3) :: da,da2
  real(PS) :: ita,ita_w, aLd, aRd
  real(PS),parameter  :: ita_1=20.0,ita_2=0.05,eps=0.01
  integer :: i,j,im1,ip1
  
!!c  open(1,file='check_cppmy.dat') 
!!c  do j=1,n3
!!c     write(1,'(15ES15.6)') (CPPM(i,j),i=1,11)
!!c  end do
!!c  close(1)
!!c  stop
  
  ! calculate the slopes
  do i=max(2,ns-1),min(n3-1,ne+1)
     im1=max(1,i-1)
     ip1=min(n3,i+1)
     da(i)=CPPM(1,i)*a(im1)+CPPM(2,i)*a(i)+CPPM(3,i)*a(ip1)
     
     if((a(ip1)-a(i))*(a(i)-a(im1))>0.0) then
        da(i)=min(abs(da(i)),2.0*abs(a(i)-a(im1)),2.0*abs(a(ip1)-a(i)))*sign(1.0_PS,da(i))
     else
        da(i)=0.0
     end if
     
     da2(i)=CPPM(4,i)*a(im1)+CPPM(5,i)*a(i)+CPPM(6,i)*a(ip1)
  end do

  ! no gradient at the boundary
  if(ns<=2) then
     da(1)=0.0
     da2(1)=0.0
  endif
  if(ne>=n3-1) then
     da(n3)=0.0
     da2(n3)=0.0
  endif
 
  ! interpolation at the left and right interfaces of cell i
  do i=max(2,ns-1),min(n3-2,ne)
     aLR(i,2)=CPPM(7,i)*a(i)+CPPM(8,i)*a(i+1)+CPPM(9,i)*da(i)+CPPM(10,i)*da(i+1)
     aLR(i+1,1)=aLR(i,2)
  end do
  if(ns<=2) then
     aLR(2,1)=a(1)
  endif
  if(ne>=n3-1) then
     aLR(n3-1,2)=a(n3)
  endif

  if(istp==1) then
     
     ! treatment of discontinuity in the cell
     do i=max(2,ns),min(n3-1,ne)
        if(-da2(i+1)*da2(i-1)>0.0.or.abs(a(i+1)-a(i-1))-eps*min(abs(a(i+1)),abs(a(i-1)))>0.0) then
           ! dimensionless
           ita_w=CPPM(11,i)*(da2(i+1)-da2(i-1))/(a(i+1)-a(i-1))

!!c           ita_w=-((da2(i+1)-da2(i-1))/(0.5*dx(i+1)+dx(i)+0.5*dx(i-1)))*&
!!c                (( (0.5*dx(i)+0.5*dx(i-1))**3.0+(0.5*dx(i+1)+0.5*dx(i))**3.0)/(a(i+1)-a(i-1)))
        else
           ita_w=0.0
        end if
        ita=max(0.0_PS,min(ita_1*(ita_w-ita_2),1.0_PS))
        
        aLd=a(i-1)+0.5*da(i-1)
        aRd=a(i+1)-0.5*da(i+1)
        
        aLR(i,1)=aLR(i,1)*(1.0-ita)+aLd*ita
        aLR(i,2)=aLR(i,2)*(1.0-ita)+aRd*ita
     end do
  end if

  ! monotonicity algorithm
  do i=max(2,ns),min(n3-1,ne)
     if((aLR(i,2)-a(i))*(a(i)-aLR(i,1))<0.0) then
        aLR(i,1)=a(i)
        aLR(i,2)=a(i)          
     else if((aLR(i,2)-aLR(i,1))*(a(i)-0.5*(aLR(i,1)+aLR(i,2)))>&
          (aLR(i,2)-aLR(i,1))*(aLR(i,2)-aLR(i,1))/6.0) then
        aLR(i,1)=3.0*a(i)-2.0*aLR(i,2)
     else if((aLR(i,2)-aLR(i,1))*(a(i)-0.5*(aLR(i,1)+aLR(i,2)))<&
          -(aLR(i,2)-aLR(i,1))*(aLR(i,2)-aLR(i,1))/6.0) then
        aLR(i,2)=3.0*a(i)-2.0*aLR(i,1)
     end if
     
  end do
  if(ns==1) then
     aLR(1,1)=a(1)
     aLR(1,2)=a(1)
  endif
  if(ne==n3) then
     aLR(n3,1)=a(n3)
     aLR(n3,2)=a(n3)
  endif
 
  ! calculation of coefficient
  do i=ns,ne
     a6(i)=6.0*(a(i)-0.5*(aLR(i,1)+aLR(i,2)))
  end do

end subroutine cal_aLRa6_y


subroutine cal_aLRa6_z(i,j,ns,ne,k1eta,n1,n1mx,n2mx,n3mx,lbin,istp,CPPM,CPPME,a,aLR,a6)
  implicit none
  ! ns and ne is the starting and ending grid cell for aLR and a6
  integer,intent(in) :: i,j,ns,ne
  integer,intent(in) :: k1eta,n1,n1mx,n2mx,n3mx,lbin,istp
  real(PS) :: a(lbin),a6(n1)
  real(PS) :: aLR(n1,2)
  real(RP) :: CPPM(11,n1mx),CPPME(11)
  real(PS),dimension(n1) :: da,da2
  real(PS) :: ita,ita_w, aLd, aRd
  real(PS),parameter  :: ita_1=20.0,ita_2=0.05,eps=0.01
  integer :: n,m,k,l
  
  
!!c  open(1,file='check_cppm_eta.dat') 
!!c    write(1,'(11ES15.6)') (CPPME(m),m=1,11)
!!c  close(1)
!!c  stop

  ! calculate the slopes
  if(ns<=k1eta+1) then
     da(k1eta)=CPPME(1)*a(k1eta-1)+CPPME(2)*a(k1eta)+CPPME(3)*a(k1eta+1)
     if((a(k1eta+1)-a(k1eta))*(a(k1eta)-a(k1eta-1))>0.0) then
        da(k1eta)=min(abs(da(k1eta)),2.0*abs(a(k1eta)-a(k1eta-1)),2.0*abs(a(k1eta+1)-a(k1eta)))*sign(1.0_PS,da(k1eta))
     else
        da(k1eta)=0.0
     end if
     da2(k1eta)=CPPME(4)*a(k1eta-1)+CPPME(5)*a(k1eta)+CPPME(6)*a(k1eta+1)
  end if

  do k=max(ns-1,k1eta+1),min(n1-1,ne+1)
     da(k)=CPPM(1,k)*a(k-1)+CPPM(2,k)*a(k)+CPPM(3,k)*a(k+1)
     
     if((a(k+1)-a(k))*(a(k)-a(k-1))>0.0) then
        da(k)=min(abs(da(k)),2.0*abs(a(k)-a(k-1)),2.0*abs(a(k+1)-a(k)))*sign(1.0_PS,da(k))
     else
        da(k)=0.0
     end if
     
     da2(k)=CPPM(4,k)*a(k-1)+CPPM(5,k)*a(k)+CPPM(6,k)*a(k+1)
  end do

  if(ns==k1eta) then
     da(k1eta-1)=0.0
     da2(k1eta-1)=0.0
  endif
  if(ne>=n1-1) then
     da(n1)=0.0
     da2(n1)=0.0
  endif 
  
  ! interpolation at the left and right interfaces of cell i
  if(ns<=k1eta+1) then
     aLR(k1eta,2)=CPPME(7)*a(k1eta)+CPPME(8)*a(k1eta+1)+CPPME(9)*da(k1eta)+CPPME(10)*da(k1eta+1)
     aLR(k1eta+1,1)=aLR(k1eta,2)
  end if
  do k=max(ns-1,k1eta+1),min(n1-2,ne)
     aLR(k,2)=CPPM(7,k)*a(k)+CPPM(8,k)*a(k+1)+CPPM(9,k)*da(k)+CPPM(10,k)*da(k+1)
     aLR(k+1,1)=aLR(k,2)
  end do
  if(ns==k1eta) then
     aLR(k1eta,1)=aLR(k1eta,2)
!!c     aLR(k1eta,1)=0.0
  endif
  if(ne>=n1-1) then
     aLR(n1-1,2)=a(n1)
  endif

  if(istp==1) then
     if(ns==k1eta) then
        if(-da2(k1eta+1)*da2(k1eta-1)>0.0.or.abs(a(k1eta+1)-a(k1eta-1))-eps*min(abs(a(k1eta+1)),abs(a(k1eta-1)))>0.0) then
           ita_w=CPPME(11)*(da2(k1eta+1)-da2(k1eta-1))/(a(k1eta+1)-a(k1eta-1))
        else
           ita_w=0.0
        end if
        ita=max(0.0_PS,min(ita_1*(ita_w-ita_2),1.0_PS))
        
        aLd=a(k1eta-1)+0.5*da(k1eta-1)
        aRd=a(k1eta+1)-0.5*da(k1eta+1)
        
        aLR(k1eta,1)=aLR(k1eta,1)*(1.0-ita)+aLd*ita
        aLR(k1eta,2)=aLR(k1eta,2)*(1.0-ita)+aRd*ita
     end if
     
     ! treatment of discontinuity in the cell
     do k=max(ns,k1eta+1),min(n1-1,ne)
        if(-da2(k+1)*da2(k-1)>0.0.or.abs(a(k+1)-a(k-1))-eps*min(abs(a(k+1)),abs(a(k-1)))>0.0) then
           ita_w=CPPM(11,k)*(da2(k+1)-da2(k-1))/(a(k+1)-a(k-1))

!!c           ita_w=-((da2(i+1)-da2(i-1))/(0.5*dx(i+1)+dx(i)+0.5*dx(i-1)))*&
!!c                (( (0.5*dx(i)+0.5*dx(i-1))**3.0+(0.5*dx(i+1)+0.5*dx(i))**3.0)/(a(i+1)-a(i-1)))
        else
           ita_w=0.0
        end if
        ita=max(0.0_PS,min(ita_1*(ita_w-ita_2),1.0_PS))
        
        aLd=a(k-1)+0.5*da(k-1)
        aRd=a(k+1)-0.5*da(k+1)
        
        aLR(k,1)=aLR(k,1)*(1.0-ita)+aLd*ita
        aLR(k,2)=aLR(k,2)*(1.0-ita)+aRd*ita
     end do
  end if
  
  ! monotonicity algorithm
  do k=ns,min(n1-1,ne)


     if((aLR(k,2)-a(k))*(a(k)-aLR(k,1))<0.0) then
       aLR(k,1)=a(k)
       aLR(k,2)=a(k)          
     else 
!tmp       if(abs(alr(k,1))>1.0e+20.or.abs(alr(k,2))>1.0e+20) then
!tmp         write(*,*) "ppm check",ns,ne,n1,k1eta
!tmp         do l=ns,min(n1-1,ne)
!tmp           write(*,*) "check aLR",l,aLR(l,1),aLR(l,2),a(l)
!tmp         enddo
!tmp
!tmp         write(*,*) CPPME(7),CPPME(8),CPPME(9),CPPME(10)
!tmp         write(*,*) a(k1eta),a(k1eta+1),da(k1eta),da(k1eta+1)
!tmp       endif

       if((aLR(k,2)-aLR(k,1))*(a(k)-0.5*(aLR(k,1)+aLR(k,2)))>&
          (aLR(k,2)-aLR(k,1))*(aLR(k,2)-aLR(k,1))/6.0) then


         aLR(k,1)=3.0*a(k)-2.0*aLR(k,2)
       else if((aLR(k,2)-aLR(k,1))*(a(k)-0.5*(aLR(k,1)+aLR(k,2)))<&
          -(aLR(k,2)-aLR(k,1))*(aLR(k,2)-aLR(k,1))/6.0) then
         aLR(k,2)=3.0*a(k)-2.0*aLR(k,1)
       endif
     end if
     
  end do
! this is for just upstream scheme
!  if(ns==k1eta) then
!     aLR(k1eta,1)=a(k1eta)
!     aLR(k1eta,2)=a(k1eta)
!  endif
  if(ne==n1) then
     aLR(n1,1)=a(n1)
     aLR(n1,2)=a(n1)
  endif 
  
  ! calculation of coefficient
  do k=ns,ne
     a6(k)=6.0*(a(k)-0.5*(aLR(k,1)+aLR(k,2)))
  end do
  
end subroutine cal_aLRa6_z

subroutine cal_flux_z(k,k1eta,n1,n1mx,lbin,zz,dz,dt,den,a,aLR,a6,zarr,w,flux,from)
  use scale_prc, only: &
     PRC_abort
  implicit none
!!c  integer,intent(in) :: k,k1eta,n1,n1mx
!!c  real,dimension(n1) :: a,a6,dz,den
!!c  real,dimension(n1mx) :: z
!!c  real,dimension(n1,2) :: aLR
!!c  real,intent(in)    :: w,zarr
!!c  real,intent(inout) :: flux
  integer,intent(in) :: k,k1eta,n1,n1mx,lbin
  real(PS) :: a(lbin),a6(n1),dz(n1)
  real(RP) :: den(lbin),zz(n1mx),zarr
  real(PS),dimension(n1,2) :: aLR
  real(PS),intent(in)    :: w
  real(PS),intent(inout) :: flux
  character*(*) :: from

  real(PS) :: dt,s1,s2,yisum,zdep
  integer :: kdep,kck
  integer :: ka,n
  
  
  ! find departure point
  ! search for departure point
  kck=0
  zdep=zarr-w*dt

  if(zdep<zz(k1eta)-dz(k1eta)) then
     zdep=zz(k1eta)-dz(k1eta)
     kck=1
  else if(zdep>zz(n1-1)) then
     ! case of periodic boundary
     zdep=zz(n1-1)
     kck=1
  end if
 
  ! find the interface departure point
!!c  kdep=int((zdep-zz(1))/dz(2))+2
  if(w<0.0) then
    do n=k+1,n1
       if(zz(n)>zdep.and.zz(n-1)<=zdep) then
          kdep=n
          exit
       end if
    end do
  else
    do n=k,k1eta,-1
!!c       if(zz(n)>zdep.and.zz(n-1)<=zdep) then
       if(zz(n)>=zdep.and.zz(n-1)<zdep) then
          kdep=n
          exit
       end if
    end do
 end if
  
  
!!c  ka=kdep-1
  ka=kdep
  
  flux=0.0
  if(w>0.0) then
     ! ka is the cell number for which the interpolation function is used.
     ka=kdep
     ! calculate integer flux
     yisum=0.0
!!c     if(kck==1) then
!!c        do n=ka+1,n1-1
!!c           flux=flux+a(n)*den(n)*dvol(n)
!!c           yisum=yisum+dvol(n)
!!c        end do
!!c        do n=2,k
!!c           flux=flux+a(n)*den(n)*dvol(n)
!!c           yisum=yisum+dvol(n)
!!c        end do
!!c     else
!!c        do n=ka+1,k
!!c           if(k>=1) then
!!c              flux=flux+a(n)*den(n)*dvol(n)
!!c              yisum=yisum+dvol(n)
!!c           else
!!c              ! periodic condition
!!c              flux=flux+a(n1-1-k)*den(n1-1-k)*dvol(n1-1-k)
!!c              yisum=yisum+dvol(n1-1-k)
!!c           end if
!!c        end do
!!c     end if
     do n=ka+1,k
        flux=flux+a(n)*den(n)*dz(n)
!!c        flux=flux+a(n)*dz(n)
        yisum=yisum+dz(n)
     end do

     ! calculate fractional flux
     s1=zz(kdep)-zdep
     s2=s1/dz(ka)
     if(yisum+s1>0.0_PS) then
        flux=(flux+&
!!c             (aLR(ka,2)-0.5*s2*(aLR(ka,2)-aLR(ka,1)-(1.0-2.0*s2/3.0)*a6(ka)))*s1&
             den(ka)*(aLR(ka,2)-0.5_PS*s2*(aLR(ka,2)-aLR(ka,1)-(1.0_PS-2.0_PS*s2/3.0_PS)*a6(ka)))*s1&
             )/(yisum+s1)
     end if

     if ( debug ) then
        if(flux>1.0e+10_PS.or.flux<-1.0e+10_PS) then
           write(*,*) from
           write(*,'("flux large in slppm 1",2I5,20ES15.6)') kdep,ka,flux,w,den(ka),aLR(ka,1),aLR(ka,2),a6(ka),s1,yisum,dz(ka),s2,zdep
        end if
     end if

  elseif(w<0.0_PS) then
     ka=kdep
     ! calculate integer flux
     yisum=0.0_PS
!!c     if(kck==1) then
!!c        do n=2,ka-1
!!c           flux=flux+a(n)*den(n)*dvol(n)
!!c           yisum=yisum+dvol(n)
!!c        end do
!!c        do n=i+1,n1-1
!!c           flux=flux+a(n)*den(n)*dvol(n)
!!c           yisum=yisum+dz(n)
!!c        end do
!!c     else
!!c        do n=i+1,ka-1
!!c           if(k<=n1) then
!!c              flux=flux+a(n)*den(n)*dz(n)
!!c              yisum=yisum+dz(n)
!!c           else
!!c              ! periodic condition
!!c              flux=flux+a(k-n1+1)*den(k-n1+1)*dz(k-n1+1)
!!c              yisum=yisum+dz(k-n1+1)
!!c           end if
!!c        end do
!!c     end if
     do n=k+1,ka-1
        flux=flux+a(n)*den(n)*dz(n)
!!c        flux=flux+a(n)*dz(n)
        yisum=yisum+dz(n)
     end do
     ! calculate fractional flux
     if(kdep-1>=1) then
        s1=zdep-zz(kdep-1)
     else 
        LOG_ERROR("cal_flux_z",*) "something wrong in vdy_scladv_slppm flux_z",kdep
        call PRC_abort
     end if
     s2=s1/dz(ka)
     if(yisum+s1>0.0) then
        flux=(flux+&
!!c             (aLR(ka,1)+0.5*s2*(aLR(ka,2)-aLR(ka,1)+(1.0-2.0*s2/3.0)*a6(ka)))*s1&
             den(ka)*(aLR(ka,1)+0.5*s2*(aLR(ka,2)-aLR(ka,1)+(1.0-2.0*s2/3.0)*a6(ka)))*s1&
             )/(yisum+s1)
     end if
     if ( debug ) then
        if(flux>1.0e+10.or.flux<-1.0e+10) then
           write(*,*) from
           write(*,'("flux large in slppm 2",2I5,20ES15.6)') kdep,ka,flux,w,den(ka),aLR(ka,1),aLR(ka,2),a6(ka),s1,yisum,dz(ka),s2,zdep
        end if
     end if
  end if
  flux=flux*w

  if ( debug ) then
     if(flux>1.0e+10.or.flux<-1.0e+10) then
        write(*,*) from
        write(*,'("flux large in slppm 3",2I5,20ES15.6)') kdep,ka,flux,w,den(ka),aLR(ka,1),aLR(ka,2),a6(ka),s1,yisum,dz(ka),s2,zdep
     end if
  end if
  
end subroutine cal_flux_z
subroutine cal_flux_y(j,ns,n3,lbin,yy,dy,dt,den,dvol,dzv,coslv,a,aLR,a6,yarr,v,flux)
  use scale_prc, only: &
     PRC_abort
  implicit none
  integer,intent(in) :: j,ns,n3,lbin
  real(PS) :: a(lbin),a6(n3),yy(lbin),dzv(lbin),den(lbin),dvol(n3),coslv(lbin)
  real(PS),dimension(n3,2) :: aLR
  real(PS),intent(in)    :: v,yarr,dy,dt
  real(PS),intent(inout) :: flux
  real(PS) :: s1,s2,yisum,ydep
  integer :: jdep,jck
  integer :: n,ja
  
!!c  open(1,file='check_yy.dat') 
!!c  do n=1,n3
!!c     write(1,*) yy(n)
!!c  end do
!!c  close(1)
!!c  stop
  ! find departure point
  ! search for departure point
  jck=0
  ydep=yarr-v*dt

  ! find the interface departure point
  jdep=int((ydep-yy(1))/dy)+2
  jdep=max(1,min(n3,jdep))    
  
  ja=jdep
  
  flux=0.0
  if(v>0.0) then
     ! ja is the cell number for which the interpolation function is used.
     ja=jdep
     ! calculate integer flux
     yisum=0.0

     if(jck==0) then
        do n=j,ja+1,-1
           if(dzv(n)>0.0) then
              flux=flux+a(n)*den(n)*dvol(n)
              yisum=yisum+dy
!!c              yisum=yisum+dvol(n)
           else
              ! case for hitting barrir
              flux=flux/yisum
              goto 100
           end if
        end do
     elseif(jck==1) then
        do n=j,2,-1
           if(dzv(n-1)>0.0) then
              flux=flux+a(n)*den(n)*dvol(n)
              yisum=yisum+dy
!!c              yisum=yisum+dvol(n)
           else
              ! case for hitting barrir
              flux=flux/yisum
              goto 100
           end if
        end do
        do n=n3-1,ja+1,-1
           if(dzv(n-1)>0.0) then
              flux=flux+a(n)*den(n)*dvol(n)
              yisum=yisum+dy
!!c              yisum=yisum+dvol(n)
           else
              ! case for hitting barrir
              flux=flux/yisum
              goto 100
           end if
        end do
     end if

     ! calculate fractional flux
     s1=yy(jdep)-ydep
     s2=s1/dy
     if(yisum+s1>0.0) then
!!c        flux=flux+&
!!c             den(ja)*dzv(ja)*coslv(ja)*(aLR(ja,2)-0.5*s2*(aLR(ja,2)-aLR(ja,1)-(1.0-2.0*s2/3.0)*a6(ja)))
        flux=(flux+&
!!c             dzv(ja)*(aLR(ja,2)-0.5*s2*(aLR(ja,2)-aLR(ja,1)-(1.0-2.0*s2/3.0)*a6(ja)))*s1&
             den(ja)*dzv(ja)*coslv(ja)*(aLR(ja,2)-0.5*s2*(aLR(ja,2)-aLR(ja,1)-(1.0-2.0*s2/3.0)*a6(ja)))*s1&
             )/(yisum+s1)
     end if
  elseif(v<0.0) then
     ja=jdep
     ! calculate integer flux
     yisum=0.0
     if(jck==0) then
        do n=j+1,ja-1
           if(dzv(n-1)>0.0) then
              flux=flux+a(n)*den(n)*dvol(n)
              yisum=yisum+dy
!!c              yisum=yisum+dvol(n)
           else
              ! case for hitting barrir
              flux=flux/yisum
              goto 100
           end if
        end do
     elseif(jck==1) then
        do n=j+1,n3-1
           if(dzv(n-1)>0.0) then
              flux=flux+a(n)*den(n)*dvol(n)
              yisum=yisum+dy
!!c              yisum=yisum+dvol(n)
           else
              ! case for hitting barrir
              flux=flux/yisum
              goto 100
           end if
        end do
        do n=2,ja-1
           if(dzv(n-1)>0.0) then
              flux=flux+a(n)*den(n)*dvol(n)
              yisum=yisum+dy
!!c              yisum=yisum+dvol(n)
           else
              ! case for hitting barrir
              flux=flux/yisum
              goto 100
           end if
        end do
     end if
     ! calculate fractional flux
     if(jdep-1>=1) then
        s1=ydep-yy(jdep-1)
     else 
        LOG_ERROR("cal_flux_y",*) "something wrong in vdy_scladv_slppm flux_y",jdep
        call PRC_abort
     end if
     s2=s1/dy
     if(yisum+s1>0.0) then
!!c        flux=flux+&
!!c             den(ja)*dzv(ja-1)*coslv(ja-1)*(aLR(ja,1)+0.5*s2*(aLR(ja,2)-aLR(ja,1)+(1.0-2.0*s2/3.0)*a6(ja)))
        flux=(flux+&
             den(ja)*dzv(ja-1)*coslv(ja-1)*(aLR(ja,1)+0.5*s2*(aLR(ja,2)-aLR(ja,1)+(1.0-2.0*s2/3.0)*a6(ja)))*s1&
!!c             dzv(ja-1)*(aLR(ja,1)+0.5*s2*(aLR(ja,2)-aLR(ja,1)+(1.0-2.0*s2/3.0)*a6(ja)))*s1&
             )/(yisum+s1)
     end if
  end if
100 continue
  flux=flux*v
  
end subroutine cal_flux_y

subroutine cal_flux_x(i,ns,n2,lbin,xx,dx,dt,den,dvol,dzu,coslu,a,aLR,a6,xarr,u,flux)
  use scale_prc, only: &
     PRC_abort
  implicit none
  integer,intent(in) :: i,ns,n2,lbin
  real(PS) :: a(lbin),a6(n2),xx(lbin),dzu(lbin),den(lbin),dvol(n2),coslu
  real(PS),dimension(n2,2) :: aLR
  real(PS),intent(in)    :: u,xarr,dx,dt
  real(PS),intent(inout) :: flux
  real(PS) :: s2,s1,xisum,xdep
  integer :: idep,ick
  integer :: n,ia
  
!!c  open(1,file='check_yy.dat') 
!!c  do n=1,n3
!!c     write(1,*) y(n)
!!c  end do
!!c  close(1)
!!c  stop
  ! find departure point
  ! search for departure point
  ick=0
!!c7a  xdep=xarr-u*dt

  ! expand u in the normal dimension of x
  ! because other coefs were calculated with dx sizing
  xdep=xarr-u/coslu*dt

  ! find the interface departure point
!!c7a  idep=int((xdep-xx(1))/(dx*coslu))+2
  idep=int((xdep-xx(1))/dx)+2
  idep=max(1,min(n2,idep))    
  
  ia=idep
  
  flux=0.0
  if(u>0.0) then
     ! ia is the cell number for which the interpolation function is used.
     ia=idep
     ! calculate integer flux
     xisum=0.0

     if(ick==0) then
        do n=i,ia+1,-1
           if(dzu(n)>0.0) then
              flux=flux+a(n)*den(n)*dvol(n)
              xisum=xisum+dx
!!c7a              xisum=xisum+dx*coslu
!!c              xisum=xisum+dvol(n)
           else
              ! case for hitting barrir
              flux=flux/xisum
              goto 100
           end if
        end do
     elseif(ick==1) then
        do n=i,2,-1
           if(dzu(n-1)>0.0) then
              flux=flux+a(n)*den(n)*dvol(n)
              xisum=xisum+dx
!!c              xisum=xisum+dvol(n)
!!c7a              xisum=xisum+dx*coslu
           else
              ! case for hitting barrir
              flux=flux/xisum
              goto 100
           end if
        end do
        do n=n2-1,ia+1,-1
           if(dzu(n-1)>0.0) then
              flux=flux+a(n)*den(n)*dvol(n)
              xisum=xisum+dx
!!c7a              xisum=xisum+dx*coslu
!!c              xisum=xisum+dvol(n)
           else
              ! case for hitting barrir
              flux=flux/xisum
              goto 100
           end if
        end do
     end if

     ! calculate fractional flux
     s1=xx(idep)-xdep
     s2=s1/dx
!!c7a     s2=s1/(dx*coslu)
     if(xisum+s1>0.0) then
!!c        flux=flux+&
!!c             den(ia)*dzu(ia)*(aLR(ia,2)-0.5*s2*(aLR(ia,2)-aLR(ia,1)-(1.0-2.0*s2/3.0)*a6(ia)))
        flux=(flux+&
!!c             dzu(ia)*(aLR(ia,2)-0.5*s2*(aLR(ia,2)-aLR(ia,1)-(1.0-2.0*s2/3.0)*a6(ia)))*s1&
             den(ia)*dzu(ia)*(aLR(ia,2)-0.5*s2*(aLR(ia,2)-aLR(ia,1)-(1.0-2.0*s2/3.0)*a6(ia)))*s1&
             )/(xisum+s1)
     end if
  elseif(u<0.0) then
     ia=idep
     ! calculate integer flux
     xisum=0.0
     if(ick==0) then
        do n=i+1,ia-1
           if(dzu(n-1)>0.0) then
              flux=flux+a(n)*den(n)*dvol(n)
              xisum=xisum+dx
!!c              xisum=xisum+dvol(n)
!!c7a              xisum=xisum+dx*coslu
           else
              ! case for hitting barrir
              flux=flux/xisum
              goto 100
           end if
        end do
     elseif(ick==1) then
        do n=i+1,n2-1
           if(dzu(n-1)>0.0) then
              flux=flux+a(n)*den(n)*dvol(n)
              xisum=xisum+dx
!!c7a              xisum=xisum+dx*coslu
!!c              xisum=xisum+dvol(n)
           else
              ! case for hitting barrir
              flux=flux/xisum
              goto 100
           end if
        end do
        do n=2,ia-1
           if(dzu(n-1)>0.0) then
              flux=flux+a(n)*den(n)*dvol(n)
              xisum=xisum+dx
!!c              xisum=xisum+dvol(n)
!!c7a              xisum=xisum+dx*coslu
           else
              ! case for hitting barrir
              flux=flux/xisum
              goto 100
           end if
        end do
     end if
     ! calculate fractional flux
     if(idep-1>=1) then
        s1=xdep-xx(idep-1)
     else
        LOG_ERROR("cal_flux_x",*) "something wrong in vdy_scladv_slppm flux_x",idep
        call PRC_abort
     end if
     s2=s1/dx
!!c7a     s2=s1/(dx*coslu)
     if(xisum+s1>0.0) then
!!c        flux=flux+&
!!c             den(ia)*dzu(ia-1)*(aLR(ia,1)+0.5*s2*(aLR(ia,2)-aLR(ia,1)+(1.0-2.0*s2/3.0)*a6(ia)))
        flux=(flux+&
             den(ia)*dzu(ia-1)*(aLR(ia,1)+0.5*s2*(aLR(ia,2)-aLR(ia,1)+(1.0-2.0*s2/3.0)*a6(ia)))*s1&
!!c             dzu(ia-1)*(aLR(ia,1)+0.5*s2*(aLR(ia,2)-aLR(ia,1)+(1.0-2.0*s2/3.0)*a6(ia)))*s1&
             )/(xisum+s1)
     end if
  end if
100 continue
  flux=flux*u
  
end subroutine cal_flux_x

end module mod_scladv_slppm
