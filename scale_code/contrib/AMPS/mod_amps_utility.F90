#include "scalelib.h"
!OCL SERIAL
module mod_amps_utility

  use scale_io
  use acc_amps
  use com_amps
  implicit none
  private

  public :: RDAPTB
  public :: RDCETB
  public :: RDSTTB
  public :: init_inherent_growth_par
  public :: init_osmo_par
  public :: init_normal_lut
  public :: init_inv_normal_lut
  public :: read_AMPSTASK
  public :: read_seed
  public :: assign_seeds
  public :: config_varindex
  public :: get_growth_mode
  public :: get_growth_mode_max
  public :: get_inact
  public :: get_inact_tropic
  public :: get_inact_tropic2
  public :: get_len_c2a
  public :: get_len_s1
  public :: get_len_s3
  public :: cal_growth_mode_inl_vec
  public :: cal_growth_mode_hex_inl_vec
  public :: cal_mmass
  public :: cal_mmass_scale
  public :: cal_linprms_vec_s
  public :: cal_lincubprms_vec
  public :: cal_transbin_vec
  public :: findcond1
  public :: findcond2
  public :: sclsedprz_original
  public :: sclsedaer
  public :: moistthermo2
  public :: diabatic
  public :: negadj
  public :: negadj_bin2
  public :: negadj_bin_ap
  public :: integ_mictend1
  public :: integ_mictend2
  public :: shift_bin_vec
  public :: fill_bin_ap
  public :: fill_bulk_ap
  public :: osm_ammsul
  public :: osm_sodchl
  public :: QSPARM2
  public :: random_genvar
  public :: srand2
  public :: getznorm2
  public :: cdfnor
  public :: acos_m
  public :: fGM
  public :: zbrent
  public :: get_cmod_inh

  ! variables necessary for sequential random number generation
  type random_genvar
    sequence
    integer :: JSEED,IFRST,isect_seed,NEXTN
  end type random_genvar

  real(RP), private, parameter :: r = 287.04_RP
  real(RP), private, parameter :: rm = 461.50_RP
  real(RP), private, parameter :: cpd = 1004.64_RP
  real(RP), private, parameter :: alvl = 2.501E+6_RP
  real(RP), private, parameter :: alvi = 2.834E+6_RP
  real(RP), private, parameter :: g = 9.80665_RP
  real(RP), private, parameter :: epsvap = R / rm

contains

!!$      subroutine set_seedsize
!!$      use com_amps
!!$      implicit none
!!$
!!$         call random_seed(size=sizeseed)
!!$         write(fid_alog,*) "size of seeds for random generator",sizeseed
!!$      return
!!$      end subroutine set_seedsize

!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     read the lookup tables
!
!     - inputs -
!     DRCETB ...... directory where lookup talbes are stored.
!                   e.g., /home/nms/tempei/exp1/looktables
!
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine RDCETB(nbr,ncr,nbi,nci)
      use com_amps
      implicit none
      integer, intent(in) :: nbr,nbi,ncr,nci
!      character*100   DRCETB
      character(len=110) :: ifname
      character(len=300) :: header
      real(PS) :: dum1,dum2,dum3
      integer :: i,j,k,l,nrow,ncol,npar,nt
      integer :: io_rdcetb
!     +++++++++++ read drop_drop.dat +++++++++++++++++++++++++++++
      io_rdcetb = IO_get_available_fid()
      open(io_rdcetb, file=trim(DRCETB)//"/drop_drop_Rey4.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_drpdrp,nc_drpdrp,xs_drpdrp,dx_drpdrp &
           ,ys_drpdrp,dy_drpdrp
      do i = 1, nr_drpdrp
         read(io_rdcetb,*) (drpdrp(i,j), j = 1, nc_drpdrp)
      end do
      close(io_rdcetb)

!     +++++++++++ read hex_drop_Nre_Ec.dat +++++++++++++++++++++++++++++
      open(io_rdcetb, file=trim(DRCETB)//"/hex_drop_Nre_Ec.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_hexdrp,nc_hexdrp,xs_hexdrp,dx_hexdrp &
           ,ys_hexdrp,dy_hexdrp
      do i = 1, nr_hexdrp
         read(io_rdcetb,*) (hexdrp(i,j), j = 1, nc_hexdrp)
      end do
      close(io_rdcetb)

!     +++++++++++ read bbc_drop_Nre_Ec.dat +++++++++++++++++++++++++++++
      open(io_rdcetb, file=trim(DRCETB)//"/bbc_drop_Nre_Ec.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_bbcdrp,nc_bbcdrp,xs_bbcdrp,dx_bbcdrp &
           ,ys_bbcdrp,dy_bbcdrp
      do i = 1, nr_bbcdrp
         read(io_rdcetb,*) (bbcdrp(i,j), j = 1, nc_bbcdrp)
      end do
      close(io_rdcetb)

!     +++++++++++ read col_drop_Nre_Ec.dat +++++++++++++++++++++++++++++
      open(io_rdcetb, file=trim(DRCETB)//"/col_drop_Nre_Ec.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_coldrp,nc_coldrp,xs_coldrp,dx_coldrp &
           ,ys_coldrp,dy_coldrp
      do i = 1, nr_coldrp
         read(io_rdcetb,*) (coldrp(i,j), j = 1, nc_coldrp)
      end do
      close(io_rdcetb)

!     +++++++++++ read grp01_ratNre_Ec.dat +++++++++++++++++++++++++++++
!ccc      ifname=DRCETB(1:LEN(DRCETB))//"/grap01_drop.dat"
      open(io_rdcetb, file=trim(DRCETB)//"/grp01_ratNre_Ec.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_gp1drp,nc_gp1drp,xs_gp1drp,dx_gp1drp &
           ,ys_gp1drp,dy_gp1drp
      do i = 1, nr_gp1drp
         read(io_rdcetb,*) (gp1drp(i,j), j = 1, nc_gp1drp)
      end do
      close(io_rdcetb)

!     +++++++++++ read grp04_ratNre_Ec.dat +++++++++++++++++++++++++++++
!ccc      ifname=DRCETB(1:LEN(DRCETB))//"/grap04_drop.dat"
      open(io_rdcetb, file=trim(DRCETB)//"/grp04_ratNre_Ec.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_gp4drp,nc_gp4drp,xs_gp4drp,dx_gp4drp &
           ,ys_gp4drp,dy_gp4drp
      do i = 1, nr_gp4drp
         read(io_rdcetb,*) (gp4drp(i,j), j = 1, nc_gp4drp)
      end do
      close(io_rdcetb)

!     +++++++++++ read grp08_ratNre_Ec.dat +++++++++++++++++++++++++++++
!ccc      ifname=DRCETB(1:LEN(DRCETB))//"/grap08_drop.dat"
      open(io_rdcetb, file=trim(DRCETB)//"/grp08_ratNre_Ec.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nr_gp8drp,nc_gp8drp,xs_gp8drp,dx_gp8drp &
           ,ys_gp8drp,dy_gp8drp
      do i = 1, nr_gp8drp
         read(io_rdcetb,*) (gp8drp(i,j), j = 1, nc_gp8drp)
      end do
      close(io_rdcetb)


!     +++++++++++++++++ read frequency of habit (T<-20) ++++++
      open(io_rdcetb, file=trim(DRCETB)//"/pol_frq.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (pol_frq(i,j), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/pla_frq.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (pla_frq(i,j), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/col_frq.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (col_frq(i,j), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/ros_frq.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (ros_frq(i,j), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/ppo_frq.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (ppo_frq(i,j), j = 1, ncol)
      end do
      close(io_rdcetb)

!     +++ normalize the frequency +++
      do i=1,nrow
         do j=1,ncol
            pol_frq(i,j)=max(0.0_RP,pol_frq(i,j))
            pla_frq(i,j)=max(0.0_RP,pla_frq(i,j))
            col_frq(i,j)=max(0.0_RP,col_frq(i,j))
            ros_frq(i,j)=max(0.0_RP,ros_frq(i,j))
            ppo_frq(i,j)=max(0.0_RP,ppo_frq(i,j))

            dum1=pol_frq(i,j)+pla_frq(i,j)+col_frq(i,j)
            pol_frq(i,j)=pol_frq(i,j)/dum1
            pla_frq(i,j)=pla_frq(i,j)/dum1
            col_frq(i,j)=col_frq(i,j)/dum1

!            dum2=pol_frq(i,j)+pla_frq(i,j)+col_frq(i,j)
!            write(fid_alog,*) "i,j,bf sum of frq, af",i,j,dum1,dum2
         end do
      end do

!     +++++++++++++++++ read map data for diagnose a and c ++++++
      open(io_rdcetb, file=trim(DRCETB)//"/tmp_map_col.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (mtac_map_col(i,j,1), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/tmd_map_col.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (mtac_map_col(i,j,2), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/lmt_mass_col.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) dum1,lmt_mass_col(i),dum2,dum3
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/tmp_map_pla.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (mtac_map_pla(i,j,1), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/tmd_map_pla.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) (mtac_map_pla(i,j,2), j = 1, ncol)
      end do
      close(io_rdcetb)
      open(io_rdcetb, file=trim(DRCETB)//"/lmt_mass_pla.dat")
      read(io_rdcetb,*) header
      read(io_rdcetb,*) nrow, ncol
      do i = 1, nrow
         read(io_rdcetb,*) dum1,lmt_mass_pla(i),dum2,dum3
      end do
      close(io_rdcetb)
      call cal_indx_lmtmass()

      if ( IsMaster .and. debug ) then
        write(fid_alog,*) "bottom of RDCETB"
      end if

      return
      end subroutine RDCETB

      subroutine RDAPTB()
      use com_amps
      implicit none
      character(len=300) :: header
      integer :: i, j, k,l,id
      integer :: io_rdaptb
      return
!     +++++++++++ read ap_act.bin.dat +++++++++++++++++++++++++++++
      io_rdaptb = IO_get_available_fid()
      open(io_rdaptb, file=trim(DRAPTB)//"/ap_act.bin.dat",form='unformatted')
      read(io_rdaptb) (apt_st(i),i=1,4),(apt_add(i),i=1,4),(apt_max(i),i=1,4) &
           ,(napt(i),i=1,4)
      do i=1,napt(1)
         do j=1,napt(2)
            do k=1,napt(3)
               do l=1,napt(4)
                  id=(j-1)*napt(3)*napt(4)+(k-1)*napt(4)+l
                  read(io_rdaptb) frac_apact(1,id),frac_apact(2,id)
               end do
            end do
         end do
         if ( IsMaster .and. debug ) then
           if(i==1+nint((ap_lnsig(1)-apt_st(1))/apt_add(1))) then
              write(fid_alog,*) "ap_lnsig(1), and chosen sigma:",ap_lnsig(1), &
                   apt_st(1)+real(i-1,PS_KIND)*apt_add(1)
           end if
           if(i==1+nint((ap_lnsig(2)-apt_st(1))/apt_add(1))) then
              write(fid_alog,*) "ap_lnsig(2), and chosen sigma:",ap_lnsig(2), &
                   apt_st(1)+real(i-1,PS_KIND)*apt_add(1)
           end if
         end if
      end do
      close(io_rdaptb)
      return
      end subroutine RDAPTB


      subroutine RDSTTB()
      use com_amps
      implicit none
      real(PS) :: dum1,dum2
      integer i,j
      integer :: io_rdsttb
!     +++++++++++ read standard normal distribution dat +++++++++++++++++++++++++++++
      io_rdsttb = IO_get_available_fid()
      open(io_rdsttb, file=trim(DRSTTB)//"/stdnorm.dat",form='formatted')
      do i=1,451
         read(io_rdsttb,*) dum1,dum2,(znorm(j,i),j=1,4)
      end do
      close(io_rdsttb)

      return
      end subroutine RDSTTB


!!$      real(PS) function getznorm(which,x)
!!$!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$!     return      x
!!$!                / z**(which-1) * exp(-z**2.0/2.0) dz
!!$!                -inf
!!$!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$      implicit none
!!$      real(PS) :: znorm
!!$      common/STTBL/znorm(4,451)
!!$      real(PS) :: x
!!$      integer :: which,ix
!!$      if(x>=0.0_RP) then
!!$         getznorm=0.5_RP+znorm(which,min(max(1,nint(x*100.0_RP)+1),451))
!!$      else
!!$         getznorm=0.5_RP-znorm(which,min(max(1,nint(-x*100.0_RP)+1),451))
!!$      end if
!!$      end function getznorm

      function getznorm2(x) result(p)
        real(DS) :: x,mean,sd,q,bound,p!,out
        integer :: status
        mean=0.0_DS
        sd=1.0_DS
        call cdfnor (1, p, q, x, mean, sd, status, bound )
      end function getznorm2

!!$      real(PS) function getdrpdrp(NreL,rrat)
!!$      use com_amps
!!$      implicit none
!!$      real(PS) :: NreL,rrat,NreL_p,rrat_p
!!$      integer :: i1,j1
!!$      real(PS) :: x1,y1,wx,wy
!!$      real(PS) :: ec_min
!!$      parameter (ec_min=0.0_RP)
!!$
!!$      NreL_p=max(NreL,1.0e-10_RP)
!!$      if(rrat<0.0_RP) then
!!$        write(fid_alog,*) "something is not right getdrpdrp",rrat,NreL_p
!!$        stop
!!$      endif
!!$      if(rrat>1.0_RP) then
!!$        rrat_p=1.0_RP/rrat
!!$      else
!!$        rrat_p=rrat
!!$      endif
!!$
!!$!      write(fid_alog,'("ck drpdrp",10ES15.6)') rrat,NreL
!!$!     *           ,rrat_p,NreL_p,xs_drpdrp
!!$!     *           ,dx_drpdrp
!!$!     *           ,ys_drpdrp,dy_drpdrp
!!$
!!$      j1=max(1,min(nc_drpdrp-1 &
!!$           ,int((rrat_p-xs_drpdrp)/dx_drpdrp)+1))
!!$      x1=real(j1-1)*dx_drpdrp+xs_drpdrp
!!$
!!$      i1=max(1,min(nr_drpdrp-1 &
!!$           ,int((log10(NreL_p)-ys_drpdrp)/dy_drpdrp)+1))
!!$      y1=real(i1-1)*dy_drpdrp+ys_drpdrp
!!$
!!$      wx=max(0.0_RP,min(1.0_RP,(rrat_p-x1)/dx_drpdrp))
!!$      wy=max(0.0_RP,min(1.0_RP,(log10(NreL_p)-y1)/dy_drpdrp))
!!$
!!$      getdrpdrp=min(1.0_RP,max(ec_min, &
!!$           (1.0_RP-wx)*(1.0_RP-wy)*drpdrp(i1,j1)+ &
!!$           (1.0_RP-wx)*wy*drpdrp(i1,j1+1)+ &
!!$           wx*(1.0-wy)*drpdrp(i1+1,j1)+ &
!!$           wx*wy*drpdrp(i1+1,j1+1)))
!!$
!!$      end function getdrpdrp

!!$      real(PS) function getgpxdrp(g_den,NreL,rrat)
!!$      use com_amps
!!$      implicit none
!!$      real(PS) :: g_den,NreL,rrat,NreL_p
!!$      integer :: i1,j1
!!$      real(PS) :: x1,y1,wx,wy
!!$      real(PS) :: ec_min
!!$      parameter (ec_min=0.0_RP)
!!$
!!$!ccc      if(NreL<=0.0_RP) then
!!$!ccc         write(fid_alog,*) "NreL num is neg or 0",NreL
!!$!ccc      end if
!!$      NreL_p=max(NreL,1.0e-10_RP)
!!$
!!$
!!$
!!$      if(g_den.le.0.2_RP) then
!!$         j1=max(1,min(nc_gp1drp-1 &
!!$              ,int((rrat-xs_gp1drp)/dx_gp1drp)+1))
!!$         x1=real(j1-1,PS_KIND)*dx_gp1drp+xs_gp1drp
!!$
!!$         i1=max(1,min(nr_gp1drp-1 &
!!$              ,int((log10(NreL_p)-ys_gp1drp)/dy_gp1drp)+1))
!!$         y1=real(i1-1,PS_KIND)*dy_gp1drp+ys_gp1drp
!!$
!!$!     do not extrapolate
!!$         wx=min(1.0_RP,max(0.0_RP,(rrat-x1)/dx_gp1drp))
!!$         wy=min(1.0_RP,max(0.0_RP,(log10(NreL_p)-y1)/dy_gp1drp))
!!$
!!$         getgpxdrp=min(1.0_RP,max(ec_min, &
!!$              (1.0_RP-wx)*(1.0_RP-wy)*gp1drp(i1,j1)+ &
!!$              (1.0_RP-wx)*wy*gp1drp(i1,j1+1)+ &
!!$              wx*(1.0_RP-wy)*gp1drp(i1+1,j1)+ &
!!$              wx*wy*gp1drp(i1+1,j1+1)))
!!$
!!$      else if( 0.2_RP .le. g_den .and. g_den .le. 0.6_RP ) then
!!$         j1=max(1,min(nc_gp4drp-1 &
!!$              ,int((rrat-xs_gp4drp)/dx_gp4drp)+1))
!!$         x1=real(j1-1,PS_KIND)*dx_gp4drp+xs_gp4drp
!!$
!!$         i1=max(1,min(nr_gp4drp-1 &
!!$              ,int((log10(NreL_p)-ys_gp4drp)/dy_gp4drp)+1))
!!$         y1=real(i1-1,PS_KIND)*dy_gp4drp+ys_gp4drp
!!$
!!$!     do not extrapolate
!!$         wx=min(1.0_RP,max(0.0_RP,(rrat-x1)/dx_gp4drp))
!!$         wy=min(1.0_RP,max(0.0_RP,(log10(NreL_p)-y1)/dy_gp4drp))
!!$
!!$         getgpxdrp=min(1.0_RP,max(ec_min, &
!!$              (1.0_RP-wx)*(1.0_RP-wy)*gp4drp(i1,j1)+ &
!!$              (1.0_RP-wx)*wy*gp4drp(i1,j1+1)+ &
!!$              wx*(1.0_RP-wy)*gp4drp(i1+1,j1)+ &
!!$              wx*wy*gp4drp(i1+1,j1+1)))
!!$
!!$      else if( 0.6_RP .le. g_den ) then
!!$         j1=max(1,min(nc_gp8drp-1 &
!!$              ,int((rrat-xs_gp8drp)/dx_gp8drp)+1))
!!$         x1=real(j1-1,PS_KIND)*dx_gp8drp+xs_gp8drp
!!$
!!$         i1=max(1,min(nr_gp8drp-1 &
!!$              ,int((log10(NreL_p)-ys_gp8drp)/dy_gp8drp)+1))
!!$         y1=real(i1-1,PS_KIND)*dy_gp8drp+ys_gp8drp
!!$
!!$!     do not extrapolate
!!$         wx=min(1.0_RP,max(0.0_RP,(rrat-x1)/dx_gp8drp))
!!$         wy=min(1.0_RP,max(0.0_RP,(log10(NreL_p)-y1)/dy_gp8drp))
!!$
!!$         getgpxdrp=min(1.0_RP,max(ec_min, &
!!$              (1.0_RP-wx)*(1.0_RP-wy)*gp8drp(i1,j1)+ &
!!$              (1.0_RP-wx)*wy*gp8drp(i1,j1+1)+ &
!!$              wx*(1.0_RP-wy)*gp8drp(i1+1,j1)+ &
!!$              wx*wy*gp8drp(i1+1,j1+1)))
!!$
!!$
!!$      end if
!!$      end function getgpxdrp

!!$      real(PS) function getcrydrp(habit,Nre_r,Nre_s)
!!$      use com_amps
!!$      implicit none
!!$      real(PS) :: Nre_r,Nre_s,Nre_rp,Nre_sp
!!$      integer :: habit,i1,j1
!!$      real(PS) :: x1,y1,wx,wy
!!$      real(PS) :: ec_min
!!$      parameter (ec_min=0.0_RP)
!!$
!!$
!!$      Nre_rp=max(Nre_r,1.0e-10_RP)
!!$      Nre_sp=max(Nre_s,1.0e-10_RP)
!!$!ccc      if(Nre_r<=0.0_RP) then
!!$!ccc         write(fid_alog,*) "Nre_r num is neg or 0",Nre_r
!!$!ccc      end if
!!$!ccc      if(Nre_s<=0.0_RP) then
!!$!ccc         write(fid_alog,*) "Nre_s num is neg or 0",Nre_s
!!$!ccc      end if
!!$
!!$
!!$      if(habit.eq.1.or.habit.eq.5.or.habit.eq.6) then
!!$!     assume that irregular crystal is close to hexagonal plate
!!$         j1=max(1,min(nc_hexdrp-1 &
!!$              ,int((log10(Nre_rp)-xs_hexdrp)/dx_hexdrp)+1))
!!$         x1=real(j1-1,PS_KIND)*dx_hexdrp+xs_hexdrp
!!$
!!$         i1=max(1,min(nr_hexdrp-1 &
!!$              ,int((log10(Nre_sp)-ys_hexdrp)/dy_hexdrp)+1))
!!$         y1=real(i1-1,PS_KIND)*dy_hexdrp+ys_hexdrp
!!$
!!$!     do not extrapolate!
!!$         wx=min(1.0_RP,max(0.0_RP,(log10(Nre_rp)-x1)/dx_hexdrp))
!!$         wy=min(1.0_RP,max(0.0_RP,(log10(Nre_sp)-y1)/dy_hexdrp))
!!$
!!$         getcrydrp=min(1.0_RP,max(ec_min, &
!!$              (1.0_RP-wx)*(1.0_RP-wy)*hexdrp(i1,j1)+ &
!!$              (1.0_RP-wx)*wy*hexdrp(i1,j1+1)+ &
!!$              wx*(1.0_RP-wy)*hexdrp(i1+1,j1)+ &
!!$              wx*wy*hexdrp(i1+1,j1+1)))
!!$
!!$      else if(habit.eq.2.or.habit.eq.4) then
!!$!     assume that rosette bullets char is close to broad branced crystal.
!!$         j1=max(1,min(nc_bbcdrp-1 &
!!$              ,int((log10(Nre_rp)-xs_bbcdrp)/dx_bbcdrp)+1))
!!$         x1=real(j1-1,PS_KIND)*dx_bbcdrp+xs_bbcdrp
!!$
!!$         i1=max(1,min(nr_bbcdrp-1 &
!!$              ,int((log10(Nre_sp)-ys_bbcdrp)/dy_bbcdrp)+1))
!!$         y1=real(i1-1,PS_KIND)*dy_bbcdrp+ys_bbcdrp
!!$
!!$!     do not extrapolate!
!!$         wx=min(1.0_RP,max(0.0_RP,(log10(Nre_rp)-x1)/dx_bbcdrp))
!!$         wy=min(1.0_RP,max(0.0_RP,(log10(Nre_sp)-y1)/dy_bbcdrp))
!!$
!!$         getcrydrp=min(1.0_RP,max(ec_min, &
!!$              (1.0_RP-wx)*(1.0_RP-wy)*bbcdrp(i1,j1)+ &
!!$              (1.0_RP-wx)*wy*bbcdrp(i1,j1+1)+ &
!!$              wx*(1.0_RP-wy)*bbcdrp(i1+1,j1)+ &
!!$              wx*wy*bbcdrp(i1+1,j1+1)))
!!$
!!$      else if(habit.eq.3) then
!!$         j1=max(1,min(nc_coldrp-1 &
!!$              ,int((log10(Nre_rp)-xs_coldrp)/dx_coldrp)+1))
!!$         x1=real(j1-1,PS_KIND)*dx_coldrp+xs_coldrp
!!$
!!$         i1=max(1,min(nr_coldrp-1 &
!!$              ,int((log10(Nre_sp)-ys_coldrp)/dy_coldrp)+1))
!!$         y1=real(i1-1,PS_KIND)*dy_coldrp+ys_coldrp
!!$
!!$!     do not extrapolate!
!!$         wx=min(1.0_RP,max(0.0_RP,(log10(Nre_rp)-x1)/dx_coldrp))
!!$         wy=min(1.0_RP,max(0.0_RP,(log10(Nre_sp)-y1)/dy_coldrp))
!!$
!!$         getcrydrp=min(1.0_RP,max(ec_min, &
!!$              (1.0_RP-wx)*(1.0_RP-wy)*coldrp(i1,j1)+ &
!!$              (1.0_RP-wx)*wy*coldrp(i1,j1+1)+ &
!!$              wx*(1.0-wy)*coldrp(i1+1,j1)+ &
!!$              wx*wy*coldrp(i1+1,j1+1)))
!!$      else
!!$         write(fid_alog,*) "no efficiency is defined for the habit",habit
!!$      end if
!!$      end function getcrydrp

!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     read the axis ratio-diameter relation data
!
!     - inputs -
!     IFPDTB ...... filename of ratio-diameter relation data
!
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!ccc      subroutine RDPDTB(IFPDTB)
!ccc      character*100   IFPDTB
!ccc      character*110   ifname
!ccc      character*300   header
!ccc      real dum
!ccc      integer i
!ccc      common/PDTBL/npd,spd,epd,dpd,phidia(96)
!cccC     +++++++++++ read spd_*.dat +++++++++++++++++++++++++++++
!ccc      open(101, file=IFPDTB)
!ccc      read(101,*) npd, spd,epd,dpd
!ccc      do i = 1, npd
!ccc         read(101,*) dum, phidia(i)
!ccc      end do
!ccc      close(101)
!ccc      end

      integer function get_growth_mode(tmp,si)
      use com_amps
      implicit none
      real(PS) :: tmp,si,frq(3)
      real(PS) :: dtmp,sttmp,edtmp,tmpc
      parameter(dtmp=1.0_RP,sttmp=-70.0_RP,edtmp=-20.0_RP)
      real(PS) :: dsi,stsi,edsi
      parameter(dsi=0.01_RP,stsi=0.0_RP,edsi=1.0_RP)
      integer :: itmpmin,itmpmax
      parameter(itmpmin=1,itmpmax=51)
      integer :: isimin,isimax
      parameter(isimin=1,isimax=101)
      integer :: itmp,isi
      integer :: i_tmp_le0,i_tmp_gtm45,i_x_lepol,i_x_ltppo,i_x_ltros &
           ,i_x_lepolpla

      real(PS) :: x1,x2
!!!      real rand2
!
!!!      integer JSEED,IFRST,isect_seed,NEXTN
!!!      COMMON /IRNGPRIVA/JSEED,IFRST,isect_seed,NEXTN
!!!!$omp threadprivate(/IRNGPRIVA/)

!
!cc      integer seed_priv,isect_seed
!cc      common/IRNGPRIV/seed_priv(256),isect_seed
!ccc$omp threadprivate(/IRNGPRIV/)
!
      integer omp_get_thread_num,omp_in_parallel,omp_get_num_threads

!

      !tmpc=tmp-273.16
      !itmp=max(itmpmin,min(nint((tmpc-sttmp)/dtmp)+1,itmpmax))
      !isi=max(isimin,min(nint(abs(si-stsi)/dsi)+1,isimax))

!
!      write(fid_alog,'("from get_growth_mode",2I4,50I8)') isect_seed
!     *          ,omp_get_thread_num(),seed_priv(1:sizeseed)


!     roll a dice
!ccc$omp critical
!cc      call random_seed(put=seed_priv)
      !call random_number(x1)
      !call random_number(x2)
!cc      call random_seed(get=seed_priv)
!ccc$omp end critical
!!!      x1=RAND2()
!!!      x2=RAND2()
!tmp      write(fid_alog,'("bf get_growth_mode",4I5,3I15,ES15.6)') isect_seed
!tmp     *          ,omp_in_parallel()
!tmp     *          ,omp_get_num_threads(),omp_get_thread_num()
!tmp     *          ,JSEED,IFRST,NEXTN,x
!tmp      write(fid_alog,*) "roll 1:",x

      !i_tmp_le0=0.5_RP*(1.0_RP+sign(1.0_RP,273.16_RP-tmp))
      !i_tmp_gtm45=0.5_RP*(1.0_RP-sign(1.0_RP,-45.0_RP-tmpc))

      !i_x_lepol=0.5_RP*(1.0+sign(1.0_RP,pol_frq(itmp,isi)-x1))

      !i_x_ltppo=0.5_RP*(1.0_RP-sign(1.0_RP,x2-ppo_frq(itmp,isi)))
      !i_x_ltros=0.5_RP*(1.0_RP-sign(1.0_RP,x2-ros_frq(itmp,isi)))

      !i_x_lepolpla=0.5_RP*(1.0+sign(1.0_RP,max(pol_frq(itmp,isi),0.0_RP)+&
      !                               max(pla_frq(itmp,isi),0.0_RP)-x1))

      ! CHIARUI
      !get_growth_mode= i_x_lepol*( &
      !        i_x_ltppo*i_tmp_gtm45 * 5 &
      !    +(1-i_x_ltppo*i_tmp_gtm45)*i_x_ltros*(1-i_tmp_gtm45)* 4 &
      !    +(1-i_x_ltppo*i_tmp_gtm45)*(1-i_x_ltros*(1-i_tmp_gtm45))* 1 &
      !                           )+ &
      !              (1-i_x_lepol)*i_x_lepolpla*  2 + &
      !              (1-i_x_lepol)*(1-i_x_lepolpla)*  3

      tmpc=tmp-273.16
      itmp=max(itmpmin,min(nint((tmpc-sttmp)/dtmp)+1,itmpmax))
      isi=max(isimin,min(nint(abs(si-stsi)/dsi)+1,isimax))

      if(pol_frq(itmp,isi)== &
            max(pol_frq(itmp,isi),pla_frq(itmp,isi),col_frq(itmp,isi))) &
            then
         get_growth_mode=1
!       assume that planar and columnar polycrystals are mutually exclusive.
         if(tmpc>-45.0.and.&
            ppo_frq(itmp,isi)== &
            max(ppo_frq(itmp,isi),1.0-ppo_frq(itmp,isi)-ros_frq(itmp,isi))) then
            get_growth_mode=5
         elseif(tmpc<=-45.0.and.&
            ros_frq(itmp,isi)== &
            max(ros_frq(itmp,isi),1.0-ppo_frq(itmp,isi)-ros_frq(itmp,isi))) then
            ! assume that rosetta starts only below T<-45
            get_growth_mode=4
         end if
      elseif(pla_frq(itmp,isi)== &
            max(pla_frq(itmp,isi),col_frq(itmp,isi))) then
         get_growth_mode=2
      else
         get_growth_mode=3
      end if
    !get_growth_mode=3

!org      if(x<=pol_frq(itmp,isi)) then
!org         get_growth_mode=1
!ccc$omp critical
!cc         call random_seed(put=seed_priv)
!!!         call random_number(x)
!cc         call random_seed(get=seed_priv)
!ccc$omp end critical
!org          x=RAND2()
!tmp          write(fid_alog,*) "roll 2:",x
!       assume that planar and columnar polycrystals are mutually exclusive.
!org         if(x<ppo_frq(itmp,isi).and.tmpc>-45.0) then
!org            get_growth_mode=5
!            write(fid_alog,'("planar polycrystals nucleated",4ES15.6)')
!     *        x,ppo_frq(itmp,isi),tmpc,si
!org
!org         elseif(x<ros_frq(itmp,isi).and.tmpc<=-45.0) then
            ! assume that rosetta starts only below T<-45
!ccc            if(tmp-273.16>-35) then
!ccc               write(fid_alog,*) itmp,isi, ros_frq(itmp,isi)
!ccc            end if
!org            get_growth_mode=4
!org         endif
!org      else if(x<=max(pol_frq(itmp,isi),0.0_RP) &
!org                 +max(pla_frq(itmp,isi),0.0_RP)) then
!org         get_growth_mode=2
!org      else
!org         get_growth_mode=3
!org      end if

!tmp note don't forget to remove the following
!tmp      get_growth_mode=1

!      if(tmpc<-20.0_RP) then
!      write(fid_alog,'("Tc,si,x,mode,pol,pla,col",3ES15.6,I3,3ES15.6)')
!     *          tmpc,si,x
!     *         ,get_growth_mode
!     *         ,pol_frq(itmp,isi),pla_frq(itmp,isi),col_frq(itmp,isi)
!
!      endif

      return
      end function get_growth_mode

subroutine cal_growth_mode_inl_vec(igm,nbin,L,ihabit_gm_random &
                                  ,tmp,si,rdsd)
 use com_amps,only: pol_frq,ppo_frq,ros_frq,pla_frq,col_frq
 implicit none
 integer,intent(in) :: nbin,L,ihabit_gm_random
 integer,dimension(*) :: igm
 real(PS),intent(in),dimension(*) :: tmp,si
 type(random_genvar),intent(inout) :: rdsd
 !
 ! new space
 !
 ! random number
 real(PS),dimension(nbin*L)  :: rnum1,rnum2

 real(PS) :: dtmp,sttmp,edtmp,tmpc
 parameter(dtmp=1.0_PS,sttmp=-70.0_PS,edtmp=-20.0_PS)
 real(PS) :: dsi,stsi,edsi
 parameter(dsi=0.01_PS,stsi=0.0_PS,edsi=1.0_PS)
 integer :: itmpmin,itmpmax
 parameter(itmpmin=1,itmpmax=51)
 integer :: isimin,isimax
 parameter(isimin=1,isimax=101)
 integer :: itmp,isi
 integer :: i,n,in

 real(PS) :: x1,x2
!

  if(ihabit_gm_random.eq.1) then

!CDIR NOVECTOR
    do in=1,nbin*L
      call rand2_ty(rnum1(in),rdsd)
      call rand2_ty(rnum2(in),rdsd)
    enddo

!CDIR NODEP
    do in=1,nbin*L
      n=(in-1)/NBIN+1
      i=in-(n-1)*NBIN

      tmpc=tmp(n)-273.16_PS
      itmp=max(itmpmin,min(nint((tmpc-sttmp)/dtmp)+1,itmpmax))
      isi=max(isimin,min(nint(abs(si(n)-stsi)/dsi)+1,isimax))

!
!      write(fid_alog,'("from get_growth_mode",2I4,50I8)') isect_seed
!     *          ,omp_get_thread_num(),seed_priv(1:sizeseed)


!     roll a dice

!cc      x1=RAND2()
      x1 = rnum1(in)



!cc      x2=RAND2()
      x2 = rnum2(in)

      if(x1<=pol_frq(itmp,isi)) then
         igm(in)=1
!tmp          write(fid_alog,*) "roll 2:",x
!       assume that planar and columnar polycrystals are mutually exclusive.
         if(x2<ppo_frq(itmp,isi).and.tmpc>-45.0) then
            igm(in)=5
!            write(fid_alog,'("planar polycrystals nucleated",4ES15.6)')
!     *        x,ppo_frq(itmp,isi),tmpc,si

         elseif(x2<ros_frq(itmp,isi).and.tmpc<=-45.0) then
            ! assume that rosetta starts only below T<-45
!ccc            if(tmp-273.16>-35) then
!ccc               write(fid_alog,*) itmp,isi, ros_frq(itmp,isi)
!ccc            end if
            igm(in)=4
         endif
      else if(x1<=max(pol_frq(itmp,isi),0.0_RP) &
                 +max(pla_frq(itmp,isi),0.0_RP)) then
         igm(in)=2
      else
         igm(in)=3
      end if
    enddo

  else
!!c    write(fid_alog,*) "in max gm"

!CDIR NODEP
    do in=1,nbin*L
      n=(in-1)/NBIN+1
      i=in-(n-1)*NBIN

      tmpc=tmp(n)-273.16
      itmp=max(itmpmin,min(nint((tmpc-sttmp)/dtmp)+1,itmpmax))
      isi=max(isimin,min(nint(abs(si(n)-stsi)/dsi)+1,isimax))

      if(pol_frq(itmp,isi)== &
            max(pol_frq(itmp,isi),pla_frq(itmp,isi),col_frq(itmp,isi))) &
            then
         igm(in)=1
!       assume that planar and columnar polycrystals are mutually exclusive.
         if(tmpc>-45.0.and.&
            ppo_frq(itmp,isi)== &
            max(ppo_frq(itmp,isi),1.0-ppo_frq(itmp,isi)-ros_frq(itmp,isi))) then
            igm(in)=5
         elseif(tmpc<=-45.0.and.&
            ros_frq(itmp,isi)== &
            max(ros_frq(itmp,isi),1.0-ppo_frq(itmp,isi)-ros_frq(itmp,isi))) then
            ! assume that rosetta starts only below T<-45
            igm(in)=4
         end if
      elseif(pla_frq(itmp,isi)== &
            max(pla_frq(itmp,isi),col_frq(itmp,isi))) then
         igm(in)=2
      else
         igm(in)=3
      end if
    enddo

  endif


  return
end subroutine cal_growth_mode_inl_vec

!!$      function get_growth_mode_inl(tmp,si,rnum1,rnum2) result(igm)
!!$      use com_amps
!!$      implicit none
!!$      integer :: igm
!!$      real(PS),intent(in) :: tmp,si,rnum1,rnum2
!!$      ! new space
!!$      real(PS) :: dtmp,sttmp,edtmp,tmpc
!!$      parameter(dtmp=1.0_RP,sttmp=-70.0_RP,edtmp=-20.0_RP)
!!$      real(PS) :: dsi,stsi,edsi
!!$      parameter(dsi=0.01_RP,stsi=0.0_RP,edsi=1.0_RP)
!!$      integer :: itmpmin,itmpmax
!!$      parameter(itmpmin=1,itmpmax=51)
!!$      integer :: isimin,isimax
!!$      parameter(isimin=1,isimax=101)
!!$      integer :: itmp,isi
!!$      integer :: i_tmp_le0,i_tmp_gtm45,i_x_lepol,i_x_ltppo,i_x_ltros &
!!$           ,i_x_lepolpla
!!$
!!$      real(PS) :: x1,x2
!!$!!!      real rand2
!!$!
!!$!
!!$!
!!$!
!!$      integer omp_get_thread_num,omp_in_parallel,omp_get_num_threads
!!$
!!$!
!!$
!!$      tmpc=tmp-273.16
!!$      itmp=max(itmpmin,min(nint((tmpc-sttmp)/dtmp)+1,itmpmax))
!!$      isi=max(isimin,min(nint(abs(si-stsi)/dsi)+1,isimax))
!!$
!!$!
!!$!      write(fid_alog,'("from get_growth_mode",2I4,50I8)') isect_seed
!!$!     *          ,omp_get_thread_num(),seed_priv(1:sizeseed)
!!$
!!$
!!$!     roll a dice
!!$
!!$!cc      x1=RAND2()
!!$      x1 = rnum1
!!$
!!$
!!$
!!$!cc      x2=RAND2()
!!$      x2 = rnum2
!!$
!!$!tmp      write(fid_alog,'("bf get_growth_mode",4I5,3I15,ES15.6)') isect_seed
!!$!tmp     *          ,omp_in_parallel()
!!$!tmp     *          ,omp_get_num_threads(),omp_get_thread_num()
!!$!tmp     *          ,JSEED,IFRST,NEXTN,x
!!$!tmp      write(fid_alog,*) "roll 1:",x1,x2
!!$
!!$!try      i_tmp_le0=0.5_RP*(1.0_RP+sign(1.0_RP,273.16_RP-tmp))
!!$!try      i_tmp_gtm45=0.5_RP*(1.0_RP-sign(1.0_RP,-45.0_RP-tmpc))
!!$
!!$!try      i_x_lepol=0.5_RP*(1.0_RP+sign(1.0_RP,pol_frq(itmp,isi)-x1))
!!$
!!$!try      i_x_ltppo=0.5_RP*(1.0_RP-sign(1.0_RP,x2-ppo_frq(itmp,isi)))
!!$!try      i_x_ltros=0.5_RP*(1.0_RP-sign(1.0_RP,x2-ros_frq(itmp,isi)))
!!$
!!$!try      i_x_lepolpla=0.5_RP*(1.0_RP+sign(1.0_RP,max(pol_frq(itmp,isi),0.0_RP)+&
!!$!try                                     max(pla_frq(itmp,isi),0.0_RP)-x1))
!!$
!!$
!!$!try      igm= i_x_lepol*( &
!!$!try              i_x_ltppo*i_tmp_gtm45 * 5 &
!!$!try          +(1-i_x_ltppo*i_tmp_gtm45)*i_x_ltros*(1-i_tmp_gtm45)* 4 &
!!$!try          +(1-i_x_ltppo*i_tmp_gtm45)*(1-i_x_ltros*(1-i_tmp_gtm45))* 1 &
!!$!try                                 )+ &
!!$!try                    (1-i_x_lepol)*i_x_lepolpla*  2 + &
!!$!try                    (1-i_x_lepol)*(1-i_x_lepolpla)*  3
!!$
!!$      if(tmp>273.16) then
!!$         igm=0
!!$      elseif(x1<=pol_frq(itmp,isi)) then
!!$         igm=1
!!$!tmp          write(fid_alog,*) "roll 2:",x
!!$!       assume that planar and columnar polycrystals are mutually exclusive.
!!$         if(x2<ppo_frq(itmp,isi).and.tmpc>-45.0) then
!!$            igm=5
!!$!            write(fid_alog,'("planar polycrystals nucleated",4ES15.6)')
!!$!     *        x,ppo_frq(itmp,isi),tmpc,si
!!$
!!$         elseif(x2<ros_frq(itmp,isi).and.tmpc<=-45.0) then
!!$            ! assume that rosetta starts only below T<-45
!!$!ccc            if(tmp-273.16>-35) then
!!$!ccc               write(fid_alog,*) itmp,isi, ros_frq(itmp,isi)
!!$!ccc            end if
!!$            igm=4
!!$         endif
!!$      else if(x1<=max(pol_frq(itmp,isi),0.0_RP) &
!!$                 +max(pla_frq(itmp,isi),0.0_RP)) then
!!$         igm=2
!!$      else
!!$         igm=3
!!$      end if
!!$
!!$
!!$!      if(tmpc<-20.0_RP) then
!!$!      write(fid_alog,'("Tc,si,x,mode,pol,pla,col",3ES15.6,I3,3ES15.6)')
!!$!     *          tmpc,si,x
!!$!     *         ,get_growth_mode
!!$!     *         ,pol_frq(itmp,isi),pla_frq(itmp,isi),col_frq(itmp,isi)
!!$!
!!$!      endif
!!$
!!$      return
!!$      end function get_growth_mode_inl


subroutine cal_growth_mode_hex_inl_vec(igm,nbin,L,ihabit_gm_random &
                                      ,tmp,si,rdsd)
 use com_amps, only: pla_frq,col_frq
 implicit none
 integer,intent(in) :: nbin,L,ihabit_gm_random
 integer,dimension(*) :: igm
 real(PS),intent(in),dimension(*) :: tmp,si
 type(random_genvar),intent(inout) :: rdsd
 !
 ! new space
 !
 real(PS),dimension(nbin*L)  :: rnum
 !
 real(PS) :: dtmp,sttmp,edtmp,tmpc
 parameter(dtmp=1.0_RP,sttmp=-70.0_RP,edtmp=-20.0_RP)
 real(PS) :: dsi,stsi,edsi
 parameter(dsi=0.01_RP,stsi=0.0_RP,edsi=1.0_RP)
 integer :: itmpmin,itmpmax
 parameter(itmpmin=1,itmpmax=51)
 integer :: isimin,isimax
 parameter(isimin=1,isimax=101)
 integer :: itmp,isi
 integer :: i,n,in
 real(PS) :: x
!!!      real rand2
!
  if(ihabit_gm_random.eq.1) then

!CDIR NOVECTOR
    do in=1,nbin*L
      call rand2_ty(rnum(in),rdsd)
    enddo

!CDIR NODEP
    do in=1,nbin*L
      n=(in-1)/NBIN+1
      i=in-(n-1)*NBIN

      tmpc=tmp(n)-273.16
      itmp=max(itmpmin,min(nint((tmpc-sttmp)/dtmp)+1,itmpmax))
      isi=max(isimin,min(nint(abs(si(n)-stsi)/dsi)+1,isimax))

      x = rnum(in)

!!!      write(fid_alog,*) "roll 3:",x

      if(x<=max(pla_frq(itmp,isi),0.0_RP)) then
        igm(in)=1
      else
        igm(in)=2
      end if
    enddo
  else

!CDIR NODEP
    do in=1,nbin*L
      n=(in-1)/NBIN+1
      i=in-(n-1)*NBIN

      tmpc=tmp(n)-273.16
      itmp=max(itmpmin,min(nint((tmpc-sttmp)/dtmp)+1,itmpmax))
      isi=max(isimin,min(nint(abs(si(n)-stsi)/dsi)+1,isimax))

      if(pla_frq(itmp,isi)>col_frq(itmp,isi)) then
!org      elseif(pla_frq(itmp,isi)== &
!org            max(pla_frq(itmp,isi),col_frq(itmp,isi))) then
         igm(in)=1
      else
         igm(in)=2
      end if
    enddo

  endif

  return
end subroutine cal_growth_mode_hex_inl_vec

!!$      function get_growth_mode_hex_inl(tmp,si,rnum) result(igm)
!!$      use com_amps, only: pla_frq,col_frq
!!$      implicit none
!!$      integer :: igm
!!$      real(PS),intent(in) :: tmp,si,rnum
!!$      ! new space
!!$      real(PS) :: dtmp,sttmp,edtmp,tmpc
!!$      parameter(dtmp=1.0_RP,sttmp=-70.0_RP,edtmp=-20.0_RP)
!!$      real(PS) :: dsi,stsi,edsi
!!$      parameter(dsi=0.01_RP,stsi=0.0_RP,edsi=1.0_RP)
!!$      integer :: itmpmin,itmpmax
!!$      parameter(itmpmin=1,itmpmax=51)
!!$      integer :: isimin,isimax
!!$      parameter(isimin=1,isimax=101)
!!$      integer :: itmp,isi
!!$      integer :: i_tmp_le0,i_x_lepla,i_testv_gt0
!!$      real(PS) :: x
!!$!!!      real rand2
!!$!
!!$      tmpc=tmp-273.16
!!$      itmp=max(itmpmin,min(nint((tmpc-sttmp)/dtmp)+1,itmpmax))
!!$      isi=max(isimin,min(nint(abs(si-stsi)/dsi)+1,isimax))
!!$
!!$!     roll a dice
!!$!!!!      call random_number(x)
!!$!!!      x=RAND2()
!!$      x = rnum
!!$
!!$!!!      write(fid_alog,*) "roll 3:",x
!!$
!!$!try      i_tmp_le0=0.5_RP*(1.0+sign(1.0_RP,273.16_RP-tmp))
!!$
!!$!try      x=x*(max(pla_frq(itmp,isi),0.0_RP)+max(col_frq(itmp,isi),0.0_RP))
!!$
!!$!try      i_x_lepla=0.5_RP*(1.0+sign(1.0_RP,max(pla_frq(itmp,isi),0.0_RP)-x))
!!$
!!$      if(tmp>273.16_RP) then
!!$        igm=0
!!$      elseif(x<=max(pla_frq(itmp,isi),0.0_RP)) then
!!$        igm=1
!!$      else
!!$        igm=2
!!$      end if
!!$
!!$!try      igm=i_tmp_le0*(&
!!$!try              i_x_lepla*1+(1-i_x_lepla)*2&
!!$!try                )
!!$
!!$      return
!!$      end function get_growth_mode_hex_inl

!!$      integer function get_growth_mode_hex(tmp,si)
!!$      use com_amps
!!$      implicit none
!!$      real(PS) :: tmp,si,frq(3)
!!$      real(PS) :: dtmp,sttmp,edtmp,tmpc
!!$      parameter(dtmp=1.0_RP,sttmp=-70.0_RP,edtmp=-20.0_RP)
!!$      real(PS) :: dsi,stsi,edsi
!!$      parameter(dsi=0.01_RP,stsi=0.0_RP,edsi=1.0_RP)
!!$      integer :: itmpmin,itmpmax
!!$      parameter(itmpmin=1,itmpmax=51)
!!$      integer :: isimin,isimax
!!$      parameter(isimin=1,isimax=101)
!!$      integer :: itmp,isi
!!$      real(PS) :: x
!!$!tmp      real rand2
!!$!
!!$      if(tmp>273.16) then
!!$         get_growth_mode_hex=0
!!$         return
!!$      endif
!!$      tmpc=tmp-273.16
!!$      itmp=max(itmpmin,min(nint((tmpc-sttmp)/dtmp)+1,itmpmax))
!!$      isi=max(isimin,min(nint(abs(si-stsi)/dsi)+1,isimax))
!!$
!!$!     roll a dice
!!$!      call random_number(x)
!!$       x=RAND2()
!!$!ccc      write(fid_alog,*) "roll 1:",x
!!$
!!$      x=x*(max(pla_frq(itmp,isi),0.0_RP)+max(col_frq(itmp,isi),0.0_RP))
!!$
!!$      if(x<=max(pla_frq(itmp,isi),0.0_RP)) then
!!$         get_growth_mode_hex=1
!!$      else
!!$         get_growth_mode_hex=2
!!$      end if
!!$
!!$!tmp      get_growth_mode_hex=1
!!$
!!$      return
!!$      end function get_growth_mode_hex

      integer function get_growth_mode_max(tmp,si)
      use com_amps
      implicit none
      real(PS) :: tmp,si,frq(3)
      real(PS) :: dtmp,sttmp,edtmp,tmpc
      parameter(dtmp=1.0_RP,sttmp=-70.0_RP,edtmp=-20.0_RP)
      real(PS) :: dsi,stsi,edsi
      parameter(dsi=0.01_RP,stsi=0.0_RP,edsi=1.0_RP)
      integer :: itmpmin,itmpmax
      parameter(itmpmin=1,itmpmax=51)
      integer :: isimin,isimax
      parameter(isimin=1,isimax=101)
      integer :: itmp,isi
      real(PS) :: x
      if(tmp>273.16) then
         get_growth_mode_max=0
         return
      endif
      tmpc=tmp-273.16
      itmp=max(itmpmin,min(nint((tmpc-sttmp)/dtmp)+1,itmpmax))
      isi=max(isimin,min(nint(abs(si-stsi)/dsi)+1,isimax))

      if(pol_frq(itmp,isi)== &
            max(pol_frq(itmp,isi),pla_frq(itmp,isi),col_frq(itmp,isi))) &
            then
         get_growth_mode_max=1
!       assume that planar and columnar polycrystals are mutually exclusive.
         if(tmpc>-45.0.and.&
            ppo_frq(itmp,isi)== &
            max(ppo_frq(itmp,isi),1.0-ppo_frq(itmp,isi)-ros_frq(itmp,isi))) then
            get_growth_mode_max=5
         elseif(tmpc<=-45.0.and.&
            ros_frq(itmp,isi)== &
            max(ros_frq(itmp,isi),1.0-ppo_frq(itmp,isi)-ros_frq(itmp,isi))) then
            ! assume that rosetta starts only below T<-45
            get_growth_mode_max=4
         end if
      elseif(pla_frq(itmp,isi)== &
            max(pla_frq(itmp,isi),col_frq(itmp,isi))) then
         get_growth_mode_max=2
      else
         get_growth_mode_max=3
      end if

      return
      end function get_growth_mode_max

!!$      function get_growth_mode_hex_max(tmp,si) result(igm)
!!$      use com_amps, only: pla_frq,col_frq
!!$      implicit none
!!$      integer :: igm
!!$      real(PS),intent(in) :: tmp,si
!!$      real(PS) :: dtmp,sttmp,edtmp,tmpc
!!$      parameter(dtmp=1.0_RP,sttmp=-70.0_RP,edtmp=-20.0_RP)
!!$      real(PS) :: dsi,stsi,edsi
!!$      parameter(dsi=0.01_RP,stsi=0.0_RP,edsi=1.0_RP)
!!$      integer :: itmpmin,itmpmax
!!$      parameter(itmpmin=1,itmpmax=51)
!!$      integer :: isimin,isimax
!!$      parameter(isimin=1,isimax=101)
!!$      integer :: itmp,isi
!!$!try      integer i_tmp_le0,i_pla_gtcol
!!$
!!$      tmpc=tmp-273.16
!!$      itmp=max(itmpmin,min(nint((tmpc-sttmp)/dtmp)+1,itmpmax))
!!$      isi=max(isimin,min(nint(abs(si-stsi)/dsi)+1,isimax))
!!$
!!$!try      i_tmp_le0=0.5_RP*(1.0_RP+sign(1.0_RP,273.16_RP-tmp))
!!$
!!$!try      i_pla_gtcol=0.5_RP*(1.0-sign(1.0_RP,col_frq(itmp,isi)-pla_frq(itmp,isi)))
!!$
!!$!try      igm=i_tmp_le0*(&
!!$!try              i_pla_gtcol*1+(1-i_pla_gtcol)*2&
!!$!try                )
!!$
!!$
!!$      if(tmp>273.16) then
!!$         igm=0
!!$      elseif(pla_frq(itmp,isi)>col_frq(itmp,isi)) then
!!$!org            max(pla_frq(itmp,isi),
!!$!org      elseif(pla_frq(itmp,isi)== &
!!$!org            max(pla_frq(itmp,isi),col_frq(itmp,isi))) then
!!$         igm=1
!!$      else
!!$         igm=2
!!$      end if
!!$
!!$      return
!!$      end function get_growth_mode_hex_max


  subroutine init_inherent_growth_par
  use com_amps
  implicit none

  if ( IsMaster .and. debug ) then
    write(fid_alog,*) "Initializing inherent growth parameterization"
  end if
  ! +++ definition +++
   nok_igp=23

   x_igp = (/-60.0,-55.0, -50.0, -45.0, -40.0, -35.0, -30.0, -27.0, -25.0, -23.0, -21.0, -20.0, -17.0, &
             -15.0, -12.0, -10.0, -8.0, -6.0, -5.0, -4.0, -3.5, -2.5, -1.5 /)


!!c    a(1,:) = (/  0.0008,  -0.0120,  0.0, 2.3000/)
!!c    a(2,:) = (/ 0.0044, -0.0339, -0.0600, 2.1000/)
!!c    a(3,:) = (/ 0.0012, -0.0018, -0.0706,  1.5000/)
!!c    a(1,:)=(/ 0.0017 , -0.0287 , 0.0 , 4.5000 /)
!!c    a(2,:)=(/ 0.0126 , -0.1036 , -0.1565 , 4.0000 /)
!!c    a(3,:)=(/ 0.0053 , -0.0145 , -0.2487 , 2.2000 /)


!!c    a(1,:)=(/ 0.0003 , -0.0074 , -0.0700 , 3.5000 /)
!!c    a(2,:)=(/ 0.0009 , -0.0120 , -0.1231 , 3.0000 /)
!!c    a(3,:)=(/ 0.0083 , -0.0445 , -0.1737 , 2.2000 /)
!!c    a(4,:) = (/    0.0,       0.0,         0.0,    1.2500/)
!!c    a(5,:) = (/ -0.0025,    0.0220,        0.0,    1.2500/)
!!c    a(6,:) = (/ -0.0073,   0.0375,   0.0640,    1.3800/)
!!c    a(7,:) = (/ -0.0433,    0.0981,    0.1269,    1.6000/)
!!c    a(8,:) = (/  0.0700,   -0.3150,      0.0,    1.9000/)
!!c    a(9,:) = (/  0.3556,   -0.4356, -0.4200,    1.2000/)
    a_igp(1,:)=(/ 0.000000e+00 , 0.000000e+00 , -3.000000e-02 , 2.300000e+00 /)
    a_igp(2,:)=(/ 0.000000e+00 , 0.000000e+00 , -3.000000e-02 , 2.150000e+00 /)
    a_igp(3,:)=(/ -4.800000e-04 , 2.400000e-03 , -3.000000e-02 , 2.000000e+00 /)
    a_igp(4,:)=(/ 1.586667e-03 , -1.353333e-02 , -4.200000e-02 , 1.850000e+00 /)
    a_igp(5,:)=(/ 1.512821e-03 , -5.897436e-03 , -5.833333e-02 , 1.500000e+00 /)
    a_igp(6,:)=(/ -9.597381e-05 , 8.490998e-04 , -3.846154e-03 , 1.250000e+00 /)
    a_igp(7,:)=(/ -1.176598e-04 , 9.293226e-05 , -2.553191e-03 , 1.240000e+00 /)
    a_igp(8,:)=(/ 2.040230e-03 , -6.494253e-03 , -5.172414e-03 , 1.230000e+00 /)
    a_igp(9,:)=(/ -1.439394e-03 , 3.712121e-03 , -6.666667e-03 , 1.210000e+00 /)
    a_igp(10,:)=(/ -1.597052e-03 , -1.726044e-02 , -9.090909e-03 , 1.200000e+00 /)
    a_igp(11,:)=(/ 4.920791e-01 , -7.947818e-01 , -9.729730e-02 , 1.100000e+00 /)

    b_igp(1,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(2,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(3,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(4,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(5,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(6,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(7,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(8,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(9,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(10,:) = (/ 0.0 , 0.0 , 0.0 , 0.7600 /)
    b_igp(11,:) = (/ 0.0431 , -0.1031 , 0.0 , 0.7600 /)

    a_igp(12,:) = (/-0.0012,0.0363,-0.2244,0.7000/)
!!c    a_igp(10,:) = (/ -0.0074, 0.0807, -0.3355, 0.8000 /)
    a_igp(13,:) = (/ -0.0005, 0.0141, -0.0511, 0.3200 /)
    a_igp(14,:) = (/  0.0035, 0.0109, -0.0011, 0.2700 /)
!    a_igp(8,:) = (/ -0.0027,  0.0293,  0.0, 0.2700 /)
    a_igp(15,:) = (/ -0.0063, 0.0427, 0.1597, 0.4600 /)
!    a_igp(9,:) = (/ -0.0169, 0.0928, 0.1021, 0.4600 /)
    a_igp(16,:) = (/  0.0212, 0.0051, 0.2552, 0.9000 /)
!    a_igp(10,:) = (/ -0.0200, 0.0798, 0.2702, 0.9000 /)
    a_igp(17,:) = (/ -0.1109, 0.1320, 0.5294, 1.6000 /)
    a_igp(18,:) = (/  0.1061, -0.5332, -0.2729, 2.3000 /)
    a_igp(19,:) = (/  0.2758, -0.2148, -1.0209, 1.6000 /)
    a_igp(20,:) = (/ -0.0118, 0.6125, -0.6233, 0.6400 /)
    a_igp(21,:) = (/ -0.2951, 0.5948, -0.0197, 0.4800 /)
    a_igp(22,:) = (/  0.0550,   -0.0995, 0.1244, 0.7600 /)
    a_igp(23,:) = (/ -0.0000, 0.0108, 0.0906, 0.8400 /)

    b_igp(12:23,1:4)=a_igp(12:23,1:4)

  end subroutine init_inherent_growth_par

  subroutine read_AMPSTASK
    use scale_prc, only: &
       PRC_abort
    use com_amps
    implicit none
      !integer,parameter :: log_fid=101,ctl_fid=102
      integer :: LOG_FID, CTL_FID
      integer :: ierr

!     +++ input variables +++
      namelist /AMPS_param/&
         level_comp, &
         debug, &
         debug_level, &
         coll_level, &
         out_type, &
         T_print_period, &
         output_format, &
         token_a, &
         token_c, &
         token_r, &
         token_s, &
         dtype_c, &
         dtype_r, &
         dtype_s, &
         dtype_a, &
         hbreak_c, &
         hbreak_r, &
         hbreak_s, &
         flagp_c, &
         flagp_r, &
         flagp_s, &
         flagp_a, &
         act_type, &
         srat_c, &
         srat_r, &
         sth_r, &
         srat_s, &
         srat_a, &
         fcon_c, &
         sadd_c, &
         sadd_r, &
         sadd_s, &
         sadd_a, &
         minmass_c, &
         minmass_r, &
         minmass_s, &
         minmass_a, &
         M_aps, &
         M_api, &
         den_aps, &
         den_api, &
         nu_aps, &
         phi_aps, &
         ap_lnsig, &
         ap_mean, &
         N_ap_ini, &
         eps_ap, &
         ap_mean_cp, &
         ap_sig_cp, &
         APSNAME, &
         DRCETB, &
         DRAPTB, &
         DRSTTB, &
         micexfg, &
         n_step_cl, &
         n_step_vp, &
         ihabit_gm_random, &
         CCNMAX, &
         frac_dust, &
         CRIC_RN_IMM

        CTL_FID = IO_get_available_fid()
        open(CTL_FID, &
          file='AMPSTASK.F', &
          form='formatted',   &
          status='old',       &
          iostat=ierr)
        if(ierr/=0) then
           LOG_ERROR("read_AMPSTASK",*) 'Cannot open AMPS PARAMETER file, please check.'
           call PRC_abort
        end if
        rewind(CTL_FID)
        read(CTL_FID,nml=AMPS_param)
        close(CTL_FID)

        if ( IsMaster .and. debug ) then
          LOG_FID = IO_get_available_fid()
          open(LOG_FID, &
            file='AMPSTASK'//'_'//string_fid_amps//'.log', &
            form='formatted',   &
            status='replace',       &
            iostat=ierr)

          write(LOG_FID,nml=AMPS_param)
          write(LOG_FID,*) "THIS IS STANDARD CODE"
          close(LOG_FID)
        end if

      end subroutine read_AMPSTASK
! AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS
! AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS
      subroutine config_varindex
      use scale_prc, only: &
         PRC_abort
      use maxdims
      use par_amps
      use com_amps
      implicit none
!
!     Configure the variables index
!
!
      integer :: ivar

      integer :: ierr

      !---------------------------------------------
      !
      !  Ice configuration
      !
      !    total mass of ice particle
      imt_q=1
      ivar=1

      !    rime
      ivar=ivar+1
      imr_q=ivar
      i1_mcp_ice=ivar

      !    aggregate
      ivar=ivar+1
      ima_q=ivar

      !    crystal
      ivar=ivar+1
      imc_q=ivar

      !    melt water
      ivar=ivar+1
      imw_q=ivar

      !    frozen water (nucleation)
      ivar=ivar+1
      imf_q=ivar

      !    total aerosol mass
      ivar=ivar+1
      imat_q=ivar

      !    soluble aerosol mass
      ivar=ivar+1
      imas_q=ivar

      i2_mcp_ice=ivar
      nvar_mcp_ice=i2_mcp_ice-i1_mcp_ice+1

      !    number concentration
      ivar=ivar+1
      icon_q=ivar

      !    circiumcsribing volume (dry)
      ivar=ivar+1
      i1_vcp_ice=ivar
      ivcs_q=ivar

      i2_vcp_ice=ivar
      nvar_vcp_ice=i2_vcp_ice-i1_vcp_ice+1


      !    a axis length
      ivar=ivar+1
      iacr_q=ivar
      i1_acp_ice=ivar

      !    c axis length
      ivar=ivar+1
      iccr_q=ivar

      !    dendritic length
      ivar=ivar+1
      idcr_q=ivar

      !    center of gravity coordinates, a
      ivar=ivar+1
      iag_q=ivar

      !    center of gravity coordinates, c
      ivar=ivar+1
      icg_q=ivar
      i2_acp_ice=ivar

      nvar_acp_ice=i2_acp_ice-i1_acp_ice+1

      !    extra crystalline structure
      ivar=ivar+1
      inex_q=ivar
      i1_ccp_ice=ivar

      !    activated IN concentration for contact parameter diagnosis
!org      ivar=ivar+1
!org      icnc_q=ivar
      i2_ccp_ice=ivar

      nvar_ccp_ice=i2_ccp_ice-i1_ccp_ice+1

      ! total non-mass component variables
      nvar_nonmcp_ice=nvar_acp_ice+nvar_vcp_ice+nvar_ccp_ice

      ! total number of Property Variables
      nvar_pv_ice=nvar_mcp_ice+nvar_nonmcp_ice

      !
      ! array specification for g%MS%mass,g%MS%dmassdt,new_M
      ivar=1
      imt=1

      ivar=ivar+1
      imr=ivar

      ivar=ivar+1
      ima=ivar

      ivar=ivar+1
      imc=ivar

      ivar=ivar+1
      imw=ivar

      ivar=ivar+1
      imf=ivar

      ivar=ivar+1
      imat=ivar

      ivar=ivar+1
      imas=ivar

      ivar=ivar+1
      imai=ivar

      ! array specification for ratio_M
      imr_m=imr-1

      ima_m=ima-1

      imc_m=imc-1

      imw_m=imw-1

      imf_m=imf-1

      imat_m=imat-1

      imas_m=imas-1

      imai_m=imai-1

      ! array specification for new_Q
      ivar=0

      ivar=ivar+1
      ivcs=ivar

      ivar=ivar+1
      iacr=ivar

      ivar=ivar+1
      iccr=ivar

      ivar=ivar+1
      idcr=ivar

      ivar=ivar+1
      iag=ivar

      ivar=ivar+1
      icg=ivar

      ivar=ivar+1
      inex=ivar

!org      ivar=ivar+1
!org      icnc=ivar

      ! array specification for g%MSvol
!org      icnc_v=1

      !---------------------------------------------
      !
      !  Liquid configuration
      !
      !---------------------------------------------
      ! array specfication for qrpv
      !    total mass of liq particle
      ivar=1
      rmt_q=ivar

      !    total aerosol mass
      ivar=ivar+1
      rmat_q=ivar
      i1_mcp_liq=ivar

      !    soluble aerosol mass
      ivar=ivar+1
      rmas_q=ivar
      i2_mcp_liq=ivar

      nvar_mcp_liq=i2_mcp_liq-i1_mcp_liq+1

      !    number concentration
      ivar=ivar+1
      rcon_q=ivar

      !    activated IN concentration for contact parameter diagnosis
!org      ivar=ivar+1
!org      rcnc_q=ivar

!org      i1_ccp_liq=ivar
!org      i2_ccp_liq=ivar
!org      nvar_ccp_liq=i2_ccp_liq-i1_ccp_liq+1

      i1_ccp_liq=1
      i2_ccp_liq=-9
      nvar_ccp_liq=0

      ! total non-mass component variables
      nvar_nonmcp_liq=nvar_ccp_liq

      ! total number of Property Variables
      nvar_pv_liq=nvar_mcp_liq+nvar_nonmcp_liq


      ! array specification for g%MS%mass,g%MS%dmassdt,new_M
      ivar=1
      rmt=1

      ivar=ivar+1
      rmat=ivar

      ivar=ivar+1
      rmas=ivar

      ivar=ivar+1
      rmai=ivar


      ! array specification for ratio_M
      rmat_m=rmat-1

      ivar=ivar+1
      rmas_m=rmas-1

      ivar=ivar+1
      rmai_m=rmai-1

      ! array specification for new_Q
!org      rcnc=1

      ! array specification for g%MS%vol,g%MS%dvoldt
!org      ivar=1
!org      rcnc_v=ivar

      !---------------------------------------------
      !
      !  Aerosol configuration
      !
      !---------------------------------------------
      ! array specfication for qapv
      ivar=1
      amt_q=1

      ivar=ivar+1
      acon_q=ivar

      ivar=ivar+1
      ams_q=ivar
      i1_mcp_aer=ivar
      i2_mcp_aer=ivar

      nvar_mcp_aer=i2_mcp_aer-i1_mcp_aer+1

      !    activated IN concentration for contact parameter diagnosis
!org      ivar=ivar+1
!org      acnc_q=ivar

!org      i1_ccp_aer=ivar
!org      i2_ccp_aer=ivar
!org      nvar_ccp_aer=i2_ccp_aer-i1_ccp_aer+1

      i1_ccp_aer=1
      i2_ccp_aer=-9
      nvar_ccp_aer=0

      ! total number of Property Variables
      nvar_pv_aer=nvar_mcp_aer+nvar_ccp_aer

      ! array specification for g%MS%mass,g%MS%dmassdt,new_M
      ivar=1
      amt=1

      ivar=ivar+1
      ams=ivar

      ivar=ivar+1
      ami=ivar

      ! array specification for new_Q
!org      ivar=1
!org      acnc=ivar

      ! array specification for g%MS%vol,g%MS%dvoldt
!org      ivar=1
!org      acnc_v=ivar


      if ( IsMaster .and. debug ) then
        write(fid_alog,*) " "
        write(fid_alog,*) "*** Ice configuration *** "

        write(fid_alog,100) imt_q,imr_q,ima_q,imc_q,imw_q,imf_q,imat_q,imas_q
100 format('Total mass:',I3,', Rim mass:',I3,', Agg mass:',I3,', Crystal mass:',I3,/ &
          ,'Melt-water mass:',I3,', Frozen mass:',I3,', Total aerosol mass:',I3,/ &
          ,'Soluble mass: ',I3)

        write(fid_alog,105) nvar_mcp_ice,i1_mcp_ice,i2_mcp_ice
105 format('Total number of mass components:',I3,', Starting index:',I3,', Ending index:',I3)

        write(fid_alog,110) icon_q,ivcs_q,iacr_q,iccr_q,idcr_q,iag_q,icg_q,inex_q
110 format('Con:',I3,', Volume:',I3,/ &
          ,'a-axis:',I3,', c-axis:',I3,', dendric:',I3,/ &
          ,'cg a-axis:',I3,', cg c-axis:',I3,', no of extra cry:',I3)


        write(fid_alog,115) nvar_vcp_ice,i1_vcp_ice,i2_vcp_ice
115 format('Total number of volme components:',I3,', Starting index:',I3,', Ending index:',I3)

        write(fid_alog,120) nvar_acp_ice,i1_acp_ice,i2_acp_ice
120 format('Total number of axis components:',I3,', Starting index:',I3,', Ending index:',I3)

        write(fid_alog,125) nvar_nonmcp_ice
125 format('Number of non-mass variables:',I3)

        write(fid_alog,130) nvar_pv_ice
130 format('Number of Property variables:',I3)


        write(fid_alog,*) " "
        write(fid_alog,*) "*** Liquid configuration ***"

        write(fid_alog,200) rmt_q,rmat_q,rmas_q
200 format('Total mass:',I3,', Total aerosol mass:',I3,', Soluble mass:',I3)

        write(fid_alog,205) nvar_mcp_liq,i1_mcp_liq,i2_mcp_liq
205 format('Total number of mass components:',I3,', Starting index:',I3,', Ending index:',I3)

        write(fid_alog,210) rcon_q
210 format('Con:',I3)

        write(fid_alog,215) nvar_nonmcp_liq
215 format('Number of non-mass variables:',I3)

        write(fid_alog,220) nvar_pv_liq
220 format('Number of Property variables:',I3)

        write(fid_alog,*) " "
        write(fid_alog,*) "*** Aerosol configuration ***"

        write(fid_alog,300) amt_q,ams_q
300 format('Total mass:',I3,', Soluble mass:',I3)

        write(fid_alog,305) nvar_mcp_aer,i1_mcp_aer,i2_mcp_aer
305 format('Total number of mass components:',I3,', Starting index:',I3,', Ending index:',I3)

        write(fid_alog,310) acon_q
310 format('Con:',I3)

        write(fid_alog,320) nvar_pv_aer
320 format('Number of Property variables:',I3)
      end if

!
!  Check the maxdims
!
      ! max number of mass components and volume components for ice particles
      if(mxnmasscomp<nvar_mcp_ice) then
        LOG_ERROR("config_varindex",*) "Error! Increase mxnmasscomp.",mxnmasscomp,nvar_mcp_ice
        call PRC_abort
      endif
      if(mxnvol<nvar_vcp_ice) then
        LOG_ERROR("config_varindex",*) "Error! Increase mxnvol.",mxnvol,nvar_vcp_ice
        call PRC_abort
      endif

      ! axis lengths, and non-mass component
      ! variables for ice particles.
      if(mxnaxis<nvar_acp_ice) then
        LOG_ERROR("config_varindex",*) "Error! Increase mxnvol.",mxnaxis,nvar_acp_ice
        call PRC_abort
      endif
      if(mxnnonmc<nvar_nonmcp_ice) then
        LOG_ERROR("config_varindex",*) "Error! Increase mxnvol.",mxnnonmc,nvar_nonmcp_ice
        call PRC_abort
      endif
      !
      ! max number of mass components for liquid particles
      if(mxnmasscomp_r<nvar_mcp_liq) then
        LOG_ERROR("config_varindex",*) "Error! Increase mxnmasscomp.",mxnmasscomp_r,nvar_mcp_liq
        call PRC_abort
      endif

      end subroutine config_varindex
! AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS
! AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS

  function get_inact(sv2) result(N_IN)
    real(PS),intent(in) :: sv2
    real(PS):: N_IN
    real(PS), parameter           :: b = 0.1296_RP
    real(PS), parameter           :: a = -0.639_RP
    ! ------------------------------------------------------------------------------------
    ! Meyers et al (1992) diagnostic value
    N_IN = exp(a+b*100.0_RP*min(sv2,1.0_RP))*1.0e-3_PS
    ! ------------------------------------------------------------------------------------
  end function get_inact

! added for KWAJEX
  function get_inact_tropic(sv2,T) result(N_IN)
    ! modified Meyer formula according to tropics condition by Phillips et al. 2007.
    real(PS),intent(in) :: sv2,T
    real(PS):: N_IN
    real(PS), parameter           :: b = 12.96_RP
    real(PS), parameter           :: a = -0.639_RP
    real(PS), parameter           :: psi=0.06_RP
    ! ------------------------------------------------------------------------------------
    if(T>=243.16) then
      ! Meyers et al (1992) diagnostic value
!org      N_IN = exp(a+b*sv2)*1.0e-3_PS*psi
      N_IN = exp(min(a+b*sv2,80.0_RP))*1.0e-3_PS*psi
    else
      ! DeMott et al (2003)'s formula
!org      N_IN = 1.0e-3_RP*(exp(b*(sv2-0.1_RP)))**0.3_RP
      N_IN = 1.0e-3_RP*(exp(min(b*(sv2-0.1_RP),80.0_RP)))**0.3_RP
    endif
    ! ------------------------------------------------------------------------------------
  end function get_inact_tropic
! end added for KWAJEX

!!$  function get_inact_cnt(ev,sv2,T,mean_mass_ap,eps_map,den_ai,Ax,miv,N_dep,dt) result(N_IN)
!!$    ! classical nucleation theory (see p410 of Khvorostyanov and Curry 2014)
!!$    ! follows Savre and Ekman (2014) for deposition freezing
!!$    real(PS),intent(in) :: ev,sv2,T
!!$    real(PS),intent(in) :: mean_mass_ap,eps_map,den_ai
!!$    ! time step
!!$    real(PS),intent(in) :: dt
!!$    ! pre-exponential factor and contact parameter, number of deposition IN
!!$    real(PS),intent(in) :: Ax,miv,N_dep
!!$    real(PS):: N_IN
!!$    ! bolzman constant [erg/K]
!!$    real(PS),parameter :: ak=1.38e-16
!!$    ! ak * nus where water molecule frequency of vibration nus=1.0e+13 [s^{-1}]
!!$    real(PS),parameter :: akn=1.38e-3
!!$    real(PS) :: fgeo
!!$    real(PS) :: sigma_iv,r_d
!!$    real(PS) :: dGdes,dGg,aJdhet,alogAdep
!!$    ! 16*pi/3.0
!!$    real(PS),parameter :: coef16pi3=16.7551608191456
!!$
!!$    ! ------------------------------------------------------------------------------------
!!$    ! surface tension between ice and air
!!$    !   Hoose and Mohler (2012)
!!$    sigma_iv=76.1-0.155*(T-T_0)+28.5+0.25*(T-T_0)
!!$
!!$    ! geometrical factor under large limit of convex aerosol
!!$    fgeo=0.25*(2.0+miv)*(1.0-miv)*(1.0-miv)
!!$
!!$    ! insoluble aerosol radius
!!$    r_d=(mean_mass_ap*(1.0-eps_map)/&
!!$         den_ai/coef4pi3)**(1.0/3.0)
!!$
!!$    ! condensation coef = 1, and denstity of water molucule in supercooled
!!$    ! region = 1 [g/cm3]
!!$    alogAdep=2.0*log(r_d)+2.0*log(ev)-log(akn*T)+0.5*(log(sigma_iv)-log(ak*T*fgeo))
!!$
!!$    dGdes=Ax*exp(0.03678*(T-T_0))
!!$    dGg=fgeo*coef16pi3*sigma_iv*sigma_iv*sigma_iv/(den_i*R_v*T*log(1.0_PS+sv2))**2
!!$
!!$    aJdhet=exp(alogAdep+max(dGdes-dGg,0.0_PS)/(ak*T))
!!$
!!$    N_IN = N_dep*(1.0-exp(-aJdhet*dt))
!!$    ! ------------------------------------------------------------------------------------
!!$  end function get_inact_cnt

! added for KWAJEX
      real(PS) function get_inact_tropic2(sv2,T)
! modified Meyer formula according to tropics condition by Phillips et al. 2007.
        real(PS) :: sv2,T
        real(PS) :: N_IN
        real(PS), parameter           :: b = 12.96
        real(PS), parameter           :: a = -0.639
        real(PS), parameter           :: psi=0.06
!------------------------------------------------------------------------------------
        if(T>=243.16) then
! Meyers et al (1992) diagnostic value
           N_IN = exp(a+b*sv2)*1.0e-3*psi
        else
! DeMott et al (2003)'s formula
           N_IN = 1.0e-3*(exp(b*(sv2-0.1)))**0.3
        endif
! ------------------------------------------------------------------------------------
        get_inact_tropic2=N_IN
      end function get_inact_tropic2
! end added for KWAJEX

!!$      real(PS) function get_inact_tropic3(sv2,T,kden)
!!$! modified Meyer formula according to tropics condition by Phillips et al. 2007.
!!$        real(PS) :: sv2,T,kden
!!$        real(PS) :: N_IN
!!$        real(PS), parameter           :: b     = 12.96
!!$        real(PS), parameter           :: a     = -0.639
!!$        real(PS), parameter           :: c1    = 1000.0
!!$        real(PS), parameter           :: psi   = 0.058707
!!$        real(PS), parameter           :: gamma = 2.0
!!$        real(PS), parameter           :: rhoc  = 0.67
!!$!------------------------------------------------------------------------------------
!!$        if(T>=243.16) then
!!$! Meyers et al (1992) diagnostic value
!!$           N_IN = psi*gamma/rhoc*c1*exp(a+b*sv2)/kden*0.000001
!!$        else
!!$! DeMott et al (2003)'s formula
!!$           N_IN = gamma/rhoc*c1*(exp(b*(sv2-0.1)))**0.3/kden*0.000001
!!$        endif
!!$! ------------------------------------------------------------------------------------
!!$        get_inact_tropic3=N_IN
!!$      end function get_inact_tropic3

!!$      real(PS) function get_inact_pamarcmip(T)
!!$        implicit none
!!$        real(PS) :: T
!!$        real(PS) :: N_IN
!!$!------------------------------------------------------------------------------------
!!$        N_IN = (782.476D0 - 141.94D0*log(T))
!!$! ------------------------------------------------------------------------------------
!!$        get_inact_pamarcmip=exp(N_IN)*1.0D-3
!!$      end function get_inact_pamarcmip

!!$      real(PS) function getfrac_apact(S_m,eta,zeta,index)
!!$      use com_amps
!!$      implicit none
!!$
!!$!     index=1: mass fraction of activated CCN
!!$!     index=2: con fraction of activated CCN
!!$!ccc      real ap_lnsig,ap_mean,apt_st,apt_add,apt_max,frac_apact
!!$!ccc     *     ,M_ap,den_ap,nu_ap,phi_ap,eps_ap,ap_mass
!!$!ccc      common/RAPTBL/ap_lnsig,ap_mean,apt_st(4),apt_add(4)
!!$!ccc     *     ,apt_max(4),frac_apact(2,8320),M_ap,den_ap,nu_ap
!!$!ccc     *     ,phi_ap,eps_ap,ap_mass(2)
!!$!ccc      integer napt
!!$!ccc      common/IAPTBL/napt(4)
!!$
!!$      real(PS) :: S_m,eta,zeta,z1,e1,s1,w1,w2,w3
!!$      integer :: i1,i2,i3,id1,id2,id3,id4,id5,id6,id7,id8 &
!!$           ,index
!!$
!!$      i1=max(1,min(int((S_m-apt_st(2))/apt_add(2)),napt(2)-1))
!!$      i2=max(1,min(int((eta-apt_st(3))/apt_add(3)),napt(3)-1))
!!$      i3=max(1,min(int((zeta-apt_st(4))/apt_add(4)),napt(4)-1))
!!$
!!$      id1=(i1-1)*napt(3)*napt(4)+(i2-1)*napt(4)+i3
!!$      id2=(i1-1)*napt(3)*napt(4)+i2*napt(4)+i3
!!$      id3=(i1-1)*napt(3)*napt(4)+(i2-1)*napt(4)+i3+1
!!$      id4=(i1-1)*napt(3)*napt(4)+i2*napt(4)+i3+1
!!$      id5=i1*napt(3)*napt(4)+(i2-1)*napt(4)+i3
!!$      id6=i1*napt(3)*napt(4)+i2*napt(4)+i3
!!$      id7=i1*napt(3)*napt(4)+(i2-1)*napt(4)+i3+1
!!$      id8=i1*napt(3)*napt(4)+i2*napt(4)+i3+1
!!$
!!$      z1=apt_st(4)+real(i3-1,PS_KIND)*apt_add(4)
!!$      e1=apt_st(3)+real(i2-1,PS_KIND)*apt_add(3)
!!$      s1=apt_st(2)+real(i1-1,PS_KIND)*apt_add(2)
!!$
!!$      w1=max(0.0_RP,min(1.0_RP,(zeta-z1)/apt_add(4)))
!!$      w2=max(0.0_RP,min(1.0_RP,(eta-e1)/apt_add(3)))
!!$      w3=max(0.0_RP,min(1.0_RP,(S_m-s1)/apt_add(2)))
!!$
!!$      getfrac_apact= &
!!$           (1.0_RP-w3)*((1.0_RP-w1)*(1.0_RP-w2)*frac_apact(index,id1)+ &
!!$           (1.0_RP-w1)*w2*frac_apact(index,id2)+ &
!!$           w1*(1.0_RP-w2)*frac_apact(index,id3)+ &
!!$           w1*w2*frac_apact(index,id4))+ &
!!$           w3*((1.0_RP-w1)*(1.0_RP-w2)*frac_apact(index,id5)+ &
!!$           (1.0_RP-w1)*w2*frac_apact(index,id6)+ &
!!$           w1*(1.0_RP-w2)*frac_apact(index,id7)+ &
!!$           w1*w2*frac_apact(index,id8))
!!$      return
!!$      end function getfrac_apact

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     print out necessary variables for radiative transfer calculation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$      subroutine print_radcal(n1,n2,n3 &
!!$           ,u,v,w &
!!$           ,den,theta,qv,qc,qr,qi &
!!$           ,thp,pp,pi01d &
!!$           ,alatt,alont &
!!$           ,k1eta,zs,zz,deltax,deltay &
!!$           ,sst,pctlnd &
!!$           ,time,ngr &
!!$           ,npa,nba,nca,qap &
!!$           ,npr,nbr,ncr,qrp &
!!$           ,npi,nbi,nci,qip &
!!$           ,pcprl,pcpri)
!!$      use maxdims
!!$      use com_amps
!!$      implicit none
!!$
!!$!      INCLUDE 'common.nms'
!!$      integer :: n1,n2,n3,k1eta(n2,*)
!!$      integer :: npa,nba,nca,npr,nbr,ncr,npi,nbi,nci,ngr
!!$      real(PS) :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3),time
!!$      real(PS) :: den(n1,n2,*),theta(n1,n2,*) &
!!$          ,qv(n1,n2,*),qc(n1,n2,*),qr(n1,n2,*),qi(n1,n2,*) &
!!$          ,thp(n1,n2,*),pp(n1,n2,*),pi01d(*) &
!!$          ,alatt(n2,*),alont(n2,*) &
!!$          ,zz(*),zs(n2,*),deltax,deltay &
!!$          ,sst(n2,*),pctlnd(n2,*) &
!!$          ,qap(npa,nba,nca,n1,n2,*),qrp(npr,nbr,ncr,n1,n2,*) &
!!$          ,qip(npi,nbi,nci,n1,n2,*) &
!!$          ,pcprl(npr,nbr,ncr,n2,*),pcpri(npi,nbi,nci,n2,*)
!!$      integer :: i,j,k,l,n,kint,j1,j2,jj1,jj2,ica,ipa,iba,ipr,icr,ibr &
!!$              ,ipi,ici,ibi
!!$
!!$      real(PS) :: thdry,pird,var,r,cp,cpr,p00
!!$      parameter (var=-999.9,r=287.0,cp=1004.0,cpr=cp/r,p00=1.0e+5)
!!$      integer :: inn
!!$      real(PS) :: tv,pv,denv,qvv,wave,uave,vave,pcprl_tot,pcpri_tot
!!$      CHARACTER(len=10) :: numfil
!!$
!!$      inn=int(time/60)
!!$      if(inn.lt.100) then
!!$         write(numfil,1011) inn,ngr
!!$      elseif(inn.lt.1000) then
!!$         write(numfil,1002) inn,ngr
!!$      else
!!$         write(numfil,1003) inn,ngr
!!$      endif
!!$ 1011 format('prof',i2,'.',i1)
!!$ 1002 format('prof',i3,'.',i1)
!!$ 1003 format('prof',i4,'.',i1)
!!$      open(4,file='Z.loc',status='unknown',form='formatted')
!!$      write(4,*) n1,n2,n3
!!$      write(4,*) deltax,deltay
!!$      do k=1,n1
!!$         write(4,'(i5,f12.4)') k,zz(k)
!!$      enddo
!!$      close(4)
!!$
!!$
!!$      open(ngr,file=numfil,status='unknown',form='formatted')
!!$      write(ngr,'(f10.2)')  time
!!$      write(ngr,'(3I5)')  n1,n2,n3
!!$      write(ngr,'(3I5)')  nbr,nbi
!!$      write(ngr,'(100ES15.6)') (binbr(i),i=1,nbr+1)
!!$      write(ngr,'(100ES15.6)') (binbi(i),i=1,nbi+1)
!!$!      open(14,file='check_qr_bm.dat',status='unknown',form='formatted')
!!$      do i=2,n2-1
!!$         do j=2,n3-1
!!$            pcprl_tot=0.0
!!$            do k=1,ncr
!!$               do l=1,nbr
!!$                  pcprl_tot=pcprl_tot+pcprl(1,l,k,i,j)
!!$               end do
!!$            end do
!!$            pcprl_tot=pcprl_tot*3600.0
!!$            pcpri_tot=0.0
!!$            do k=1,nci
!!$               do l=1,nbi
!!$                  pcpri_tot=pcpri_tot+pcpri(1,l,k,i,j)
!!$               end do
!!$            end do
!!$            pcpri_tot=pcpri_tot*3600.0
!!$
!!$            write(ngr,1000) i,j,k1eta(i,j) &
!!$                    ,alatt(i,j),alont(i,j),zs(i,j),sst(i,j) &
!!$                    ,pctlnd(i,j),pcprl_tot,pcpri_tot
!!$            do k=1,n1-1
!!$               if(k.ge.k1eta(i,j))then
!!$                  thdry=THETA(k,i,j)/(1.+.61*qv(k,i,j))
!!$                  PIRD=(PP(k,i,j)+PI01D(K))/CP
!!$                  tv=thdry*PIRD
!!$                  pv=PIRD**CPR*P00*0.01
!!$                  denv=den(k,i,j)
!!$                  qvv=qv(k,i,j)
!!$!                  qupc(k)=
!!$!                  qrpc(k)=
!!$!                  qagg(k)=
!!$!                  qrag(k)=
!!$!                  qgra(k)=
!!$                  wave=(W(k,i,j)+W(max(1,k-1),i,j))*.5
!!$                  uave=(u(k,i,j)+u(k,max(1,i-1),j))*.5
!!$                  vave=(v(k,i,j)+v(k,i,max(1,j-1)))*.5
!!$               else
!!$                  tv=var
!!$                  pv=var
!!$                  denv=var
!!$                  qvv=var
!!$                  wave=var
!!$                  uave=var
!!$                  vave=var
!!$               endif
!!$               write(ngr,1001) tv,pv,denv,uave,vave,wave,qvv
!!$!               write(14,'(3I5,2ES15.6)') k,i,j
!!$!     *               ,qrp(rmt_q,1,1,k,i,j),qrp(rcon_q,1,1,k,i,j)
!!$            end do
!!$         end do
!!$      end do
!!$ 1000 format(3I5,7f12.4)
!!$ 1001 format(12E16.6)
!!$      close(ngr)
!!$!      close(14)
!!$
!!$      return
!!$      end subroutine print_radcal

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     print out necessary variables for radiative transfer calculation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$      subroutine print_radcal_bulk(n1,n2,n3 &
!!$           ,u,v,w &
!!$           ,den,theta,qv,qc,qr,qi &
!!$           ,thp,pp,pi01d &
!!$           ,alatt,alont &
!!$           ,k1eta,zs,zz,deltax,deltay &
!!$           ,sst,pctlnd &
!!$           ,time,ngr &
!!$           ,rdistp,pdistp,sdistp,adistp,gdistp &
!!$           ,npa,nba,nca,qap &
!!$           ,npr,nbr,ncr,qrp &
!!$           ,npi,nbi,nci,qip &
!!$           ,ircnfl,ipcnfl,iscnfl,iacnfl,igcnfl &
!!$           ,pcprl,pcpri)
!!$      use maxdims
!!$      use com_amps
!!$      implicit none
!!$
!!$      integer :: n1,n2,n3,k1eta(n2,*)
!!$      integer :: npa,nba,nca,npr,nbr,ncr,npi,nbi,nci,ngr &
!!$           ,ircnfl,ipcnfl,iscnfl,iacnfl,igcnfl
!!$!
!!$      real(PS) :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3),time
!!$      real(PS) :: den(n1,n2,*),theta(n1,n2,*) &
!!$          ,qv(n1,n2,*),qc(n1,n2,*),qr(n1,n2,*),qi(n1,n2,*) &
!!$          ,thp(n1,n2,*),pp(n1,n2,*),pi01d(*) &
!!$          ,alatt(n2,*),alont(n2,*) &
!!$          ,zz(*),zs(n2,*),deltax,deltay &
!!$          ,sst(n2,*),pctlnd(n2,*) &
!!$          ,qap(npa,nba,nca,n1,n2,*),qrp(npr,nbr,ncr,n1,n2,*) &
!!$          ,qip(npi,nbi,nci,n1,n2,*) &
!!$          ,pcprl(npr,nbr,ncr,n2,*),pcpri(npi,nbi,nci,n2,*) &
!!$          ,rdistp,pdistp,sdistp,adistp,gdistp
!!$!
!!$      integer :: i,j,k,l,n,kint,j1,j2,jj1,jj2,ica,ipa,iba,ipr,icr,ibr &
!!$              ,ipi,ici,ibi
!!$      real(PS) :: thdry,pird,var,cpr  &
!!$          ,CLWC,RLWC,GIWC,PIWC,AIWC,SIWC
!!$      integer :: inn
!!$      real(PS) :: tv,pv,denv,qvv,wave,uave,vave,pcprl_tot,pcpri_tot
!!$      CHARACTER(len=16) :: numfil
!!$
!!$      var=-99.0
!!$      cpr=cpd/r
!!$
!!$      inn=int(time/60)
!!$      if(inn.lt.100) then
!!$         write(numfil,1011) inn,ngr
!!$      elseif(inn.lt.1000) then
!!$         write(numfil,1002) inn,ngr
!!$      else
!!$         write(numfil,1003) inn,ngr
!!$      endif
!!$ 1011 format('prof_bmp_',i2,'m.',i1)
!!$ 1002 format('prof_bmp_',i3,'m.',i1)
!!$ 1003 format('prof_bmp_',i4,'m.',i1)
!!$
!!$      open(ngr,file=numfil,status='unknown',form='formatted')
!!$
!!$      write(ngr,'("%UW-NMS file name:",a20)') numfil
!!$      write(ngr,*) "%The microphysical outputs are organized for" &
!!$      //" increasing grid point number. The first grid point is at" &
!!$      //" the south-west corner of the domain"
!!$      write(ngr,*) "%Domain Information:"
!!$      write(ngr,'(I12,a)') n2-2 &
!!$                ,"    !number of grid points west-east direction"
!!$      write(ngr,'(I12,a)') n3-2 &
!!$                ,"    !number of grid points south-north direction"
!!$      write(ngr,'(F11.7,a)') deltax/1000.0 &
!!$                ,"     !horizontal spacing (km)"
!!$      write(ngr,'(I12,a)') n1-1 &
!!$                ,"    !number of vertical levels"
!!$      write(ngr,*) "%Vertical level of the atmosphere above the" &
!!$      //" considered grid point (m):"
!!$      do k=1,n1-1
!!$         write(ngr,'(F11.4)') zz(k)
!!$      end do
!!$
!!$      write(ngr,*) "%Microphysical Information:"
!!$      write(ngr,*) "          6     !Number of Hydrometeors"
!!$      write(ngr,*) "%PSD: 1=MONO-DISPERSED,2=CONSTANT-INTERCEPT," &
!!$      //" 3=CONSTANT-SLOPE"
!!$      write(ngr,*) "%PSD  RMIN(cm) RMAX(cm) SLOPE(cm-1)/INTERCEPT(cm-4)" &
!!$      //" DENSITY(gr/cm3)"
!!$      write(ngr,*) "1    0.0009  0.0011      0.000    1.000"
!!$      write(ngr,*) "2    0.0011  0.3000      0.080    1.000"
!!$      write(ngr,*) "2    0.0130  0.3000      0.071    0.900"
!!$      write(ngr,*) "1    0.0110  0.0130      0.000    0.100"
!!$      write(ngr,*) "2    0.0130  0.5000      0.014    0.200"
!!$      write(ngr,*) "3    0.0130  0.6000      3.000    0.015  0.600"
!!$      write(ngr,*) "%For each grid point the following information" &
!!$      //" is provided:"
!!$      write(ngr,*) "%grid point number(index)-land fraction(pcland)-" &
!!$      //"first" &
!!$      //" atmospheric level(iztop)-latitude(lat)-longitude(lon)-" &
!!$      //" surface Horizontal Wind (m/sec)-surface rain rate(rr)"
!!$      write(ngr,*) "%For each vertical level above the considered grid point is provided:"
!!$      write(ngr,*) "%Temp(k)-Pressure(mb)-Watervapor Mixing(g/Kg)- Cloud "//&
!!$      " LWC(g/m^3)-Rain LWC(g/m^3)-Graupel IWC(g/m^3)-Pristine ice IWC(g/"// &
!!$      "m^3)-Snow IWC(g/m^3)-Aggregates IWC(g/m^3)- Vertical Wind (m/sec)"
!!$
!!$
!!$
!!$!      open(14,file='check_qr_bm.dat',status='unknown',form='formatted')
!!$      do i=2,n2-1
!!$         do j=2,n3-1
!!$            pcprl_tot=0.0
!!$            do k=1,ncr
!!$               do l=1,nbr
!!$                  pcprl_tot=pcprl_tot+pcprl(1,l,k,i,j)
!!$               end do
!!$            end do
!!$            pcprl_tot=pcprl_tot*3600.0
!!$            pcpri_tot=0.0
!!$            do k=1,nci
!!$               do l=1,nbi
!!$                  pcpri_tot=pcpri_tot+pcpri(1,l,k,i,j)
!!$               end do
!!$            end do
!!$            pcpri_tot=pcpri_tot*3600.0
!!$
!!$            k=k1eta(i,j)
!!$            uave=(u(k,i,j)+u(k,max(1,i-1),j))*.5
!!$            vave=(v(k,i,j)+v(k,i,max(1,j-1)))*.5
!!$
!!$            write(ngr,1000) j-1,pctlnd(i,j),k1eta(i,j) &
!!$                 ,alatt(i,j),alont(i,j),sqrt(uave**2.0+vave**2.0) &
!!$                 ,pcprl_tot+pcpri_tot
!!$
!!$            do k=1,n1-1
!!$               if(k.ge.k1eta(i,j))then
!!$                  thdry=THETA(k,i,j)/(1.+.61*qv(k,i,j))
!!$                  PIRD=(PP(k,i,j)+PI01D(K))/CPD
!!$                  tv=thdry*PIRD
!!$                  pv=PIRD**CPR*P00*0.01
!!$                  denv=den(k,i,j)
!!$                  qvv=qv(k,i,j)*1.0e+3
!!$!                  qupc(k)=
!!$!                  qrpc(k)=
!!$!                  qagg(k)=
!!$!                  qrag(k)=
!!$!                  qgra(k)=
!!$                  wave=(W(k,i,j)+W(max(1,k-1),i,j))*.5
!!$                  uave=(u(k,i,j)+u(k,max(1,i-1),j))*.5
!!$                  vave=(v(k,i,j)+v(k,i,max(1,j-1)))*.5
!!$
!!$!                 in g/m^3
!!$                  CLWC=qc(k,i,j)*denv*1.0e+3
!!$                  if(ipcnfl.gt.0) then
!!$                     RLWC=qrp(1,1,1,k,i,j)*denv*1.0e+3
!!$                  else
!!$                     RLWC=0.0
!!$                  end if
!!$                  if(ipcnfl.gt.0) then
!!$                     PIWC=qip(1,1,1,k,i,j)*denv*1.0e+3
!!$                  else
!!$                     PIWC=0.0
!!$                  end if
!!$                  if(iscnfl.gt.0) then
!!$                     SIWC=qip(1,1,2,k,i,j)*denv*1.0e+3
!!$                  else
!!$                     SIWC=0.0
!!$                  end if
!!$                  if(iacnfl.gt.0) then
!!$                     AIWC=qip(1,1,3,k,i,j)*denv*1.0e+3
!!$                  else
!!$                     AIWC=0.0
!!$                  end if
!!$                  if(igcnfl.gt.0) then
!!$                     GIWC=qip(1,1,4,k,i,j)*denv*1.0e+3
!!$                  else
!!$                     GIWC=0.0
!!$                  end if
!!$               else
!!$                  tv=var
!!$                  pv=var
!!$                  qvv=var
!!$                  wave=var
!!$                  CLWC=var
!!$                  RLWC=var
!!$                  GIWC=var
!!$                  PIWC=var
!!$                  SIWC=var
!!$                  AIWC=var
!!$               endif
!!$               write(ngr,1001) tv,pv,qvv &
!!$                   ,CLWC,RLWC &
!!$                   ,GIWC,PIWC,SIWC,AIWC &
!!$                   ,wave
!!$            end do
!!$         end do
!!$      end do
!!$ 1000 format(I5,f9.4,I4,4f10.3)
!!$ 1001 format(10F13.4)
!!$      close(ngr)
!!$
!!$      return
!!$      end subroutine print_radcal_bulk

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     print out necessary variables for trajectory analysis
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$      subroutine print_trajana(n1,n2,n3 &
!!$           ,k1eta &
!!$           ,u,v,w &
!!$           ,time,ngr &
!!$           ,npr,nbr,ncr,qrp &
!!$           ,npi,nbi,nci,qip)
!!$      use par_amps
!!$      use maxdims
!!$      use com_amps
!!$      implicit none
!!$
!!$      integer :: n1,n2,n3,k1eta(n2,*)
!!$      integer :: npr,nbr,ncr,npi,nbi,nci,ngr
!!$      real(PS) :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3),time
!!$      real(PS) :: qrp(npr,nbr,ncr,n1,n2,*)  &
!!$          ,qip(npi,nbi,nci,n1,n2,*)
!!$      integer :: i,j,k,l,n,kint,j1,j2,jj1,jj2,ipr,icr,ibr &
!!$              ,ipi,ici,ibi
!!$
!!$!     new space
!!$      real(PS) :: vtmr(n1,n2,n3),vtmi(n1,n2,n3),vtm(n1,n2,n3) &
!!$          ,vtcr(n1,n2,n3),vtci(n1,n2,n3),vtc(n1,n2,n3) &
!!$          ,uave(n1,n2,n3),vave(n1,n2,n3),wave(n1,n2,n3)
!!$      real(PS) :: totmr,totcr,totmi,totci,totm,totc
!!$      real(PS) :: thdry,pird,var,cpr
!!$      integer :: inn
!!$      CHARACTER(len=20) :: numfil
!!$
!!$      var=-999.9
!!$      cpr=cpd/r
!!$
!!$      write(fid_alog,*) "I am writing trajectory file"
!!$
!!$
!!$!     calculate mass and con weighted terminal velocities
!!$      vtmr=0.0
!!$      vtmi=0.0
!!$      vtm=0.0
!!$      vtcr=0.0
!!$      vtci=0.0
!!$      vtc=0.0
!!$      uave=var
!!$      vave=var
!!$      wave=var
!!$      do i=2,n2-1
!!$         do j=2,n3-1
!!$            do k=2,n1-1
!!$               if(k.ge.k1eta(i,j))then
!!$                  wave(k,i,j)=(W(k,i,j)+W(max(1,k-1),i,j))*.5
!!$                  uave(k,i,j)=(u(k,i,j)+u(k,max(1,i-1),j))*.5
!!$                  vave(k,i,j)=(v(k,i,j)+v(k,i,max(1,j-1)))*.5
!!$
!!$                  totmr=0.0
!!$                  totcr=0.0
!!$                  do icr=1,ncr
!!$                     do ibr=1,nbr
!!$                        totmr=totmr+qrp(rmt_q,ibr,icr,k,i,j)
!!$                        totcr=totcr+qrp(rcon_q,ibr,icr,k,i,j)
!!$                        vtmr(k,i,j)=vtmr(k,i,j) &
!!$                    +qrp(npr-1,ibr,icr,k,i,j)*qrp(rmt_q,ibr,icr,k,i,j)
!!$                        vtcr(k,i,j)=vtcr(k,i,j) &
!!$                    +qrp(npr,ibr,icr,k,i,j)*qrp(rcon_q,ibr,icr,k,i,j)
!!$                     end do
!!$                  end do
!!$
!!$                  totmi=0.0
!!$                  totci=0.0
!!$                  do ici=1,nci
!!$                     do ibi=1,nbi
!!$                        totmi=totmi+qip(imt_q,ibi,ici,k,i,j)
!!$                        totci=totci+qip(icon_q,ibi,ici,k,i,j)
!!$                        vtmi(k,i,j)=vtmi(k,i,j) &
!!$                    +qip(npi-1,ibi,ici,k,i,j)*qip(imt_q,ibi,ici,k,i,j)
!!$                        vtci(k,i,j)=vtci(k,i,j) &
!!$                    +qip(npi,ibi,ici,k,i,j)*qip(icon_q,ibi,ici,k,i,j)
!!$                     end do
!!$                  end do
!!$
!!$                  totm=totmr+totmi
!!$                  totc=totcr+totci
!!$                  if(totm>1.0e-30)  &
!!$                          vtm(k,i,j)=(vtmr(k,i,j)+vtmi(k,i,j))/totm
!!$                  if(totc>1.0e-30) &
!!$                          vtc(k,i,j)=(vtcr(k,i,j)+vtci(k,i,j))/totc
!!$
!!$                  if(totmr>1.0e-30) vtmr(k,i,j)=vtmr(k,i,j)/totmr
!!$                  if(totcr>1.0e-30) vtcr(k,i,j)=vtcr(k,i,j)/totcr
!!$
!!$                  if(totmi>1.0e-30) vtmi(k,i,j)=vtmi(k,i,j)/totmi
!!$                  if(totci>1.0e-30) vtci(k,i,j)=vtci(k,i,j)/totci
!!$
!!$               else
!!$                  wave(k,i,j)=var
!!$                  uave(k,i,j)=var
!!$                  vave(k,i,j)=var
!!$               endif
!!$            end do
!!$         end do
!!$      end do
!!$
!!$!     write output file.
!!$      inn=int(time+0.001)
!!$      if(inn.lt.100) then
!!$         write(numfil,1011) inn,ngr
!!$      elseif(inn.lt.1000) then
!!$         write(numfil,1002) inn,ngr
!!$      elseif(inn.lt.10000) then
!!$         write(numfil,1003) inn,ngr
!!$      elseif(inn.lt.100000) then
!!$         write(numfil,1004) inn,ngr
!!$      else
!!$         write(numfil,1005) inn,ngr
!!$      endif
!!$ 1011 format('traj',i2,'.',i1)
!!$ 1002 format('traj',i3,'.',i1)
!!$ 1003 format('traj',i4,'.',i1)
!!$ 1004 format('traj',i5,'.',i1)
!!$ 1005 format('traj',i6,'.',i1)
!!$!      open(ngr,file=numfil,status='unknown',form='formatted')
!!$      open(ngr,file=trim(numfil),status='unknown',form='unformatted')
!!$      write(ngr)  time,real(n1,PS_KIND),real(n2,PS_KIND),real(n3,PS_KIND) &
!!$                  ,((v(k,2,j),k=1,n1),j=1,n3) &
!!$                  ,((w(k,2,j),k=1,n1),j=1,n3) &
!!$                  ,((vtmr(k,2,j),k=1,n1),j=1,n3) &
!!$                  ,((vtcr(k,2,j),k=1,n1),j=1,n3) &
!!$                  ,((vtmi(k,2,j),k=1,n1),j=1,n3) &
!!$                  ,((vtci(k,2,j),k=1,n1),j=1,n3) &
!!$                  ,((vtm(k,2,j),k=1,n1),j=1,n3) &
!!$                  ,((vtc(k,2,j),k=1,n1),j=1,n3)
!!$
!!$      close(ngr)
!!$      return
!!$      end subroutine print_trajana



      real(PS) function get_len_c2a(amass)
!     mass-dimension relationship of bullet roseete
      implicit none
      real(PS) :: amass
!     parameters for mass-dimension relationship from Chiruta and Wang (2003)
      real(PS) :: a_ros,b_ros,N_ros,den_i
      parameter( a_ros=0.3257,b_ros=0.5206,N_ros=4.0,den_i = 0.91668) ! sense test ros original N_ros=4.0, changed to 2.0
      get_len_c2a=(amass/(den_i*a_ros*N_ros**b_ros))**(1.0/3.0)
      end function get_len_c2a

      real(PS) function get_len_s1(amass)
!     mass-dimension relationship of S1
      implicit none
      real(PS) :: amass
!     side plane crystal mass-length relation from Mitchel et al (1990)
      real(PS) :: a,b
      parameter(a=0.0041900509,b=2.3)
!     assume that the shape is cylinder that has the aspect ratio of 0.25 for all sizes, so
!     when it exceed max density, radiuas has to be increased.
!      real den_max,c1
!      parameter(den_max=1.5*0.25/(1.0+0.25**2)**1.5*0.91668
!     *     ,c1=(1.0/(den_max*4.0*3.14159265358979/3.0))**(1.0/3.0)
!     *     /(1.0+0.25**2)**0.5 )
      real(PS) :: c1
      parameter(c1=0.885565704066706)
      get_len_s1=max(0.5*(amass/a)**(1.0/b),c1*amass**(1.0/3.0))
      end function get_len_s1

      real(PS) function get_len_s3(amass)
!     mass-dimension relationship of S3
      implicit none
      real(PS) :: amass
!     S3 crystal mass-length relation from Mitchel et al (1990)
      real(PS) :: a,b
      parameter(a=0.0027696359,b=2.1)
!     assume that the shape is cylinder that has the aspect ratio of 0.25 for all sizes, so
!     when it exceed max density, radiuas has to be increased.
!      real den_max,c1
!      parameter(den_max=1.5*0.25/(1.0+0.25**2)**1.5*0.91668
!     *     ,c1=(1.0/(den_max*4.0*3.14159265358979/3.0))**(1.0/3.0)
!     *     /(1.0+0.25**2)**0.5)
      real(PS) :: c1
      parameter(c1=0.885565704066706)
      get_len_s3=max(0.5*(amass/a)**(1.0/b),c1*amass**(1.0/3.0))
      end function get_len_s3

!!$      subroutine cal_mtac(growth_mode_hex,tmp,mass,alen,clen)
!!$      use com_amps
!!$      implicit none
!!$      real(PS) :: tmp,mass,alen,clen
!!$      integer :: growth_mode_hex,iswitch
!!$      real(PS) :: dtmp,sttmp,edtmp,tmpc,mass_log
!!$      parameter(dtmp=1.0_RP,sttmp=-50.0_RP,edtmp=-1.0_RP)
!!$      real(PS) :: dmass,stmass,edmass
!!$      parameter(dmass=0.1_RP,stmass=-12.0_RP,edmass=-3.0_RP)
!!$      integer :: itmpmin,itmpmax
!!$      parameter(itmpmin=1,itmpmax=50)
!!$      integer :: imassmin,imassmax
!!$      parameter(imassmin=1,imassmax=91)
!!$      integer :: itmp,imass
!!$      real(PS) :: t1,m1,w1,w2,den,phi
!!$      tmpc=tmp-273.16
!!$      mass_log=log10(mass)
!!$
!!$      itmp=max(itmpmin,min(int((tmpc-sttmp)/dtmp)+1,itmpmax-1))
!!$
!!$
!!$      if(growth_mode_hex==1) then
!!$         mass_log=min(min(lmt_mass_pla(itmp),lmt_mass_pla(itmp+1)) &
!!$                 ,mass_log)
!!$
!!$         imass=max(imassmin,min(int(abs(mass_log-stmass)/dmass)+1 &
!!$                ,min(i_lmt_mass_pla(itmp),i_lmt_mass_pla(itmp+1))-1))
!!$
!!$         t1=sttmp+real(itmp-1,PS_KIND)*dtmp
!!$         m1=stmass+real(imass-1,PS_KIND)*dmass
!!$
!!$         w1=max(0.0_RP,min(1.0_RP,(tmpc-t1)/dtmp))
!!$         w2=max(0.0_RP,min(1.0_RP,(mass_log-m1)/dmass))
!!$
!!$         iswitch=1
!!$         phi=10.0_RP**( &
!!$              (1.0_RP-w1)*(1.0_RP-w2)*mtac_map_pla(itmp,imass,iswitch)+ &
!!$              (1.0_RP-w1)*w2*mtac_map_pla(itmp,imass+1,iswitch)+ &
!!$              w1*(1.0_RP-w2)*mtac_map_pla(itmp+1,imass,iswitch)+ &
!!$              w1*w2*mtac_map_pla(itmp+1,imass+1,iswitch) &
!!$              )
!!$         iswitch=2
!!$         den=10.0_RP**( &
!!$              (1.0_RP-w1)*(1.0_RP-w2)*mtac_map_pla(itmp,imass,iswitch)+ &
!!$              (1.0_RP-w1)*w2*mtac_map_pla(itmp,imass+1,iswitch)+ &
!!$              w1*(1.0_RP-w2)*mtac_map_pla(itmp+1,imass,iswitch)+ &
!!$              w1*w2*mtac_map_pla(itmp+1,imass+1,iswitch) &
!!$              )
!!$      elseif(growth_mode_hex==2) then
!!$         mass_log=min(min(lmt_mass_col(itmp),lmt_mass_col(itmp+1)) &
!!$                 ,mass_log)
!!$
!!$         imass=max(imassmin,min(int(abs(mass_log-stmass)/dmass)+1 &
!!$                ,min(i_lmt_mass_col(itmp),i_lmt_mass_col(itmp+1))-1))
!!$
!!$!         imass=max(imassmin,min(int(abs(mass_log-stmass)/dmass)+1
!!$!     *          ,imassmax-1))
!!$
!!$         t1=sttmp+real(itmp-1,PS_KIND)*dtmp
!!$         m1=stmass+real(imass-1,PS_KIND)*dmass
!!$
!!$         w1=max(0.0_RP,min(1.0_RP,(tmpc-t1)/dtmp))
!!$         w2=max(0.0_RP,min(1.0_RP,(mass_log-m1)/dmass))
!!$
!!$         iswitch=1
!!$         phi=10.0_RP**( &
!!$              (1.0_RP-w1)*(1.0_RP-w2)*mtac_map_col(itmp,imass,iswitch)+ &
!!$              (1.0_RP-w1)*w2*mtac_map_col(itmp,imass+1,iswitch)+ &
!!$              w1*(1.0_RP-w2)*mtac_map_col(itmp+1,imass,iswitch)+ &
!!$              w1*w2*mtac_map_col(itmp+1,imass+1,iswitch) &
!!$              )
!!$         iswitch=2
!!$         den=10.0_RP**( &
!!$              (1.0_RP-w1)*(1.0_RP-w2)*mtac_map_col(itmp,imass,iswitch)+ &
!!$              (1.0_RP-w1)*w2*mtac_map_col(itmp,imass+1,iswitch)+ &
!!$              w1*(1.0_RP-w2)*mtac_map_col(itmp+1,imass,iswitch)+ &
!!$              w1*w2*mtac_map_col(itmp+1,imass+1,iswitch) &
!!$              )
!!$      end if
!!$
!!$
!!$      alen=(mass/(den*4.188790205_RP*((1.0+phi**2)**1.5_RP)))**0.33333333_RP
!!$      clen=alen*phi
!!$
!!$!!c      write(fid_alog,'("check mtac:gmh,tmpc,mass,mlog,den,phi,plog,alen,clen"
!!$!!c     *        ,I5,f6.1,9ES15.6)')
!!$!!c     *growth_mode_hex,tmpc,mass,mass_log,den,phi,log10(phi),alen,clen
!!$!!c      write(fid_alog,'("phi",2I5,8ES15.6)') itmp,imass
!!$!!c     *     ,mtac_map_pla(itmp,imass,1),mtac_map_pla(itmp+1,imass,1)
!!$!!c     *     ,mtac_map_pla(itmp,imass+1,1),mtac_map_pla(itmp+1,imass+1,1)
!!$!!c     *     ,mtac_map_col(itmp,imass,1),mtac_map_col(itmp+1,imass,1)
!!$!!c     *     ,mtac_map_col(itmp,imass+1,1),mtac_map_col(itmp+1,imass+1,1)
!!$!!c      write(fid_alog,'("den",2I5,8ES15.6)') itmp,imass
!!$!!c     *     ,mtac_map_pla(itmp,imass,2),mtac_map_pla(itmp+1,imass,2)
!!$!!c     *     ,mtac_map_pla(itmp,imass+1,2),mtac_map_pla(itmp+1,imass+1,2)
!!$!!c     *     ,mtac_map_col(itmp,imass,2),mtac_map_col(itmp+1,imass,2)
!!$!!c     *     ,mtac_map_col(itmp,imass+1,2),mtac_map_col(itmp+1,imass+1,2)
!!$      return
!!$      end subroutine cal_mtac

      subroutine cal_indx_lmtmass()
      use com_amps
      implicit none
      real(PS) :: dtmp,sttmp,edtmp,tmpc,mass_log
      parameter(dtmp=1.0_RP,sttmp=-50.0_RP,edtmp=-1.0_RP)
      real(PS) :: dmass,stmass,edmass
      parameter(dmass=0.1_RP,stmass=-12.0_RP,edmass=-3.0_RP)
      integer :: itmpmin,itmpmax
      parameter(itmpmin=1,itmpmax=50)
      integer :: imassmin,imassmax
      parameter(imassmin=1,imassmax=91)
      integer :: itmp,imass
      real(PS) :: t1,m1,w1,w2,den,phi

      if ( IsMaster .and. debug ) then
        write(fid_alog,*) "pla"
      end if
      do itmp=itmpmin,itmpmax
         i_lmt_mass_pla(itmp)=max(imassmin &
               ,min(int(abs(lmt_mass_pla(itmp)-stmass)/dmass)+1 &
               ,imassmax-1))
         lmt_mass_pla(itmp)=stmass+real(i_lmt_mass_pla(itmp)-1,PS_KIND)*dmass
         if ( IsMaster .and. debug ) then
           write(fid_alog,*) itmp,i_lmt_mass_pla(itmp),lmt_mass_pla(itmp)
         end if
      end do
      if ( IsMaster .and. debug ) then
        write(fid_alog,*) "col"
      end if
      do itmp=itmpmin,itmpmax
         i_lmt_mass_col(itmp)=max(imassmin &
               ,min(int(abs(lmt_mass_col(itmp)-stmass)/dmass)+1 &
               ,imassmax-1))
         lmt_mass_col(itmp)=stmass+real(i_lmt_mass_col(itmp)-1,PS_KIND)*dmass
         if ( IsMaster .and. debug ) then
           write(fid_alog,*) itmp,i_lmt_mass_col(itmp),lmt_mass_col(itmp)
         end if
      end do
      return
      end subroutine cal_indx_lmtmass

      subroutine read_seed(time)
      use com_amps
      implicit none
      real(MP_KIND) :: time,rv(2)
      integer :: i,inn,issz
      CHARACTER(len=13) :: numfil

!!!      if(time<0.0_RP) then
!!!         call random_seed(size=issz)
!!!      write(*,*) "size of seed",issz,size(seed)

!!!         call random_number(rv)
!!!         seed=int(rv*100000);
!!!         open(31,file='seed_0m')
!!!         write(*,*) "writing seeds for random generator to: ",'seed_0m'
!!!         do i=1,sizeseed
!!!            write(31,*) seed(i)
!!!            write(*,*) seed(i)
!!!         end do
!!!         close(31)
!!!         stop
!!      else
         inn=int(time/60)
         if(inn.lt.10) then
            write(numfil,1010) inn
         elseif(inn.lt.100) then
            write(numfil,1011) inn
         elseif(inn.lt.1000) then
            write(numfil,1002) inn
         else
            write(numfil,1003) inn
         endif
 1010 format('seedini_',i1,'m')
 1011 format('seedini_',i2,'m')
 1002 format('seedini_',i3,'m')
 1003 format('seedini_',i4,'m')
         open(31,file=numfil,action='READ')
! CHIARUI >>>>>>>>>> 8/4/2019
         !write(*,*) "number of seeds for random generator (before): ", sizeseed
         sizeseed = 2
         !write(*,*) "number of seeds for random generator (after): ", sizeseed
! CHIARUI <<<<<<<<<<
         !write(*,*) "reading seeds for random generator from: ",numfil
         do i=1,sizeseed
            read(31,*) seed(i)
            !write(*,*) seed(i)
         end do
         close(31)
!tmp      end if
      call random_seed(size=issz)
      !write(*,*) "size of seed",issz,size(seed)
      call random_seed(put=seed)

      end subroutine read_seed

!!$      subroutine write_seed(time)
!!$      use com_amps
!!$      implicit none
!!$      real(PS) :: time
!!$      integer :: i,inn
!!$      CHARACTER(len=10) :: numfil
!!$      inn=int(time/60.0_RP)
!!$      if(inn.lt.100) then
!!$         write(numfil,1011) inn
!!$      elseif(inn.lt.1000) then
!!$         write(numfil,1002) inn
!!$      else
!!$         write(numfil,1003) inn
!!$      endif
!!$ 1011 format('seed_',i2,'m')
!!$ 1002 format('seed_',i3,'m')
!!$ 1003 format('seed_',i4,'m')
!!$
!!$      call random_seed(get=seed_old)
!!$      open(31,file=numfil)
!!$      write(*,*) "writing seeds for random generator to: ",numfil
!!$      do i=1,sizeseed
!!$         write(31,*) seed_old(i)
!!$         write(*,*) seed_old(i)
!!$      end do
!!$      close(31)
!!$      write(*,*) "starting seeds were"
!!$      do i=1,sizeseed
!!$         write(*,*) seed(i)
!!$      end do
!!$      end subroutine write_seed


      subroutine assign_seeds(nsect,seed_sec)
      use com_amps
      implicit none
      integer :: nsect
      integer :: seed_sec(*)
      real(PS),allocatable :: rv(:)
      integer :: i,j,var_status

!      write(fid_alog,*) "sizeseed,nsect",sizeseed,nsect
!      allocate(seed_sec(sizeseed,nsect),stat=var_status)
!      if(var_status.ne.0) then
!        write(fid_alog,*) "seed_sec not allocated",var_status
!        stop
!      endif

      allocate(rv(nsect))

!!c      write(fid_alog,*) "asign_seeds"
!!c      do i=1,sizeseed
!!c        write(fid_alog,*) seed(i)
!!c      enddo
      call random_seed(put=seed)

      call random_number(rv)
      do i=1,nsect
        seed_sec(i)=int(rv(i)*1.0e+9)
      enddo
      deallocate(rv)

!  newly assign seed

      allocate(rv(sizeseed))

      call random_number(rv)
      do i=1,sizeseed
        seed(i)=int(rv(i)*1.0e+9)
      enddo

!!c      write(fid_alog,*) "next asign_seeds"
!!c      do i=1,sizeseed
!!c        write(fid_alog,*) seed(i)
!!c      enddo

      deallocate(rv)

      end subroutine assign_seeds

!!$      subroutine delete_seeds(nsect,seed_sec)
!!$      use com_amps
!!$      implicit none
!!$      integer nsect
!!$      integer seed_sec(sizeseed,*)
!!$
!!$      seed(1:sizeseed)=seed_sec(1:sizeseed,nsect)
!!$!tmp      write(fid_alog,'("delete_seeds",50I15)') seed(1:sizeseed)
!!$
!!$      end subroutine delete_seeds

!!$      subroutine set_seedpriv(seed_priv,isect_seed,isect,seed_sec,is)
!!$      use com_amps
!!$      implicit none
!!$      integer seed_sec(sizeseed,*),isect_seed
!!$      integer seed_priv(*),isect,is
!!$
!!$      if(is.eq.0) then
!!$        isect_seed=isect
!!$        seed_priv(1:sizeseed)=seed_sec(1:sizeseed,isect)
!!$      else
!!$        seed_sec(1:sizeseed,isect)=seed_priv(1:sizeseed)
!!$      endif
!!$      end subroutine set_seedpriv

! added for sheba
!!$      subroutine cal_tvol(tvol,nzp,nxp,nyp,nfpt,nfpth &
!!$                         ,k1eta,zs,zz,xx,yy)
!!$      implicit none
!!$      integer :: nzp,nxp,nyp,k1eta(nxp,*),nfpt,nfpth
!!$      real(PS) :: tvol,zs(nxp,*),zz(*),xx(*),yy(*)
!!$      integer :: k,i,j
!!$      tvol=0.0
!!$!      do j=nfpth+1,nyp-nfpth
!!$      do j=4,nyp-4
!!$        do i=2,nxp-1
!!$          do k=k1eta(i,j),nzp-nfpt
!!$            tvol=tvol+(xx(i)-xx(i-1))*(yy(j)-yy(j-1)) &
!!$                       *(zz(k)-max(zs(i,j),zz(k-1)))
!!$          enddo
!!$        enddo
!!$      enddo
!!$      write(*,*) "Tvol:",tvol,nfpt,nfpth
!!$      return
!!$      end subroutine cal_tvol
!uwnms      subroutine cal_tvol(tvol,nzp,nxp,nyp,nfpt,nfpth &
!uwnms                         ,k1eta,zs,zz,xx,yy)
!uwnms      use module_mpinms, only: &
!uwnms                  mpirank_e,mpirank_w,mpirank_s,mpirank_n &
!uwnms                 ,mpi_get_sum,mpi_get_max0
!uwnms      use com_amps, only: fid_alog
!uwnms      implicit none
!uwnms      real tvol,zs(nxp,*),zz(*),xx(*),yy(*)
!uwnms      integer nzp,nxp,nyp,k1eta(nxp,*),nfpt,nfpth
!uwnms      integer k,i,j,i1str,i2end,j1str,j2end
!uwnms      real,dimension(1) :: x1d
!uwnms!
!uwnms      tvol=0.0
!uwnms      if(mpirank_e.eq.-9) then
!uwnms        i2end=nxp-4
!uwnms      else
!uwnms        i2end=nxp-1
!uwnms      endif
!uwnms      if(mpirank_w.eq.-9) then
!uwnms        i1str=4
!uwnms      else
!uwnms        i1str=2
!uwnms      endif
!uwnms      if(mpirank_n.eq.-9) then
!uwnms        j2end=nyp-4
!uwnms      else
!uwnms        j2end=nyp-1
!uwnms      endif
!uwnms      if(mpirank_s.eq.-9) then
!uwnms        j1str=4
!uwnms      else
!uwnms        j1str=2
!uwnms      endif
!uwnms      if(nxp==3) then
!uwnms        i1str=2
!uwnms        i2end=2
!uwnms      endif
!uwnms      do j=j1str,j2end
!uwnms        do i=i1str,i2end
!uwnms          do k=k1eta(i,j),nzp-nfpt
!uwnms            tvol=tvol+(xx(i)-xx(i-1))*(yy(j)-yy(j-1)) &
!uwnms                       *(zz(k)-max(zs(i,j),zz(k-1)))
!uwnms          enddo
!uwnms        enddo
!uwnms      enddo
!uwnms      x1d(1)=tvol
!uwnms      call mpi_get_sum(x1d,1)
!uwnms      tvol=x1d(1)
!uwnms      write(fid_alog,*) "Tvol:",tvol,nfpt,nfpth
!uwnms      return
!uwnms      end subroutine cal_tvol
! added for end sheba
!
!     The following codes are taken from Park and Miller (1988)
!           www.cisl.ucar.edu/zine/96/spring/articles/3.random-6.html
!
!!!      SUBROUTINE SRAND2(rdsd,ISEED,isect)
      SUBROUTINE SRAND2(jseed,ifrst,isect_seed,ISEED,isect)
!
!  This subroutine sets the integer seed to be used with the
!  companion RAND function to the value of ISEED.  A flag is
!  set to indicate that the sequence of pseudo-random numbers
!  for the specified seed should start from the beginning.
!
!!!      type(random_genvar),intent(inout) :: rdsd
      integer,intent(in) :: iseed,isect
      integer JSEED,IFRST,isect_seed
!!!      integer JSEED,IFRST,isect_seed,NEXTN
!!!      COMMON /IRNGPRIVA/JSEED,IFRST,isect_seed,NEXTN
!!!!$omp threadprivate(/IRNGPRIVA/)
!
      JSEED = ISEED
      IFRST = 0
      isect_seed=isect
!      rdsd%JSEED = ISEED
!      rdsd%IFRST = 0
!      rdsd%isect_seed=isect
!
      RETURN
      END SUBROUTINE SRAND2

!!$      REAL(PS) FUNCTION RAND2()
!!$!
!!$!  This function returns a pseudo-random number for each invocation.
!!$!  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal
!!$!  standard number generator whose Pascal code appears in the article:
!!$!
!!$!     Park, Steven K. and Miller, Keith W., "Random Number Generators:
!!$!     Good Ones are Hard to Find", Communications of the ACM,
!!$!     October, 1988.
!!$!
!!$      integer :: mplier,modlus,mobymp,momdmp
!!$      PARAMETER (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773, &
!!$                 MOMDMP=2836)
!!$!
!!$      integer :: JSEED,IFRST,isect_seed,NEXTN
!!$!!!      COMMON /IRNGPRIVA/JSEED,IFRST,isect_seed,NEXTN
!!$!!!!$omp threadprivate(/IRNGPRIVA/)
!!$      INTEGER :: HVLUE, LVLUE, TESTV
!!$!
!!$      IF (IFRST .EQ. 0) THEN
!!$        NEXTN = JSEED
!!$        IFRST = 1
!!$      ENDIF
!!$!
!!$      HVLUE = NEXTN / MOBYMP
!!$      LVLUE = MOD(NEXTN, MOBYMP)
!!$      TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
!!$      IF (TESTV .GT. 0) THEN
!!$        NEXTN = TESTV
!!$      ELSE
!!$        NEXTN = TESTV + MODLUS
!!$      ENDIF
!!$      RAND2 = REAL(NEXTN,PS_KIND)/REAL(MODLUS,PS_KIND)
!!$!
!!$      RETURN
!!$      END function rand2

      subroutine RAND2_TY(rnv,rdsd)
!
!  This function returns a pseudo-random number for each invocation.
!  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal
!  standard number generator whose Pascal code appears in the article:
!
!     Park, Steven K. and Miller, Keith W., "Random Number Generators:
!     Good Ones are Hard to Find", Communications of the ACM,
!     October, 1988.
!
      real(PS),intent(inout) :: rnv
      type(random_genvar),intent(inout) :: rdsd
!
      integer :: mplier,modlus,mobymp,momdmp
      PARAMETER (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773, &
                 MOMDMP=2836)
!
      INTEGER :: HVLUE, LVLUE, TESTV
!
      IF (rdsd%IFRST .EQ. 0) THEN
        rdsd%NEXTN = rdsd%JSEED
        rdsd%IFRST = 1
      ENDIF
!
      HVLUE = rdsd%NEXTN / MOBYMP
      LVLUE = MOD(rdsd%NEXTN, MOBYMP)
      TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
!
      IF (TESTV .GT. 0) THEN
        rdsd%NEXTN = TESTV
      ELSE
        rdsd%NEXTN = TESTV + MODLUS
      ENDIF
      rnv = REAL(rdsd%NEXTN,PS_KIND)/REAL(MODLUS,PS_KIND)
!
      RETURN
      END subroutine rand2_TY

      SUBROUTINE QSPARM2(ESTBAR,ESITBAR)
!        based on Murphy and Koop (2005)
      implicit none
      real(DS), DIMENSION(150) ::  ESTBAR
      real(DS), DIMENSION(111) ::  ESITBAR
      real(DS) :: T
      integer :: K, JD
      T=163.
      DO K=1,111
        T=T+1.
        ESITBAR(K)=exp(9.550426-5723.265/T+3.53068*log(T)-0.00728332*T)
      ENDDO

      T=163.
      DO JD=1,150
        T=T+1.
        ESTBAR(JD)=exp(54.842763-6763.22/T-4.210*log(T)+0.000367*T &
                  + tanh(0.0415*(T-218.8))*(53.878-1331.22/T &
                  - 9.44523*log(T)+0.014025*T))
      end do
      RETURN
      END subroutine QSPARM2

!!$      subroutine snowamt(snow,prcp,t,tsoil,isntyp)
!!$      real(PS),intent(inout) :: snow
!!$      real(PS),intent(in) :: prcp,t,tsoil
!!$      real(PS) :: ratio
!!$      character *(*) :: isntyp
!!$!
!!$!      write(fid_alog,*) "snowamt>isntyp",isntyp
!!$!
!!$      if(prcp.eq.0.0_RP) then
!!$        ratio=0.
!!$      else
!!$        ratio=max(0.0_RP,min(30.0_RP,5.0_RP*sqrt(1.2_RP*(max(0.0_RP,(276.16_RP-t))))))
!!$!       ratio=max(8.0_RP,min(30.0_RP,5.0_RP*sqrt(1.2_RP*(max(0.0_RP,(276.16_RP-t))))))
!!$      endif
!!$      if(isntyp.eq.'RAIN') then
!!$        ratio=0.
!!$      elseif(isntyp.eq.'GRAUP') then
!!$        ratio=1.8
!!$      elseif(isntyp.eq.'AGG'.and.t.lt.273.16_RP) then
!!$        ratio=max(ratio,20.0_RP)
!!$      elseif(isntyp.eq.'SNOW'.and.t.lt.273.16_RP) then
!!$        ratio=max(ratio,10.0_RP)
!!$      elseif(isntyp.eq.'PRIS'.and.t.lt.273.16_RP) then
!!$        ratio=max(ratio,10.0_RP)
!!$      endif
!!$      ratio=ratio/25.4
!!$      snow=max(0.0_RP,prcp*ratio)
!!$!     if(tsoil.gt.273.16) then
!!$!       snow=snow*min(1.0_RP,max(0.0_RP,(274_RP-tsoil)/0.84_RP))
!!$!     endif
!!$      return
!!$      end subroutine snowamt

      SUBROUTINE VELFILL(VT,N1,Z,dn01d)
      integer,intent(in) :: n1
      real(PS),DIMENSION(*) :: VT
      real(RP),DIMENSION(*) :: Z,dn01d
      integer :: k,kk,kup
      real(PS) :: dz1
      k=1
5     continue
      IF(VT(K).GE.0.) THEN
        DO KK=K+1,N1
          IF(VT(KK).LT.0.) THEN
            KUP=KK
            GO TO 20
          ENDIF
        ENDDO
        if(k.eq.1) then
          do kk=k+1,n1
            vt(kk)=0.
          enddo
        else
          do kk=k,n1
            vt(kk)=vt(k-1)*sqrt(dn01d(k-1)/dn01d(kk))
          enddo
        endif
        return
20      CONTINUE
        IF(k.eq.1) THEN
          DO KK=1,KUP-1
            VT(KK)=VT(KUP)*sqrt(dn01d(kUP)/dn01d(kk))
          ENDDO
        ELSE
          DO KK=k,KUP-1
            DZ1=(Z(KK)-Z(K-1))/(Z(KUP)-Z(K-1))
            VT(KK)=(1.-DZ1)*VT(K-1) &
                        +DZ1*VT(KUP)
          ENDDO
        ENDIF
        K=KUP+1
      else
        k=k+1
      ENDIF
      IF(K.le.N1) GO TO 5
      RETURN
      END subroutine velfill

      subroutine findcond1(k1m,k2m,k1c,k2c,k1r,k2r,k1i,k2i,micptr,qr,qi &
                ,n1,qc,k1eta,ncr,nci,i,j,level,kbnd,idir,from)
      implicit none
! arguments
      integer :: micptr(*),k1eta(*),n1,level,kbnd,idir &
               ,k1c,k2c,k1r,k2r,k1i,k2i,ncr,nci &
               ,k1m,k2m,i,j
      real(MP_KIND) :: qr(*),qi(*),qc(*)
      character *(*) :: from
! new space
      integer :: npts,ncycle,l,ks,ke,k
!
!      write(fid_alog,*) "findcond1 (1)",i,j,omp_in_parallel()
!     *          ,omp_get_num_threads()
!     *          ,omp_get_thread_num(),loc(k),loc(qr)

      k1c=10000
      k2c=-10000
      k1r=10000
      k2r=-10000
      k1i=10000
      k2i=-10000

      if(kbnd.eq.4) then
        npts=3
        ncycle=npts*2+1
        l=n1-ncycle
!        ks=k1eta+npts-1
        ks=2+npts-1
        ke=n1-1-npts
!        write(fid_alog,*) "I am here",kbnd,ks,ke,npts,ncycle,l
      else
        ks=1
        if(idir.eq.3) ks=k1eta(1)
!ccc        ks=k1eta
!ccc        ke=n1-1
        ke=n1
      end if
      if(level.ge.3) then
        do k=ks,ke

!tmp          if(qi(k).ge.1.0e-1) then
!tmp             write(fid_alog,'("fd1>something wrong qi:",a,3I5,10ES15.6)')  &
!tmp               from,k,i,j,qc(k),qr(k),qi(k)
!tmp          endif
!tmp          if(qr(k).ge.1.0e-1) then
!tmp             write(fid_alog,'("fd1>something wrong qr:",a,3I5,10ES15.6)')  &
!tmp               from,k,i,j,qc(k),qr(k),qi(k)
!tmp          endif

          if(idir.ne.3) then
            if(j.lt.k1eta(k)) then
              micptr(k)=0
              cycle
            endif
          end if
          micptr(k)=0

          if(qr(k).ne.0..and.ncr.gt.0) then
            k1r=min(k1r,k)
            k2r=max(k2r,k)
            micptr(k)=1
          endif
          if(qi(k).ne.0..and.nci.gt.0) then
            k1i=min(k1i,k)
            k2i=max(k2i,k)
            micptr(k)=1
          endif
          if(qc(k).gt.0.) then
            micptr(k)=1
            k1c=min(k1c,k)
            k2c=max(k2c,k)
            if(level.eq.5) then
               k1r=min(k1r,k)
               k2r=max(k2r,k)
            end if
          endif
        enddo
        if(kbnd.eq.4) then
          do k=2,ks-1
            if(idir.ne.3) then
              if(j.lt.k1eta(k)) then
                micptr(k)=0
                cycle
              endif
            end if
            micptr(k)=0
            if(qr(k+l).ne.0..and.ncr.gt.0) then
              k1r=min(k1r,k)
              k2r=max(k2r,k)
              micptr(k)=1
            endif
            if(qi(k+l).ne.0..and.nci.gt.0) then
              k1i=min(k1i,k)
              k2i=max(k2i,k)
              micptr(k)=1
            endif
            if(qc(k+l).gt.0.) then
              micptr(k)=1
              k1c=min(k1c,k)
              k2c=max(k2c,k)
              if(level.eq.5) then
                 k1r=min(k1r,k)
                 k2r=max(k2r,k)
              end if
            endif
          end do
          do k=ke+1,n1-1
            if(idir.ne.3) then
              if(j.lt.k1eta(k)) then
                micptr(k)=0
                cycle
              endif
            end if
            micptr(k)=0
            if(qr(k-l-npts).ne.0..and.ncr.gt.0) then
              k1r=min(k1r,k)
              k2r=max(k2r,k)
              micptr(k)=1
            endif
            if(qi(k-l-npts).ne.0..and.nci.gt.0) then
              k1i=min(k1i,k)
              k2i=max(k2i,k)
              micptr(k)=1
            endif
            if(qc(k-l-npts).gt.0.) then
              micptr(k)=1
              k1c=min(k1c,k)
              k2c=max(k2c,k)
              if(level.eq.5) then
                 k1r=min(k1r,k)
                 k2r=max(k2r,k)
              end if
             endif
           end do
         end if
         k1m=min(k1c,k1r,k1i)
         k2m=max(k2c,k2r,k2i)
      elseif(level.eq.2) then
        do k=ks,ke
          if(idir.ne.3) then
            if(j.lt.k1eta(k)) then
!7b              micptr(k)=0
              cycle
            endif
          end if
!7b          micptr(k)=0

          if(qc(k).gt.0.) then
!7b            micptr(k)=1
            k1c=min(k1c,k)
            k2c=max(k2c,k)
          endif
        enddo
        if(kbnd.eq.4) then
          do k=2,ks-1
            if(idir.ne.3) then
              if(j.lt.k1eta(k)) then
!7b                micptr(k)=0
                cycle
              endif
            end if
!7b            micptr(k)=0
            if(qc(k+l).gt.0.) then
!7b              micptr(k)=1
              k1c=min(k1c,k)
              k2c=max(k2c,k)
            endif
          end do
          do k=ke+1,n1-1
            if(idir.ne.3) then
              if(j.lt.k1eta(k)) then
!7b                micptr(k)=0
                cycle
                endif
              end if
!7b              micptr(k)=0
            if(qc(k-l-npts).gt.0.) then
!7b              micptr(k)=1
              k1c=min(k1c,k)
              k2c=max(k2c,k)
            endif
          end do
        end if
        k1m=min(k1c,k1r,k1i)
        k2m=max(k2c,k2r,k2i)
      end if
!     if(i.eq.ival.and.j.eq.jval) then
!     print 1,k1m,k2m,k1c,k2c,k1r,k2r,k1i,k2i
1     format('k1m,k2m,k1c,k2c,k1r,k2r,k1i,k2i',4(2i4,2x))
!     endif

      return
      end subroutine findcond1

      subroutine findcond2(k1m,k2m,k1c,k2c,k1r,k2r,k1i,k2i,micptr &
                ,k1br,k2br,k1bi,k2bi &
                ,qr,qi &
                ,n1,npr,nbr,ncr,npi,nbi,nci,level,nbhzcl &
                ,qc,qrp,qip,qtp,qro,qio &
                ,k1eta,ival,jval,from,idir)
      use par_amps
      use com_amps, only: fid_alog
      implicit none
! arguments
      integer :: micptr(*),k1eta(*),n1,level &
               ,npr,nbr,ncr,npi,nbi,nci,nbhzcl &
               ,k1c,k2c,k1r,k2r,k1i,k2i &
               ,k1m,k2m,ival,jval,idir &
               ,k1br(nbr,*),k2br(nbr,*),k1bi(nbi,*),k2bi(nbi,*)
      real(MP_KIND) :: qr(*),qi(*),qc(*),qrp(npr,nbr,ncr,*),qip(npi,nbi,nci,*) &
           ,qtp(*),qro(nbr,ncr,*),qio(nbi,nci,*)
      character *(*) :: from
! new space
      integer :: k,ipr,ibr,icr,ipi,ibi,ici &
             ,k1,ibr_st,iii
!
! NOTE: these limit parameters are very important. This makes the model
!       slower or faster, and also affect predicted fields!
!
      real(PS) :: RRLMT,RILMT,RRLMTB,RILMTB
      parameter(RRLMT=1.0d-10,RILMT=1.0d-15)
      parameter(RRLMTB=1.0d-22,RILMTB=1.0d-22)
!      parameter (qrmin=1.e-8,qimin=1.e-13)
!

!     Find limits of the microphysics and also remove
!     residual mnd insignificant icrophysics
      k1c=10000
      k2c=-10000
      k1r=10000
      k2r=-10000
      k1i=10000
      k2i=-10000

      k1m=10000
      k2m=-10000

      do icr=1,ncr
        do ibr=1,nbr
          k1br(ibr,icr)=10000
          k2br(ibr,icr)=-10000
        end do
      end do
      do ici=1,nci
        do ibi=1,nbi
          k1bi(ibi,ici)=10000
          k2bi(ibi,ici)=-10000
        end do
      end do

!ccc      k1=2
      k1=1
      if(idir.eq.3) k1=k1eta(1)
      if(level.eq.3) then
         ibr_st=1

!ccc         do k=k1,n1-1
         do k=k1,n1
           qr(k)=0.
           do ibr=ibr_st,nbr
             do icr=1,ncr
               if(qrp(rmt_q,ibr,icr,k).eq.0.0_RP.and. &
                  qro(ibr,icr,k).eq.0.0_RP) then
!
!                 do nothing (the qrpv is not copied back to shared memory)
!
                 cycle
               elseif((qrp(rmt_q,ibr,icr,k).lt.RRLMT).or. &
                   (idir.lt.3.and.jval.lt.k1eta(k)).or. &
                   (qtp(k).le.0.0_RP) &
                  ) then
                 do ipr=1,npr
                   qrp(ipr,ibr,icr,k)=0.0_RP
                 enddo

                 k1r=min(k1r,k)
                 k2r=max(k2r,k)
                 k1br(ibr,icr)=min(k1br(ibr,icr),k)
                 k2br(ibr,icr)=max(k2br(ibr,icr),k)
               elseif(qrp(rmt_q,ibr,icr,k).ge.RRLMT) then
                 k1r=min(k1r,k)
                 k2r=max(k2r,k)
                 qr(k)=qr(k)+qrp(rmt_q,ibr,icr,k)
                 micptr(k)=1

                 k1br(ibr,icr)=min(k1br(ibr,icr),k)
                 k2br(ibr,icr)=max(k2br(ibr,icr),k)
               endif
             enddo
!             if(qrp(rcon_q,ibr,1,k)>1.0e+8) then
!                write(fid_alog,*) "find2 from qrp2",from,qrp(rcon_q,ibr,1,k)
!             end if
           enddo
           if(debug.and.qr(k).gt.1.e-1) then
             print 1,from,nbr,ncr,ival,jval,k,k1 &
                     ,qr(k),qi(k),qtp(k),qc(k)
1            format('large qr in findcond2 called from: ',a15, &
                      ' nbr,ncr,ival,jval,k,k1eta',6i5 &
                      ,' qr(k),qi(k),qtp(k),qc(k) ',3p4f10.3)
             write(fid_alog,*) "qrp1",((qrp(rmt_q,ibr,icr,k),ibr=1,nbr),icr=1,ncr)
!             stop
           endif

           qi(k)=0.
           do ibi=1,nbi
             do ici=1,nci
               if(qip(1,ibi,ici,k).eq.0.0_RP.and. &
                  qio(ibi,ici,k).eq.0.0_RP) then
!
!                 do nothing (the qipv is not copied back to shared memory)
!
                 cycle
               elseif((qip(1,ibi,ici,k).lt.RILMT).or. &
                  (idir.lt.3.and.jval.lt.k1eta(k)).or. &
                  (qtp(k).le.0.0_RP) &
                 ) then
                 do ipi=1,npi
                   qip(ipi,ibi,ici,k)=0.0_RP
                 enddo
                 k1i=min(k1i,k)
                 k2i=max(k2i,k)
                 k1bi(ibi,ici)=min(k1bi(ibi,ici),k)
                 k2bi(ibi,ici)=max(k2bi(ibi,ici),k)
               elseif(qip(1,ibi,ici,k).ge.RILMT) then
                 k1i=min(k1i,k)
                 k2i=max(k2i,k)
                 qi(k)=qi(k)+qip(1,ibi,ici,k)
                 micptr(k)=1

                 k1bi(ibi,ici)=min(k1bi(ibi,ici),k)
                 k2bi(ibi,ici)=max(k2bi(ibi,ici),k)

               endif
             enddo
           enddo

           if(debug.and.qi(k).gt.1.e-1) then
              print 2,from,nbi,nci,ival,jval,k &
                     ,qr(k),qi(k),qtp(k),qc(k)
2             format('large qi in findcond2 called from: ',a15, &
                          ' nbi,nci,ival,jval,k',5i5 &
                      ,' qr(k),qi(k),qtp(k),qc(k) ',3p4f10.3)
!              stop
              write(fid_alog,*) "qip1",((qip(1,ibi,ici,k),ibi=1,nbi),ici=1,nci)
           endif

!           write(fid_alog,'("in fd2:k,i,j,k1r,k2r,k1i,k2i",7I7,10ES15.6)')
!     *              k,ival,jval
!     *             ,k1r,k2r,k1i,k2i,qr(k),qi(k),qtp(k)
!     *             ,((qrp(rmt_q,ibr,icr,k),ibr=1,nbr),icr=1,ncr)
!     *             ,((qip(1,ibi,ici,k),ibi=1,nbi),ici=1,nci)
!
!           if(qi(k).ne.0.0_RP) then
!             write(fid_alog,*) "fd2>something wrong qi: ",from,k,qtp(k)
!     *         ,qc(k),qr(k),qi(k)
!             do ibi=1,nbi
!               do ici=1,nci
!                 write(fid_alog,'(3I5,ES15.6)') k,ibi,ici,qip(1,ibi,ici,k)
!               enddo
!             enddo
!             stop
!           endif
!
!           if(qr(k).ne.0.0_RP) then
!             write(fid_alog,*) "fd2>something wrong qr: ",from,k,qtp(k)
!     *         ,qc(k),qr(k),qi(k)
!             do ibr=1,nbr
!               do icr=1,ncr
!                 write(fid_alog,'(3I5,ES15.6)') k,ibr,icr,qrp(rmt_q,ibr,icr,k)
!               enddo
!             enddo
!             stop
!           endif

           if(qc(k).gt.0.) then
             micptr(k)=1
             k1c=min(k1c,k)
             k2c=max(k2c,k)
           endif
!
!          This is necessary to turn on the saturated boxes without hydrometeors
!
           if(micptr(k).eq.1) then
             k1m=min(k1m,k)
             k2m=max(k2m,k)
           endif
         enddo
      elseif(level.ge.5) then
!        qr starts from nbhzcl+1 bin of qrp for bin microphysics (SHIPS)
!        liquid bin spectrum that uses nbhzcl bins as haze and cloud droplet.
!         write(fid_alog,*) "ck:nbhzcl",from,level,nbhzcl
         ibr_st=nbhzcl+1

!ccc         do k=k1,n1-1
         do k=k1,n1
           qc(k)=0.0_RP
           qr(k)=0.0_RP
           do icr=1,ncr
             do ibr=1,ibr_st-1
               if(qrp(rmt_q,ibr,icr,k).eq.0.0_RP.and. &
                  qro(ibr,icr,k).eq.0.0_RP) then
!
!                 do nothing (the qrpv is not copied back to shared memory)
!
                 cycle
               elseif((qrp(rmt_q,ibr,icr,k).lt.RRLMTB).or. &
                   (idir.lt.3.and.jval.lt.k1eta(k)).or. &
                   (qtp(k).le.0.0_RP) &
                  ) then
                 do ipr=1,npr
                   qrp(ipr,ibr,icr,k)=0.0_RP
                 enddo
!
!            need to update variables that are deleted.
!
                 k1c=min(k1c,k)
                 k2c=max(k2c,k)
                 k1r=min(k1r,k)
                 k2r=max(k2r,k)
                 k1br(ibr,icr)=min(k1br(ibr,icr),k)
                 k2br(ibr,icr)=max(k2br(ibr,icr),k)

               elseif(qrp(rmt_q,ibr,icr,k).ge.RRLMTB) then
                 k1c=min(k1c,k)
                 k2c=max(k2c,k)
                 k1r=min(k1r,k)
                 k2r=max(k2r,k)
                 qc(k)=qc(k) &
                  +max(qrp(rmt_q,ibr,icr,k)-qrp(rmat_q,ibr,icr,k),0.0_RP)
                 micptr(k)=1
                 k1br(ibr,icr)=min(k1br(ibr,icr),k)
                 k2br(ibr,icr)=max(k2br(ibr,icr),k)
               endif
             enddo
             do ibr=ibr_st,nbr
               if(qrp(rmt_q,ibr,icr,k).eq.0.0_RP.and. &
                  qro(ibr,icr,k).eq.0.0_RP) then
!
!                 do nothing (the qrpv is not copied back to shared memory)
!
                 cycle
               elseif((qrp(rmt_q,ibr,icr,k).lt.RRLMTB).or. &
                   (idir.lt.3.and.jval.lt.k1eta(k)).or. &
                   (qtp(k).le.0.0_RP) &
                  ) then
                 do ipr=1,npr
                   qrp(ipr,ibr,icr,k)=0.0_RP
                 enddo
!
!            need to update variables that are deleted.
!
                 k1r=min(k1r,k)
                 k2r=max(k2r,k)
                 k1br(ibr,icr)=min(k1br(ibr,icr),k)
                 k2br(ibr,icr)=max(k2br(ibr,icr),k)


               elseif(qrp(rmt_q,ibr,icr,k).ge.RRLMTB) then
                 k1r=min(k1r,k)
                 k2r=max(k2r,k)
                 qr(k)=qr(k) &
                   +max(qrp(rmt_q,ibr,icr,k)-qrp(rmat_q,ibr,icr,k),0.0_RP)
                 micptr(k)=1
                 k1br(ibr,icr)=min(k1br(ibr,icr),k)
                 k2br(ibr,icr)=max(k2br(ibr,icr),k)
               endif
             enddo
           enddo
           if(debug.and.qr(k).gt.1.e-2) then
             print 1,from,ibr,icr,ival,jval,k,k1 &
                     ,qr(k),qi(k),qtp(k),qc(k)
             write(fid_alog,'("qrp_m",30ES15.6)') (qrp(rmt_q,ibr,1,k),ibr=1,nbr)
             write(fid_alog,'("qrp_c",30ES15.6)') (qrp(rcon_q,ibr,1,k),ibr=1,nbr)
             write(fid_alog,'("qrp_vtm",30ES15.6)') (qrp(npr,ibr,1,k),ibr=1,nbr)
!         stop
           endif
           qi(k)=0.
           do ibi=1,nbi
             do ici=1,nci
               if(qip(imt_q,ibi,ici,k).eq.0.0_RP.and. &
                  qio(ibi,ici,k).eq.0.0_RP) then
!
!                 do nothing (the qipv is not copied back to shared memory)
!
                 cycle
               elseif((qip(imt_q,ibi,ici,k).lt.RILMT).or. &
                   (idir.lt.3.and.jval.lt.k1eta(k)).or. &
                   (qtp(k).le.0.0_RP) &
                 ) then
                 do ipi=1,npi
                   qip(ipi,ibi,ici,k)=0.0_RP
                  enddo
!
!            need to update variables that are deleted.
!
                 k1i=min(k1i,k)
                 k2i=max(k2i,k)
                 k1bi(ibi,ici)=min(k1bi(ibi,ici),k)
                 k2bi(ibi,ici)=max(k2bi(ibi,ici),k)

               elseif(qip(imt_q,ibi,ici,k).ge.RILMT) then
                 k1i=min(k1i,k)
                 k2i=max(k2i,k)
! <<< 2015/03 T. Hashino changed for KiD
!org                 qi(k)=qi(k) &
!org                  +max(0.0,qip(imt_q,ibi,ici,k)-qip(imat_q,ibi,ici,k))

                 qi(k)=qi(k) &
                  +max(0.0_RP,qip(imt_q,ibi,ici,k)-qip(imat_q,ibi,ici,k)-qip(imw_q,ibi,ici,k))
                 qr(k)=qr(k) &
                  +max(0.0_RP,qip(imw_q,ibi,ici,k))
! >>> 2015/03 T. Hashino changed for KiD

                 micptr(k)=1

                 k1bi(ibi,ici)=min(k1bi(ibi,ici),k)
                 k2bi(ibi,ici)=max(k2bi(ibi,ici),k)

                 if(debug.and.qi(k).gt.1.e-2) then
                   print 2,from,ibi,ici,ival,jval,k &
                         ,qr(k),qi(k),qtp(k),qc(k)
!             stop
                 endif
               endif
             enddo
           enddo
!
!          This is necessary to turn on the saturated boxes without hydrometeors
!
           if(micptr(k).eq.1) then
             k1m=min(k1m,k)
             k2m=max(k2m,k)
           endif

         enddo
      elseif(level.eq.2) then
!ccc         do k=k1,n1-1
         do k=k1,n1
           if(qc(k).gt.0.) then
!7b             micptr(k)=1
             k1c=min(k1c,k)
             k2c=max(k2c,k)
           endif
         enddo
      end if
      k1m=min(k1m,k1c,k1r,k1i)
      k2m=max(k2m,k2c,k2r,k2i)
      return
      end subroutine findcond2


      ! <<< 2015/03 T. Hashino changed below
!org      subroutine sclsedprz(fqrp,ql,qtp &
!org           ,fqldiab &
      subroutine sclsedprz_original(ql,qmlt,qc,qtp &
! <<< 2015/03 T. Hashino changed below
           ,ispray &
           ,n1,n2,n3,n4,npr,nbr,ncr,k1r,k2r,k1m,k2m &
           ,k1br,k2br &
           ,qrp,den,wb,wacc &
           ,spdsfc &
           ,theta,qv,thskin,pgnd &
           ,k1eta,z,zz,dzzm,dzvm,dz1 &         ! z center grid, zz half grid
           ,dtl,i,j,CPPM,CPPME,isnow,iadvv &
!           ,ictpc,ictsw,ictag,ictgr &
!           ,ircnfl,ipcnfl,iscnfl,iacnfl,igcnfl &
           ,level,nbhzcl,kmicv,imicv,jmicv &
           ,U &       ! [IN]    u velocity,  cell center
           ,V &       ! [IN]    v velocity,  cell center
           ,TEMP &    ! [IN]    temperature, cell center
           ,dens &    ! [IN]    dry+moist density, cell center
           ,momz &    ! [IN]    z momentum, cell face
           ,CV_WATER &! [IN]    specific heat water
           ,CV_ICE &  ! [IN]    specific heat ice
           ,den_t  &  ! [INOUT] density tendency, cell center
           ,momz_t &  ! [INOUT] z momentum tendency, cell face (added for scale momentum flux output)
           ,rhou_t &  ! [INOUT] rho u velocity tendency, cell center
           ,rhov_t &  ! [INOUT] rho v velocity tendency, cell center
           ,rhoe_t &  ! [INOUT] rho e velocity tendency, cell center
           ,sflx)     ! [INOUT] surface flux
        use scale_prc, only: &
           PRC_abort
        use par_amps
        use maxdims
        use mod_scladv_slppm, only: &
           cal_aLRa6_z, &
           cal_flux_z
        implicit none

      integer,intent(in) :: n1,n2,n3,n4 &
           ,npr,nbr,ncr &
           ,k1br(nbr,*),k2br(nbr,*) &
           ,i,j,k1eta &
           ,isnow,ispray,iadvv,level,nbhzcl
!           ,ictpc,ictsw,ictag,ictgr &
!           ,ircnfl,ipcnfl,iscnfl,iacnfl,igcnfl

      integer,intent(inout) :: k1r,k2r,k1m,k2m
      integer,intent(in) :: imicv(*),jmicv(*),kmicv(*)

      real(RP),intent(in) :: theta(*),qv(*),thskin &
           ,pgnd,den(*),wb(*),wacc(*) &
           ,z(*),zz(*),dzzm(*),dzvm(*) &
           ,spdsfc &
           ,CPPM(11,*),CPPME(11,n2mx,*)
      real(DS),intent(in) :: dtl
      real(RP),intent(inout) :: qrp(npr,nbr,ncr,*)
      real(RP),intent(inout) :: &
           ql(*),qmlt(*),qc(*),qtp(*)

      real(RP), intent(in)     :: dz1

      real(RP), intent(in)    :: momz(n1), U(n1), V(n1), TEMP(n1), dens(n1)
      real(RP), intent(in)    :: CV_WATER, CV_ICE
      real(RP), intent(inout) :: den_t(n1), momz_t(n1), rhou_t(n1), rhov_t(n1), rhoe_t(n1), sflx

! new space
      integer :: k,ka,kb,k1,k2 &
           ,ipr,ibr,icr,niter,iter &
           ,inc,ninc &
           ,isfcset,ibr2,ibr1,kk,nspray
      real(PS) :: aLR(n1,2),a6(n1),dzzz(n1) &
           ,vel(n1,2),vt(n1,2) &
           ,rold(n1),rold_1(n1) &
           ,vel_d(n1),flxsp(n1),fluxf(n1)
      real(PS) :: fluxf_scale(n1), flux(n1), eflx(n1)
      real(PS) :: fluxf_2, rold_2(n1), aLR_2(n1,2),a6_2(n1) ! for scale flux in PPM
      real(PS) :: mass_convert
      integer :: iwv,nwv
      character(len=10) :: isntyp

      real(PS) :: snow &
             ,diab &
             ,tair,tsoil,pignd &
             ,dtnew &
             ,cfl,wt &
             ,dz &
             ,am0,am0p1,am0m1,rm1,rp1,r0,r1,r2,r00,am00,dfmsdr0 &
             ,rinc, s2
      real(PS) :: cpdrd, gcp, coef, denl, rmax, rcld
      integer,dimension(n1) :: ierror1
      integer :: nprx(nprmx)
      integer :: n1mxd,n2mxd,n3mxd
      integer :: istp,ntn
      integer :: jmt_q,jcon_q

! local variables for semi-lagrangian mid-point rule
      real(RP) :: rfdz2(n1), dvel_d(n1), semilag_dist(n1), Z_src
      integer  :: k_dst, k_src(n1)

      cpdrd=cpd/r
      gcp=g/cpd
      coef=4./3.*3.14
      denl=1.e3
      rmax=500.
      rcld=10.

      fluxf_scale = 0.0_RP
      flux = 0.0_RP
      eflx = 0.0_RP
!
!-------------------------------------------------------------------------------------
!     k1 and k2 are the vertical limits of integration determined by the points with
!     non zero total rain (to get started)
!---------------------------------------------------------------------------------------
!
!      return
      do icr=1,ncr
        if(level.eq.3.and.isnow.eq.1) then
!          if(icr.eq.ictpc) then
!            if(ipcnfl.eq.5) then
!              nprx(icr)=3
!            else
!              nprx(icr)=2
!            endif
!          elseif(icr.eq.ictsw) then
!            if(iscnfl.eq.5) then
!              nprx(icr)=3
!            else
!              nprx(icr)=2
!            endif
!          elseif(icr.eq.ictag) then
!            if(iacnfl.eq.5) then
!              nprx(icr)=3
!            else
!              nprx(icr)=2
!            endif
!          elseif(icr.eq.ictgr) then
!            if(igcnfl.eq.5) then
!              nprx(icr)=3
!            else
!              nprx(icr)=2
!            endif
!          endif
!        elseif(level.eq.3) then
!          if(ircnfl.eq.5) then
!            nprx(icr)=3
!          else
!            nprx(icr)=2
!          endif
        else
          nprx(icr)=npr
        endif
      enddo

      n1mxd=n1mx
      n2mxd=n2mx
      n3mxd=n3mx
      istp=0

      if(level.eq.3) then
        nwv=1
        jmt_q=1
        jcon_q=2
      elseif(level.eq.4) then
        nwv=2
        jmt_q=1
        jcon_q=2
      elseif(level.eq.5) then
        nwv=2
        if(isnow.eq.1) then
          jmt_q=imt_q
          jcon_q=icon_q
        else
          jmt_q=rmt_q
          jcon_q=rcon_q
        endif
      end if

      dzzz(k1eta)=dz1
      do k=k1eta+1,n1
        dzzz(k)=1.0/dzzm(k)
      end do

      do k=1,n1
        rold(k)=0.
        rold_1(k)=0.
        flxsp(k)=0.
        vt(k,1)=0.
        vt(k,2)=0.
        vel(k,1)=0.
        vel(k,2)=0.
!tmp         flm(k)=0.
!tmp         fluxa(k)=0.
        fluxf(k)=0.
!
!
!tmp         w(k)=wb(k)
!tmp         if(k.gt.k1eta) then
!tmp            w(k)=wb(k)
!tmp         elseif(k.eq.k1eta) then
!tmp            dz=(zz(k)-zz(k-1))
!tmp            w(k)=max(min(w(k)*dtl/dz,1.),-1.)*dz/dtl
!tmp         else
!tmp            w(k)=0.
!tmp         endif
      enddo

      isfcset=0
      if(ispray.eq.1.and.ibr.eq.1.and.icr.eq.1.and.spdsfc.gt.5.) then
         do ibr=1,nbr
            ibr2=min(nbr,ibr+1)
            ibr1=max(1,ibr-1)
            am0=rainmass(ibr)
            am0p1=rainmass(ibr+1)
            am0m1=rainmass(ibr-1)
            rm1=(3.*am0m1/(4.*3.1415*denl))**.33333
            rp1=(3.*am0p1/(4.*3.1415*denl))**.33333
            r0=(3.*am0/(4.*3.1415*denl))**.33333
            if(ibr.eq.1) then
               r1=0
            else
               r1=0.5*(r0+rm1)
            endif
            if(ibr.eq.nbr) then
               r2=2*r0-r1
            else
               r2=0.5*(r0+rp1)
            endif
            ninc=10
            rinc=(r2-r1)/10.
            flxsp(ibr)=0
            do inc=1,ninc
               r00=(float(inc)-0.5)*rinc
               am00=4./3*3.1415*r0**3*denl
               call sprayms(dfmsdr0,r00,spdsfc)
               flxsp(ibr)=flxsp(ibr)+dfmsdr0*rinc*coef*denl/den(k1eta) &
                    *(r00*1.e-6)**3
            enddo
            if(flxsp(ibr)/(zz(k1eta)-zz(k1eta-1))*dtl.gt.1.e-6)  &
                 then
               nspray=1
               k1r=k1eta
               k2r=max(k1r,k2r)
               k1m=min(k1m,k1r)
               k2m=max(k2r,k2m)
            endif
         enddo
      endif
!
!
      if(k1r.gt.k2r) return

      k1=max(k1r-1,k1eta)
      k2=min(k2r+1,n1-1)
      if(iadvv.eq.6) then
         k1=max(k1eta,k1-2)
         k2=min(n1,k2+2)
      endif

!----------------------------------------------------------------------------------------
!     calculate mass flux
!----------------------------------------------------------------------------------------
      do icr=1,ncr
        do ibr=1,nbr
!
!-------------------------------------------------------------------------------------
!         k1 and k2 are the vertical limits of integration determined by the points with
!         non zero terminal verlocity for this specific bin and category
!---------------------------------------------------------------------------------------
!
!
          k1=max(k1br(ibr,icr)-1,k1eta)
          k2=min(k2br(ibr,icr)+1,n1-1)
          if(k1.gt.k2) cycle

          if(iadvv.eq.6) then
            k1=max(k1eta,k1-2)
            k2=min(n1,k2+2)
          endif
!
!         determine the primary (mass) terminal velocity
!
          do iwv=1,nwv
            do k=k1,k2
              vt(k,iwv)=qrp(npr-nwv+iwv,ibr,icr,k)
              if(qrp(1,ibr,icr,k).ne.0.0.and.vt(k,iwv).le.0.) then
                vt(k,iwv)=sign(min(abs(vt(k,iwv)*wacc(k)),50.0_RP),wacc(k))
              else
                vt(k,iwv)=0.
              endif
            enddo
          enddo
!
!----------------------------------------------------------------------------------
!                 Find  velocity at vertical grid interface points and then
!                 look for cfl violations and how many iterations are necessary to
!                 keep it stable.  This is primarily for mesoscale applications
!                 where vertical resolution is high and we want rain to fall several
!                 vertical grid cells in a single dtl length timestep
!----------------------------------------------------------------------------
!

          cfl=0.
          do k=max(k1,k1eta+1),k2
            cfl=max(real(cfl,DS), &
                max(abs(vt(k,1)),abs(vt(k-1,1)) &
                   ,abs(vt(k,2)),abs(vt(k-1,2)))*dtl &
             /(zz(k)-zz(k-1)))
          enddo
!
!----------------------------------------------------------------------------------------
!     Find how many iterations (small timesteps)  necessary for stability
!----------------------------------------------------------------------------------------
!
          niter=1.+cfl/0.95
!
!----------------------------------------------------------------------------------
!          Now begin loop through all parameters except last one, which is the terminal velocity
!---------------------------------------------------------------------------------
!
!
!           ka and kb are the limits of the vertical integration expanded
!           for multiple iterations
!
          ka=max(k1eta,k1-niter+1)
          kb=min(n1-1,k2+niter-1)

!!!            write(*,'("ck den:",I5,10ES15.6)') k,z(ka-1),den(ka-1)
!tmp            write(*,*) "niter ",niter,ka,kb

          do iwv=1,nwv
            ! this seems to adjust grids with positive (upwards) terminal velocity
            call velfill(vt(ka-1,iwv),kb-ka+3,z(ka-1),den(ka-1))
!
!           now get terminal velocity at "w" point
!
!            do k=ka-1,kb+1
            do k=ka,kb ! CHIARUI BUG
              if(iadvv.eq.12) then
                vel(k,iwv)=0.5*(vt(k,iwv)+vt(k+1,iwv))
              else
                vel(k,iwv)=0.5*(vt(k,iwv)+vt(k+1,iwv))
              endif
            enddo
            if(ka.eq.k1eta) then
              if(iadvv.eq.12) then
                vel(ka-1,iwv)=vt(ka,iwv)
              else
                vel(ka-1,iwv)=vt(ka,iwv)
              endif
            else
              if(iadvv.eq.12) then
                vel(ka-1,iwv)=0.5*(vt(ka-1,iwv)+vt(ka,iwv))
              else
                vel(ka-1,iwv)=0.5*(vt(ka-1,iwv)+vt(ka,iwv))
              endif
            endif
            if(kb.eq.n1-1) then
              if(iadvv.eq.12) then
                vel(kb+1,iwv)=vt(kb+1,iwv)
              else
                vel(kb+1,iwv)=vt(kb+1,iwv)
              endif
            else
              if(iadvv.eq.12) then
                vel(kb+1,iwv)=0.5*(vt(kb+1,iwv)+vt(kb+2,iwv))
              else
                vel(kb+1,iwv)=0.5*(vt(kb+1,iwv)+vt(kb+2,iwv))
              endif
            endif
          end do
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Parameter Loop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
          do ipr=1,nprx(icr)-nwv   ! begin ipr "parameter" loop

!
!           initialize temporary prediction arrays
!
            do k=1,n1
              rold(k)=0.
              vel_d(k)=0.
            enddo

            ! seems not so cheap
            !if ((isnow == 1 .and. ipr <= 8) .or. (isnow /= 1 .and. ipr <= 3)) then
            !  do k = 1, n1
            !    den1(k) = den(k)
            !  enddo
            !else
            !  do k = 1, n1
            !    den1(k) = den(k)*0.001D0
            !  enddo
            !endif

            if(level==3.or. &
              (level==4.and.ipr==1).or.&
              (level==5.and.(&
                  (isnow==1.and. &
                  (ipr==imt_q.or.ipr==imr_q.or.ipr==ima_q.or. &
                   ipr==imc_q.or.ipr==imw_q.or.ipr==imf_q.or. &
                   ipr==imat_q.or.ipr==imas_q)) &
                  .or. &
                  (isnow==0.and. &
                  (ipr==rmt_q.or.ipr==rmat_q.or.ipr==rmas_q)))) ) &
               then
!                 case of mass weighted velocity
              iwv=1
            else
!                 case of con weighted velocity
              iwv=2
            end if

            if(level.eq.3) then
              do k=ka-1,kb+1
                rold(k)=qrp(ipr,ibr,icr,k)
                vel_d(k)=vel(k,1)
                if((ipr.eq.2.and.isnow.eq.0).or. &
                   (ipr.eq.2.and.isnow.eq.1.and.icr.ne.1)) then
                  vel_d(k)=vel_d(k)/2.1
                endif
!             if(vel_d(k).lt.-1..and.isnow.eq.0) then
!               print 116,k,ipr,nprx(icr),iwv,nwv,vel_d(k),rold(k)
!116             format('k,ipr,nx,iwv,nwv,vel ',5i5,f10.1,e15.5)
!             endif
!             print 45,k,rold(k),vel_d(k)
!45            format('k,rold,vel_d',i5,2e15.5)
              enddo
            else
              do k=ka-1,kb+1
                rold(k)=qrp(ipr,ibr,icr,k)
                vel_d(k)=vel(k,iwv)
              enddo
            endif

            advscheme_if: if(iadvv.eq.12) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Semi-Lagrangian PPM
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              dtnew=dtl/float(niter)
!              if(niter.gt.1) then
!                write(*,'("Warning:multiple step in sedimentation",8I5,10ES15.6)')  &
!                     kmicv(k),imicv(k),jmicv(k),ibr,icr,niter &
!                    ,ka,kb,cfl
!                do k=max(k1eta,ka-1),min(n1-1,kb+1)
!                        write(*,'("k,i,j,vel_d",3I5,10ES15.6)')  &
!                         kmicv(k),imicv(k),jmicv(k),vel_d(k) &
!                        ,qrp(jmt_q,ibr,icr,k),qrp(jcon_q,ibr,icr,k) &
!                        ,qrp(jcon_q,ibr,icr,k)/max(qrp(jcon_q,ibr,icr,k),1.0e-25)
!                enddo
!              endif

              do iter=1,niter

!     semi-Lagrangian PPM method (modo Ong Chia Rui)
!     1. calculation of points at left and right interfaces and coefficient.
                call cal_aLRa6_z(i,j,max(k1eta,ka-1),min(n1,kb+1) &
                       ,k1eta,n1,n1mxd,n2mxd,n3mxd &
                       ,n1,istp &
                       ,CPPM,CPPME &
                       ,rold(1),aLR,a6)

                if (isnow == 1 .and. (ipr==imt_q)) then
                  do k=1,n1
                    rold_2(k) = max(0.0_RP,rold(k) - qrp(imw_q,ibr,icr,k) - qrp(imat_q,ibr,icr,k))
                  enddo
                  call cal_aLRa6_z(i,j,max(k1eta,ka-1),min(n1,kb+1) &
                         ,k1eta,n1,n1mxd,n2mxd,n3mxd &
                         ,n1,istp &
                         ,CPPM,CPPME &
                         ,rold_2(1),aLR_2,a6_2)
                elseif (isnow == 1 .and. (ipr==imw_q)) then
                  rold_2 = rold
                  aLR_2 = aLR
                  a6_2 = a6
                elseif (isnow /= 1 .and. (ipr==rmt_q)) then
                  do k=1,n1
                    rold_2(k) = max(0.0_RP,rold(k) - qrp(rmat_q,ibr,icr,k))
                  enddo
                  call cal_aLRa6_z(i,j,max(k1eta,ka-1),min(n1,kb+1) &
                         ,k1eta,n1,n1mxd,n2mxd,n3mxd &
                         ,n1,istp &
                         ,CPPM,CPPME &
                         ,rold_2(1),aLR_2,a6_2)
                endif

!     2. calculation of flux
                if (isnow == 1) then
                  if (ipr <= 8) then
                    mass_convert = 1.0D0
                  else
                    mass_convert = 0.001D0
                  endif
                else
                  if (ipr <= 3) then
                    mass_convert = 1.0D0
                  else
                    mass_convert = 0.001D0
                  endif
                endif
                do k=max(k1eta,ka-1),min(n1-1,kb)
!
!                 k is the interface position
!
                  if(vel_d(k).gt.0.0_RP) then
                    s2=(zz(k)-(zz(k)-vel_d(k)*dtnew))/dzzz(k)
                    fluxf(k)= &
                      den(k)*(aLR(k,2)-0.5*s2*(aLR(k,2)-aLR(k,1) &
                         -(1.0-2.0*s2/3.0)*a6(k))) &
                     *vel_d(k)*mass_convert

                    if (isnow == 1 .and. (ipr==imt_q)) then
                      fluxf_2= &
                        den(k)*(aLR_2(k,2)-0.5*s2*(aLR_2(k,2)-aLR_2(k,1) &
                           -(1.0-2.0*s2/3.0)*a6_2(k)))/float(niter) &
                       *vel_d(k)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_ICE &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    elseif (isnow == 1 .and. (ipr==imw_q)) then
                      fluxf_2= &
                        den(k)*(aLR_2(k,2)-0.5*s2*(aLR_2(k,2)-aLR_2(k,1) &
                           -(1.0-2.0*s2/3.0)*a6_2(k)))/float(niter) &
                       *vel_d(k)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_WATER &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    elseif (isnow /= 1 .and. ipr==rmt_q) then
                      fluxf_2= &
                        den(k)*(aLR_2(k,2)-0.5*s2*(aLR_2(k,2)-aLR_2(k,1) &
                           -(1.0-2.0*s2/3.0)*a6_2(k)))/float(niter) &
                       *vel_d(k)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_WATER &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    endif
                    LOG_ERROR("sclsedprz_original",*) "TVEL>0"
                    call PRC_abort
                  else
                    s2=((zz(k)-vel_d(k)*dtnew)-zz(k))/dzzz(k+1)
                    fluxf(k)= &
                      den(k+1)*(aLR(k+1,1)+0.5*s2*(aLR(k+1,2)-aLR(k+1,1) &
                          +(1.0-2.0*s2/3.0)*a6(k+1))) &
                     *vel_d(k)*mass_convert

                    if (isnow == 1 .and. (ipr==imt_q)) then
                      fluxf_2= &
                        den(k+1)*(aLR_2(k+1,2)-0.5*s2*(aLR_2(k+1,2)-aLR_2(k+1,1) &
                           -(1.0-2.0*s2/3.0)*a6_2(k+1)))/float(niter) &
                       *vel_d(k)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_ICE &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    elseif (isnow == 1 .and. (ipr==imw_q)) then
                      fluxf_2= &
                        den(k+1)*(aLR_2(k+1,2)-0.5*s2*(aLR_2(k+1,2)-aLR_2(k+1,1) &
                           -(1.0-2.0*s2/3.0)*a6_2(k+1)))/float(niter) &
                       *vel_d(k)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_WATER &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    elseif (isnow /= 1 .and. ipr==rmt_q) then
                      fluxf_2= &
                        den(k+1)*(aLR_2(k+1,2)-0.5*s2*(aLR_2(k+1,2)-aLR_2(k+1,1) &
                           -(1.0-2.0*s2/3.0)*a6_2(k)))/float(niter) &
                       *vel_d(k)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_WATER &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    endif
                  endif
                end do
                fluxf(k1eta-1)=fluxf(k1eta)
                if (ka-1 < k1eta) then
                  k = ka - 1
                  s2=((zz(k+1)-vel_d(k+1)*dtnew)-zz(k+1))/dzzz(k+2)
                  if (isnow == 1 .and. (ipr==imt_q)) then
                      fluxf_2= &
                        den(k+2)*(aLR_2(k+2,2)-0.5*s2*(aLR_2(k+2,2)-aLR_2(k+2,1) &
                           -(1.0-2.0*s2/3.0)*a6_2(k+2)))/float(niter) &
                       *vel_d(k+1)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_ICE &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    elseif (isnow == 1 .and. (ipr==imw_q)) then
                      fluxf_2= &
                        den(k+2)*(aLR_2(k+2,2)-0.5*s2*(aLR_2(k+2,2)-aLR_2(k+2,1) &
                           -(1.0-2.0*s2/3.0)*a6_2(k+2)))/float(niter) &
                       *vel_d(k+1)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_WATER &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    elseif (isnow /= 1 .and. ipr==rmt_q) then
                      fluxf_2= &
                        den(k+2)*(aLR_2(k+2,2)-0.5*s2*(aLR_2(k+2,2)-aLR_2(k+2,1) &
                           -(1.0-2.0*s2/3.0)*a6_2(k+1)))/float(niter) &
                       *vel_d(k+1)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_WATER &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    endif
                endif


!
!------------------------------------------------------------------------
!                 Precip Rate and accumulation Rate for no iteration case
!------------------------------------------------------------------------
!
                if(ka.eq.k1eta) then
                  if(ipr.eq.1) then
                     ! scale surface flux
                     sflx = sflx + fluxf(k1eta-1)/float(niter)
                  endif
!
!     surface spray flux
!
                  if(ipr.eq.1.and.icr.eq.1.and.ispray.eq.1) then
                    fluxf(k1eta-1)=fluxf(k1eta-1)+flxsp(ibr)
                  endif
                endif
!
!     predict settling now!
!
                do k=ka,kb
                  if(ka.eq.k1eta) then
                    rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                      /((z(k)-z(k-1))*den(k)*mass_convert)
                  else
                    rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                      /((zz(k)-zz(k-1))*den(k)*mass_convert)
                  endif
                enddo

              enddo

            else if(iadvv.eq.22) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Semi-Lagrangian mid-point rule
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Multiple loops
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

              dtnew=dtl/float(niter)
              if(niter.gt.1) then
                !write(*,*) "niter is more than 1",niter,dtnew
                do iter=1,niter

                  ! ka-1 is KS-1, for boundary value, v(ka-1) = v(ka)
                  ! kb is KE
                  ! CDZ = 1 / dzzm
                  do k=ka,kb-1
                    rfdz2(k) = 1.0_RP / (1.0_RP/dzzm(k) + 1.0_RP/dzzm(k+1))
                    dvel_d(k) = vel_d(k) - vel_d(k-1)
                  enddo
                  dvel_d(kb) = vel_d(kb) - vel_d(kb-1)
                  !if (ka-1 == k1eta) dvel_d(k) = 0.0_RP

                  ! find backward distance
                  do k=ka,kb-1
                    semilag_dist(k) = - vel_d(k)*dtnew &
                                      + vel_d(k)*dtnew**2/2.0_RP  * (dvel_d(k) + dvel_d(k+1))*rfdz2(k) &
                                      - vel_d(k)*dtnew**3/4.0_RP  * (dvel_d(k+1)*dzvm(k+1) - dvel_d(k)*dzvm(k))*rfdz2(k)
                    semilag_dist(k) = max(semilag_dist(k),0.0_RP)
                  enddo
                  semilag_dist(ka-1) = - vel_d(ka-1)*dtnew &
                                       + vel_d(ka-1)*dtnew**2/2.0_RP  * dvel_d(ka)*dzzm(ka)
                  semilag_dist(ka-1) = max(semilag_dist(ka-1),0.0_RP)

                  ! not allow overlap of backward distance
                  do k=kb-2,ka,-1
                    semilag_dist(k) = min(semilag_dist(k),semilag_dist(k+1)+1.0_RP/dzzm(k+1))
                  enddo

                  ! search backward cell number
                  do k_dst=ka-1,kb-1
                    Z_src = zz(k_dst) + semilag_dist(k_dst)

                    k_src(k_dst) = k_dst
                    do k=k_dst,kb-1
                      if (Z_src > zz(k) .and. Z_src <= zz(k+1)) then
                        k_src(k_dst) = k
                      endif
                    enddo
                    if (Z_src > zz(kb)) k_src(k_dst) = kb
                  enddo

                  if (isnow == 1) then
                    if (ipr <= 8) then
                      mass_convert = 1.0D0
                    else
                      mass_convert = 0.001D0
                    endif
                  else
                    if (ipr <= 3) then
                      mass_convert = 1.0D0
                    else
                      mass_convert = 0.001D0
                    endif
                  endif

                  fluxf = 0.0_RP

                  do k_dst=ka-1,kb-1
                    if(vel_d(k_dst).gt.0.0_RP) then
                       LOG_ERROR("sclsedprz_original",*) "TVEL>0"
                       call PRC_abort
                    endif

                    do k=k_dst,k_src(k_dst)-1
                      fluxf(k_dst)=fluxf(k_dst)-rold(k+1)*den(k+1)*mass_convert/dzzm(k+1)/dtnew
                      semilag_dist(k) = semilag_dist(k) - 1.0_RP/ dzzm(k)

                      if (isnow == 1 .and. (ipr==imt_q)) then
                        fluxf_2 = -max(0.0_RP,rold(k+1) - qrp(imw_q,ibr,icr,k+1) - qrp(imat_q,ibr,icr,k+1))*den(k+1)/float(niter)/dzzm(k+1)/dtnew
                        fluxf_scale(k_dst) = fluxf_scale(k_dst) + fluxf_2
                        eflx(k_dst) = eflx(k_dst) + fluxf_2*TEMP(k+1)*CV_ICE &
                                + fluxf_2/dzvm(k)*g
                      elseif (isnow == 1 .and. ipr==imw_q) then
                        fluxf_2 = -rold(k+1)*den(k+1)/float(niter)/dzzm(k+1)/dtnew
                        fluxf_scale(k_dst) = fluxf_scale(k_dst) + fluxf_2
                        eflx(k_dst) = eflx(k_dst) + fluxf_2*TEMP(k+1)*CV_WATER &
                                + fluxf_2/dzvm(k)*g
                      elseif (isnow /= 1 .and. ipr==rmt_q) then
                        fluxf_2 = -max(0.0_RP,rold(k+1) - qrp(rmat_q,ibr,icr,k+1))*den(k+1)/float(niter)/dzzm(k+1)/dtnew
                        fluxf_scale(k_dst) = fluxf_scale(k_dst) + fluxf_2
                        eflx(k_dst) = eflx(k_dst) + fluxf_2*TEMP(k+1)*CV_WATER &
                                + fluxf_2/dzvm(k)*g
                      endif
                    enddo

                    if (k_src(k_dst) < kb) then
                      fluxf(k_dst)=fluxf(k_dst)-rold(k_src(k_dst)+1)*den(k_src(k_dst)+1)*mass_convert*semilag_dist(k_dst)/dtnew
                      if (isnow == 1 .and. (ipr==imt_q)) then
                        fluxf_2 = -max(0.0_RP,rold(k_src(k_dst)+1) - qrp(imw_q,ibr,icr,k_src(k_dst)+1) - qrp(imat_q,ibr,icr,k_src(k_dst)+1))*den(k_src(k_dst)+1)/float(niter)*semilag_dist(k_dst)/dtnew
                        fluxf_scale(k_dst) = fluxf_scale(k_dst) + fluxf_2
                        eflx(k_dst) = eflx(k_dst) + fluxf_2*TEMP(k_src(k_dst)+1)*CV_ICE &
                                + fluxf_2/dzvm(k_src(k_dst))*g
                      elseif (isnow == 1 .and. ipr==imw_q) then
                        fluxf_2 = -rold(k_src(k_dst)+1)*den(k_src(k_dst)+1)/float(niter)*semilag_dist(k_dst)/dtnew
                        fluxf_scale(k_dst) = fluxf_scale(k_dst) + fluxf_2
                        eflx(k_dst) = eflx(k_dst) + fluxf_2*TEMP(k_src(k_dst)+1)*CV_WATER &
                                + fluxf_2/dzvm(k_src(k_dst))*g
                      elseif (isnow /= 1 .and. ipr==rmt_q) then
                        fluxf_2 = -max(0.0_RP,rold(k_src(k_dst)+1) - qrp(rmat_q,ibr,icr,k_src(k_dst)+1))*den(k_src(k_dst)+1)/float(niter)*semilag_dist(k_dst)/dtnew
                        fluxf_scale(k_dst) = fluxf_scale(k_dst) + fluxf_2
                        eflx(k_dst) = eflx(k_dst) + fluxf_2*TEMP(k_src(k_dst)+1)*CV_WATER &
                                + fluxf_2/dzvm(k_src(k_dst))*g
                      endif
                    endif
                  enddo
!
!
!
!
!------------------------------------------------------------------------
!                 Precip Rate and accumulation Rate for iteration case
!------------------------------------------------------------------------
!
                  if(ka.eq.k1eta) then
                    if(ipr.eq.1) then
                        ! scale surface flux
                        sflx = sflx + fluxf(k1eta-1)/float(niter)
                      endif
!
!     surface spray flux
!
                    if(ipr.eq.1.and.icr.eq.1.and.ispray.eq.1) then
                      fluxf(k1eta-1)=fluxf(k1eta-1)+flxsp(ibr)
                    endif
                  endif
!
!     predict settling now!
!
                  do k=ka,kb
                    if(k.eq.k1eta) then
                      rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                        /((z(k)-z(k-1))*den(k)*mass_convert)
                      !rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                      !  /((z(k)-z(k-1))*den(k))
                    else
                      rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                        /((zz(k)-zz(k-1))*den(k)*mass_convert)
                      !rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                      !  /((zz(k)-zz(k-1))*den(k))
                    endif
                  enddo
                enddo

              else

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Single loops
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


                dvel_d = 0.0_RP
                rfdz2 = 0.0_RP
                semilag_dist = 0.0_RP
                k_src = 0
                ! ka-1 is KS-1, for boundary value, v(ka-1) = v(ka)
                ! kb is KE
                ! CDZ = 1 / dzzm
                do k=ka,kb-1
                  rfdz2(k) = 1.0_RP / (1.0_RP/dzzm(k) + 1.0_RP/dzzm(k+1))
                  dvel_d(k) = vel_d(k) - vel_d(k-1)
                enddo
                dvel_d(kb) = vel_d(kb) - vel_d(kb-1)
                !if (ka-1 == k1eta) dvel_d(k) = 0.0_RP

                ! find backward distance
                do k=ka,kb-1
                  semilag_dist(k) = - vel_d(k)*dtl &
                                    + vel_d(k)*dtl**2/2.0_RP  * (dvel_d(k) + dvel_d(k+1))*rfdz2(k) &
                                    - vel_d(k)*dtl**3/4.0_RP  * (dvel_d(k+1)*dzvm(k+1) - dvel_d(k)*dzvm(k))*rfdz2(k)
                  semilag_dist(k) = max(semilag_dist(k),0.0_RP)
                enddo
                semilag_dist(ka-1) = - vel_d(ka-1)*dtl &
                                     + vel_d(ka-1)*dtl**2/2.0_RP  * dvel_d(ka)*dzzm(ka)
                semilag_dist(ka-1) = max(semilag_dist(ka-1),0.0_RP)

                ! not allow overlap of backward distance
                do k=kb-2,ka,-1
                  semilag_dist(k) = min(semilag_dist(k),semilag_dist(k+1)+1.0_RP/dzzm(k+1))
                enddo

                ! search backward cell number
                do k_dst=ka-1,kb-1
                  Z_src = zz(k_dst) + semilag_dist(k_dst)

                  k_src(k_dst) = k_dst
                  do k=k_dst,kb-1
                    if (Z_src > zz(k) .and. Z_src <= zz(k+1)) then
                      k_src(k_dst) = k
                    endif
                  enddo
                  if (Z_src > zz(kb)) k_src(k_dst) = kb
                enddo

                if (isnow == 1) then
                  if (ipr <= 8) then
                    mass_convert = 1.0D0
                  else
                    mass_convert = 0.001D0
                  endif
                else
                  if (ipr <= 3) then
                    mass_convert = 1.0D0
                  else
                    mass_convert = 0.001D0
                  endif
                endif

                fluxf = 0.0_RP

                do k_dst=ka-1,kb-1
                  if(vel_d(k_dst).gt.0.0_RP) then
                     LOG_ERROR("sclsedprz_original",*) "TVEL>0"
                     call PRC_abort
                  endif

                  do k=k_dst,k_src(k_dst)-1
                    fluxf(k_dst)=fluxf(k_dst)-rold(k+1)*den(k+1)*mass_convert/dzzm(k+1)/dtl
                    semilag_dist(k) = semilag_dist(k) - 1.0_RP/ dzzm(k)

                    if (isnow == 1 .and. (ipr==imt_q)) then
                      fluxf_2 = -max(0.0_RP,rold(k+1) - qrp(imw_q,ibr,icr,k+1) - qrp(imat_q,ibr,icr,k+1))*den(k+1)/dzzm(k+1)/dtl
                      fluxf_scale(k_dst) = fluxf_scale(k_dst) + fluxf_2
                      eflx(k_dst) = eflx(k_dst) + fluxf_2*TEMP(k+1)*CV_ICE &
                                  + fluxf_2/dzvm(k)*g
                    elseif (isnow == 1 .and. ipr==imw_q) then
                      fluxf_2 = -rold(k+1)*den(k+1)/dzzm(k+1)/dtl
                      fluxf_scale(k_dst) = fluxf_scale(k_dst) + fluxf_2
                      eflx(k_dst) = eflx(k_dst) + fluxf_2*TEMP(k+1)*CV_WATER &
                                  + fluxf_2/dzvm(k)*g
                    elseif (isnow /= 1 .and. ipr==rmt_q) then
                      fluxf_2 = -max(0.0_RP,rold(k+1) - qrp(rmat_q,ibr,icr,k+1))*den(k+1)/dzzm(k+1)/dtl
                      fluxf_scale(k_dst) = fluxf_scale(k_dst) + fluxf_2
                      eflx(k_dst) = eflx(k_dst) + fluxf_2*TEMP(k+1)*CV_WATER &
                                  + fluxf_2/dzvm(k)*g
                    endif
                  enddo

                  if (k_src(k_dst) < kb) then
                    fluxf(k_dst)=fluxf(k_dst)-rold(k_src(k_dst)+1)*den(k_src(k_dst)+1)*mass_convert*semilag_dist(k_dst)/dtl
                    if (isnow == 1 .and. (ipr==imt_q)) then
                      fluxf_2 = -max(0.0_RP,rold(k_src(k_dst)+1) - qrp(imw_q,ibr,icr,k_src(k_dst)+1) - qrp(imat_q,ibr,icr,k_src(k_dst)+1))*den(k_src(k_dst)+1)*semilag_dist(k_dst)/dtl
                      fluxf_scale(k_dst) = fluxf_scale(k_dst) + fluxf_2
                      eflx(k_dst) = eflx(k_dst) + fluxf_2*TEMP(k_src(k_dst)+1)*CV_ICE &
                                  + fluxf_2/dzvm(k_src(k_dst))*g
                    elseif (isnow == 1 .and. ipr==imw_q) then
                      fluxf_2 = -rold(k_src(k_dst)+1)*den(k_src(k_dst)+1)*semilag_dist(k_dst)/dtl
                      fluxf_scale(k_dst) = fluxf_scale(k_dst) + fluxf_2
                      eflx(k_dst) = eflx(k_dst) + fluxf_2*TEMP(k_src(k_dst)+1)*CV_WATER &
                                  + fluxf_2/dzvm(k_src(k_dst))*g
                    elseif (isnow /= 1 .and. ipr==rmt_q) then
                      fluxf_2 = -max(0.0_RP,rold(k_src(k_dst)+1) - qrp(rmat_q,ibr,icr,k_src(k_dst)+1))*den(k_src(k_dst)+1)*semilag_dist(k_dst)/dtl
                      fluxf_scale(k_dst) = fluxf_scale(k_dst) + fluxf_2
                      eflx(k_dst) = eflx(k_dst) + fluxf_2*TEMP(k_src(k_dst)+1)*CV_WATER &
                                  + fluxf_2/dzvm(k_src(k_dst))*g
                    endif
                  endif
                enddo

                if(ka.eq.k1eta) then
!
!------------------------------------------------------------------------
!                 Precip Rate and accumulation Rate for no iteration case
!------------------------------------------------------------------------
!
                  if(ipr.eq.1) then
                      ! scale surface flux
                      sflx = sflx + fluxf(k1eta-1)
                  endif
!
!     surface spray flux
!
                  if(ipr.eq.1.and.icr.eq.1.and.ispray.eq.1) then
                    fluxf(k1eta-1)=fluxf(k1eta-1)+flxsp(ibr)
                  endif
                endif
!
!     predict settling now!
!
                do k=ka,kb
                  if(k.eq.k1eta) then
                    rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                      /((z(k)-z(k-1))*den(k)*mass_convert)
                    !rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                    !  /((z(k)-z(k-1))*den(k))
                  else
                    rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                      /((zz(k)-zz(k-1))*den(k)*mass_convert)
                    !rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                    !  /((zz(k)-zz(k-1))*den(k))
                  endif
                enddo

              endif

            else
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Upstream linear
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Multiple iterations
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
              dtnew=dtl/float(niter)
              if(niter.gt.1) then
                !write(*,*) "niter is more than 1",niter,dtnew
                do iter=1,niter

!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!                 Forward - Upstream precip setteling over niter iterations
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
                  if (isnow == 1) then
                    if (ipr <= 8) then
                      mass_convert = 1.0D0
                    else
                      mass_convert = 0.001D0
                    endif
                  else
                    if (ipr <= 3) then
                      mass_convert = 1.0D0
                    else
                      mass_convert = 0.001D0
                    endif
                  endif
                  do k=ka-1,kb
                    if(vel_d(k).gt.0.0_RP) then
                      fluxf(k)=vel_d(k)*rold(k)*den(k)*mass_convert

                      if (isnow == 1 .and. (ipr==imt_q)) then
                        fluxf_2 = vel_d(k)*max(0.0_RP,rold(k) - qrp(imw_q,ibr,icr,k) - qrp(imat_q,ibr,icr,k))*den(k)/float(niter)
                        fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                        eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_ICE &
                                + fluxf_2/dzvm(k)*g
                      elseif (isnow == 1 .and. (ipr==imw_q)) then
                        fluxf_2 = vel_d(k)*rold(k)*den(k)/float(niter)
                        fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                        eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_WATER &
                                + fluxf_2/dzvm(k)*g
                      elseif (isnow /= 1 .and. ipr==rmt_q) then
                        fluxf_2 = vel_d(k)*max(0.0_RP,rold(k) - qrp(rmat_q,ibr,icr,k))*den(k)/float(niter)
                        fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                        eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_WATER &
                                + fluxf_2/dzvm(k)*g
                      endif
                      LOG_ERROR("sclsedprz_original",*) "TVEL>0"
                      call PRC_abort
                    else
                      fluxf(k)=vel_d(k)*rold(k+1)*den(k+1)*mass_convert

                      if (isnow == 1 .and. (ipr==imt_q)) then
                        fluxf_2 = vel_d(k)*max(0.0_RP,rold(k+1) - qrp(imw_q,ibr,icr,k+1) - qrp(imat_q,ibr,icr,k+1))*den(k+1)/float(niter)
                        fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                        eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_ICE &
                                + fluxf_2/dzvm(k)*g
                      elseif (isnow == 1 .and. ipr==imw_q) then
                        fluxf_2 = vel_d(k)*rold(k+1)*den(k+1)/float(niter)
                        fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                        eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_WATER &
                                + fluxf_2/dzvm(k)*g
                      elseif (isnow /= 1 .and. ipr==rmt_q) then
                        fluxf_2 = vel_d(k)*max(0.0_RP,rold(k+1) - qrp(rmat_q,ibr,icr,k+1))*den(k+1)/float(niter)
                        fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                        eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_WATER &
                                + fluxf_2/dzvm(k)*g
                      endif
                    endif
                  enddo
!
!
!
!
!------------------------------------------------------------------------
!                 Precip Rate and accumulation Rate for iteration case
!------------------------------------------------------------------------
!
                  if(ka.eq.k1eta) then
                    if(ipr.eq.1) then
                        ! scale surface flux
                        sflx = sflx + fluxf(k1eta-1)/float(niter)
                    endif
!
!     surface spray flux
!
                    if(ipr.eq.1.and.icr.eq.1.and.ispray.eq.1) then
                      fluxf(k1eta-1)=fluxf(k1eta-1)+flxsp(ibr)
                    endif
                  endif
!
!     predict settling now!
!
                  do k=ka,kb
                    if(k.eq.k1eta) then
                      rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                        /((z(k)-z(k-1))*den(k)*mass_convert)
                      !rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                      !  /((z(k)-z(k-1))*den(k))
                    else
                      rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                        /((zz(k)-zz(k-1))*den(k)*mass_convert)
                      !rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                      !  /((zz(k)-zz(k-1))*den(k))
                    endif
                  enddo
                enddo
              else
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Single Iteration
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!                 Forward - Upstream vertical advection just 1 iteration
!------------------------------------------------------------------------
!------------------------------------------------------------------------
                if (isnow == 1) then
                  if (ipr <= 8) then
                    mass_convert = 1.0D0
                  else
                    mass_convert = 0.001D0
                  endif
                else
                  if (ipr <= 3) then
                    mass_convert = 1.0D0
                  else
                    mass_convert = 0.001D0
                  endif
                endif
                do k=ka-1,kb
                  if(vel_d(k).gt.0.0_RP) then
                    fluxf(k)=vel_d(k)*rold(k)*den(k)*mass_convert

                    if (isnow == 1 .and. (ipr==imt_q)) then
                      fluxf_2 = vel_d(k)*max(0.0_RP,rold(k) - qrp(imw_q,ibr,icr,k) - qrp(imat_q,ibr,icr,k))*den(k)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_ICE &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    elseif (isnow == 1 .and. (ipr==imw_q)) then
                      fluxf_2 = vel_d(k)*rold(k)*den(k)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_WATER &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    elseif (isnow /= 1 .and. ipr==rmt_q) then
                      fluxf_2 = vel_d(k)*max(0.0_RP,rold(k) - qrp(rmat_q,ibr,icr,k))*den(k)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_WATER &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    endif
                    LOG_ERROR("sclsedprz_original",*) "TVEL>0"
                    call PRC_abort
                  else
                    fluxf(k)=vel_d(k)*rold(k+1)*den(k+1)*mass_convert

                    if (isnow == 1 .and. (ipr==imt_q)) then
                      fluxf_2 = vel_d(k)*max(0.0_RP,rold(k+1) - qrp(imw_q,ibr,icr,k+1) - qrp(imat_q,ibr,icr,k+1))*den(k+1)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_ICE &
                              + fluxf_2/dzvm(k)*g           ! potential energy
                    elseif (isnow == 1 .and. (ipr==imw_q)) then
                      fluxf_2 = vel_d(k)*rold(k+1)*den(k+1)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_WATER &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    elseif (isnow /= 1 .and. ipr==rmt_q) then
                      fluxf_2 = vel_d(k)*max(0.0_RP,rold(k+1) - qrp(rmat_q,ibr,icr,k+1))*den(k+1)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_WATER &
                              + fluxf_2/dzvm(k)*g           ! potential energy
                    endif
                  endif
                enddo
                if(ka.eq.k1eta) then
!
!------------------------------------------------------------------------
!                 Precip Rate and accumulation Rate for no iteration case
!------------------------------------------------------------------------
!
                  if(ipr.eq.1) then
                      ! scale surface flux
                      sflx = sflx + fluxf(k1eta-1)
                  endif
!
!     surface spray flux
!
                  if(ipr.eq.1.and.icr.eq.1.and.ispray.eq.1) then
                    fluxf(k1eta-1)=fluxf(k1eta-1)+flxsp(ibr)
                  endif
                endif
!
!     predict settling now!
!
                do k=ka,kb
                  if(k.eq.k1eta) then
                    rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                      /((z(k)-z(k-1))*den(k)*mass_convert)
                    !rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                    !  /((z(k)-z(k-1))*den(k))
                  else
                    rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                      /((zz(k)-zz(k-1))*den(k)*mass_convert)
                    !rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                    !  /((zz(k)-zz(k-1))*den(k))
                  endif
                enddo
              endif
            endif advscheme_if
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     End of iteration/noiteration branch
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

            if(level==3.or.level==4) then
              if(ipr.eq.1) then
                do k=ka,kb
                  diab=rold(k)-qrp(ipr,ibr,icr,k)
                  ql(k)=ql(k)+diab
                  qtp(k)=qtp(k)+diab
!     write(*,*) "k,diab,rold,qrp"
!     *                 ,k,diab,rold(k),qrp(ipr,ibr,icr,k)
                end do
              endif
            elseif(level>=5) then
              if(isnow==0) then
                if(ipr==rmt_q) then
                  do k=ka-1,kb+1
                    rold_1(k)=rold(k)- &
                    max(0.0_RP,qrp(rmt_q,ibr,icr,k) &
                             -qrp(rmat_q,ibr,icr,k))
                  enddo
                elseif(ipr==rmat_q) then
                  if(ibr.le.nbhzcl) then
                    do k=ka,kb
                      diab=rold_1(k)-rold(k)
                      qc(k)=qc(k)+diab
                      qtp(k)=qtp(k)+diab
!                           if(qtp(k)<=0.0_RP) then
!                              write(*,37) kmicv(k),imicv(k),jmicv(k) &
!                                         ,k,isnow,qtp(k)
!                              qtp(k)=1.0e-30_RP
!                           end if
                    enddo
                  else
                    do k=ka,kb
                      diab=rold_1(k)-rold(k)
                      ql(k)=ql(k)+diab
                      qtp(k)=qtp(k)+diab
!                           if(qtp(k)<=0.0_RP) then
!                              write(*,37) kmicv(k),imicv(k),jmicv(k) &
!                                         ,k,isnow,qtp(k)
!                              qtp(k)=1.0e-30_RP
!                           end if
                    end do
                  endif
                end if
              elseif(isnow==1) then
                if(ipr==imt_q) then
                  do k=ka-1,kb+1
! <<< 2015/03 T. Hashino changed for KiD
!org                           rold_1(k)=rold(k)- &
!org                                max(0.0_RP,qrp(imt_q,ibr,icr,k)- &
!org                                        qrp(imat_q,ibr,icr,k))
                    rold_1(k)=rold(k)- &
                            max(0.0_RP,qrp(imt_q,ibr,icr,k)- &
                                    qrp(imat_q,ibr,icr,k)-&
                                    qrp(imw_q,ibr,icr,k))
                  enddo
                elseif(ipr==imw_q) then
                  do k=ka-1,kb+1
                    rold_1(k)=rold_1(k)-rold(k)
                    diab=rold(k)-qrp(imw_q,ibr,icr,k)
                    qmlt(k)=qmlt(k)+diab
                    qtp(k)=qtp(k)+diab
                  enddo
! >>> 2015/03 T. Hashino changed for KiD
                elseif(ipr==imat_q) then
                  ierror1(ka:kb)=0
                  do k=ka,kb
                    diab=rold_1(k)-rold(k)
                    ql(k)=ql(k)+diab
                    qtp(k)=qtp(k)+diab
                    if(qtp(k)<=0.0_RP) then
                      ierror1(k)=1
                    end if
                  end do
                  if(any(ierror1(ka:kb)>0)) then
                    do k=ka,kb
                      if(ierror1(k)>0) then
                         if(debug) write(*,37) kmicv(k),imicv(k),jmicv(k) &
                               ,k,isnow,qtp(k)
37       format("scladvprz>qtp(k) is neg at k,i,j,kk,isnow,qtp" &
               ,5I5,ES15.6)
                        qtp(k)=1.0e-30
                      endif
                    enddo
                  endif

                end if
              endif
            endif

            do k=ka,kb
              qrp(ipr,ibr,icr,k)=rold(k)
            enddo

          enddo               !ipr loop
        enddo                  !ibr loop
      enddo                     !icr loop

!
!  momentum flux for scale
!
      ! z center grid, zz half grid
      do k = k1eta, n1 ! k1eta = 2
         if(k .eq. k1eta) then
            den_t(k) = den_t(k)-(fluxf_scale(k)-fluxf_scale(k-1))/(z(k)-z(k-1))
         else
            den_t(k) = den_t(k)-(fluxf_scale(k)-fluxf_scale(k-1))/(zz(k)-zz(k-1))
         endif
      enddo
      flux(n1) = 0.0_RP
      ! AMPS z-grid index: 1~n1, 1 is underground, n1=KE-KS+2. SCALE z-grid index: KE-1.
      ! please make reference to scale_atmos_phy_mp_common.f90/ATMOS_PHY_MP_precipitation_mom
      do k = 2, n1-1
        flux(k) = (fluxf_scale(k) + fluxf_scale(k-1))*momz(k)/(dens(k) + dens(k+1))
      enddo
      !flux(n1-1) = 0.0_RP

      ! z momentum
      do k = 2, n1-1
        momz_t(k) = momz_t(k) - (flux(k+1) - flux(k))*dzvm(k)
      enddo
      momz_t(n1) = 0.0_RP

      ! x momentum
      do k = 1, n1-1
        flux(k) = fluxf_scale(k)*U(k+1)
      enddo
      flux(n1) = 0.0_RP
      do k = 2, n1
        rhou_t(k) = rhou_t(k) - (flux(k) - flux(k-1))*dzzm(k)
      enddo

      ! y momentum
      do k = 1, n1-1
        flux(k) = fluxf_scale(k)*V(k+1)
      enddo
      flux(n1) = 0.0_RP
      do k = 2, n1
        rhov_t(k) = rhov_t(k) - (flux(k) - flux(k-1))*dzzm(k)
      enddo

      ! internal energy
      do k = 2, n1
        rhoe_t(k) = rhoe_t(k) - (eflx(k) - eflx(k-1))*dzzm(k)
      enddo
      return
      end subroutine sclsedprz_original


! <<< 2015/06 T. Hashino added
      subroutine sclsedaer(ica_sed1,ica_sed2 &
           ,n1,n2,n3,npa,nba,nca,k1a,k2a &
           ,qap,den,wacc &
           ,k1eta,z,zz,dzzm,dz1,dtl,i,j,CPPM,CPPME &
           ,iadvv &
           ,level,kmicv,imicv,jmicv)
      use par_amps
      use maxdims
      use mod_scladv_slppm, only : cal_aLRa6_z,cal_flux_z
      implicit none

      integer,intent(in) :: ica_sed1,ica_sed2 &
           ,n1,n2,n3 &
           ,npa,nba,nca &
           ,i,j,k1eta &
           ,iadvv,level
      integer,intent(inout) :: k1a,k2a
      integer,intent(in) :: imicv(*),jmicv(*),kmicv(*)

      real(RP),intent(in) :: wacc(*),den(*) &
           ,z(*),zz(*),dzzm(*) &
           ,CPPM(11,*),CPPME(11,n2,*)
      real(DS),intent(in) :: dtl
      real(RP),intent(inout) :: qap(npa,nba,nca,*)

! new space
      integer :: k,ka,kb,k1,k2 &
           ,ipa,iba,ica,niter,iter &
           ,inc,ninc &
           ,ibr2,ibr1
      real(RP) :: aLR(n1,2),a6(n1),dz1,dzzz(n1) &
           ,vel(n1,2),vt(n1,2) &
           ,rold(n1),rold_1(n1) &
           ,vel_d(n1),fluxf(n1)
      integer :: iwv,nwv

      real(RP) :: dtnew &
             ,cfl,wt &
             ,dz
!      real, parameter :: cpd=1004. &
!                        ,rd=287. &
!                        ,cpdrd=cpd/rd &
!                        ,p00=1.e5 &
!                        ,g=9.8 &
!                        ,gcp=g/cpd &
!                        ,coef=4./3.*3.14 &
!                        ,denl=1.e3 &
!                        ,rmax=500. &
!                        ,rcld=10.

      integer :: n1mxd,n2mxd,n3mxd
      integer :: istp,ntn
      integer :: jmt_q,jcon_q
!
!-------------------------------------------------------------------------------------
!     k1 and k2 are the vertical limits of integration determined by the points with
!     non zero total rain (to get started)
!---------------------------------------------------------------------------------------
!
!      return

      n1mxd=n1
      n2mxd=n2
      n3mxd=n3
      istp=0

      nwv=2
      jmt_q=1
      jcon_q=2

      dzzz(k1eta)=dz1
      do k=k1eta+1,n1
         dzzz(k)=1.0/dzzm(k)
      end do

      do k=1,n1
         rold(k)=0.
         rold_1(k)=0.
         vt(k,1:2)=0.
         vel(k,1:2)=0.
         fluxf(k)=0.
      enddo
!
!
      if(k1a.gt.k2a) return

      k1=max(k1a-1,k1eta)
      k2=min(k2a+1,n1-1)
      if(iadvv.eq.6) then
         k1=max(k1eta,k1-2)
         k2=min(n1,k2+2)
      endif

!----------------------------------------------------------------------------------------
!     calculate mass flux
!----------------------------------------------------------------------------------------
!
      do ica=ica_sed1,ica_sed2
         do iba=1,nba
!
!-------------------------------------------------------------------------------------
!         k1 and k2 are the vertical limits of integration determined by the points with
!         non zero terminal verlocity for this specific bin and category
!---------------------------------------------------------------------------------------
!
!
!           determine the primary (mass) terminal velocity
!
            do k=k1,k2
               do iwv=1,nwv
                  if(qap(jmt_q,iba,ica,k).gt.0.0_RP.and. &
                     qap(npa-nwv+iwv,iba,ica,k).lt.0.0_RP) then
                     vt(k,iwv)=abs(qap(npa-nwv+iwv,iba,ica,k))*wacc(k)
                  elseif(qap(jmt_q,iba,ica,k).lt.1.e-4_RP.and. &
                         qap(npa-nwv+iwv,iba,ica,k).ge.0.0_RP) then
                     vt(k,iwv)=0.
                  else
                     if(qap(npa-nwv+iwv,iba,ica,k).gt.0.0_RP) then
             print 444,k,i,j,iba,ica,qap(npa-nwv+iwv,iba,ica,k) &
                 ,qap(jmt_q,iba,ica,k)
444           format('Terminal velocity not set properly  k,i,j' &
                      ,6i5,'  vt=',3pf12.6,' qr=',0pe12.4)
!             stop
                     endif
                     vt(k,iwv)=0.
                  endif
!                  if(k.eq.k1) then
!                    write(*,'("ck vt",2I5,10ES15.6)') k1,iwv,vt(k,iwv),vt(k+1,iwv) &
!                          ,qap(npa-nwv+iwv,iba,ica,k)
!                  endif
                  if(vt(k,iwv).gt.0..or.vt(k,iwv).lt.-500.0_RP) then
                     print *,'terminal velocity is wrong here'
                     print *,'iba,ica',iba,ica
                     print *,'kmic,imic,jmic',kmicv(k) &
                                 ,imicv(k),jmicv(k)
                     print *,'k,iwv,nwv,vt(k,iwv)=',k,iwv,nwv,vt(k,iwv)
                     print 1122,(k,ipa,wacc(k),qap(ipa,iba,ica,k) &
                      ,ipa=1,npa)
 1122              format('k,ipa,wacc,qap(ipa,iba,ica)',2i5,2e15.5)
                     stop
                  endif
               end do
            enddo
!
!----------------------------------------------------------------------------------
!                 Find  velocity at vertical grid interface points and then
!                 look for cfl violations and how many iterations are necessary to
!                 keep it stable.  This is primarily for mesoscale applications
!                 where vertical resolution is high and we want rain to fall several
!                 vertical grid cells in a single dtl length timestep
!----------------------------------------------------------------------------
!

            cfl=0.0_RP
            do k=max(k1,k1eta+1),k2
               cfl=max(real(cfl,DS), &
                   max(abs(vt(k,1)),abs(vt(k-1,1)) &
                      ,abs(vt(k,2)),abs(vt(k-1,2)))*dtl &
                 /(zz(k)-zz(k-1)))
            enddo
!
!----------------------------------------------------------------------------------------
!     Find how many iterations (small timesteps)  necessary for stability
!----------------------------------------------------------------------------------------
!
            niter=1.+cfl/0.95
!
!----------------------------------------------------------------------------------
!          Now begin loop through all parameters except last one, which is the terminal velocity
!---------------------------------------------------------------------------------
!
!
!           ka and kb are the limits of the vertical integration expanded
!           for multiple iterations
!
            ka=max(k1eta,k1-niter+1)
            kb=min(n1-1,k2+niter-1)

!!!            write(*,'("ck den:",I5,10ES15.6)') k,z(ka-1),den(ka-1)
!tmp            write(*,*) "niter ",niter,ka,kb

            do iwv=1,nwv
               call velfill(vt(ka-1,iwv),min(kb-ka+3,n1),z(ka-1),den(ka-1))
               do k=ka-1,kb+1
                  vel(k,iwv)=0.5*(vt(k,iwv)+vt(min(n1,k+1),iwv))
!                  write(*,'("ck vel:",2I5,10ES15.6)') iwv,k,vel(k,iwv) &
!                         ,vt(k,iwv),vt(k+1,iwv)
               enddo
            end do
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Parameter Loop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
            ! assuming only mass mixing ratio is predicted.
            iwv=1
            do ipa=1,npa-nwv   ! begin ipa "parameter" loop
               do k=ka-1,kb+1
                  rold(k)=qap(ipa,iba,ica,k)
                  vel_d(k)=vel(k,iwv)
               enddo

               if(iadvv.eq.12) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Semi-Lagrangian PPM
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  dtnew=dtl/float(niter)
                  if(debug .and. niter.gt.1) then
                     write(*,'("Warning:multiple step in sedimentation",8I5,10ES15.6)') &
                        kmicv(k),imicv(k),jmicv(k),iba,ica,niter &
                       ,ka,kb,cfl
                     do k=max(k1eta,ka-1),min(n1-1,kb+1)
                        write(*,'("k,i,j,vel_d",3I5,10ES15.6)')  &
                         kmicv(k),imicv(k),jmicv(k),vel_d(k) &
                        ,qap(jmt_q,iba,ica,k),qap(jcon_q,iba,ica,k) &
                        ,qap(jcon_q,iba,ica,k)/max(qap(jcon_q,iba,ica,k),1.0e-25_RP)
                     enddo
                  endif

                  do iter=1,niter

!     semi-Lagrangian PPM method
!     1. calculation of points at left and right interfaces and coefficient.
                     call cal_aLRa6_z(i,j,max(k1eta,ka-1),min(n1,kb+1) &
                       ,k1eta,n1,n1mxd,n2mxd,n3mxd &
                       ,n1,istp &
                       ,CPPM,CPPME &
                       ,rold(1),aLR,a6)

!     2. calculation of flux
                     do k=max(k1eta,ka-1),min(n1-1,kb)
!     i is the interface position
                        call cal_flux_z(k,k1eta,n1,n1mxd,n1 &
                          ,zz,dzzz,dtnew &
                          ,den,rold(1),aLR,a6 &
                          ,zz(k),vel_d(k),fluxf(k),'sed 1')
                     end do

                     if(ka.eq.k1eta) then
                        fluxf(k1eta-1)=vt(k1eta,iwv) &
                                     *rold(k1eta)*den(k1eta)
                        wt=min(1.0_RP,max(dzcrit,(zz(k1eta)-zz(k1eta-1))) &
                             /(zz(k1eta+1)-zz(k1eta)))
                        fluxf(k1eta-1)=fluxf(k1eta) &
                             +(fluxf(k1eta-1)-fluxf(k1eta))*wt
                     endif
!
!     predict settling now!
!
                     do k=ka,kb
                        rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                          /((zz(k)-zz(k-1))*den(k))

!tmp                       if(ipr==1) then
!tmp                         write(*,'("ck sed4:",I5,10ES15.6)') k,rold(k),fluxf(k-1:k),den(k),dtnew,z(k-1:k)
!tmp                       endif
                     enddo
                  enddo


               else
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Upstream linear
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Multiple iterations
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
                  if(niter.gt.1) then
                     dtnew=dtl/float(niter)
                     if(debug) write(*,*) "niter is more than 1",niter,dtnew
                     do iter=1,niter

!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!                 Forward - Upstream precip setteling over niter iterations
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
                        do k=ka-1,kb
                           if(vel_d(k).gt.0.0_RP) then
                              fluxf(k)=vel_d(k) &
                                   *0.5_RP*(den(k)+den(k+1))*rold(k)
                           else
                              fluxf(k)=vel_d(k) &
                                   *0.5_RP*(den(k)+den(k+1))*rold(k+1)
                           endif
                        enddo
!
!
                        if(ka.eq.k1eta) then
                           fluxf(k1eta-1)= &
                                vt(k1eta,iwv)*rold(k1eta)*den(k1eta)
                           wt=min(1.0_RP,max(dzcrit,(zz(k1eta)-zz(k1eta-1))) &
                                /(zz(k1eta+1)-zz(k1eta)))
                           fluxf(k1eta-1)=fluxf(k1eta) &
                                +(fluxf(k1eta-1)-fluxf(k1eta))*wt
                        endif
!
!
!     predict settling now!
!
                        do k=ka,kb
                           rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                                /((zz(k)-zz(k-1))*den(k))
                        enddo
                     enddo
                  else
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Single Iteration
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!                 Forward - Upstream vertical advection just 1 iteration
!------------------------------------------------------------------------
!------------------------------------------------------------------------

                     do k=ka-1,kb
                        if(vel_d(k).gt.0.0_RP) then
                           fluxf(k)=vel_d(k)*0.5 &
                                *(den(k)+den(k+1))*rold(k)
                        else
                           fluxf(k)=vel_d(k)*0.5 &
                                *(den(k)+den(k+1))*rold(k+1)
                        endif
                     enddo
                     if(ka.eq.k1eta) then
                        fluxf(k1eta-1)= &
                            vt(k1eta,iwv)*rold(k1eta)*den(k1eta)
                        wt=min(1.0_RP,max(dzcrit,(zz(k1eta)-zz(k1eta-1))) &
                             /(zz(k1eta+1)-zz(k1eta)))
                        fluxf(k1eta-1)=fluxf(k1eta) &
                             +(fluxf(k1eta-1)-fluxf(k1eta))*wt
                     endif
!
!     predict settling now!
!
                     do k=ka,kb
                        rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                             /((zz(k)-zz(k-1))*den(k))
!                        write(*,'("ck_sed",5I5,10ES15.6)') ipr,ibr,icr,kmicv(k),imicv(k),rold(k) &
!                              ,fluxf(k),fluxf(k-1),den(k),zz(k),zz(k-1)

                     enddo
                  endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     End of iteration/noiteration branch
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
               end if
               do k=ka,kb
                  qap(ipa,iba,ica,k)=rold(k)
               enddo

            enddo               !ipa loop
         enddo                  !iba loop
      enddo                     !ica loop

      return
    end subroutine sclsedaer

    subroutine moistthermo2(theta, &
                            den, &
                            t, &
                            ptot, &
                            qtp, &
                            qv, &
                            qc, &
                            qr, &
                            qi, &
                            micptr,  &
                            npr,nbr,ncr,npi,nbi,nci, &
                            level,nbhzcl, &
                            thp, &
                            qrp, &
                            qip, &
                            pb, &
                            estbar, &
                            esitbar, &
                            pi0, &
                            nmic,imic,jmic,kmic,ivis, &
                            from,isect)
!mod              ,from,isect,tke,qvcmx,qvcmn,iptcldy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     theta(k) : virtual potential temperature
!     th: potential temperature
!     thp: etheta il
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use scale_prc, only: &
         PRC_abort
      use par_amps
      use com_amps, only: fid_alog
      implicit none
!     arguments
      integer :: micptr(*),level,nbhzcl  &
             ,nmic,imic(*),jmic(*),kmic(*)  &
             ,npr,nbr,ncr,npi,nbi,nci       &
             ,ivis,isect
!             ,ivis,isect,iptcldy
      real(MP_KIND) :: qv(*),qr(*),qi(*),qc(*),t(*),ptot(*)  &
             ,thp(*),qtp(*),pb(*),theta(*),den(*)  &
             ,qrp(npr,nbr,ncr,*),qip(npi,nbi,nci,*) &
             ,pi0(*)
      real(DS) :: estbar(*),esitbar(*)
! ,tke(*),qvcmx(*),qvcmn(*)

      character *(*) :: from
!     newspace
      integer :: k,i,iter  &
             ,ipr,ibr,icr,ipi,ibi,ici,niter1,J &
             ,ibr_st
      real(PS) :: error,tol  &
          ,pi,fct   &
          ,qvc,th,th2,wt,es,qs,qc2,qvis  &
          ,esi,qex(nmic),svw_lmt,ov,sw,til
      real(PS) :: cv, ep, gama, gamai, cpr, rcp, aklv, akiv
!        number of weighted velocity
      integer :: nwv,kk
!!!      external nan

      real(PS) :: RRLMT,RILMT,RRLMTB,RILMTB
!
!     NOTE: these limit parameters are very important. This makes the model
!           slower or faster, and also affect predicted fields!
!
      parameter(RRLMT=1.0d-10,RILMT=1.0d-15)
      parameter(RRLMTB=1.0d-22,RILMTB=1.0d-22)

      integer :: ierr

      cv=cpd-r
      ep=r/rm
      gama=cpd/cv
      gamai=1.0_RP/gama
      cpr=cpd/r
      rcp=r/cpd
      tol=0.001_RP
      aklv=alvl/cpd
      akiv=alvi/cpd
      svw_lmt=0.0_RP


      if (level.ge.5) then
!
!        qr starts from nbhzcl+1 bin of qrp for bin microphysics (SHIPS)
!        liquid bin spectrum that uses nbhzcl bins as haze and cloud droplet.
!         write(fid_alog,*) "ck:nbhzcl",from,level,nbhzcl
         ibr_st=nbhzcl+1
         nwv=2
         do k=1,nmic
            qr(k)=0.
            qi(k)=0.
            qc(k)=0.
            micptr(k)=0
            do icr=1,ncr
               do ibr=1,ibr_st-1
!7b                  if(qrp(rmt_q,ibr,icr,k).gt.0.0_RP) then
                  if(qrp(rmt_q,ibr,icr,k).ge.RRLMTB) then
                     qc(k)=qc(k) &
                  +max(0.0_RP,qrp(rmt_q,ibr,icr,k)-qrp(rmat_q,ibr,icr,k))
                     micptr(k)=1
                  else
                     do ipr=1,npr
                        qrp(ipr,ibr,icr,k)=0.0_RP
                     end do
                  end if
               end do
               do ibr=ibr_st,nbr
!7b                  if(qrp(rmt_q,ibr,icr,k).gt.0.0_RP) then
                  if(qrp(rmt_q,ibr,icr,k).ge.RRLMTB) then
                     qr(k)=qr(k) &
                   +max(0.0_RP,qrp(rmt_q,ibr,icr,k)-qrp(rmat_q,ibr,icr,k))
                     micptr(k)=1
                  else
                     do ipr=1,npr
                        qrp(ipr,ibr,icr,k)=0.0_RP
                     end do
                  endif
               end do
            enddo
            do ibi=1,nbi
               do ici=1,nci
!7b                  if(qip(imt_q,ibi,ici,k).gt.0.0_RP) then
                  if(qip(imt_q,ibi,ici,k).ge.RILMTB) then
! <<< 2015/03 T. Hashino modified for KiD
!org                     qi(k)=qi(k) &
!org                     +max(0.0_RP,qip(imt_q,ibi,ici,k)-qip(imat_q,ibi,ici,k))
                     qi(k)=qi(k) &
                     +max(0.0_RP,qip(imt_q,ibi,ici,k)-qip(imat_q,ibi,ici,k)-qip(imw_q,ibi,ici,k))

                     qr(k)=qr(k) &
                     +max(0.0_RP,qip(imw_q,ibi,ici,k))
! >>> 2015/03 T. Hashino modified for KiD
                     micptr(k)=1
                  else
                     do ipi=1,npi
                        qip(ipi,ibi,ici,k)=0.0_RP
                     end do
                  endif
               enddo
            enddo
!
!           Integrity check of predicted precipitation
!
            if(qr(k)+qi(k).gt.qtp(k)) then
               qtp(k)=max(qtp(k),0.0_RP)
               if( qc(k)+qr(k)+qi(k) .gt. 1.0e-20_RP ) then
                  fct=max(0.0_RP,qtp(k)/(qc(k)+qr(k)+qi(k)))
               else
                  qtp(k)=0.0_RP
                  fct=0.0_RP
                  nwv=0
               end if

               qc(k)=qc(k)*fct
               qr(k)=qr(k)*fct
               qi(k)=qi(k)*fct
               do icr=1,ncr
                 do ibr=1,nbr
                    do ipr=1,npr-nwv
                       qrp(ipr,ibr,icr,k)=qrp(ipr,ibr,icr,k)*fct
                    end do
                 end do
               end do
               do ici=1,nci
                  do ibi=1,nbi
                     do ipi=1,npi-nwv
                       qip(ipi,ibi,ici,k)=qip(ipi,ibi,ici,k)*fct
                     end do
                  end do
               end do
            endif

!---------- This is initial guess ------
!
          qv(k)=max(0.0_RP,qtp(k)-qr(k)-qi(k)-qc(k))
          if(ivis.eq.0) then
            pi=(pi0(k)+pb(k))/cpd
            ptot(k)=pi**cpr*1.e5_RP
            t(k)=thp(k)*pi
            th=thp(k)*(1.0_RP+(aklv*(qr(k)+qc(k))  &
                         +akiv*qi(k))/max(t(k),253.0_RP))
            ov=qv(k)

            if(ivis.le.1) then
               !qv(k)=max(0.0_RP,qtp(k)-qr(k)-qi(k)-qc(k))

               til=thp(k)*pi
               t(k)=til*(1.0_RP+ &
                  (aklv*(qr(k)+qc(k))+akiv*qi(k))/253.0_RP)
               if(t(k)>253.0_RP) then
                  t(k)=0.5_RP*(til+sqrt(til**2+4.0_RP*til* &
                  (aklv*(qr(k)+qc(k))+akiv*qi(k))))
                  if(t(k)<253.0_RP) then
                     LOG_ERROR("moistthermo2",*) "mo2 t(k) is not right",t(k)
                     call PRC_abort
                  endif
               endif
               th=thp(k)*(1.0_RP+(aklv*(qr(k)+qc(k)) &
                         +akiv*qi(k))/max(t(k),253.0_RP))

               !theta(k)=th*(1.0_RP+0.61_RP*qv(k))
               theta(k)=th*(epsvap+qv(k))/epsvap/(1.0_RP+qv(k))

               !den(k)=ptot(k)/(r*theta(k)*pi) ! CHIARUI: DENSITY STAYS CONSTANT IN SCALE

               i=max(1,min(int(t(k))-163,149))
               wt=max(min(t(k)-float(i+163),1.0_RP),0.0_RP)
               es=estbar(i)*(1.-wt)+estbar(i+1)*wt
               qs=0.622_RP*es/max(ptot(k)-0.378_RP*es,es)
!              sw=ptot(k)*qv(k)/(0.622+qv(k))/es-1.0
               sw=ptot(k)*qv(k)/(0.622_RP*(1.0+0.61_RP*qv(k)))/es-1.0_RP

!              if(sw.gt.0.05) then
!                write(fid_alog,'("mo2 sup:",a,3I5,20ES15.6)')
!    *               from,kmic(k),imic(k),jmic(k)
!    *              ,ptot(k),qv(k),qr(k)+qc(k),qi(k),t(k),sw
!              endif

               !if(ivis.eq.1) then
               !   i=max(1,min(int(t(k))-163,149))
               !   wt=max(min(t(k)-float(i+163),1.0D0),0.0D0)
               !   es=estbar(i)*(1.-wt)+estbar(i+1)*wt
               !   qs=0.622D0*es/max(ptot(k)-0.378D0*es,es)
               !   if(qv(k).gt.qs) micptr(k)=1
               !end if
            elseif(ivis.eq.2) then
               qvc=max(0.0_RP,qtp(k)-qr(k)-qi(k))
               do iter=1,100
!
                  t(k)=th*pi
                  i=max(1,min(int(t(k))-163,149))
                  wt=max(min(t(k)-float(i+163),1.0_RP),0.0_RP)
                  es=estbar(i)*(1.-wt)+estbar(i+1)*wt
                  qs=0.622_RP*es/max(ptot(k)-0.378_RP*es,es)
                  qc(k)=max(qvc-qs,0.0_RP)
                  qv(k)=qvc-qc(k)
                  th2=thp(k)*(1.+(aklv*(qr(k)+qc(k))+akiv*qi(k)) &
                       /max(t(k),253.0_RP))
                  error=th2-th
                  th=th+0.3_RP*error
                  if(abs(error).le.tol) go to 50
               enddo
50             continue
               theta(k)=th*(1.+0.61_RP*qv(k))
               t(k)=th*pi
               den(k)=ptot(k)/(r*theta(k)*pi)
               i=max(1,min(int(t(k))-163,149))
               wt=max(min(t(k)-float(i+163),1.0_RP),0.0_RP)
               es=estbar(i)*(1.-wt)+estbar(i+1)*wt
               qs=0.622_RP*es/max(ptot(k)-0.378_RP*es,es)
               qc(k)=max(qvc-qs,0.0_RP)
               qv(k)=qvc-qc(k)
               if(qc(k).gt.0.0_RP) micptr(k)=1
            end if

            if(debug .and. den(k).lt.0.0_RP) then
               write(fid_alog,*) from
               write(fid_alog,'("moisthermo2:density is negative k,den(k),ptot(k),theta(k):",I5,3ES15.6)') &
                    k,den(k) ,ptot(k),theta(k)
               write(fid_alog,*) "qv,th,thp,t",qv(k),th,thp(k),t(k)
               write(fid_alog,*) "qr,qc,qi",qr(k),qc(k),qi(k)
               write(fid_alog,*) "qtp",qtp(k)
             end if
!            if(qtp(k)-(qc(k)+qr(k)+qi(k)+qv(k))<0.0_RP) then
!               write(fid_alog,'("mo2>qt,qc,qr,qi,qv",6ES15.6)')
!     *                qtp(k),qc(k),qr(k),qi(k),qv(k)
!            end if

            i=max(1,min(int(t(k))-163,149))
            wt=max(min(t(k)-(i+163),1.0_RP),0.0_RP)
            es=estbar(i)*(1.-wt)+estbar(i+1)*wt
            qs=.622*es/max(ptot(k)-0.378_RP*es,es)


!        turn on the microphysics pointer if supersaturated over ice,
!        but it is called only before microphysics scheme.
!
            if(ivis>=1.and.t(k).lt.273.16_RP) then
               J=MAX(1,MIN(int(T(K))-163,110))
               wt=MAX(MIN(T(K)-(J+163),1.0_RP),0.0_RP)
               esi=esitbar(J)*(1.-wt)+esitbar(J+1)*wt
               qvis=0.622_RP*esi/max((ptot(k)-0.378_RP*esi),esi)
               if(qv(k).gt.qvis) micptr(k)=1
            end if

          endif ! chiarui test

         enddo


      end if
      return
      end subroutine moistthermo2

      subroutine diabatic(fthp,fqtp,qv,L,fqldiab,fqidiab,thp,qtp &
                         ,theta,qc,qr,qi,dtl,pb,pi0,k1m,k2m &
                         ,kmic,imic,jmic,k1etat,idir,from)
      use com_amps, only: fid_alog
      implicit none
      integer :: L,kmic(*),imic(*),jmic(*),k1etat(*) &
             ,k1m,k2m,i,j,k,k1,kk,idir
      real(PS) :: fqldiab(*),fqidiab(*),fthp(*),fqtp(*),qv(*),theta(*) &
          ,pi0(*),pb(*),th,t,thp(*),qtp(*) &
          ,qi(*),qc(*),qr(*) &
          ,fthdiab,ql2,qi2,thil1,thil2,ql
      real(DS),intent(in) :: dtl
      real(RP) :: cv, ep, gama, gamai, cpr, rcp, aklv, akiv
      character *(*) :: from

      cv=cpd-r
      ep=r/rm
      gama=cpd/cv
      gamai=1.0_RP/gama
      cpr=cpd/r
      rcp=r/cpd
      aklv=alvl/cpd
      akiv=alvi/cpd

      if(k2m.lt.0) return
      if(idir.eq.1) then
        k1=k1etat(1)
      else
        k1=1
      endif
      do k=k1,L
        if(idir.ne.1.and.kmic(k).lt.k1etat(k)) cycle
        if(fqldiab(k).ne.0..or.fqidiab(k).ne.0.0_RP) then
!         fqtp(k)=fqtp(k)+fqldiab(k)+fqidiab(k)
          qv(k)=qtp(k)-qc(k)-qr(k)-qi(k)
          th=theta(k)/(1.0_RP+0.61_RP*qv(k))
          t=th*(pb(k)+pi0(k))/cpd
!
!         ql=qc(k)+qr(k)
!         ql2=ql+fqldiab(k)*dtl
!         qi2=qi(k)+fqidiab(k)*dtl
!         thil1=th/(1.+(aklv*ql+akiv*qi(k))/max(t,253.0_RP))
!         thil2=th/(1.+(aklv*ql2+akiv*qi2)/max(t,253.0_RP))
!         fthdiab=(thil2-thil1)/dtl
!         fthp(k)=fthp(k)+fthdiab
!
          if(th.lt.100..or.th.gt.1000..or.t.lt.100..or.t.gt.320.0_RP) then
            print 10,from,k1,kmic(k),imic(k),jmic(k) &
                 ,t,th,theta(k),qv(k),pb(k),pi0(k)
10          format('bad input to diabatic1 call from ',a16,/ &
            'k1,k,i,j,t,th,theta(k),qv,pb,pi0 ',4i5,6e15.5)
            print 11,(kk,thp(kk),theta(kk),qv(kk),qc(kk),qr(kk),qi(kk) &
                  ,pb(kk),pi0(kk),kk=1,l)
11         format( &
           'k,thp(k),theta(k),qv(k),qc(k),qr(k),qi(k),pb(k),pi0(k)' &
           ,/,(i5,0p2f10.3,3p4f10.3,0p2f10.1))
!            stop
          endif

!
!         Keep theta constant as precip falls
!
          thp(k)=th/(1.0D0 &
                +(aklv*(qc(k)+qr(k))+akiv*qi(k))/max(t,253.0_PS))

!          thp(k)=thp(k)*(1.+(aklv*(qr(k)+qc(k)-fqldiab(k)) &
!                            +akiv*(qi(k)      -fqidiab(k))) &
!                            /max(t,253.0_PS)) &
!                       /(1.+(aklv*(qr(k)+qc(k)) &
!                            +akiv*(qi(k))) &
!                            /max(t,253.0_PS))
!
        endif
      enddo

      return
      end subroutine diabatic

      real(PS) function rainmass(i)
        use scale_prc, only: &
           PRC_abort
        use com_amps, only: &
           fid_alog
        integer :: i
        LOG_ERROR("rainmass",*) "ERROR! rainmass is not ready yet"
        call PRC_abort

        rainmass=i
        return
      end function rainmass

      subroutine negadj(n1,npr,nbr,ncr,npi,nbi,nci,i,kk &
                       ,qtp,thp,pc,theta,qv &
                       ,qr,qi,DEN &
                       ,pi01d &
                       ,iccnfl,ircnfl,ipcnfl,iscnfl,iacnfl,igcnfl &
                       ,ictcd,ictrn,ictpc,ictsw,ictag,ictgr)
      use maxdims
      implicit none
      integer :: k,ka,kb,k1,n1,i,kk &
             ,npr,nbr,ncr,npi,nbi,nci &
             ,iccnfl,ircnfl,ipcnfl,iscnfl,iacnfl,igcnfl &
             ,ictcd,ictrn,ictpc,ictsw,ictag,ictgr
      real(PS) :: qtp(*),thp(*) &
          ,pc(*),theta(*),qv(*) &
          ,qr(npr,nbr,ncr,*),qi(npi,nbi,nci,*) &
          ,DEN(*) &
          ,rltot,qtmin,f,vl,vi,rcd,dt,th,dlq,dic,t,dth,rtnew &
          ,di,dl &
          ,pi01d(*) &
          ,rc(n1),rr(n1),rp_chiarui(n1),rs(n1),ra(n1) &
          ,rg(n1),rlqic(n1),rlqicfx(n1) &
          ,ric(n1),rlq(n1) &
          ,cp,r,cv,rm,g,ep,gama,gamai,cpr,rcp,p00,p00i,alvl &
          ,alvi,aklv,akiv
      parameter(cp=1004.0_RP,r=287.0_RP,cv=cp-r,rm=461.5,g=9.8 &
               ,ep=r/rm,gama=cp/cv,gamai=1./gama,cpr=cp/r,rcp=r/cp &
               ,p00=1.e5,p00i=1./p00 &
               ,alvl=2.50e6,alvi=2.834e6,aklv=alvl/cp,akiv=alvi/cp)
!
!      This routine assures all positive definites remain positive
!

      ! limit for mass mixing ratio
      real(PS),parameter :: RCCRIT=0.1e-15 & ! 1 1/m^3 x 0.5 micron
                       ,RRCRIT=1.0E-12 & ! 0.001 #/m^3  x 100 micron
                       ,RPCRIT=0.1E-15 & ! 1 1/m^3 x 0.5 micron
                       ,RSCRIT=1.E-12 & ! 0.001 #/m^3 x 100 micron
                       ,RACRIT=1.E-12 &
                       ,RGCRIT=1.E-12
      ! particle mass limit to limit number concentration
      ! cloud droples: 0.5 to 100 micron   in radius
      real(PS),parameter :: AMC0=0.5E-15,AMCMX=4.0E-9
      ! rain drops: 10 micron to 5 mm
      real(PS),parameter :: AMR0=4.E-12,AMRMX=5.0E-4
      ! pristine crystals: 13 micron to 2.3 cm (assuming density 0.1 g/cm^3)
      real(PS),parameter :: AMP0=1.E-12,AMPMX=50.E-4
      ! snow flakes: 13 micron to 4.9 cm (assuming density 0.1 g/cm^3)
      real(PS),parameter :: AMS0=1.E-12,AMSMX=5.E-2
      ! aggregates: 0.28 mm to 2.3 cm (assuming density 0.1 g/cm^3)
      real(PS),parameter :: AMA0=1.E-8 ,AMAMX=50.E-4
      ! graupel: 0.29 mm to 14 cm (assuming density 0.1 g/cm^3)
      real(PS),parameter :: AMG0=4.E-8,AMGMX=5.0


          k1=1
          ka=1
          kb=n1
          IF(ICCNFL.GT.0) THEN
            DO K=k1,N1
              rc(K)=qr(1,1,ictcd,k)
              if(iccnfl.eq.5.and.rc(k).eq.0.0_RP) then
                qr(2,1,ictcd,k)=0.0_RP
              endif
            ENDDO
          ELSE
            DO K=k1,N1
              rc(K)=0.0_RP
            ENDDO
          ENDIF
          IF(IRCNFL.GT.0) THEN
            DO K=k1,N1
              rr(K)=qr(1,1,ictrn,k)
              if(ircnfl.eq.5.and.rr(k).eq.0.0_RP) then
                qr(2,1,ictrn,k)=0.0_RP
              endif
            ENDDO
          ELSE
            DO K=k1,N1
              rr(K)=0.0_RP
            ENDDO
          ENDIF
          IF(IPCNFL.GT.0) THEN
            DO K=k1,N1
              rp_chiarui(K)=qi(1,1,ictpc,k)
              if(ipcnfl.eq.5.and.rp_chiarui(k).eq.0.0_RP) then
                qi(2,1,ictpc,k)=0.0_RP
              endif
            ENDDO
          ELSE
            DO K=k1,N1
              rp_chiarui(K)=0.0_RP
            ENDDO
          ENDIF
          IF(ISCNFL.GT.0) THEN
            DO K=k1,N1
              rs(K)=qi(1,1,ictsw,k)
              if(iscnfl.eq.5.and.rs(k).eq.0.0_RP) then
                qi(2,1,ictsw,k)=0.
              endif
            ENDDO
          ELSE
            DO K=k1,N1
              rs(K)=0.0_RP
            ENDDO
          ENDIF
          IF(IACNFL.GT.0) THEN
            DO K=k1,N1
              ra(K)=qi(1,1,ictag,k)
              if(iacnfl.eq.5.and.ra(k).eq.0.0_RP) then
                qi(2,1,ictag,K)=0.
              endif
            ENDDO
          ELSE
            DO K=k1,N1
              ra(K)=0.0_RP
            ENDDO
          ENDIF
          IF(IGCNFL.GT.0) THEN
            DO K=k1,N1
              rg(K)=qi(1,1,ictgr,k)
              if(igcnfl.eq.5.and.rg(k).eq.0.0_RP) then
                qi(2,1,ictgr,k)=0.0_RP
              endif
            ENDDO
          ELSE
            DO K=k1,N1
              rg(K)=0.
            ENDDO
          ENDIF

          rltot=0.
          qtmin=100.
          DO K=k1,N1
            ric(K)=rp_chiarui(K)+rs(K)+ra(K)+rg(K)
            rlq(K)=rc(k)+rr(K)
            rlqic(K)=ric(K)+rlq(K)
            rlqicfx(K)=MAX(rlqic(K),0.0_RP)
            rltot=rltot+rlqic(k)
            qtmin=min(qtmin,qtp(k))
          ENDDO
          IF(QTMIN.LE.0.OR.RLTOT.NE.0) THEN
            DO K=k1,N1
              if(rlqic(k).ne.0.0_RP.or.qtp(k).lt.0.0_RP) then
                F=den(k)
                IF(rc(K)*F.LT.RCCRIT) rc(K)=0.0_RP
                IF(rr(K)*F.LT.RRCRIT) rr(K)=0.0_RP
                IF(rp_chiarui(K)*F.LT.RPCRIT) rp_chiarui(K)=0.0_RP
                IF(rs(K)*F.LT.RSCRIT) rs(K)=0.0_RP
                IF(ra(K)*F.LT.RACRIT) ra(K)=0.0_RP
                IF(rg(K)*F.LT.RGCRIT) rg(K)=0.0_RP

                rc(K)=MIN(rc(K),rlqicfx(K)-rr(k)-rp_chiarui(K)-rs(K) &
                     -ra(K)-rg(K))
                IF(rc(K).LT.0.0_RP) THEN
                  rc(K)=0.0_RP
                  rr(K)=MIN(rr(K),rlqicfx(K)-rp_chiarui(K)-rs(K) &
                       -ra(K)-rg(K))
                  IF(rr(K).LT.0.0_RP) THEN
                    rr(K)=0.0_RP
                    rp_chiarui(K)=MIN(rp_chiarui(K),rlqicfx(K)-rs(K)-ra(K) &
                         -rg(K))
                    IF(rp_chiarui(K).LT.0.0_RP) THEN
                      rp_chiarui(K)=0.0_RP
                      rs(K)=MIN(rs(K),rlqicfx(K)-ra(K)-rg(K))
                      IF(rs(K).LT.0.0_RP) THEN
                        rs(K)=0.0_RP
                        ra(K)=MIN(ra(K),rlqicfx(K)-rg(K))
                        IF(ra(K).LT.0.0_RP) THEN
                          ra(K)=0.0_RP
                          rg(K)=MIN(rlqicfx(K),rg(K))
                        endif
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
                VI=rp_chiarui(K)+rs(K)+ra(K)+rg(K)
                VL=rc(k)+rr(k)
                rcd=max(0.0_RP,VL+VI)
                DI=VI-ric(K)
                DL=VL-rlq(K)
                TH=theta(k)/(1.0_RP+0.61_RP*qv(k))
                DLQ=MAX(DL+MIN(DI,0.0_RP),0.0_RP)
                DIC=MAX(DI+MIN(DL,0.0_RP),0.0_RP)
                T=TH*(PI01D(K)+pc(k))/CP
                dth=-thp(k)**2/TH &
                     *(AKLV*DLQ+AKIV*DIC)/MAX(T,253.0_RP)
                rtnew=MAX(0.0_RP,qtp(k)-rlqic(K))+rcd
                thp(k)=thp(k)+dth
                qtp(k)=rtnew
                qv(k)=max(0.0_RP,min(qv(k),rtnew-rcd))
                if ( debug ) then
                   if(qtp(k).lt.0.0_RP.or.qv(k).gt.qtp(k))then
                      write(*,*) "something wrong",qtp(k),qv(k)
                   endif
                end if

                IF(ICCNFL.GT.0) THEN
                  IF(ICCNFL.EQ.5) THEN
                    qr(2,1,ictcd,k)=qr(2,1,ictcd,k)*rc(K) &
                              /MAX(qr(1,1,ictcd,k),RCCRIT)
                    qr(2,1,ictcd,k)=rc(K)/MAX(AMC0,MIN(AMCMX, &
                    rc(K)/MAX(qr(2,1,ictcd,k),1.E-20_RP) ))
                  ENDIF
                  qr(1,1,ictcd,k)=rc(K)
                ENDIF
                IF(IRCNFL.GT.0) THEN
                  IF(IRCNFL.EQ.5) THEN
                    qr(2,1,ictrn,k)=qr(2,1,ictrn,k)*rr(K) &
                              /MAX(qr(1,1,ictrn,k),RRCRIT)
                    qr(2,1,ictrn,k)=rr(K)/MAX(AMR0,MIN(AMRMX, &
                    rr(K)/MAX(qr(2,1,ictrn,k),1.E-20_RP) ))
                  ENDIF
                  qr(1,1,ictrn,k)=rr(K)
                ENDIF
                IF(IPCNFL.GT.0) THEN
                  IF(IPCNFL.EQ.5) THEN
                    qi(2,1,ictpc,k)=qi(2,1,ictpc,k)*rp_chiarui(K) &
                              /MAX(qi(1,1,ictpc,k),RPCRIT)
                    qi(2,1,ictpc,k)=rp_chiarui(K)/MAX(AMP0,MIN(AMPMX, &
                    rp_chiarui(K)/MAX(qi(2,1,ictpc,k),1.E-20_RP) ))
                  ENDIF
                  qi(1,1,ictpc,k)=rp_chiarui(K)
                ENDIF
                IF(ISCNFL.GT.0) THEN
                  IF(ISCNFL.EQ.5) THEN
                    qi(2,1,ictsw,k)=qi(2,1,ictsw,k)*rs(K) &
                              /MAX(qi(1,1,ictsw,k),RSCRIT)
                    qi(2,1,ictsw,k)=rs(K)/MAX(AMS0,MIN(AMSMX, &
                    rs(K)/MAX(qi(2,1,ictsw,k),1.E-20_RP) ))
                  ENDIF
                  qi(1,1,ictsw,k)=rs(K)
                ENDIF
                IF(IACNFL.GT.0) THEN
                  IF(IACNFL.EQ.5) THEN
                    qi(2,1,ictag,k)=qi(2,1,ictag,k)*ra(K) &
                              /MAX(qi(1,1,ictag,k),RACRIT)
                    qi(2,1,ictag,k)=ra(K)/MAX(AMA0,MIN(AMAMX, &
                    ra(K)/MAX(qi(2,1,ictag,k),1.E-20_RP) ))
                  ENDIF
                  qi(1,1,ictag,k)=ra(K)
                ENDIF
                IF(IGCNFL.GT.0) THEN
                  IF(IGCNFL.EQ.5) THEN
                    qi(2,1,ictgr,k)=qi(2,1,ictgr,k)*rg(K) &
                              /MAX(qi(1,1,ictgr,k),RGCRIT)
                    qi(2,1,ictgr,k)=rg(K)/MAX(AMG0,MIN(AMGMX, &
                    rg(K)/MAX(qi(2,1,ictgr,k),1.E-20_RP) ))
                  ENDIF
                  qi(1,1,ictgr,k)=rg(K)
                ENDIF
              ENDIF
            ENDDO
          ENDIF
      RETURN
      END subroutine negadj
!
      subroutine negadj_bin2(llbin,n1,npr,nbr,ncr,npi,nbi,nci,i,kk &
           ,qtp,thp,pc,theta,qv &
           ,qc,qr,qi &
           ,pi01d &
           ,ircnfl,ipcnfl &
           ,imicv,jmicv,kmicv)
      use com_amps
      use par_amps
      use maxdims
      implicit none
!      include 'com_micro'
      integer :: k,ka,kb,k1,llbin,n1,i,kk,iq &
             ,npr,nbr,ncr,npi,nbi,nci &
             ,ircnfl,ipcnfl
      integer :: imicv(llbin),jmicv(llbin),kmicv(llbin)
      real(MP_KIND) :: qtp(llbin),thp(llbin) &
          ,pc(llbin),theta(llbin),qv(llbin),qc(llbin) &
          ,pi01d(*) &
          ,qr(npr,nbr,ncr,llbin),qi(npi,nbi,nci,llbin)
      real(PS) :: rltot,qtmin,f,vi,vm,rcd,dt,th,dlq,dic,t,dth,rtnew &
          ,rlqicfx_left &
          ,di,dl &
          ,amr0,amrmx,amrmn,amp0,ampmx,ampmn &
          ,ams0,amsmx,amsmn,ama0,amamx,amamn &
          ,amg0,amgmx,amgmn &
          ,rr(n1),ri(n1),rml(n1),rlqic(n1),rlqicfx(n1) &
          ,ric(n1),rlq(n1),rmelt(n1) &
!          ,rr(lbnmx),ri(lbnmx),rlqic(lbnmx),rlqicfx(lbnmx) &
!          ,ric(lbnmx),rlq(lbnmx) &
          ,cp,r,cv,rm,g,ep,gama,gamai,cpr,rcp,p00,p00i,alvl &
          ,alvi,aklv,akiv
      parameter(cp=1004.0_RP,r=287.0_RP,cv=cp-r,rm=461.5,g=9.8 &
               ,ep=r/rm,gama=cp/cv,gamai=1./gama,cpr=cp/r,rcp=r/cp &
               ,p00=1.e5,p00i=1./p00 &
               ,alvl=2.50e6,alvi=2.834e6,aklv=alvl/cp,akiv=alvi/cp)
!
!      This routine assures all positive definites remain positive
!
      real(PS) :: RRLMT,RILMT,RRLMTB,RILMTB
      parameter(RRLMT=1.0d-10,RILMT=1.0d-15)
      parameter(RRLMTB=1.0d-22,RILMTB=1.0d-22)
      parameter( &
           AMR0=3.E-6,AMRMX=5.0E-4 &
           ,AMRMN=5.E-4,AMP0=1.E-12,AMPMX=50.E-4,AMPMN=1.E-8 &
           ,AMS0=1.E-12,AMSMX=5.E-2,AMSMN=1.e-8,AMA0=1.E-8 &
           ,AMAMX=50.E-4,AMAMN=37.6E-3,AMG0=4.E-8,AMGMX=5.0 &
           ,AMGMN=28.3E-3)
!      integer mark_cmi(n1),mark_cmr(n1)
!     *     ,marki(nbi,N1),markr(nbr,N1)
      real(PS) :: oqr(nbr,n1),oqi(nbi,n1),mod_rat1,mod_rat2,test_sum
      real(PS) :: nqr(nbr,n1),nqi(nbi,n1)
!      real binbi(nbi),binbr(nbr)
!       write(fid_alog,'(30ES15.7)') (binbi(ka),ka=1,6)
!       write(fid_alog,'(30ES15.7)') (binbr(ka),ka=1,11)

!tmp      return

!tmp      write(fid_alog,*) "negadj: rmt,rmat,rcon",rmt_q,rmat_q,rcon_q

      k1=1
      IF(IRCNFL.GT.0) THEN
         DO K=k1,N1
            rr(K)=0.0
            do ka=1,nbr
               oqr(ka,k)=qr(rmt_q,ka,1,k)-qr(rmat_q,ka,1,k)
               nqr(ka,k)=oqr(ka,k)
!7b               if(qr(rmt_q,ka,1,k).lt.RRLMT*qtp(k)) then
               if(qr(rmt_q,ka,1,k).lt.RRLMTB) then
                  do kb=1,npr
                     qr(kb,ka,1,k)=0.0_RP
                  end do
               else
                  if(oqr(ka,k)<0.01*qr(rmat_q,ka,1,k)) then
                     if(qr(rmat_q,ka,1,k)>0.0_RP) then
                        qr(rmt_q,ka,1,k)=1.01_RP*qr(rmat_q,ka,1,k)
                        oqr(ka,k)=qr(rmt_q,ka,1,k)-qr(rmat_q,ka,1,k)
                        nqr(ka,k)=oqr(ka,k)
                        rr(K)=rr(K)+oqr(ka,k)
                     else
                        do kb=1,npr
                           qr(kb,ka,1,k)=0.0_RP
                        end do
                     end if
                  else
                     rr(K)=rr(K)+oqr(ka,k)
!ccc                  if(qr(rmt_q,ka,1,k).eq.0.0_RP) then
!ccc                     do kb=2,npr
!ccc                        qr(kb,ka,1,k)=0.0_RP
!ccc                     end do
!ccc                  end if
                  end if
               end if
            end do
         ENDDO
      ELSE
         DO K=k1,N1
            rr(K)=0.
         ENDDO
      ENDIF
      IF(IPCNFL.GT.0) THEN
         DO K=k1,N1
            ri(K)=0.0
! <<< 2015/03 T. Hashino modified for KiD
            rml(K)=0.0
            do ka=1,nbi
               oqi(ka,k)=qi(imt_q,ka,1,k)-qi(imat_q,ka,1,k)
!               oqi(ka,k)=qi(imt_q,ka,1,k)-qi(imat_q,ka,1,k)-qi(imw_q,ka,1,k)
!               oqmlt(ka,k)=qi(imw_q,ka,1,k)
               nqi(ka,k)=oqi(ka,k)
!               nqmlt(ka,k)=oqmlt(ka,k)
! << 2015/03 T. Hashino modified for KiD
!7b               if(qi(1,ka,1,k).lt.RILMT*qtp(k)) then
               if(qi(imt_q,ka,1,k).lt.RILMTB) then
                  do kb=1,npi
                     qi(kb,ka,1,k)=0.0
                  end do
               else
!!c                  mod_rat1=qi(imt_q,ka,1,k)-qi(imr_q,ka,1,k)-qi(imc_q,ka,1,k)
!!c                  if(qi(imc_q,ka,1,k)<mod_rat1) then
!!c                     write(fid_alog,*) "negadj 1:k,bin,massag produced",k,ka
!!c     *                   ,mod_rat1,mod_rat1,qi(1,ka,1,k),qi(imc_q,ka,1,k)
!!c                  end if

                  if(oqi(ka,k)<0.01*qi(imat_q,ka,1,k)) then
                     if(qi(imat_q,ka,1,k)>0.0_RP) then
                        qi(imt_q,ka,1,k)=1.01*qi(imat_q,ka,1,k)
                        qi(imc_q,ka,1,k)=qi(imt_q,ka,1,k)
                        qi(imr_q,ka,1,k)=0.0_RP
                        qi(ima_q,ka,1,k)=0.0_RP
                        qi(imw_q,ka,1,k)=0.0_RP
                        oqi(ka,k)=qi(imt_q,ka,1,k)-qi(imat_q,ka,1,k)
                        nqi(ka,k)=oqi(ka,k)
                        ri(K)=ri(K)+oqi(ka,k)
                     else
                        do kb=1,npi
                           qi(kb,ka,1,k)=0.0_RP
                        end do
                     end if
                  else
                     ri(K)=ri(K)+oqi(ka,k)
! <<< 2015/03 T. Hashino mod
                     rml(K)=rml(K)+qi(imw_q,ka,1,k)
! >>> 2015/03 T. Hashino mod
!ccc               if(qi(1,ka,1,k).eq.0.0_RP) then
!ccc                  do kb=2,npi
!ccc                     qi(kb,ka,1,k)=0.0_RP
!ccc                  end do
!ccc               end if
                  end if
!!c                  mod_rat1=qi(imt_q,ka,1,k)-qi(imr_q,ka,1,k)-qi(imc_q,ka,1,k)
!!c                  if(qi(10,ka,1,k)<mod_rat1) then
!!c                     write(fid_alog,*) "negadj 3:k,bin,massag produced",k,ka
!!c     *                   ,mod_rat1,mod_rat1,qi(imt_q,ka,1,k),qi(imc_q,ka,1,k)
!!c                  end if
               end if
            end do
         ENDDO
      ELSE
         DO K=k1,N1
            ri(K)=0.
            rml(K)=0.
         ENDDO
      ENDIF

      rltot=0.0
      qtmin=100.
      DO K=k1,N1
         ric(K)=ri(K)
         rlq(K)=rr(K)
         rmelt(K)=rml(K)
         rlqic(K)=ric(K)+rlq(K)
         rlqicfx(K)=MAX(rlqic(K),0.0_RP)
         rltot=rltot+rlqic(k)
         qtmin=min(qtmin,qtp(k))
      ENDDO
!         write(fid_alog,*) "here 1"
      IF(QTMIN.LE.0.OR.RLTOT.NE.0) THEN
         DO K=k1,N1
            if(rlqic(k).ne.0.0_RP .or. qtp(k).lt.0.0_RP) then
               if(qtp(k)<0.0_RP) then
                  if ( debug ) then
                     write(fid_alog,*) "qtp is neg at k",k,qtp(k)
                     write(fid_alog,*) "qv,rlqicf ",qv(k),rlqicfx(k)
                  end if
                  qtp(k)=max(0.0_RP,min(qv(k),rlqicfx(k)))
!ccc                  qtp(k)=max(qtp(k),qv(k),0.0_RP)
               end if
!ccc               rlqicfx_left=qtp(k)-qv(k)-rlqicfx(k)
               rlqicfx_left=qtp(k)-rlqicfx(k)
!ccc               do ka = 1, nbr
!ccc                  if(qr(rmt_q,ka,1,k).lt.RRLMT*qtp(k)) then
!ccc                     do kb=1,npr
!ccc                        qr(kb,ka,1,k)=0.0_RP
!ccc                     end do
!ccc                  end if
!ccc               end do
!ccc               do ka = 1, nbi
!ccc                  oqi(ka,1)=qi(imt_q,ka,1,k)
!ccc                  if(qi(1,ka,1,k).lt.RILMT*qtp(k)) then
!ccc                     do kb=1,npi
!ccc                        qi(kb,ka,1,k)=0.0_RP
!ccc                     end do
!ccc                  end if
!ccc               end do

               if(rlqicfx_left>=0.0_RP) then
                  goto 777
!ccc                  write(fid_alog,*) "rlqicf_left is neg:",rlqicfx_left
!ccc                  write(fid_alog,*) "qtp(k),qv(k):",qtp(k),qv(k)
!ccc                  do ka = 1, nbr
!ccc                     do kb=1,npr
!ccc                        qr(kb,ka,1,k)=0.0_RP
!ccc                     end do
!ccc                  end do
!ccc                  do ka = 1, nbi
!ccc                     do kb=1,npi
!ccc                        qi(kb,ka,1,k)=0.0_RP
!ccc                     end do
!ccc                  end do
               else

                  rlqicfx_left=rlqicfx_left+nqr(1,k)
                  do ka=1,nbr-1
                     nqr(ka,k)=min(nqr(ka,k),rlqicfx_left)
                     if(nqr(ka,k).lt.0.0_RP) then
                        nqr(ka,k)=0.0_RP
                        rlqicfx_left=rlqicfx_left+nqr(ka+1,k)
                     else
                        goto 777
                     end if
                  end do
                  nqr(nbr,k)=min(nqr(nbr,k),rlqicfx_left)
                  if(nqr(nbr,k).lt.0.0_RP) then
                     nqr(nbr,k)=0.0_RP

                     rlqicfx_left=rlqicfx_left+nqi(1,k)
                     do ka = 1, nbi-1
                        nqi(ka,k)=min(nqi(ka,k),rlqicfx_left)
                        if(nqi(ka,k).lt.0.0_RP) then
                           nqi(ka,k)=0.0
                           rlqicfx_left=rlqicfx_left+nqi(ka+1,k)
                        else
                           goto 777
                        end if
                     end do
                  end if
                  nqi(nbi,k)=MIN(rlqicfx_left,nqi(nbi,k))
                  if(nqi(nbi,k)<0.0_RP) then
                     if ( debug ) then
                        write(fid_alog,*) "something wrong in neg at",k
                        write(fid_alog,*) "qtp-qv,qv",qtp(k)-qv(k),qv(k)
                        write(fid_alog,*) "qtp,rlqicfx,left",qtp(k),rlqicfx(k) &
                             ,rlqicfx_left
                        write(fid_alog,*) "nqi(nbi,k)",nqi(nbi,k)
                        write(fid_alog,*) "rr,ri",rr(k),ri(k)
                     end if
                     nqi(nbi,k)=0.0_RP
!                     stop
                  end if
               end if
 777           continue
               rr(k)=0.0
               do ka = 1, nbr
                  rr(k)=rr(k)+nqr(ka,k)
!     modify only total mass and mass components
                  if(oqr(ka,k).gt.0.0.and.oqr(ka,k).ne.nqr(ka,k)) &
                     then
                     mod_rat1=nqr(ka,k)/oqr(ka,k)
                     qr(rmt_q,ka,1,k)=qr(rmt_q,ka,1,k)*mod_rat1
                     do iq=i1_mcp_liq,i2_mcp_liq
                       qr(iq,ka,1,k)=qr(iq,ka,1,k)*mod_rat1
                     enddo
                  end if
               end do

               VI=0.0
               do ka = 1, nbi
                  VI=VI+nqi(ka,k)
!     modify only total mass and mass components
                  if(oqi(ka,k).gt.0.0.and.oqi(ka,k).ne.nqi(ka,k)) &
                     then
!cc                     mod_rat=qi(imt_q,ka,1,k)/oqi(ka,k)
!cc                     qi(imr_q,ka,1,k)=qi(imr_q,ka,1,k)*mod_rat
!cc                     qi(imc_q,ka,1,k)=qi(imc_q,ka,1,k)*mod_rat
                     mod_rat1=nqi(ka,k)/oqi(ka,k)
!                 total mass
                     qi(imt_q,ka,1,k)=qi(imt_q,ka,1,k)*mod_rat1

                     do iq=i1_mcp_ice,i2_mcp_ice
                       qi(iq,ka,1,k)=qi(iq,ka,1,k)*mod_rat1
                     enddo


!ccc                     mod_rat1=qi(imr_q,ka,1,k)/qi(imt_q,ka,1,k)
!ccc                     mod_rat2=qi(imc_q,ka,1,k)/qi(imt_q,ka,1,k)
!ccc
!ccc                     if(mod_rat2>=0.999999.and.mod_rat1<1.0e-6) then
!ccc                        qi(imc_q,ka,1,k)=qi(imt_q,ka,1,k)
!ccc                        qi(imr_q,ka,1,k)=0.0
!ccc                        qi(imw_q,ka,1,k)=0.0
!ccc                        qi(ima_q,ka,1,k)=0.0
!ccc                        write(fid_alog,*) "negadj:k,bin",k,ka
!ccc                     end if

!ccc                     mod_rat1=qi(imt_q,ka,1,k)-qi(imr_q,ka,1,k)-qi(imc_q,ka,1,k)
!ccc                     if(qi(imc_q,ka,1,k)<mod_rat1) then
!ccc                        write(fid_alog,*) "negadj 2:k,bin,massag produced",k,ka
!ccc     *                   ,mod_rat1,mod_rat1,qi(imt_q,ka,1,k),qi(imc_q,ka,1,k)
!ccc                     end if
                  end if
               end do

               vm=0.0
               do ka = 1, nbi
                 vm=vm+qi(imw_q,ka,1,k)
               end do

               rcd=max(0.0_RP,rr(k)+VI) ! K -> k
! <<< 2015/03 T. Hashino mod
!org               DI=VI-ric(K)
!org               DL=rr(K)-rlq(K)
               DI=(VI-vm)-(ric(K)-rmelt(k))
               DL=(rr(K)+vm)-(rlq(K)+rmelt(k))
! >>> 2015/03 T. Hashino mod
               TH=theta(k)/(1.0_RP+0.61_RP*qv(k))
               DLQ=MAX(DL+MIN(DI,0.0_RP),0.0_RP)
               DIC=MAX(DI+MIN(DL,0.0_RP),0.0_RP)
               T=TH*(PI01D(K)+pc(k))/CP
               dth=-thp(k)**2/TH &
                     *(AKLV*DLQ+AKIV*DIC)/MAX(T,253.0_RP)
               rtnew=MAX(0.0_RP,qtp(k)-rlqic(K))+rcd
               thp(k)=thp(k)+dth
               qtp(k)=rtnew
!               qv(k)=max(0.0_RP,min(qv(k),rtnew-rcd))

               qv(k)=max(0.0_RP,rtnew-rcd)

!               test_sum=0.0
!               do ka = 1, nbr
!                  test_sum=test_sum+qr(rmt_q,ka,1,k)
!               end do
!               if( test_sum .gt. qtp(k) ) then
!                  write(fid_alog,*) "error at neg 3", k,test_sum,qtp(k)
!               end if

!               write(fid_alog,*) "bf check_con",npr,nbr,ncr,n1,k
!     *                    ,(binbr(ka),ka=1,nbr+1)
!ccc               call check_con(npr,nbr,ncr,n1,k,qr
!ccc     *                 ,binbr,mark_cmr,markr
!ccc     *                 ,1)
!               if(qr(rmt_q,1,1,k).gt.0.0_RP) then
!                  write(fid_alog,222) imicv(k),jmicv(k),kmicv(k)
!     *                       ,qr(rcon_q,1,1,k),qr(rmt_q,1,1,k)
!222   format("con and mass of 1 bin at",3I5,2ES15.6)
!               end if
!ccc               call check_con(npi,nbi,nci,n1,k,qi
!ccc     *                 ,binbi,mark_cmi,marki
!ccc     *                 ,2)

               test_sum=0.0
               do ka = 1, nbr
                  test_sum=test_sum+qr(rmt_q,ka,1,k)-qr(rmat_q,ka,1,k)
               end do
               if( test_sum .gt. qtp(k) ) then
                  if(debug) write(fid_alog,*) "error at neg 4", k,test_sum,qtp(k)
               end if
!ccc               if( abs(qtp(k)-(test_sum+qv(k))).gt.1.0e-5) then
!ccc                  write(fid_alog,*) "error at neg 5", k,test_sum,qv(k),qtp(k)
!ccc               end if
            ENDIF
         ENDDO
      end if
      RETURN
      END subroutine negadj_bin2


      subroutine integ_mictend1(dmtend &
                           ,nmic,zz,zstv,denvm,dt &
                           ,imicv,jmicv,kmicv)

      implicit none
      integer,intent(in) :: nmic
      integer,intent(in),dimension(*) :: imicv,jmicv,kmicv
      real(PS),intent(in),dimension(*) :: zstv,zz,denvm
      real(PS),intent(inout) ::   dmtend(11,9,*)

      real(PS),intent(in) :: dt

      integer :: k,kk,i,j

      ! integrate vertically and put it in the first index
      do k=1,nmic
        kk=kmicv(k)
        do j=1,9
          do i=1,9
            dmtend(i,j,k)=dmtend(i,j,k)*(zz(kk)-zz(kk-1))*1.0e+3 ! [g/cm3/s] to [kg/m2/s]
          enddo
        enddo
      enddo
      do k=2,nmic
        kk=kmicv(k)
        do j=1,9
          do i=1,9
            dmtend(i,j,1)=dmtend(i,j,1)+dmtend(i,j,k)
          enddo
        enddo
      enddo

      end subroutine integ_mictend1

      subroutine integ_mictend2(dMcol_auto_liq,dMcol_accr_liq &
                          ,dMcol_auto_ice,dMcol_accr_ice &
                          ,dMcol_auto_rim,dMcol_accr_rim &
                          ,nmic,zz,zstv,denvm,dt &
                          ,imicv,jmicv,kmicv &
                          ,dM_auto_liq,dM_accr_liq &
                          ,dM_auto_ice,dM_accr_ice &
                          ,dM_auto_rim,dM_accr_rim )
      implicit none
      integer,intent(in) :: nmic
      integer,intent(in),dimension(*) :: imicv,jmicv,kmicv
      real(PS),intent(in),dimension(*) :: zstv,zz,denvm
      real(PS),intent(out) ::  &
                           dMcol_auto_liq,dMcol_accr_liq &
                          ,dMcol_auto_ice,dMcol_accr_ice &
                          ,dMcol_auto_rim,dMcol_accr_rim
      real(PS),intent(in),dimension(*) ::  &
                           dM_auto_liq,dM_accr_liq &
                          ,dM_auto_ice,dM_accr_ice &
                          ,dM_auto_rim,dM_accr_rim
      real(PS),intent(in) :: dt

      integer :: k,kk

      dMcol_auto_liq=0.0
      dMcol_accr_liq=0.0
      dMcol_auto_ice=0.0
      dMcol_accr_ice=0.0
      dMcol_auto_rim=0.0
      dMcol_accr_rim=0.0

      do k=1,nmic
        kk=kmicv(k)
!        write(fid_alog,*) "ck dz",kk,zz(kk)-zz(kk-1)
        dMcol_auto_liq=dMcol_auto_liq+dM_auto_liq(k)*(zz(kk)-zz(kk-1))
        dMcol_accr_liq=dMcol_accr_liq+dM_accr_liq(k)*(zz(kk)-zz(kk-1))

        dMcol_auto_ice=dMcol_auto_ice+dM_auto_ice(k)*(zz(kk)-zz(kk-1))
        dMcol_accr_ice=dMcol_accr_ice+dM_accr_ice(k)*(zz(kk)-zz(kk-1))

        dMcol_auto_rim=dMcol_auto_rim+dM_auto_rim(k)*(zz(kk)-zz(kk-1))
        dMcol_accr_rim=dMcol_accr_rim+dM_accr_rim(k)*(zz(kk)-zz(kk-1))
      enddo
      dMcol_auto_liq=dMcol_auto_liq*1.0e+3/dt   ! kg/m^2/s
      dMcol_accr_liq=dMcol_accr_liq*1.0e+3/dt   ! kg/m^2/s

      dMcol_auto_ice=dMcol_auto_ice*1.0e+3/dt   ! kg/m^2/s
      dMcol_accr_ice=dMcol_accr_ice*1.0e+3/dt   ! kg/m^2/s

      dMcol_auto_rim=dMcol_auto_rim*1.0e+3/dt   ! kg/m^2/s
      dMcol_accr_rim=dMcol_accr_rim*1.0e+3/dt   ! kg/m^2/s

      end subroutine integ_mictend2

!!$      subroutine check_con(np,nb,nc,n1,j,X,binb,mark_cm,mark,token)
!!$      use par_amps
!!$      implicit none
!!$      integer :: np,nb,nc,n1
!!$      real(PS) :: X(np,nb,nc,*),binb(*)
!!$      integer :: j
!!$      integer :: mark_cm(*),mark(nb,*)
!!$!     mean mass and middle mass
!!$      real(PS) :: mean_mass,mid_mass,old_con
!!$      integer :: message,i,k,token
!!$      integer :: jmt_q,jcon_q
!!$      real(PS) :: mod_ratio
!!$
!!$!      write(fid_alog,*) "in check_con",np,nb,nc,n1,j,(binb(i),i=1,nb+1)
!!$      if(token.eq.1) then
!!$        jmt_q=rmt_q
!!$        jcon_q=rcon_q
!!$      elseif(token.eq.2) then
!!$        jmt_q=imt_q
!!$        jcon_q=icon_q
!!$      endif
!!$
!!$      do i=1,nb
!!$         message=0
!!$         old_con=X(jcon_q,i,1,j)
!!$         if(X(jcon_q,i,1,j)>0.0_RP.and.X(jmt_q,i,1,j)>0.0_RP) then
!!$            mean_mass=X(jmt_q,i,1,j)/X(jcon_q,i,1,j)
!!$!            write(fid_alog,*) "mean_mass",mean_mass
!!$            if(mean_mass<binb(i).or.mean_mass>binb(i+1) ) then
!!$               mid_mass=(binb(i)+binb(i+1))/2.0
!!$               if( token .eq. 2 ) then
!!$!                  write(fid_alog,529) (X(k,i,1,j),k=1,10)
!!$529       format("XS",10ES15.6)
!!$!                  write(fid_alog,'("binb",11ES15.6)') (binb(k),k=1,6)
!!$               end if
!!$               X(jcon_q,i,1,j)=X(jmt_q,i,1,j)/mid_mass
!!$               mod_ratio=X(jcon_q,i,1,j)/max(old_con,1.0e-20_RP)
!!$! NOTE
!!$! only concentration weighted variables modified except mass components
!!$! since the mass compoents are already set up to sum up to total
!!$               if( token .eq. 2 ) then
!!$                  do k=i1_vcp_ice,i2_vcp_ice
!!$                    X(k,i,1,j)=X(k,i,1,j)*mod_ratio
!!$                  enddo
!!$                  do k=i1_acp_ice,i2_acp_ice
!!$                    X(k,i,1,j)=X(k,i,1,j)*mod_ratio
!!$                  enddo
!!$                  do k=i1_ccp_ice,i2_ccp_ice
!!$                    X(k,i,1,j)=X(k,i,1,j)*mod_ratio
!!$                  enddo
!!$               end if
!!$               message=2
!!$!               write(fid_alog,*)
!!$!     *            "check_con_77 > concentration was fixed (1) at"
!!$!     *                   , i,j
!!$!               write(fid_alog,'(A18,E15.6,A4, E15.6)')
!!$!     *          "           > From ", old_con, " to ", X(2,i,1,j)
!!$!               write(fid_alog,*) "           > mass,mod_ratio"
!!$!     *                    , X(1,i,1,j),mod_ratio
!!$!               write(fid_alog,*) "           > The token is ", token
!!$            end if
!!$         else if(X(jcon_q,i,1,j).ne.0.0_RP.and.X(jmt_q,i,1,j).eq.0.0_RP) then
!!$            do k=2,np
!!$               X(k,i,1,j)=0.0_RP
!!$            end do
!!$            message = 2
!!$!            write(fid_alog,*)
!!$!     *         " check_con_77 > concentration was fixed (2) at"
!!$!     *                , i,j
!!$!            write(fid_alog,'(A18,E15.6,A4, E15.6)')
!!$!     *       "           > From ", old_con, " to ", X(2,i,1,j)
!!$!            write(fid_alog,*) "           > mass", X(1,i,1,j)
!!$!
!!$!            write(fid_alog,*) "           > The token is ", token
!!$         else if(X(jcon_q,i,1,j)<=0.0_RP.and.X(jmt_q,i,1,j)>0.0_RP) then
!!$!            write(fid_alog,*)
!!$!     *         " check_con_77 > concentration was fixed (3) at"
!!$!     *                , i,j
!!$!            write(fid_alog,*) "              > The token is ", token
!!$!            write(fid_alog,'(A18, E15.6,A4, E15.6)')
!!$!     *       "           > From con", old_con, " to ", 0.0
!!$!            write(fid_alog,'(A18, E15.6,A4, E15.6)')
!!$!     *       "           > From mass", X(1,i,1,j), " to ", 0.0
!!$            do k=1,np
!!$              X(k,i,1,j)=0.0_RP
!!$            end do
!!$         else if(X(jcon_q,i,1,j)<0.0_RP) then
!!$!            write(fid_alog,*)
!!$!     *         " check_con_77 > concentration was fixed (4) at"
!!$!     *                , i,j
!!$!            write(fid_alog,*) X(1,i,1,j),X(2,i,1,j)
!!$            do k=1,np
!!$              X(k,i,1,j)=0.0_RP
!!$            end do
!!$
!!$         end if
!!$         mark(i,j)=message
!!$      end do
!!$      mark_cm(j)=message
!!$      return
!!$      end subroutine check_con

!!$      subroutine qiprint(qipv,npi,nbi,nci,nmic,from)
!!$      use com_amps, only: fid_alog
!!$      integer,intent(in) :: npi,nbi,nci,nmic
!!$      real(PS),intent(in) :: qipv(npi,nbi,nci,*)
!!$!      integer iv(*),jv(*),kv(*)
!!$      character *(*) :: from
!!$      integer :: k,ibi,kk
!!$
!!$      do k=1,nmic
!!$         do ibi=1,nbi
!!$            if(qipv(1,ibi,1,k).gt.1.e-5) then
!!$               write(fid_alog,10) from,k &
!!$                         ,(qipv(1,kk,1,k),kk=1,nbi) &
!!$                         ,(qipv(2,kk,1,k),kk=1,nbi)
!!$10        format(' qiprint ',a20,i5,/,(10ES14.6))
!!$               exit
!!$            endif
!!$         enddo
!!$      enddo
!!$      return
!!$      end subroutine qiprint

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine shift_bin_vec(n1,npr,nbr,ncr,k1r,k2r,binbr &
          ,qrpv,mmass,den,iphase &
          ,kmic,imic,jmic)
!     this subroutine has to be consistent with how to shift bins in
!     class_group
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      use maxdims
      use com_amps, only: fid_alog
!      use com_amps
      use par_amps
      implicit none
      integer :: n1,npr,nbr,ncr,k1r,k2r,iphase
      real(RP) :: qrpv(npr,nbr,ncr,*),mmass(nbr,ncr,*),den(*)
      real(PS),dimension(*),intent(in) :: binbr
      integer :: kmic(*),imic(*),jmic(*)
      integer :: naxis,nmass,ncon
      real(RP) :: den_i
      parameter(den_i= 0.91668)


      real(8),dimension(mxnbin,n1)               :: new_N
      real(8),dimension(mxnbin,n1,mxnmasscomp+1) :: new_M
      real(8),dimension(mxnbin,n1,mxnnonmc)      :: new_Q

      real(RP),dimension(mxnbin,n1,mxnnonmc) :: Qp


      ! total concentration in the shifted bin
      real(8), dimension(mxnbin+1,n1)           :: Np
      ! total mass in the shifted bin
      real(8), dimension(mxnbin+1,n1)           :: Mp

      ! averaged mass tendency after the time step for the bin
      real(8), dimension(mxnbin,n1)             :: new_mtend

      ! parameter of distribution in each bin
      real(8), dimension(mxnbin+1,n1,4)           :: a2d

!     shifted bounds of a bin, number of parameters for subdistribution
      real(RP),dimension(mxnbin+1,n1,2) :: binb3d
      ! mass tendency of shifted bin
      real(PS), dimension(mxnbin+1,n1)                     :: mtend
      real(RP),dimension(mxnbin+1,n1,mxnmasscomp) :: ratio_Mp

      ! axis ratio of an ice crystal in a shifted bin
      ! 1: c/a, 2: d/a, 3: r/a, 4: e/a
      real(PS),dimension(mxnbin+1,n1,mxnaxis-1) :: axr_p
      ! bulk sphere density of dry ice particle, and bulk crystal density in the shifted bin.
      real(PS), dimension(mxnbin+1,n1)           :: den_ip_p,den_ic_p
      ! habit in shifted bin
      integer, dimension(mxnbin+1,n1)          :: habit_p
      ! aspect ratio of circumscribing cylinder (not used)
      real(PS)              :: spx_p

      ! aspect ratio of ice particle
      real(PS), dimension(mxnbin+1,n1)              :: asr_p
      ! type in shifted bin
      integer, dimension(mxnbin+1,n1)          :: type_p


      ! ratio of ag^3 to a^3
      real(PS), dimension(mxnbin+1,n1)           :: rag_p,rcg_p
      ! number of extra ice crystals
      real(RP), dimension(mxnbin+1,n1)           :: n_exice_p

      ! activated IN fraction for contact parameter diagnosis
      real(PS), dimension(mxnbin+1,n1)           :: actINF_p


      integer,dimension(mxnbin+1,n1)           :: error_number

      integer,dimension(mxnbin+1,n1)           :: ierror1

      real(PS),dimension(mxnbin,n1) :: mmass_p
      real(PS),dimension(mxnbin,n1) :: rmod
      real(PS) :: v_sp
      real(8),dimension(n1) :: ap_dN,ap_dM
!       4*pi/3
      real(PS) :: coef3
      parameter(coef3=4.18879020478639)
!     error message
      integer :: em,icem,ncem
      integer :: ia,im,iq,k,ibr,icr,ik,j,jk
      real(PS),dimension(n1) :: den1
!     minimum fraction of aerosols in hydrometeors
      real(RP) :: min_fapt_r,min_faps_r,min_fapt_s,min_faps_s
      parameter(min_fapt_r=1.0e-18_RP,min_faps_r=1.0e-6_RP &
               ,min_fapt_s=1.0e-18_RP,min_faps_s=0.0_RP)

!     possible ratio of mean mass to bin boundaries
      real(RP) :: brat1,brat2
      parameter(brat1=1.001_RP,brat2=0.999_RP)

!     factor
      real(RP) :: frac,fct1,fct2
      parameter(frac=0.1_RP,fct1=0.1_RP,fct2=0.99_RP)

      real(PS) :: dum,max_rad,rat

!     minimum possible volume of ice crystalc
      real(RP),parameter :: V_csmin=1.1847688e-11_RP

      real(PS) :: sum_m1(max(nprmx,npimx,npamx),n1) &
             ,sum_m2(max(nprmx,npimx,npamx),n1)

      integer,dimension(n1) :: icond1
      integer,dimension(mxnbin*n1) :: icond2,icond3

      real(PS),parameter :: amin_mmass=1.0e-25
      ! bin boundary modification factor.
      real(PS),parameter    :: bbmf=0.2
!      write(fid_alog,*) "I am in shifted bin",npr,nbr,ncr,iphase

!!!      write(fid_alog,*) "k1r,k2r",k1r,k2r

!!!      return

      if(k1r>k2r) return

!   initialization

      icem=0
      ncem=0


      do k=1,n1
        den1(k)=den(k)*1.0e-3_PS

        !  these are not used.
        ap_dN(k)=0.0d+0
        ap_dM(k)=0.0d+0

        if(k>=k1r.and.k<=k2r) then
          icond1(k)=1
        else
          icond1(k)=0
        endif
      enddo

      if(iphase==1) then
!
!   in case of liquid spectrum
!
!

        do icr=1,ncr
!
!     initialization
!
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr
            new_N(ibr,k)=0.0_DS
            new_mtend(ibr,k)=0.0_DS
            rmod(ibr,k)=1.0
          enddo
          do j=1,nvar_mcp_liq+1
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              new_M(ibr,k,j)=0.0_DS
            end do
          enddo
          do j=1,nvar_nonmcp_liq
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              new_Q(ibr,k,j)=0.0_DS
            end do
          enddo

          do jk=1,(npr-2)*n1
            k=(jk-1)/(npr-2)+1
            j=jk-(k-1)*(npr-2)
            sum_m1(j,k)=0.0_PS
            sum_m2(j,k)=0.0_PS
          end do

          do ibr=1,nbr
!CDIR NODEP
            do jk=1,(npr-2)*n1
              k=(jk-1)/(npr-2)+1
              j=jk-(k-1)*(npr-2)
              sum_m1(j,k)=sum_m1(j,k)+qrpv(j,ibr,icr,k)*den1(k)
            end do
          end do

          !
          ! initialize
          !
!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond1(k)==1.and.&
               qrpv(rmt_q,ibr,icr,k)>1.0e-30_RP.and. &
               qrpv(rcon_q,ibr,icr,k)>1.0e-30_RP) then
              icond3(ik)=1
              mmass_p(ibr,k)=qrpv(rmt_q,ibr,icr,k)/qrpv(rcon_q,ibr,icr,k)
            else
              icond3(ik)=0
              mmass_p(ibr,k)=0.0_RP
            endif
          enddo

!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond3(ik)==1.and.mmass(ibr,icr,k)>=amin_mmass) then


!     calculate mass components ratio
!         for aerosols
              ratio_Mp(ibr,k,rmat_m)=min( &
                 max(qrpv(rmat_q,ibr,icr,k)/qrpv(rmt_q,ibr,icr,k) &
                        ,min_fapt_r),1.0_RP)

              ratio_Mp(ibr,k,rmas_m)=min( &
                 max(qrpv(rmas_q,ibr,icr,k)/qrpv(rmt_q,ibr,icr,k) &
                        ,min_faps_r),ratio_Mp(ibr,k,rmat_m))

              actINF_p(ibr,k)=0.0_RP

              binb3d(ibr,k,1)=max(0.0_RP,binbr(ibr) &
                                 +mmass_p(ibr,k)-mmass(ibr,icr,k))
              binb3d(ibr,k,2)=max(0.0_RP,binbr(ibr+1) &
                                 +mmass_p(ibr,k)-mmass(ibr,icr,k))


              Np(ibr,k)=qrpv(rcon_q,ibr,icr,k)*den1(k)
              Mp(ibr,k)=qrpv(rmt_q,ibr,icr,k)*den1(k)
              mtend(ibr,k)=0.0_RP
            else
              ratio_Mp(ibr,k,rmat_m)=0.0_RP
              ratio_Mp(ibr,k,rmas_m)=0.0_RP
              actINF_p(ibr,k)=0.0_RP

              binb3d(ibr,k,1)=binbr(ibr)
              binb3d(ibr,k,2)=binbr(ibr+1)
!org              Np(ibr,k)=0.0_RP
!org              Mp(ibr,k)=0.0_RP
              Np(ibr,k)=max(0.0_RP,qrpv(rcon_q,ibr,icr,k)*den1(k))
              Mp(ibr,k)=max(0.0_RP,qrpv(rmt_q,ibr,icr,k)*den1(k))
              mtend(ibr,k)=0.0_RP
            endif
          enddo

          call cal_lincubprms_vec(mxnbin+1,nbr,n1,Np,Mp,binb3d  &
                              ,a2d,error_number,"shift_liq")

!!          do k=1,n1
!!            write(fid_alog,*) "error_number af linprms",kmic(k),imic(k),jmic(k)
!!            write(fid_alog,*) (error_number(ibr,k),ibr=1,nbr)
!!          enddo

          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

!tmp            if(icond1(k)==1) then
!tmp              if(binb3d(ibr,k,1)<1.0e-20_PS.or.binb3d(ibr,k,2)<1.0e-20_PS) then
!tmp              write(fid_alog,*) "error_number",k,ibr,error_number(ibr,k),&
!tmp                   qrpv(rcon_q,ibr,icr,k)*den1(k),qrpv(rmt_q,ibr,icr,k)*den1(k),&
!tmp                   mmass_p(ibr,k),mmass(ibr,icr,k),a2d(ibr,k,1:3),binb3d(ibr,k,1:2)
!tmp            endif
!tmp            endif
            if(1<=error_number(ibr,k).and.error_number(ibr,k)<=4) then
              icem=icem+1
              ncem=ncem+1
            elseif(error_number(ibr,k)==0) then
              ncem=ncem+1
            endif
          enddo

          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond3(ik)==1.and.&
               ( (1<=error_number(ibr,k).and.error_number(ibr,k)<=4).or.&
                 mmass(ibr,icr,k)<amin_mmass) ) then
              icond2(ik)=1
            else
              icond2(ik)=0
            endif
          enddo

!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond2(ik)==1) then

              new_M(ibr,k,rmt)=new_M(ibr,k,rmt) &
                                +qrpv(rmt_q,ibr,icr,k)*den1(k)


              if(mmass_p(ibr,k)>binbr(ibr+1) &
                                .or.mmass_p(ibr,k)<binbr(ibr))then
                rmod(ibr,k)=qrpv(rmt_q,ibr,icr,k) &
                      /max(brat1*binbr(ibr) &
                      ,min(brat2*binbr(ibr+1),mmass_p(ibr,k))) &
                      /qrpv(rcon_q,ibr,icr,k)
              else
                rmod(ibr,k)=1.0
              end if

              new_N(ibr,k)=new_N(ibr,k) &
                             +qrpv(rcon_q,ibr,icr,k)*den1(k)*rmod(ibr,k)


              ! this make these cases to skip in trans_bin
              error_number(ibr,k)=10
            endif
          enddo

!     for mass components
          do iq=i1_mcp_liq,i2_mcp_liq
            im=rmt+1+iq-i1_mcp_liq
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond2(ik)==1) then

                new_M(ibr,k,im)=new_M(ibr,k,im) &
                             +qrpv(iq,ibr,icr,k)*den1(k)
              endif
            enddo
          enddo

!     for concentration component
          do iq=i1_ccp_liq,i2_ccp_liq
            im=iq-i1_ccp_liq+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond2(ik)==1) then
                new_Q(ibr,k,im)=new_Q(ibr,k,im) &
                         +qrpv(iq,ibr,icr,k)*den1(k)*rmod(ibr,k)
              endif
            end do
          end do

          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! calculation of transferred concentration and mass into original bins
          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          call cal_transbin_vec(iphase &
                           ,n1,nvar_mcp_liq &
                           ,nbr,nbr &
                           ,binbr &
                           ,error_number &
                           ,a2d,binb3d,mtend &
                           ,new_N,new_M,new_Q &
                           ,new_mtend &
                           ,ratio_Mp,den_ip_p,axr_p,spx_p &
                           ,habit_p,den_ic_p &
                           ,rag_p,rcg_p,n_exice_p &
                           ,actINF_p &
                           ,0)
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!     put new variables back into the originals
!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond1(k)==1) then
              qrpv(rmt_q,ibr,icr,k)=max(new_M(ibr,k,rmt)/den1(k),0.0_DS)
              qrpv(rcon_q,ibr,icr,k)=max(new_N(ibr,k)/den1(k),0.0_DS)
!dbg              write(fid_alog,'("af cal_transbin",2I5,2ES15.6)') &
!dbg                     k,ibr,new_N(ibr,k),new_M(ibr,k,rmt)
            endif
          enddo

! mass component
          do iq=i1_mcp_liq,i2_mcp_liq
            im=rmt+1+iq-i1_mcp_liq
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              if(icond1(k)==1) then
!tmp                write(fid_alog,*) "ck new_M",k,ibr,iq,den(k),new_M(ibr,qrpv,im)
                qrpv(iq,ibr,icr,k)=max(new_M(ibr,k,im)/den1(k),0.0_DS)
              endif
            enddo
          enddo

!   for concentration component
          do iq=i1_ccp_liq,i2_ccp_liq
            im=iq-i1_ccp_liq+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond1(k)==1) then
                qrpv(iq,ibr,icr,k)= &
                      max(new_Q(ibr,k,im)/den1(k),0.0_DS)
              endif
            enddo
          enddo

          do ibr=1,nbr
!CDIR NODEP
            do jk=1,(npr-2)*n1
              k=(jk-1)/(npr-2)+1
              j=jk-(k-1)*(npr-2)
              sum_m2(j,k)=sum_m2(j,k)+qrpv(j,ibr,icr,k)*den1(k)
            end do
          end do

          if(debug) then
            do k=1,n1
              if(sum_m1(1,k)>1.0e-30) then
                if(abs((sum_m2(1,k)-sum_m1(1,k))/sum_m1(1,k))>1.0e-4) then
                  write(fid_alog,*) "bf liq:sum", &
                    kmic(k),imic(k),jmic(k),sum_m1(1:npr-2,k)
                  write(fid_alog,*) "af liq:sum", &
                    kmic(k),imic(k),jmic(k),sum_m2(1:npr-2,k)
                  write(fid_alog,*) "np",(np(ibr,k),ibr=1,nbr)
                  write(fid_alog,*) "mp",(mp(ibr,k),ibr=1,nbr)
                  write(fid_alog,*) "error_number",(error_number(ibr,k),ibr=1,nbr)
                  write(fid_alog,*) "mmass",(mmass(ibr,icr,k),ibr=1,nbr)
                  write(fid_alog,*) "mmass_p",(mmass_p(ibr,k),ibr=1,nbr)
                  write(fid_alog,*) "rmod",(rmod(ibr,k),ibr=1,nbr)
                  write(fid_alog,*) "a2d(1)",(a2d(ibr,k,1),ibr=1,nbr)
                  write(fid_alog,*) "a2d(2)",(a2d(ibr,k,2),ibr=1,nbr)
                  write(fid_alog,*) "a2d(3)",(a2d(ibr,k,3),ibr=1,nbr)
                end if
              endif
            enddo
          endif

!     calculate mean mass for cloud droplet bin
!CDIR NODEP
          do k=1,n1
            mmass(1,icr,k)= &
                qrpv(rmt_q,1,icr,k)/max(1.0e-30_RP,qrpv(rcon_q,1,icr,k))
          enddo
        end do

        if(debug) then
          if(ncem>0) then
            write(fid_alog,301) 100.0_PS*real(icem,PS_KIND)/real(ncem,PS_KIND)
301   format("error % in liq shift_bin:",ES15.6)
          end if
        end if



      elseif(iphase==2) then
!
!  in case of ice spectrum
!

        do icr=1,ncr
!
!     initialization
!
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr
            new_N(ibr,k)=0.0_PS
            new_mtend(ibr,k)=0.0_PS
            rmod(ibr,k)=1.0
          enddo
          do j=1,nvar_mcp_ice+1  ! including total mass
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              new_M(ibr,k,j)=0.0_PS
            end do
          enddo
          do j=1,nvar_nonmcp_ice
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              new_Q(ibr,k,j)=0.0_PS
            end do
          enddo

          do jk=1,(npr-2)*n1
            k=(jk-1)/(npr-2)+1
            j=jk-(k-1)*(npr-2)
            sum_m1(j,k)=0.0
            sum_m2(j,k)=0.0
          end do

          do ibr=1,nbr
!CDIR NODEP
            do jk=1,(npr-2)*n1
              k=(jk-1)/(npr-2)+1
              j=jk-(k-1)*(npr-2)
              sum_m1(j,k)=sum_m1(j,k)+qrpv(j,ibr,icr,k)*den1(k)
            end do
          end do

          !
          ! initialize
          !
!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            mmass_p(ibr,k)=0.0
            error_number(ibr,k)=0

            if(icond1(k)==1.and.&
               qrpv(imt_q,ibr,icr,k)>1.0e-30.and. &
               qrpv(icon_q,ibr,icr,k)>1.0e-30) then

              icond3(ik)=1
              mmass_p(ibr,k)=qrpv(imt_q,ibr,icr,k)/qrpv(icon_q,ibr,icr,k)
            else
              icond3(ik)=0
              mmass_p(ibr,k)=0.0

            endif
          enddo

!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond3(ik)==1.and.mmass(ibr,icr,k)>=amin_mmass) then
!
!     calculate mass components ratio
!
!         for ice crystal mass
              ratio_Mp(ibr,k,imc_m)=min( &
                      max(max( &
                     qrpv(icon_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k)*4.763209003e-12 &
                    ,qrpv(imc_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k)) &
                      ,0.0_RP),1.0_RP)
!         for rimed mass
              ratio_Mp(ibr,k,imr_m)=min( &
                      max(qrpv(imr_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k) &
                     ,0.0_RP),1.0_RP)

!         for melt water mass
              ratio_Mp(ibr,k,imw_m)=min( &
                      max(qrpv(imw_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k) &
                     ,0.0_RP),1.0_RP)

!         for aggregation mass
              ratio_Mp(ibr,k,ima_m)=min( &
                      max(qrpv(ima_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k) &
                     ,0.0_RP),1.0_RP)

!         for freezing nucleation mass
              ratio_Mp(ibr,k,imf_m)=min( &
                      max(qrpv(imf_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k) &
                     ,0.0_RP),ratio_Mp(ibr,k,imc_m))

!         normalization
              dum = &
                ratio_Mp(ibr,k,imr_m)+ratio_Mp(ibr,k,ima_m)+ratio_Mp(ibr,k,imw_m)

              if(dum>1.0e-20_RP) then
                ratio_Mp(ibr,k,imr_m)= &
                      (1.0_RP-ratio_Mp(ibr,k,imc_m))*ratio_Mp(ibr,k,imr_m)/dum

                ratio_Mp(ibr,k,ima_m)= &
                      (1.0_RP-ratio_Mp(ibr,k,imc_m))*ratio_Mp(ibr,k,ima_m)/dum

                ratio_Mp(ibr,k,imw_m)= &
                      (1.0_RP-ratio_Mp(ibr,k,imc_m))*ratio_Mp(ibr,k,imw_m)/dum
              else
                ratio_Mp(ibr,k,imc_m)=1.0_RP
                ratio_Mp(ibr,k,imr_m)=0.0_RP
                ratio_Mp(ibr,k,ima_m)=0.0_RP
                ratio_Mp(ibr,k,imw_m)=0.0_RP
              end if

              ratio_Mp(ibr,k,imf_m)=min(ratio_Mp(ibr,k,imf_m) &
                                       ,ratio_Mp(ibr,k,imc_m))

!                     if(ratio_Mp(imc_m)<0.5) then
!                        write(fid_alog,*) "ratio_Mp(imc_m) is less than 0.5 at"
!                        write(fid_alog,*) ibr,k
!                        write(fid_alog,*) ratio_Mp(imr_m),ratio_Mp(ima_m),ratio_Mp(imc_m)
!                     end if

!         for aerosol total mass
              ratio_Mp(ibr,k,imat_m)=min( &
                      max(qrpv(imat_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k) &
                     ,min_fapt_s),1.0_RP)

!         for soluble aerosol mass
              ratio_Mp(ibr,k,imas_m)=min( &
                      max(qrpv(imas_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k) &
                     ,min_faps_s),ratio_Mp(ibr,k,imat_m))

!         for insoluble aerosol mass
!ccc                        ratio_Mp(imai_m)=max(0.0_RP,min(1.0_RP
!ccc     *                             ,ratio_Mp(imat_m)-ratio_Mp(imas_m)))

!ccc                     write(fid_alog,'("shift_bin > ratio_Mp",2I5,10ES15.6)')
!ccc     *                      ibr,k,ratio_Mp

!         for circumscribing volume
              Qp(ibr,k,ivcs)= &
                      max(qrpv(ivcs_q,ibr,icr,k)/qrpv(icon_q,ibr,icr,k) &
                                 ,mmass_p(ibr,k)/den_i)

            else
!         for crystal mass
              ratio_Mp(ibr,k,imc_m)=0.0_RP
!         for rimed mass
              ratio_Mp(ibr,k,imr_m)=0.0_RP
!         for melt water mass
              ratio_Mp(ibr,k,imw_m)=0.0_RP
!         for aggregation mass
              ratio_Mp(ibr,k,ima_m)=0.0_RP
!         for freezing nucleation mass
              ratio_Mp(ibr,k,imf_m)=0.0_RP

!         for aerosol total mass
              ratio_Mp(ibr,k,imat_m)=0.0_RP
!         for soluble aerosol mass
              ratio_Mp(ibr,k,imas_m)=0.0_RP

!         for circumscribing volume
              Qp(ibr,k,ivcs)=0.0_RP

            endif
          enddo


!         for axis component
          do iq=i1_acp_ice,i2_acp_ice
            im=nvar_vcp_ice+iq-i1_acp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond3(ik)==1.and.mmass(ibr,icr,k)>=amin_mmass) then
                Qp(ibr,k,im)=(max(qrpv(iq,ibr,icr,k)/ &
                        qrpv(icon_q,ibr,icr,k),0.0_RP))**(1.0_RP/3.0_RP)
              else
                Qp(ibr,k,im)=0.0_RP
              endif
            enddo
          enddo

!         for concentration component
          do iq=i1_ccp_ice,i2_ccp_ice
            im=nvar_vcp_ice+nvar_acp_ice+iq-i1_ccp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond3(ik)==1.and.mmass(ibr,icr,k)>=amin_mmass) then
                Qp(ibr,k,im)=qrpv(iq,ibr,icr,k)/ &
                       qrpv(icon_q,ibr,icr,k)
              else
                Qp(ibr,k,im)=0.0_RP
              endif
            enddo
          enddo

          ierror1=0
!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond3(ik)==1.and.mmass(ibr,icr,k)>=amin_mmass) then

!         Realitiy check
!!c              if(Qp(ibr,k,ivcs)<1.1847688E-11.or. &
!!c                 Qp(ibr,k,iacr)<1.0e-4.or. &
!!c                 Qp(ibr,k,iccr)<1.0e-4) then
!!c
!!c                ! case of less than 1 micron hex ice crystal
!!c
!!c                Qp(ibr,k,ivcs)=1.1847688E-11
!!c                Qp(ibr,k,iacr)=1.0e-4
!!c                Qp(ibr,k,iccr)=1.0e-4
!!c
!!c              end if
              if(Qp(ibr,k,ivcs)<=1.1847688E-11) then
                Qp(ibr,k,ivcs)=1.1847688E-11
              end if
              if(Qp(ibr,k,iacr)<1.0e-4) then
                Qp(ibr,k,iacr)=1.0e-4
              end if
              if(Qp(ibr,k,iccr)<1.0e-4) then
                Qp(ibr,k,iccr)=1.0e-4
              end if

              Qp(ibr,k,idcr)=min(Qp(ibr,k,iacr),max(0.0_RP,Qp(ibr,k,idcr)))
              if(Qp(ibr,k,iacr)<=1.0e-4_RP) then
                Qp(ibr,k,idcr)=0.0_RP
              end if
              Qp(ibr,k,iag)=max(0.0_RP,Qp(ibr,k,iag))
              Qp(ibr,k,icg)=max(0.0_RP,Qp(ibr,k,icg))
              Qp(ibr,k,inex)=min(1.0_RP,max(0.0_RP,Qp(ibr,k,inex)))

              ! c and d length ratio to the a length
              ia=1
              axr_p(ibr,k,ia)=Qp(ibr,k,iacr+ia)/Qp(ibr,k,iacr)
              ia=2
              axr_p(ibr,k,ia)=Qp(ibr,k,iacr+ia)/Qp(ibr,k,iacr)

!                        if(axr_p(2)>0.5) then
!                           write(fid_alog,*) k,ibr,axr_p(1),axr_p(2)
!                        end if

!     pristine crystal
!         ratio for center of gravities
              rag_p(ibr,k)=(Qp(ibr,k,iag)/Qp(ibr,k,iacr))**3
              rcg_p(ibr,k)=(Qp(ibr,k,icg)/Qp(ibr,k,iccr))**3

!         number of extra crystals
              n_exice_p(ibr,k)=Qp(ibr,k,inex)
!         activated IN concentration fraction
              actINF_p(ibr,k)=0.0_RP

!                          hexagonal crystals
              v_sp=coef3*((1.0+axr_p(ibr,k,1)**2)**1.5_RP)  &
                                     *Qp(ibr,k,iacr)**3
              habit_p(ibr,k)=1
              spx_p=axr_p(ibr,k,1)

              den_ic_p(ibr,k)=max(1.0e-4_RP,&
                       mmass_p(ibr,k)*ratio_Mp(ibr,k,imc_m)/v_sp)

              axr_p(ibr,k,3)=0.0_RP
              axr_p(ibr,k,4)=0.0_RP

              if(den_ic_p(ibr,k)>den_i) then
                v_sp=mmass_p(ibr,k)*ratio_Mp(ibr,k,imc_m)/den_i

                den_ic_p(ibr,k)=den_i

                Qp(ibr,k,iacr)=(v_sp/ &
                     (coef3*((1.0_RP+spx_p**2)**1.5_RP)))**(1.0_RP/3.0_RP)
                Qp(ibr,k,iccr)=Qp(ibr,k,iacr)*spx_p
              end if

              Qp(ibr,k,ivcs)=max(Qp(ibr,k,ivcs),V_csmin)

              den_ip_p(ibr,k)=mmass_p(ibr,k)/Qp(ibr,k,ivcs)

              binb3d(ibr,k,1)=max(0.0_RP,binbr(ibr) &
                             +mmass_p(ibr,k)-mmass(ibr,icr,k))
              binb3d(ibr,k,2)=max(0.0_RP,binbr(ibr+1) &
                             +mmass_p(ibr,k)-mmass(ibr,icr,k))

              Np(ibr,k)=qrpv(icon_q,ibr,icr,k)*den1(k)
              Mp(ibr,k)=qrpv(imt_q,ibr,icr,k)*den1(k)
              mtend(ibr,k)=0.0_RP

            else

              ! c and d length ratio to the a length
              ia=1
              axr_p(ibr,k,ia)=1.0_RP
              ia=2
              axr_p(ibr,k,ia)=1.0_RP

!     pristine crystal
!         ratio for center of gravities
              rag_p(ibr,k)=0.0_RP
              rcg_p(ibr,k)=0.0_RP

!         number of extra crystals
              n_exice_p(ibr,k)=0.0_RP
!         activated IN concentration fraction
              actINF_p(ibr,k)=0.0_RP

!                          hexagonal crystals
              habit_p(ibr,k)=0
              spx_p=0.0

              den_ic_p(ibr,k)=1.0_RP

              axr_p(ibr,k,3)=0.0_RP
              axr_p(ibr,k,4)=0.0_RP

              den_ip_p(ibr,k)=1.0_RP


              binb3d(ibr,k,1)=binbr(ibr)
              binb3d(ibr,k,2)=binbr(ibr+1)
              Np(ibr,k)=max(0.0_RP,qrpv(icon_q,ibr,icr,k)*den1(k))
              Mp(ibr,k)=max(0.0_RP,qrpv(imt_q,ibr,icr,k)*den1(k))
              mtend(ibr,k)=0.0_RP

            endif

!!!            if(den_ic_p(ibr,k)<1.0e-30_RP.and.np(ibr,k)>1.0e-30_RP.and.&
!!!               mp(ibr,k)>1.0e-30_RP) then
            if(den_ic_p(ibr,k)<1.0e-30_RP) then
              ierror1(ibr,k)=1
            endif
          enddo
          if(any(ierror1>0)) then
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              if(debug .and. ierror1(ibr,k)==1) then
                write(fid_alog,*) "Shift_bin_vec: Error1:",k,ibr,icond3(ik),den_ic_p(ibr,k),np(ibr,k),mp(ibr,k),&
                      mmass(ibr,icr,k),mmass_p(ibr,k),qrpv(imt_q,ibr,icr,k),qrpv(icon_q,ibr,icr,k)
              endif
            enddo
            stop
          endif

          call cal_lincubprms_vec(mxnbin+1,nbr,n1,Np,Mp,binb3d  &
                              ,a2d,error_number,"shift_ice")


          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(1<=error_number(ibr,k).and.error_number(ibr,k)<=4) then
              icem=icem+1
              ncem=ncem+1
            elseif(error_number(ibr,k)==0) then
              ncem=ncem+1
            endif
          enddo

          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond3(ik)==1.and.&
               ( (1<=error_number(ibr,k).and.error_number(ibr,k)<=4).or.&
                 mmass(ibr,icr,k)<amin_mmass) ) then
              icond2(ik)=1
            else
              icond2(ik)=0
            endif
          enddo

!          do ik=1,nbr*n1
!            k=(ik-1)/nbr+1
!            ibr=ik-(k-1)*nbr
!            write(fid_alog,*) "Shift_bin_vec: Error1-1:",k,ibr,icond3(ik),icond2(ik),den_ic_p(ibr,k),den_ip_p(ibr,k),&
!                      mmass(ibr,icr,k),mmass_p(ibr,k),&
!                      qrpv(imt_q,ibr,icr,k),qrpv(icon_q,ibr,icr,k),error_number(ibr,k),a2d(ibr,k,1:3),binb3d(ibr,k,1:2)
!          enddo


!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond2(ik)==1) then

              new_M(ibr,k,imt)=new_M(ibr,k,imt) &
                           +qrpv(imt_q,ibr,icr,k)*den1(k)

              if(mmass_p(ibr,k)>binbr(ibr+1) &
                 .or.mmass_p(ibr,k)<binbr(ibr))then

                rmod(ibr,k)=qrpv(imt_q,ibr,icr,k) &
                     /max(brat1*binbr(ibr) &
                     ,min(brat2*binbr(ibr+1),mmass_p(ibr,k))) &
                     /qrpv(icon_q,ibr,icr,k)
              else
                rmod(ibr,k)=1.0_RP
              end if

              new_N(ibr,k)=new_N(ibr,k) &
                       +qrpv(icon_q,ibr,icr,k)*den1(k)*rmod(ibr,k)

              ! this make these cases to skip in trans_bin
              error_number(ibr,k)=10
            endif
          enddo

!     for volume component
          do iq=i1_vcp_ice,i2_vcp_ice
            im=iq-i1_vcp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond2(ik)==1) then
                new_Q(ibr,k,im)=new_Q(ibr,k,im) &
                      +qrpv(iq,ibr,icr,k)*den1(k)*rmod(ibr,k)
              endif
            end do
          enddo

!     for axis component
          do iq=i1_acp_ice,i2_acp_ice
            im=iq-i1_vcp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond2(ik)==1) then
                new_Q(ibr,k,im)=new_Q(ibr,k,im) &
                      +qrpv(iq,ibr,icr,k)*den1(k)*rmod(ibr,k)
              endif
            enddo
          enddo

!     for concentration component
          do iq=i1_ccp_ice,i2_ccp_ice
            im=iq-i1_vcp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond2(ik)==1) then
                new_Q(ibr,k,im)=new_Q(ibr,k,im) &
                         +qrpv(iq,ibr,icr,k)*den1(k)*rmod(ibr,k)
              endif
            end do
          end do

          do iq=i1_mcp_ice,i2_mcp_ice
            im=imt+1+iq-i1_mcp_ice
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond2(ik)==1) then
                new_M(ibr,k,im)=new_M(ibr,k,im) &
                                  +qrpv(iq,ibr,icr,k)*den1(k)
              endif
            enddo
          enddo

!          do ik=1,nbr*n1
!            k=(ik-1)/nbr+1
!            ibr=ik-(k-1)*nbr
!
!            ierror1(ibr,k)=0
!!!!            if(icond2(ik)==0.and.error_number(ibr,k)/=10) then
!            if(error_number(ibr,k)/=10) then
!!              if(den_ic_p(ibr,k)<1.0e-30) then
!!                write(fid_alog,*) "Shift_bin_vec: Error2:",k,ibr,icond2(ik),icond3(ik),mmass(ibr,icr,k),mmass_p(ibr,k),&
!                           den_ic_p(ibr,k),np(ibr,k),mp(ibr,k),&
!!                           error_number(ibr,k)
!!!                ierror1(ibr,k)=1
!!              endif
!              if(den_ip_p(ibr,k)<1.0e-30.or.den_ip_p(ibr,k)>1.0e-30) then
!                ierror1(ibr,k)=1
!              endif
!            endif
!          enddo
!          if(any(ierror1>0)) then
!            do ik=1,nbr*n1
!              k=(ik-1)/nbr+1
!              ibr=ik-(k-1)*nbr
!              if(ierror1(ibr,k)==1) then
!!                write(fid_alog,*) "Shift_bin_vec: Error2:",k,ibr,icond2(ik),icond3(ik),mmass(ibr,icr,k),mmass_p(ibr,k),&
!!                           den_ic_p(ibr,k),np(ibr,k),mp(ibr,k),&
!!                           error_number(ibr,k)
!                write(fid_alog,*) "Shift_bin_vec: Error2:",k,ibr,icond2(ik),icond3(ik),mmass(ibr,icr,k),mmass_p(ibr,k),&
!                           den_ip_p(ibr,k),np(ibr,k),mp(ibr,k),&
!                           error_number(ibr,k)
!              endif
!            enddo
!            stop
!          endif

          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! calculation of transferred concentration and mass into original bins
          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          call cal_transbin_vec(iphase &
                           ,n1,nvar_mcp_ice &
                           ,nbr,nbr &
                           ,binbr &
                           ,error_number &
                           ,a2d,binb3d,mtend &
                           ,new_N,new_M,new_Q &
                           ,new_mtend &
                           ,ratio_Mp,den_ip_p,axr_p,spx_p &
                           ,habit_p,den_ic_p &
                           ,rag_p,rcg_p,n_exice_p &
                           ,actINF_p &
                           ,0)
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            ierror1(ibr,k)=0
!!!            if(icond2(ik)==0.and.error_number(ibr,k)/=10) then
            if(error_number(ibr,k)/=10) then
              if(den_ic_p(ibr,k)<1.0e-30) then
                ierror1(ibr,k)=1
              endif
            endif
          enddo
          if(any(ierror1>0)) then
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              if(debug .and. ierror1(ibr,k)==1) then
                write(fid_alog,*) "Shift_bin_vec: Error3:",k,ibr,icond2(ik),icond3(ik),mmass(ibr,icr,k),mmass_p(ibr,k),&
                           den_ic_p(ibr,k),np(ibr,k),mp(ibr,k),&
                           error_number(ibr,k)
              endif
            enddo
            stop
          endif

!     put new variables back into the originals
!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond1(k)==1) then
! total mass
              qrpv(imt_q,ibr,icr,k)=max(new_M(ibr,k,imt)/den1(k),0.0_DS)
! concentration
              qrpv(icon_q,ibr,icr,k)=max(new_N(ibr,k)/den1(k),0.0_DS)

            endif
          enddo

! mass component
          do iq=i1_mcp_ice,i2_mcp_ice
            im=imt+1+iq-i1_mcp_ice
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond1(k)==1) then
                qrpv(iq,ibr,icr,k)=max(new_M(ibr,k,im)/den1(k),0.0_DS)
              endif
            enddo
          enddo


!  for volume component
          do iq=i1_vcp_ice,i2_vcp_ice
            im=iq-i1_vcp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond1(k)==1) then
                qrpv(iq,ibr,icr,k)= &
                      max(new_Q(ibr,k,im)/den1(k),0.0_DS)
              endif
            end do
          end do

!   for axis component
          do iq=i1_acp_ice,i2_acp_ice
            im=iq-i1_vcp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond1(k)==1) then
                qrpv(iq,ibr,icr,k)= &
                      max(new_Q(ibr,k,im)/den1(k),0.0_DS)
              endif
            enddo
          enddo

!   for concentration component
          do iq=i1_ccp_ice,i2_ccp_ice
            im=iq-i1_vcp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond1(k)==1) then
                qrpv(iq,ibr,icr,k)= &
                      max(new_Q(ibr,k,im)/den1(k),0.0_DS)
              endif
            enddo
          enddo

          do ibr=1,nbr
!CDIR NODEP
            do jk=1,(npr-2)*n1
              k=(jk-1)/(npr-2)+1
              j=jk-(k-1)*(npr-2)
              sum_m2(j,k)=sum_m2(j,k)+qrpv(j,ibr,icr,k)*den1(k)
            end do
          end do

!cc          if(sum_m1(1)>1.0e-30) then
!cc            write(fid_alog,'("af qrpv",I5,30ES15.6)') k
!ccc     *               ,qrpv(imt_q,1:nbr,icr,k)
!ccc               end if

!!c          do k=1,n1
!!c            if(sum_m1(1,k)>1.0e-30) then
!!c              if(abs((sum_m2(1,k)-sum_m1(1,k))/sum_m1(1,k))>1.0e-4) then
!!c                write(fid_alog,*) "bf ice:sum", &
!!c                  kmic(k),imic(k),jmic(k),sum_m1(1:npr-2,k)
!!c                write(fid_alog,*) "af ice:sum", &
!!c                  kmic(k),imic(k),jmic(k),sum_m2(1:npr-2,k)
!!c              end if
!!c            end if
!!c          enddo

        end do
!        if(ncem>0) then
!           write(fid_alog,302) 100.0*real(icem)/real(ncem)
!302   format("error % in ice shift_bin:",ES15.6)
!        end if

      end if
      end subroutine shift_bin_vec

      subroutine negadj_bin_ap(n1,npa,nba,nca &
           ,npr,nbr,ncr,npi,nbi,nci &
           ,i,kk &
           ,qa,qr,qi,DEN &
           ,ircnfl,ipcnfl &
           ,imicv,jmicv,kmicv)
        use scale_prc, only: &
           PRC_abort
        use maxdims
        use com_amps
        use par_amps
        implicit none
      integer :: npa,nba,nca,npr,nbr,ncr,npi,nbi,nci &
             ,ircnfl,ipcnfl
      integer :: k,ka,kb,k1,n1,i,kk
      integer :: imicv(*),jmicv(*),kmicv(*),ica
      real(RP) :: qa(npa,nba,nca,*),qr(npr,nbr,ncr,*),qi(npi,nbi,nci,*)
      real(RP) :: DEN(*)
!     This routine assures all positive definites remain positive
!
!     The minimum possible background concentration for large particle
!     is assumed to be 1.0e-5 cm^-3
      real(PS) :: tconlmt
      parameter(tconlmt=1.0e-5)
!     The minimum radius possible for large particles
      real(PS) :: r3_lmt
      parameter(r3_lmt=1.0e-18)
!     4.0*pi/3.0
      real(PS) :: coef3
      parameter(coef3=4.18879020478639)
!     minimum possible fractions
      real(PS) :: min_fapt_r,min_faps_r,min_fapt_s,min_faps_s,sep_faps
      parameter(min_fapt_r=1.0e-18,min_faps_r=1.0e-5 &
!      parameter(min_fapt_r=1.0e-18,min_faps_r=0.0 &
               ,min_fapt_s=1.0e-18,min_faps_s=0.0 &
               ,sep_faps=1.0e-5)
      real(RP) :: mx_aprat
      parameter(mx_aprat=0.99999_RP)
      real(RP) :: m_lmt,deps_ap,den_ap,rmod1,rmod2
      real(PS) :: oapt_a(n1),oapt_r(n1),oapt_i(n1),den1
      real(PS) :: mapt_a(n1),mapt_r(n1),mapt_i(n1)
!     calcualte total AP mass

      k1=1
      IF(IRCNFL.GT.0) THEN
         DO K=k1,N1
            den1=den(k)*1.0e-3

            mapt_r(K)=0.0
            oapt_r(k)=0.0
            do ka=1,nbr
               oapt_r(k)=oapt_r(k)+qr(rmat_q,ka,1,k)
!               m_lmt=coef3*r3_lmt*den_aps(1)*qr(rcon_q,ka,1,k)
!ccc               write(fid_alog,*) den_aps(1),den_api(2)
               if(qr(rmt_q,ka,1,k)>0.0_RP) then
                  m_lmt=min_fapt_r*qr(rmt_q,ka,1,k)
                  qr(rmat_q,ka,1,k)=max( &
                           m_lmt,min(qr(rmat_q,ka,1,k) &
                          ,mx_aprat*qr(rmt_q,ka,1,k)))

                  m_lmt=min_faps_r*qr(rmat_q,ka,1,k)
                  qr(rmas_q,ka,1,k)=max( &
                       m_lmt,min(qr(rmas_q,ka,1,k),qr(rmat_q,ka,1,k)))
!                  if(qr(rmas_q,ka,1,k).le.m_lmt) then
!                     qr(rmas_q,ka,1,k)=eps_ap(1)*qr(rmat_q,ka,1,k)
!                  else
!                     qr(rmas_q,ka,1,k)=max(
!     *                   m_lmt,min(qr(rmas_q,ka,1,k),qr(rmat_q,ka,1,k)))
!                  end if
!                  if(qr(rmas_q,ka,1,k).lt.m_lmt) then
!                     qr(rmas_q,ka,1,k)=qr(rmat_q,ka,1,k)
!                  end if
               end if
               mapt_r(K)=mapt_r(K)+qr(rmat_q,ka,1,k)
            end do
            mapt_r(k)=mapt_r(k)*den1
            oapt_r(k)=oapt_r(k)*den1
         ENDDO
      ELSE
         DO K=k1,N1
            mapt_r(K)=0.
         ENDDO
      ENDIF
      IF(IPCNFL.GT.0) THEN
         DO K=k1,N1
            den1=den(k)*1.0e-3

            mapt_i(K)=0.0
            oapt_i(k)=0.0
            do ka=1,nbi
               oapt_i(k)=oapt_i(k)+qi(imat_q,ka,1,k)
!               m_lmt=coef3*r3_lmt*den_api(2)*qi(icon_q,ka,1,k)
               if(qi(imt_q,ka,1,k)>0.0_RP) then
                  m_lmt=min_fapt_s*qi(imt_q,ka,1,k)
                  qi(imat_q,ka,1,k)=max( &
                           m_lmt,min(qi(imat_q,ka,1,k) &
                          ,mx_aprat*qi(imt_q,ka,1,k)))

                  m_lmt=min_faps_s*qi(imat_q,ka,1,k)
                  qi(imas_q,ka,1,k)=max( &
                     m_lmt,min(qi(imas_q,ka,1,k),qi(imat_q,ka,1,k)))

!                  if(qi(imas_q,ka,1,k).lt.m_lmt) then
!                     qi(imas_q,ka,1,k)=eps_ap(2)*qi(imat_q,ka,1,k)
!                  else
!                     qi(imas_q,ka,1,k)=max(
!     *                     m_lmt,min(qi(imas_q,ka,1,k),qi(imat_q,ka,1,k)))
!                  end if
!                  if(qi(imas_q,ka,1,k).lt.m_lmt) then
!                     qi(imas_q,ka,1,k)=0.0
!                  end if
               end if
               mapt_i(K)=mapt_i(K)+qi(imat_q,ka,1,k)
            end do
            mapt_i(k)=mapt_i(k)*den1
            oapt_i(k)=oapt_i(k)*den1
         ENDDO
      ELSE
         DO K=k1,N1
            mapt_i(K)=0.
         ENDDO
      ENDIF
      k1=1
      do k=k1,N1
         den1=den(k)*1.0e-3

         mapt_a(k)=0.0
         oapt_a(k)=0.0
!         min_mapt=0.0

         do ica=1,nca

            oapt_a(k)=oapt_a(k)+qa(amt_q,1,ica,k)
            m_lmt=coef3*r3_lmt*den_apt(ica)

            if(qa(acon_q,1,ica,k)<tconlmt/den1) then
               qa(acon_q,1,ica,k)=tconlmt/den1
               qa(amt_q,1,ica,k)=qa(acon_q,1,ica,k)*m_lmt
               qa(ams_q,1,ica,k)=qa(amt_q,1,ica,k)*eps_ap(ica)
            end if

            if(qa(amt_q,1,ica,k)<m_lmt*qa(acon_q,1,ica,k)) then
!               write(fid_alog,*)
!     *               "negadj_bin_ap > mean mass is smaller than lmt."
!               write(fid_alog,*) "m,m_lmt",qa(amt_q,1,ica,k)/qa(acon_q,1,ica,k),m_lmt
!               write(fid_alog,*) "eps_ap,ica",eps_ap(ica),qa(amt_q,1,ica,k)
               qa(amt_q,1,ica,k)=qa(acon_q,1,ica,k)*m_lmt
               qa(ams_q,1,ica,k)=eps_ap(ica)*qa(amt_q,1,ica,k)
            end if

            if(ica==1) then
               if(qa(ams_q,1,ica,k)<sep_faps*qa(amt_q,1,ica,k)) then
                  qa(ams_q,1,ica,k)=sep_faps*qa(amt_q,1,ica,k)
               end if
            elseif(ica==2) then
               if(qa(ams_q,1,ica,k)>=sep_faps*qa(amt_q,1,ica,k)) then
                  qa(ams_q,1,ica,k)=0.99*sep_faps*qa(amt_q,1,ica,k)
               end if
            elseif(ica==3) then
               if(qa(ams_q,1,ica,k)<sep_faps*qa(amt_q,1,ica,k)) then
                  qa(ams_q,1,ica,k)=sep_faps*qa(amt_q,1,ica,k)
               end if
            elseif(ica==4) then
               if(qa(ams_q,1,ica,k)<sep_faps*qa(amt_q,1,ica,k)) then
                  qa(ams_q,1,ica,k)=sep_faps*qa(amt_q,1,ica,k)
               end if
            else
               LOG_ERROR("negadj_bin_ap",*) "negadj_ap>This cat is not defined"
               call PRC_abort
            end if
            qa(ams_q,1,ica,k)=min(qa(ams_q,1,ica,k),qa(amt_q,1,ica,k))

            mapt_a(k)=mapt_a(k)+qa(amt_q,1,ica,k)

!            min_mapt=min_mapt+m_lmt*qa(acon_q,1,ica,k)
         end do
         mapt_a(k)=mapt_a(k)*den1
         oapt_a(k)=oapt_a(k)*den1
!         min_mapt=min_mapt*den1

!         diff=max(0.0_RP,mapt_r(k)-oapt_r(k)+mapt_i(k)-oapt_i(k))
!         if(diff<mapt_a(k)-min_mapt) then
!            fac=(mapt_a(k)-diff)/mapt_a(k)
!            do ica=1,nca
!               qa(amt_q,1,ica,k)=qa(amt_q,1,ica,k)*fac
!               qa(acon_q,1,ica,k)=qa(acon_q,1,ica,k)*fac
!               qa(ams_q,1,ica,k)=qa(ams_q,1,ica,k)*fac
!            end do
!         else
!         if(mapt_r(k)+mapt_i(k)+mapt_a(k)>
!     *                      oapt_r(k)+oapt_i(k)+oapt_a(k)) then
!            write(fid_alog,*)
!     *   "negadj_bin_ap > k, mapt_r,mapt_i,mapt_a,oapt_r,oapt_i,oapt_a"
!            write(fid_alog,'(I5,6ES15.6)')
!     *          k, mapt_r(k),mapt_i(k),mapt_a(k)
!     *           ,oapt_r(k),oapt_i(k),oapt_a(k)
!            write(fid_alog,*) "k,total new,old map",k
!     *     ,mapt_r(k)+mapt_i(k)+mapt_a(k),oapt_r(k)+oapt_i(k)+oapt_a(k)
!         end if


      end do

      return
      end subroutine negadj_bin_ap

      subroutine fill_bin_ap(n3,npa,nba,nca &
           ,npr,nbr,ncr,npi,nbi,nci &
           ,qap,qrp,qip,qap01dv,denv,den01v &
           ,j1a,j2a,imicv,jmicv,kmicv)
!       This subroutine set the cloud-free aerosols back to the prescribed one.
      use maxdims, only : npamx,nbamx,ncamx
      use par_amps
      implicit none
      integer :: npa,nba,nca,npr,nbr,ncr,npi,nbi,nci,j1a,j2a
      integer :: j,n3,ipx,ibx,icx
      integer :: imicv(*),jmicv(*),kmicv(*)
      real(MP_KIND) :: qap(npa,nba,nca,*),qrp(npr,nbr,ncr,*),qip(npi,nbi,nci,*) &
          ,qap01dv(npamx,nbamx,ncamx,*),den01v(*)  ! tempei added den01v for Kid
      real(MP_KIND) :: denv(*)
      real(PS) :: den1,cr,ci,cd
! added for kwajex
      real(PS) :: fct_ap
      parameter(fct_ap=0.01)
! end added for kwajex

      do j=max(1,j1a-5),min(n3,j2a+5)
         ! density is in #/cm^3
         den1=denv(j)*1.0e-3
         cr=0.0
         cd=0.0
         do icx=1,ncr
            do ibx=1,nbr
               cr=cr+qrp(rcon_q,ibx,icx,j)
               cd=cd+qrp(rmt_q,ibx,icx,j)
            enddo
         enddo
         if(cr>1.0e-20) then
            cd=(6.0/3.14*cd/cr)**(1.0/3.0)
         else
            cd=0.0
         endif
         cr=cr*den1
         if(cr<1.0.or.cd<1.0e-4) then
!          condition for cloud free is N< 1 cm^-3 or D< 1 mic
            do icx=1,nca
               if(icx.eq.2) cycle
               do ibx=1,nba
                  do ipx=1,npa
                     qap(ipx,ibx,icx,j)=qap01dv(ipx,ibx,icx,j) &
                         *den01v(j)/denv(j)
                  end do
               end do
            end do
         endif
         ci=0.0
         do icx=1,nci
            do ibx=1,nbi
               ci=ci+qip(icon_q,ibx,icx,j)
            enddo
         enddo
         ci=ci*den1
         if(ci<1.0e-5) then
!          condition for cloud free is 0.01 liter^-1
            do icx=2,2
               do ibx=1,nba
                  do ipx=1,npa
                     qap(ipx,ibx,icx,j)=qap01dv(ipx,ibx,icx,j) &
                         *den01v(j)/denv(j)
                  end do
               end do
            end do
         endif
! added for kwajex
!         do icx=1,nca
!            do ibx=1,nba
!               if(qap(acon_q,ibx,icx,j)<fct_ap*
!     *                qap01dv(acon_q,ibx,icx,1)) then
!                  qap(acon_q,ibx,icx,j)=fct_ap*qap01dv(acon_q,ibx,icx,j) &
!                       *den01v(j)/denv(j)
!                  qap(amt_q,ibx,icx,j)=fct_ap*qap01dv(amt_q,ibx,icx,j) &
!                       *den01v(j)/denv(j)
!                  qap(ams_q,ibx,icx,j)=fct_ap*qap01dv(ams_q,ibx,icx,j) &
!                       *den01v(j)/denv(j)
!               endif
!            enddo
!         enddo
! end added for kwajex
      end do

      return
      end subroutine fill_bin_ap

      subroutine fill_bulk_ap(n3,npa,nba,nca &
           ,npr,nbr,ncr,npi,nbi,nci &
           ,qap,qrp,qip,qap01dv,denv,den01v &
           ,j1a,j2a,imicv,jmicv,kmicv)
!       This subroutine set the cloud-free aerosols back to the prescribed one.
      use maxdims, only : npamx,nbamx,ncamx
      use par_amps
      implicit none
      integer :: npa,nba,nca,npr,nbr,ncr,npi,nbi,nci,j1a,j2a
      integer :: j,n3,ipx,ibx,icx
      integer :: imicv(*),jmicv(*),kmicv(*)
      real(PS) :: qap(npa,nba,nca,*),qrp(npr,nbr,ncr,*),qip(npi,nbi,nci,*) &
          ,qap01dv(npamx,nbamx,ncamx,*),den01v(*)  ! tempei added den01v for Kid
      real(PS) :: denv(*)
      real(PS) :: den1,cr,ci,cd
      ! density of water [kg/m3]
      real(PS),parameter :: denw=1.0e+3

      do j=max(1,j1a-5),min(n3,j2a+5)
         ! density is in kg/m^3
         den1=denv(j)
         cr=0.0
         cd=0.0
         do icx=1,ncr
            do ibx=1,nbr
               cr=cr+qrp(2,ibx,icx,j)  ! [1/kg]
               cd=cd+qrp(1,ibx,icx,j)  ! [kg/kg]
            enddo
         enddo
         if(cr>1.0e-20) then
            cd=(6.0/3.14*cd/cr/denw)**(1.0/3.0)  ! [m]
         else
            cd=0.0
         endif
         cr=cr*den1
         ! cr: number concentration  [#/m^3]
         ! cd: mean radius      [m]
         if(cr<1.0e+6.or.cd<1.0e-6) then
!          condition for cloud free is N< 1 cm^-3 or D< 1 mic
            do icx=1,1
!              category 1 for CLR is dry aerosol (CCN), so only this is restored.
               do ibx=1,nba
                  do ipx=1,npa
                     qap(ipx,ibx,icx,j)=qap01dv(ipx,ibx,icx,j) &
                         *den01v(j)/denv(j)
                  end do
               end do
            end do
         endif
         ci=0.0
         do icx=1,nci
            do ibx=1,nbi
               ci=ci+qip(2,ibx,icx,j)  ! [1/kg]
            enddo
         enddo
         ci=ci*den1
         if(ci<1.0e+1) then
!          condition for cloud free is 0.01 liter^-1
            do icx=nca,nca  ! IN is assumed to be nca.
               do ibx=1,nba
                  do ipx=1,npa
                     qap(ipx,ibx,icx,j)=qap01dv(ipx,ibx,icx,j) &
                         *den01v(j)/denv(j)
                  end do
               end do
            end do
         endif
      end do

      return
      end subroutine fill_bulk_ap

      subroutine sprayms(dfmsdr0,r0,u14)
!
!     inputs:
!
!        r0     = radius of droplet (micrometer)
!        u14    = wind speed at 14 m agl (m/s)
!
!     output
!
!        dFmsdr0  = #/m**2/s/micrometer
!

      implicit none
      real(PS) :: dFsdr80,dFmsdr0,r0,r80,flx1 &
          ,dr80dr0
      real(RP) :: u14
      r80=0.518*r0**0.976
      call sea93(dFsdr80,min(r80,10._RP),u14)
      flx1=dfsdr80
      if(r80.gt.10.0_RP) then
        flx1=flx1/min(r80,37.5_RP)
      elseif(r80.gt.37.5) then
        flx1=flx1*min(r80,100._RP)**(-2.8)
      elseif(r80.gt.100) then
        flx1=flx1/r80**8
      endif
      dr80dr0=0.506*r0**(-0.024)
      dFmsdr0=3.5*flx1*dr80dr0
      return
      end subroutine sprayms

      subroutine sea93(dFsdr80,r80,u14)
!
!     inputs:
!
!        r80    = radius of droplet (micrometer)
!        u14    = wind speed at 14 m agl (m/s)
!
!     output
!
!        dFsdr80  = #/m**2/s/micrometer
!

      implicit none
      real(PS) :: dFsdr80,r80,r1,r2,f1,f2,A1,A2
      real(RP) :: u14
      parameter (r1=2.1,r2=9.2,f1=3.1,f2=3.3)
      a1=10.**(0.0676*u14+2.43)
      a2=10.**(0.959*sqrt(u14)-1.476)
      dFsdr80=a1*exp(-f1*log(r80/r1)**2)+a2*exp(-f2*log(r80/r2)**2)
      return
      end subroutine sea93

!!$      real(PS) FUNCTION TWT(RVP,THET,P,ifrom)
!!$      real(PS),intent(in) :: RVP,THET,P
!!$      integer,intent(in) :: ifrom
!!$      real(PS) :: press,rvap,piter,temper,x,aos
!!$!tmp,tsa
!!$      integer :: id
!!$!     ABS IS ABSOLUTE VALUE
!!$!     ALL ARGUMENTS AND TW (KELVIN)
!!$      PRESS=P*1.E-2
!!$      RVAP=RVP*1.E3
!!$      PITER =  PRESS
!!$!      write(fid_alog,*) "twt > pressure",PRESS
!!$      DO ID=1,10
!!$         TEMPER=THET*(PITER*1.E-3)**.286
!!$         X  =  .02*( TMR(RVAP,PITER) - TEMPER)
!!$         IF ( ABS(X).LT. 0.01  ) exit
!!$         PITER = PITER* ( 2.**(X)  )
!!$      end DO
!!$      TEMPER=THET*(PITER*1.E-3)**.286
!!$!
!!$      AOS  =   OS(TEMPER,PITER)
!!$      TWT   =  TSA( AOS,PRESS)
!!$      RETURN
!!$      END function twt

!!$      real(PS) FUNCTION TMR(W,P)
!!$      use com_amps, only: fid_alog
!!$      real(PS),intent(in) :: W,P
!!$      real(PS) :: X
!!$!     TMR(KELVIN),W(GRAMS WATER VAPOR/KILOGRAM DRY AIR),P(MILLIBAR)
!!$!     ALOG10  15   LOG TO THE BASE  TEN.
!!$      if(W*P/(622._RP + 0.378_RP*W)<=0.0_RP)then
!!$          write(fid_alog,*) "inside of log10 is negativel:",W*P/(622.0_RP + 0.378_RP*W)
!!$          write(fid_alog,'(2ES15.6)') W,P
!!$        stop
!!$        end if
!!$      X =  LOG10(   W*P/(622.0_PS+ 0.378_PS*W)  ) ! precision test, original DLOG10
!!$      TMR=10.0_RP**(0.0498646455_RP*X+2.4082965_RP)-7.07475_RP+38.9114_RP*((10.0_RP**( &
!!$        0.0915_RP*X ) - 1.2035_RP )**2 )
!!$      RETURN
!!$      END function tmr

!!$      real(PS) FUNCTION OS(T,P)
!!$      real(PS),intent(in) :: T,P
!!$!     OS AND T (KELVIN) , P (MILLIBARS )
!!$      OS=T*((1000.0_RP/P)**0.286_RP)/(EXP(-2.6518986_RP*W(T,P)/T))
!!$      RETURN
!!$      END function os

!!$      real(PS) FUNCTION TSA(OS,P)
!!$      real(PS),intent(in) :: OS,P
!!$      real(PS) :: A, TQ, D, X
!!$      integer :: ID
!!$!     TSA AND OS(KELVIN),P(MILLIBARS)
!!$!     SIGN(A,B) RREPLACES THE ALGEBRETIC SIGN OF A WITH THE SIGN OF B
!!$      A  =  OS
!!$      TQ  =  253.16
!!$      D =  120
!!$      DO 1  ID= 1,12
!!$        D = D/2.
!!$!     IF THE TEMPERATURE DIFFERENCE,X, IS SMALL,EXIT THIS LOOP
!!$        X=A*EXP(-2.6518986*W(TQ,P)/TQ)-TQ*((1000./P)**.286)
!!$        IF(ABS(X).LT.0.01)GO TO 2
!!$          TQ = TQ + SIGN(D,X)
!!$  1       CONTINUE
!!$2     TSA=TQ
!!$      RETURN
!!$      END function tsa

!!$      real(PS) FUNCTION W(T,P)
!!$      real(PS),intent(in) :: T,P
!!$      real(PS) :: X
!!$!     W(GRAMS WATER VAPOR/KILOGRAM DRY AIR ), P(MILLIBAR )
!!$      IF(T.GE.999.0_RP)GO TO 10
!!$        X  =       ESAT(T)
!!$        W  =  622.0_RP*X/(P-0.378*X)
!!$      RETURN
!!$   10 W=0.0
!!$      RETURN
!!$      END function w

!!$      real(PS) FUNCTION ESAT(T)
!!$!     ESAT(MILLIBARS),T(KELVIN)
!!$      real(PS),intent(in) :: T
!!$      real(PS),parameter :: ABZ=273.16
!!$      real(PS) :: TC
!!$      TC=T-ABZ
!!$      ESAT=6.1078*EXP((17.2693882*TC)/(TC+237.3))
!!$      RETURN
!!$      END function esat

!!$      subroutine diag_icecat(isntyp,level,ici,ictpc,ictsw,ictag,ictgr &
!!$                            ,qipv)
!!$      use par_amps
!!$      implicit none
!!$      character *(*) :: isntyp
!!$      integer :: level,ici,ictpc,ictsw,ictag,ictgr &
!!$             ,itype
!!$!tmp,get_sh_type4
!!$      real(PS) :: am1,am2,am3,am4,am5,qipv(*)
!!$
!!$      if(level.eq.3) then
!!$         if(ici.eq.ictpc) then
!!$            isntyp='PRIS'
!!$         elseif(ici.eq.ictsw) then
!!$            isntyp='SNOW'
!!$         elseif(ici.eq.ictag) then
!!$            isntyp='AGG'
!!$         elseif(ici.eq.ictgr) then
!!$            isntyp='GRAUP'
!!$         endif
!!$      elseif(level.eq.4) then
!!$         if(ici.eq.ictpc) then
!!$            isntyp='PRIS'
!!$         elseif(ici.eq.ictsw) then
!!$            isntyp='SNOW'
!!$         elseif(ici.eq.ictgr) then
!!$            isntyp='GRAUP'
!!$         endif
!!$
!!$      elseif(level.eq.5) then
!!$         am1=qipv(imt_q)
!!$         am2=qipv(imr_q)
!!$         am4=qipv(imc_q)
!!$         am5=max(qipv(imw_q),0.0_RP)
!!$         am3=max(am1-am2-am4-am5,0.0_RP)
!!$         itype=get_sh_type4(am1,am2,am3,am4,am5)
!!$         if(itype.eq.1) then
!!$            isntyp='PRIS'
!!$         elseif(itype.eq.2.or.itype.eq.4) then
!!$            isntyp='SNOW'
!!$         elseif(itype.eq.3) then
!!$            isntyp='AGG'
!!$         elseif(itype.eq.5) then
!!$            isntyp='GRAUP'
!!$         endif
!!$
!!$      end if
!!$!ccc      write(fid_alog,*) "ici,ictpc,ictsw,ictag,ictgr"
!!$!ccc     *               ,ici,ictpc,ictsw,ictag,ictgr,isntyp
!!$
!!$      return
!!$      end subroutine diag_icecat

!!$      integer function get_sh_type4(am1,am2,am3,am4,am5)
!!$      implicit none
!!$!
!!$!     am1: total mass
!!$!     am2: riming mass component
!!$!     am3: aggregation mass component
!!$!     am4: ice crystal mass component
!!$!     am5: melt-water mass component
!!$!
!!$      real(PS) :: am1,am2,am3,am4,am5,xsum,x
!!$      integer :: itype
!!$!tmp,get_sh_type4
!!$      real(PS) :: frac
!!$      parameter(frac=0.1)
!!$      if(am2.gt.(am1-am5)*frac) then
!!$         if(am3.le.(am1-am5)*frac ) then
!!$            if(am2.le.am4) then
!!$!              rimed crystal
!!$               itype=2
!!$            else
!!$!              graupel
!!$               itype=5
!!$            end if
!!$         else
!!$            if(am2.le.am3) then
!!$!              rimed aggregate
!!$               itype=4
!!$            else
!!$!              graupel
!!$               itype=5
!!$            end if
!!$         end if
!!$      else if(am3.gt.(am1-am5)*frac) then
!!$!        aggregate
!!$         itype=3
!!$      else
!!$!        pristine crystal
!!$         itype=1
!!$      end if
!!$      get_sh_type4=itype
!!$      end function get_sh_type4

!!$      subroutine azero(n,a)
!!$      implicit none
!!$      integer :: n
!!$      real(PS) :: a(n)
!!$      a(1:n)=0.0
!!$      end subroutine azero

!!$subroutine cal_transbin(iphase,clctr,igrid&
!!$     ,nmass,naxis&
!!$     ,nbr,ncr&
!!$     ,binbr,mbinb,srat,sadd &
!!$     ,par,s_bd &
!!$     ,new_N,new_M,new_Q&
!!$     ,ratio_Mp,den_ip_p,axr_p,spx_p &
!!$     ,habit_p,den_ic_p &
!!$     ,rag_p,rcg_p,n_exice_p &
!!$     ,ap_dN,ap_dM)
!!$  use maxdims
!!$  use par_amps
!!$  implicit none
!!$  ! **********************************************************************
!!$  !  calculate concentration and mass transferred to the original bins
!!$  !  and add it to the total change of those in each bin
!!$  !  This is the specialized version for ice phase.
!!$  ! **********************************************************************
!!$  integer :: iphase,clctr,igrid,nmass,naxis,nbr,ncr
!!$  ! bin boundaries
!!$  real(PS) :: binbr(*)
!!$  !     minimum bin mass, ratio and addition to define original bin boundaries
!!$  real(PS) :: mbinb,srat,sadd
!!$  ! parameter of distribution in a bin
!!$  real(8) :: par(*)
!!$  ! shifted bounds of a bin
!!$  real(8) :: s_bd(*)
!!$  ! new total concentration mapped to original bins.
!!$  real(PS) :: new_N(*)
!!$  ! new total mass mapped to original bins.
!!$  real(PS) :: new_M(mxnmasscomp+1,*)
!!$  ! new quality prognostic variables mapped to original bins
!!$  real(PS) :: new_Q(mxnnonmc,*)
!!$  ! ratio of each mass component to total mass in the shifted bin.
!!$  real(PS) :: ratio_Mp(*)
!!$  ! bulk sphere density of dry ice particle and bulk crystal density in the shifted bin.
!!$  real(RP) :: den_ip_p,den_ic_p
!!$  ! axis ratio for ice crystals and volume of cylinder
!!$  ! 1: c/a
!!$  ! 2: d/a
!!$  ! 3: r/a
!!$  ! 4: e/a
!!$  real(RP) :: axr_p(*)
!!$  ! mass tendency of the shifted mass
!!$  ! habit in shifted bin
!!$  integer :: habit_p
!!$  ! aspect ratio of circumscribing cylinder
!!$  real(RP) :: spx_p
!!$
!!$  ! ratio of ag^3 to a^3
!!$  real(PS) :: rag_p,rcg_p
!!$  ! number of extra ice crystals
!!$  real(RP) :: n_exice_p
!!$
!!$  ! concentration and mass to be moved into aerosol group
!!$  real(PS), optional  :: ap_dN, ap_dM
!!$
!!$  ! transferred total concentration and mass from bin to bin
!!$  real(RP) :: trans_dN, trans_dM
!!$  !
!!$  ! transferred volume of circumscribing spheroid
!!$  real(RP) :: trans_dvcs
!!$
!!$  ! transferred length
!!$  real(RP) :: trans_dL(mxnaxis)
!!$  ! mid mass of two boundaries
!!$  real(PS) :: mid_mass
!!$
!!$!!c     dummy parameters
!!$!!c      real a(3)
!!$!!c   real, shifted boundaries
!!$  real(8) :: left_bd, right_bd
!!$
!!$  ! dummy concentration and mass
!!$  real(PS) :: d_N, d_M
!!$
!!$  ! dummy numbers
!!$  real(PS) :: n_left, n_right, n_top, x_max
!!$  real(PS) :: n_top2, x_max2
!!$  ! error message
!!$  integer :: em
!!$  integer :: ibx, k
!!$  integer :: item(2)
!!$
!!$  ! type of subdistribution function
!!$  integer :: dstype
!!$
!!$  ! hexagonal length
!!$  integer, parameter :: npreaxis=3
!!$
!!$
!!$  ! initialize the error message
!!$  em=0
!!$  item=0
!!$  dstype=1
!!$
!!$  ! +++ check if any concentration and mass exist +++
!!$  if( par(2) <= 0.0_RP) then
!!$     return
!!$  end if
!!$
!!$  !     +++ set up the left and right boundaries of new bins +++
!!$  if( par(1) > 0.0 ) then
!!$     !     --- case of non-negative bin ---
!!$     left_bd = s_bd(1)
!!$     right_bd = s_bd(2)
!!$  else if( par(1) == -1.0 ) then
!!$     !     case of n(x'_1) < 0.0
!!$     par(1) = 0.0
!!$     left_bd = par(2)
!!$     right_bd = s_bd(2)
!!$  else if( par(1) == -2.0 ) then
!!$     !     case of n(x'_2) < 0.0
!!$     par(1) = 0.0
!!$     left_bd = s_bd(1)
!!$     right_bd = par(2)
!!$  end if
!!$
!!$  if( iphase == 1 ) then
!!$!!c     write(fid_alog,'("iphase,ibin,ngrid,nmass,naxis,nbr,ncr",8I5)') iphase,clctr,igrid,nmass,naxis,nbr,ncr
!!$!!c     write(fid_alog,'("binbr",30ES15.6)') binbr
!!$!!c     write(fid_alog,'("ratio_Mp",15ES15.6)') ratio_Mp
!!$!!c     write(fid_alog,'("axr_p",20ES15.6)') axr_p
!!$!!c     write(fid_alog,'("den,denip,spx_p,rag,rcg,nex",20ES15.6)') den_ip_p,den_ic_p,spx_p,rag_p,rcg_p,n_exice_p
!!$!!c     write(fid_alog,'("habit",5I5)') habit_p
!!$
!!$     !     --- case where the colletor is a liquid hydromteor ---
!!$     if( right_bd > binbr(nbr+1) ) then
!!$        call cal_transmom( dstype, par&
!!$             ,max(real(binbr(nbr+1),8),left_bd)&
!!$             ,right_bd,trans_dN,trans_dM)
!!$        call ck_errtrans(clctr,igrid,-9,trans_dN,trans_dM,em)
!!$        if(em==0) then
!!$           new_M(rmt,nbr)=new_M(rmt,nbr)+trans_dM
!!$           do k = 1,nmass
!!$              new_M(k+1,nbr)=&
!!$                 new_M(k+1,nbr)+ratio_Mp(k)*trans_dM
!!$           end do
!!$           !     +++ determine new concentration by mean mass +++
!!$           mid_mass=(binbr(nbr+1)+binbr(nbr))/2.0
!!$           trans_dN = trans_dM/mid_mass
!!$           new_N(nbr) = new_N(nbr) + trans_dN
!!$        end if
!!$     end if
!!$     if( left_bd < binbr(1) ) then
!!$        call cal_transmom( dstype, par&
!!$             ,left_bd &
!!$             ,min(real(binbr(1),8), right_bd)&
!!$             ,trans_dN, trans_dM)
!!$        call ck_errtrans(clctr,igrid,-9,trans_dN,trans_dM,em)
!!$     end if
!!$
!!$     !     +++ consider the collector bin and larger bins than the collector bin +++
!!$     if(clctr<=nbr) then
!!$        if( right_bd > binbr(clctr) ) then
!!$           do ibx=clctr, nbr
!!$              if( left_bd > binbr(ibx+1) ) cycle
!!$              em=0
!!$              call cal_transmom( dstype, par&
!!$                   ,max(real(binbr(ibx),8), left_bd)&
!!$                   ,min(real(binbr(ibx+1),8), right_bd) &
!!$                   ,trans_dN, trans_dM)
!!$              call ck_errtrans(clctr,igrid,ibx&
!!$                   ,trans_dN,trans_dM, em)
!!$              if( em == 0 ) then
!!$                 new_N(ibx) = new_N(ibx) + trans_dN
!!$                 new_M(rmt,ibx) = new_M(rmt,ibx) + trans_dM
!!$                 do k = 1,nmass
!!$                    new_M(k+1,ibx)=&
!!$                    new_M(k+1,ibx)+ratio_Mp(k)*trans_dM
!!$                 end do
!!$              end if
!!$              if( right_bd <= binbr(ibx+1) ) exit
!!$           end do
!!$           !     --- if the left shifted boundary is greater than the left
!!$           !     original bounary, the transfer has been done.
!!$           if( left_bd >= binbr(clctr) ) then
!!$              return
!!$           end if
!!$        end if
!!$     end if
!!$
!!$     !     +++ consider the smaller bins than the collector bin +++
!!$     do ibx=clctr-1, 1, -1
!!$        if(ibx>nbr) cycle
!!$        if( right_bd < binbr(ibx) ) cycle
!!$        em=0
!!$        call cal_transmom( dstype, par&
!!$             ,max(real(binbr(ibx),8), left_bd)&
!!$             ,min(real(binbr(ibx+1),8), right_bd)&
!!$             ,trans_dN, trans_dM)
!!$        call ck_errtrans(clctr,igrid,ibx,trans_dN,trans_dM, em)
!!$        if( em == 0 ) then
!!$           new_N(ibx) = new_N(ibx) + trans_dN
!!$           new_M(rmt,ibx) = new_M(rmt,ibx) + trans_dM
!!$           do k = 1,nmass
!!$              new_M(k+1,ibx)=&
!!$              new_M(k+1,ibx)+ratio_Mp(k)*trans_dM
!!$           end do
!!$        end if
!!$        if( left_bd >= binbr(ibx) ) then
!!$           return
!!$        end if
!!$     end do
!!$  else if(iphase==2) then
!!$     !     --- case where the colletor is a solid hydromteor ---
!!$     !     +++ case where the large shifted bin limit reaches the maximum bin limit
!!$!!c     write(fid_alog,'("iphase,ibin,ngrid,nmass,naxis,nbr,ncr",8I5)') iphase,clctr,igrid,nmass,naxis,nbr,ncr
!!$!!c     write(fid_alog,'("binbr",30ES15.6)') binbr
!!$!!c     write(fid_alog,'("ratio_Mp",15ES15.6)') ratio_Mp
!!$!!c     write(fid_alog,'("axr_p",20ES15.6)') axr_p
!!$!!c     write(fid_alog,'("den,denip,spx_p,rag,rcg,nex",20ES15.6)') den_ip_p,den_ic_p,spx_p,rag_p,rcg_p,n_exice_p
!!$!!c     write(fid_alog,'("habit",5I5)') habit_p
!!$
!!$
!!$!!c     if(habit_p>=4) then
!!$!!c        write(fid_alog,*) "habit is not right"
!!$!!c        stop
!!$!!c     end if
!!$
!!$     if( left_bd < binbr(1) ) then
!!$        call cal_transmom( dstype, par&
!!$             ,left_bd&
!!$             ,min( real(binbr(1),8), right_bd)&
!!$             ,trans_dN, trans_dM)
!!$        call ck_errtrans(clctr,igrid,-9,trans_dN,trans_dM, em)
!!$     end if
!!$     if( right_bd > binbr(nbr+1) ) then
!!$!!c        write(fid_alog,*) "Warning: The shifted limit is above the upper limit. &
!!$!!c             The extra mass will be put into the maximum bin", clctr,igrid
!!$        call cal_transmom(dstype, par, &
!!$             max(real(binbr(nbr+1),8), left_bd), &
!!$             right_bd, &
!!$             trans_dN, trans_dM)
!!$        call ck_errtrans(clctr,igrid, -9, trans_dN, trans_dM, em)
!!$        if( em == 0 ) then
!!$           new_M(imt,nbr) = new_M(imt,nbr) + trans_dM
!!$           do k = 1, nmass
!!$              new_M(k+1,nbr) = new_M(k+1,nbr) + ratio_Mp(k)*trans_dM
!!$           end do
!!$           ! +++ determine new concentration by mean mass +++
!!$           mid_mass = ( binbr(nbr+1) + binbr(nbr) )/2.0
!!$           trans_dN = trans_dM/mid_mass
!!$           new_N(nbr) = new_N(nbr) + trans_dN
!!$
!!$           call cal_transvolcs( dstype, par, den_ip_p, &
!!$                max(real(binbr(nbr+1),8), left_bd), &
!!$                right_bd, &
!!$                trans_dM*(1.0-ratio_Mp(imw_m)),trans_dvcs)
!!$           new_Q(ivcs,nbr) = new_Q(ivcs,nbr) + trans_dvcs
!!$
!!$           call cal_translen3( dstype, par, den_ic_p,axr_p,spx_p,habit_p, &
!!$                max(real(binbr(nbr+1),8), left_bd), &
!!$                right_bd, &
!!$                trans_dM*ratio_Mp(imc_m),n_exice_p,trans_dL)
!!$!!c           new_Q(2:1+naxis,nbr) = new_Q(2:1+naxis,nbr)+trans_dL(1:naxis)
!!$           new_Q(iacr:iacr+npreaxis-1,nbr) = new_Q(iacr:iacr+npreaxis-1,nbr)+trans_dL(1:npreaxis)
!!$           new_Q(iag,nbr)=new_Q(iag,nbr)+trans_dL(1)*rag_p
!!$           new_Q(icg,nbr)=new_Q(icg,nbr)+trans_dL(2)*rcg_p
!!$           new_Q(inex,nbr)=new_Q(inex,nbr)+n_exice_p*trans_dN
!!$        end if
!!$     end if
!!$
!!$     !     +++ consider the collector bin and larger bins than the collector bin +++
!!$     if(clctr<=nbr) then
!!$        if( right_bd > binbr(clctr) ) then
!!$           do ibx=clctr, nbr
!!$              if( left_bd > binbr(ibx+1) ) cycle
!!$              em=0
!!$              call cal_transmom( dstype, par&
!!$                   ,max(real(binbr(ibx),8), left_bd)&
!!$                   ,min(real(binbr(ibx+1),8), right_bd)&
!!$                   ,trans_dN, trans_dM)
!!$              call ck_errtrans(clctr,igrid,ibx&
!!$                   ,trans_dN,trans_dM,em)
!!$              if( em == 0 ) then
!!$
!!$                 new_N(ibx) = new_N(ibx) + trans_dN
!!$                 new_M(imt,ibx) = new_M(imt,ibx) + trans_dM
!!$
!!$                 if(ratio_Mp(imc_m)==1.0) then
!!$                    new_M(imc,ibx) = new_M(imc,ibx)+trans_dM
!!$
!!$                    do k=1,nmass
!!$                       if(k==imat_m.or.k==imas_m.or.k==imf_m) then
!!$                          new_M(k+1,ibx)=new_M(k+1,ibx)+ratio_Mp(k)*trans_dM
!!$                       endif
!!$                    end do
!!$
!!$
!!$                 else
!!$                    do k = 1,nmass
!!$                       new_M(k+1,ibx)=&
!!$                            new_M(k+1,ibx)+ratio_Mp(k)*trans_dM
!!$                    end do
!!$                 endif
!!$
!!$                 call cal_transvolcs( dstype, par, den_ip_p&
!!$                          ,max(real(binbr(ibx),8), left_bd)&
!!$                          ,min(real(binbr(ibx+1),8), right_bd)&
!!$                          ,trans_dM*(1.0-ratio_Mp(imw_m)),trans_dvcs)
!!$                 new_Q(ivcs,ibx) = new_Q(ivcs,ibx) + trans_dvcs
!!$
!!$                 call cal_translen3( dstype, par,den_ic_p&
!!$                      ,axr_p,spx_p,habit_p&
!!$                      ,max(real(binbr(ibx),8), left_bd)&
!!$                      ,min(real(binbr(ibx+1),8), right_bd)&
!!$                      ,trans_dM*ratio_Mp(imc_m),n_exice_p,trans_dL)
!!$!!c                 new_Q(2:1+naxis,ibx)=new_Q(2:1+naxis,ibx)&
!!$!!c                      +trans_dL(1:naxis)
!!$
!!$                 new_Q(iacr:iacr+npreaxis-1,ibx)=new_Q(iacr:iacr+npreaxis-1,ibx)+trans_dL(1:npreaxis)
!!$                 new_Q(iag,ibx)=new_Q(iag,ibx)+trans_dL(1)*rag_p
!!$                 new_Q(icg,ibx)=new_Q(icg,ibx)+trans_dL(2)*rcg_p
!!$                 new_Q(inex,ibx)=new_Q(inex,ibx)+n_exice_p*trans_dN
!!$
!!$              end if
!!$              if( right_bd <= binbr(ibx+1) ) exit
!!$           end do
!!$           !     --- if the left shifted boundary is greater than the left
!!$           !     original bounary, the transfer has been done.
!!$           if( left_bd >= binbr(clctr) ) then
!!$              return
!!$           end if
!!$        end if
!!$     end if
!!$
!!$     !     +++ consider the smaller bins than the collector bin +++
!!$     do ibx=clctr-1, 1, -1
!!$        if(ibx>nbr) cycle
!!$        if( right_bd < binbr(ibx) ) cycle
!!$        em=0
!!$        call cal_transmom( dstype, par&
!!$             ,max(real(binbr(ibx),8), left_bd)&
!!$             ,min(real(binbr(ibx+1),8), right_bd)&
!!$             ,trans_dN, trans_dM)
!!$
!!$        call ck_errtrans(clctr,igrid,ibx,trans_dN,trans_dM,em)
!!$        if( em == 0 ) then
!!$           new_N(ibx) = new_N(ibx) + trans_dN
!!$
!!$           new_M(imt,ibx) = new_M(imt,ibx) + trans_dM
!!$           if(ratio_Mp(imc_m)==1.0) then
!!$              new_M(imc,ibx)=new_M(imc,ibx)+trans_dM
!!$              do k=1,nmass
!!$                 if(k==imat_m.or.k==imas_m.or.k==imf_m) then
!!$                    new_M(k+1,ibx)=new_M(k+1,ibx)+ratio_Mp(k)*trans_dM
!!$                 endif
!!$              end do
!!$
!!$           else
!!$              do k=1,nmass
!!$                 new_M(k+1,ibx)=new_M(k+1,ibx)+ratio_Mp(k)*trans_dM
!!$              end do
!!$           end if
!!$
!!$           call cal_transvolcs( dstype, par, den_ip_p&
!!$                ,max(real(binbr(ibx),8), left_bd)&
!!$                ,min(real(binbr(ibx+1),8), right_bd) &
!!$                ,trans_dM*(1.0-ratio_Mp(imw_m)),trans_dvcs)
!!$           new_Q(ivcs,ibx) = new_Q(ivcs,ibx) + trans_dvcs
!!$
!!$           call cal_translen3( dstype, par, den_ic_p,axr_p,spx_p,habit_p&
!!$                ,max(real(binbr(ibx),8), left_bd)&
!!$                ,min(real(binbr(ibx+1),8), right_bd)&
!!$                ,trans_dM*ratio_Mp(imc_m),n_exice_p,trans_dL)
!!$!!c           new_Q(2:1+naxis,ibx)=new_Q(2:1+naxis,ibx)+trans_dL(1:naxis)
!!$           new_Q(iacr:iacr+npreaxis-1,ibx) = new_Q(iacr:iacr+npreaxis-1,ibx)+trans_dL(1:npreaxis)
!!$           new_Q(iag,ibx)=new_Q(iag,ibx)+trans_dL(1)*rag_p
!!$           new_Q(icg,ibx)=new_Q(icg,ibx)+trans_dL(2)*rcg_p
!!$           new_Q(inex,ibx)=new_Q(inex,ibx)+n_exice_p*trans_dN
!!$        end if
!!$        if( left_bd >= binbr(ibx) ) then
!!$           return
!!$        end if
!!$     end do
!!$  end if
!!$end subroutine cal_transbin

subroutine cal_transbin_vec(iphase &
     ,L,nmass &
     ,nbin,nbin_s &
     ,binb &
     ,error_number &
     ,a2d,s_bd,mtend &
     ,new_N,new_M,new_Q &
     ,new_mtend &
     ,ratio_Mp,den_ip_p,axr_p,spx_p &
     ,habit_p,den_ic_p &
     ,rag_p,rcg_p,n_exice_p &
     ,actINF_p &
     ,iaer_src,ap_dN,ap_dM,ap_dMS,ap_dNI,ap_dMV,inevp)
!!!  use acc_amps
  use maxdims
  use par_amps
  use com_amps, only: fid_alog
  implicit none
  ! **********************************************************************
  !  calculate concentration and mass transferred to the original bins
  !  and add it to the total change of those in each bin
  !  This is the specialized version for ice phase.
  !
  !  Note:
  !    the distribution parameters a2d will be modified to be
  !    in the form of
  !     n(m)=a0+a1*m+a2*m*m+a3*m*m*m
  ! **********************************************************************
  integer,intent(in) :: iphase,nmass,L
  ! numbers of original bins
  integer,intent(in) :: nbin
  ! number of bins for shifted ones
  integer,intent(in) :: nbin_s
  ! original bin boundaries
  real(PS),dimension(*),intent(in) :: binb ! precision test (PS)
  ! parameter of distribution in a bin
  real(8),dimension(mxnbin+1,L,4), intent(inout) :: a2d
  ! shifted bounds of a bin
  real(RP), dimension(mxnbin+1,L,2), intent(in) ::  s_bd
  ! mass tendency of shifted bin
  real(PS), dimension(mxnbin+1,*), intent(in)     :: mtend

  ! new total concentration mapped to original bins.
  real(8), dimension(mxnbin,L), intent(inout) :: new_N
  ! new total mass mapped to original bins.
  real(8), dimension(mxnbin,L,1+mxnmasscomp), intent(inout) :: new_M
  ! new quality prognostic variables mapped to original bins
  real(8), dimension(mxnbin,L,mxnnonmc), intent(inout) :: new_Q
  ! averaged mass tendency after the time step for the bin
  real(8), dimension(mxnbin,*), intent(inout)         :: new_mtend

  ! ratio of each mass component to total mass in the shifted bin.
  real(RP), dimension(mxnbin+1,L,mxnmasscomp), intent(in) :: ratio_Mp
  ! bulk sphere density of dry ice particle and bulk crystal density in the shifted bin.
  real(PS), dimension(mxnbin+1,*), intent(in) :: den_ip_p,den_ic_p
  ! axis ratio for ice crystals and volume of cylinder
  ! 1: c/a
  ! 2: d/a
  ! 3: r/a
  ! 4: e/a
  real(PS), dimension(mxnbin+1,L,mxnaxis-1), intent(in) :: axr_p
  ! mass tendency of the shifted mass
  ! habit in shifted bin
  integer, dimension(mxnbin+1,*), intent(in) :: habit_p
  ! aspect ratio of circumscribing cylinder (currently not used)
  real(PS), intent(in) :: spx_p

  ! ratio of ag^3 to a^3
  real(PS), dimension(mxnbin+1,*), intent(in) :: rag_p,rcg_p
  ! number of extra ice crystals
  real(RP), dimension(mxnbin+1,*), intent(in) :: n_exice_p
  ! activated IN fraction for contact parameter diagnosis
  real(PS), dimension(mxnbin+1,*), intent(in) :: actINF_p

  ! error number from cal_linprms_vec
  integer,dimension(mxnbin+1,*) :: error_number

  ! source of aerosol particles
  integer,intent(in) :: iaer_src
  ! concentration and mass to be moved into aerosol group
  real(8),dimension(LMAX,*),intent(inout),optional  :: ap_dN, ap_dM
  !   soluble mass, activated IN con, vapor
  real(8),dimension(LMAX,*),intent(inout),optional  :: ap_dMS, ap_dNI, ap_dMV
  integer,dimension(mxnbin,*),intent(inout),optional :: inevp

  ! transferred total concentration and mass from bin to bin
  real(8) :: trans_dN, trans_dM
  real(8) :: trans_dM_dum
  !
  ! transferred volume of circumscribing spheroid
  real(8) :: trans_dvcs

  ! transferred length
  real(8) :: trans_dL(mxnaxis)
  ! mid mass of two boundaries
  real(PS) :: mid_mass

!!c     dummy parameters
!!c      real a(3)
!!c   real, shifted boundaries
  real(8),dimension(mxnbin+1,LMAX) :: left_bd, right_bd
  real(8) :: aleft_in,aright_in,aboth_over,bd1,bd2

  ! dummy concentration and mass
  real(PS) :: d_N, d_M

  ! temporal variables
  real(8) :: xA,xB,xC,xD

  ! error message
  integer :: em
  integer :: ibx,k,i,n,in,icat

  ! type of subdistribution function
  integer dstype

  ! hexagonal length
  integer, parameter :: npreaxis=3

  integer :: i_epslt

  integer, dimension(mxnbin,LMAX) :: icond1
  real(8), dimension(mxnbin,LMAX) :: trans_dM_all

  real(8) :: C1,axr1,axr2
  real(8),dimension(mxnbin+1) :: binb8

  ! soluble mass fraction
  real(PS) :: eps_map
  ! minimum fraction of soluble mass for CCN
  real(PS),parameter :: lmt_frac=0.99999e-5


  ! distribution type
  dstype=1

  do ibx=1,nbin+1
    binb8(ibx)=binb(ibx)
  enddo
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! 1. set up the left and right boundaries of new bins
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!CDIR NODEP
  do in=1,nbin_s*L
    n=(in-1)/NBIN_s+1
    i=in-(n-1)*NBIN_s

    if(error_number(i,n)/=10.and.a2d(i,n,4)<=-9.98d+100) then
      ! linear distribution
      if( a2d(i,n,1) > 0.0d+0 ) then
        ! --- case of non-negative bin ---
        left_bd(i,n)=s_bd(i,n,1)
        right_bd(i,n)=s_bd(i,n,2)

      elseif( a2d(i,n,1) == -1.0d+0 ) then
        ! case of n(x'_1) < 0.0
        left_bd(i,n) = a2d(i,n,2)
        right_bd(i,n) = s_bd(i,n,2)

      else if( a2d(i,n,1) == -2.0d+0 ) then
        ! case of n(x'_2) < 0.0
        left_bd(i,n) = s_bd(i,n,1)
        right_bd(i,n) = a2d(i,n,2)

      else
        left_bd(i,n)=0.0_PS
        right_bd(i,n)=0.0_PS
      end if

      ! re-write the parameters as n(m)=a0+a1*m
      a2d(i,n,1)=max(0.0d+0,a2d(i,n,1))-a2d(i,n,2)*a2d(i,n,3)
      a2d(i,n,2)=a2d(i,n,3)
      a2d(i,n,3)=0.0d+0
      a2d(i,n,4)=0.0d+0

    elseif(error_number(i,n)/=10.and.a2d(i,n,4)>-9.98d+100) then
      ! cubic distribution
      left_bd(i,n)=s_bd(i,n,1)
      right_bd(i,n)=s_bd(i,n,2)

    else
      left_bd(i,n)=0.0_PS
      right_bd(i,n)=0.0_PS
    endif
  enddo

  if( iphase == 1 ) then
    !     --- case where the colletor is a liquid hydromteor ---

    shift_bin_loop1: do i=1,nbin_s

      icond1(:,:) = 0
      trans_dm_all(:,:) = 0


!      do in=1,nbin*L
!        n=(in-1)/nbin+1
!        ibx=in-(n-1)*nbin
      do n = 1, L
      do ibx = 1, nbin

        aleft_in=(left_bd(i,n)-binb8(ibx))*(binb8(ibx+1)-left_bd(i,n))
        aright_in=(right_bd(i,n)-binb8(ibx))*(binb8(ibx+1)-right_bd(i,n))
        aboth_over=(binb8(ibx)-left_bd(i,n))*(right_bd(i,n)-binb8(ibx+1))

        if(error_number(i,n)/=10.and.&
           (aleft_in>0.0_RP.or.aright_in>0.0_RP.or.aboth_over>=0.0_RP)) then
!!!           (aleft_in>-1.0e-30.or.aright_in>-1.0e-30.or.aboth_over>-1.0e-30)) then

          icond1(ibx,n)=1

          bd1=max(binb8(ibx), left_bd(i,n))
          bd2=min(binb8(ibx+1), right_bd(i,n))
          em=0
!!!          call cal_transmom( dstype, a2d(i,n,1) &
!!!                ,bd1,bd2 &
!!!                ,trans_dN, trans_dM)

!org          trans_dN = bd2*a2d(i,n,1)-bd2*a2d(i,n,3)*a2d(i,n,2)+bd2*a2d(i,n,3)*0.5d0*(bd1+bd2)&
!org            -bd1*a2d(i,n,1)+bd1*a2d(i,n,3)*a2d(i,n,2)-bd1*a2d(i,n,3)*0.5d0*(bd1+bd2)

!org2          trans_dN=a2d(i,n,1)*(bd2-bd1) &
!org2                  +0.5d+0*a2d(i,n,2)*(bd2*bd2-bd1*bd1) &
!org2                  +a2d(i,n,3)*(bd2*bd2*bd2-bd1*bd1*bd1)/3.0d+0 &
!org2                  +0.25d+0*a2d(i,n,4)*(bd2*bd2*bd2*bd2-bd1*bd1*bd1*bd1)

          xA=bd2+bd1
          xB=bd2-bd1
          xC=bd2*bd1
          xD=bd2*bd2+bd1*bd1
          trans_dN=a2d(i,n,1)*xB &
                  +0.5_DP*a2d(i,n,2)*xA*xB &
                  +a2d(i,n,3)*xB*(xD+xC)/3.0_DP &
                  +0.25_DP*a2d(i,n,4)*xD*xA*xB

!org          trans_dM = a2d(i,n,1)*0.5d0*bd2*bd2-a2d(i,n,1)*0.5d0*bd1*bd1 &
!org           + a2d(i,n,3)/3.0d0*bd2*bd2*bd2-a2d(i,n,3)/3.0d0*bd1*bd1*bd1 &
!org           - a2d(i,n,3)*a2d(i,n,2)*0.5d0*bd2*bd2+a2d(i,n,3)*a2d(i,n,2)*0.5d0*bd1*bd1

!org2          trans_dM=0.5d+0*a2d(i,n,1)*(bd2*bd2-bd1*bd1) &
!org2                  +a2d(i,n,2)*(bd2*bd2*bd2-bd1*bd1*bd1)/3.0d+0 &
!org2                  +0.25d+0*a2d(i,n,3)*(bd2*bd2*bd2*bd2-bd1*bd1*bd1*bd1) &
!org2                  +0.2d+0*a2d(i,n,4)*(bd2*bd2*bd2*bd2*bd2-bd1*bd1*bd1*bd1*bd1)

          trans_dM=0.5_DP*a2d(i,n,1)*xA*xB &
                  +a2d(i,n,2)*xB*(xD+xC)/3.0_DP &
                  +0.25_DP*a2d(i,n,3)*xD*xA*xB &
                  +0.2_DP*a2d(i,n,4)*xB*(xA*xA*(xD-xC)+xC*xC)

!!!          call ck_errtrans(i,n,ibx&
!!!                ,trans_dN,trans_dM, em)

          trans_dN=max(0.0e+0_DS,trans_dN)
          trans_dM=max(0.0e+0_DS,trans_dM)

          new_N(ibx,n) = new_N(ibx,n) + trans_dN
          new_M(ibx,n,rmt) = new_M(ibx,n,rmt) + trans_dM
          trans_dM_all(ibx,n)=trans_dM_all(ibx,n)+trans_dM
          new_mtend(ibx,n) = new_mtend(ibx,n) + mtend(i,n)*trans_dN

!dbg          if(ibx.eq.nbin) then
!dbg          write(fid_alog,'("ck transbin1:",3I5,20ES15.6)') n,i,ibx,trans_dN,trans_DM &
!dbg             ,a2d(i,n,1:4),left_bd(i,n),right_bd(i,n) &
!dbg             ,aleft_in,aright_in,aboth_over
!dbg          endif

        endif
      enddo
      enddo

      ! case of mass transfering above the max bin boundary,
      !   put them into the largest bin
!CDIR NODEP
      do n=1,L
        if(error_number(i,n)/=10.and.&
           right_bd(i,n)>binb8(nbin+1) ) then
          ibx=nbin

          icond1(ibx,n)=1

          bd1=max(binb8(ibx+1), left_bd(i,n))
          bd2=right_bd(i,n)
          em=0

          xA=bd2+bd1
          xB=bd2-bd1
          xC=bd2*bd1
          xD=bd2*bd2+bd1*bd1
          trans_dN=a2d(i,n,1)*xB &
                  +0.5d+0*a2d(i,n,2)*xA*xB &
                  +a2d(i,n,3)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,4)*xD*xA*xB


          trans_dM=0.5d+0*a2d(i,n,1)*xA*xB &
                  +a2d(i,n,2)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,3)*xD*xA*xB &
                  +0.2d+0*a2d(i,n,4)*xB*(xA*xA*(xD-xC)+xC*xC)

          trans_dN=max(0.0e+0_DS,trans_dN)
          trans_dM=max(0.0e+0_DS,trans_dM)

          new_N(ibx,n) = new_N(ibx,n) + trans_dN
          new_M(ibx,n,rmt) = new_M(ibx,n,rmt) + trans_dM
          trans_dM_all(ibx,n)=trans_dM_all(ibx,n)+trans_dM
          new_mtend(ibx,n) = new_mtend(ibx,n) + mtend(i,n)*trans_dN

!dbg          write(fid_alog,'("ck transbin2:",3I5,20ES15.6)') n,i,ibx,trans_dN,trans_DM &
!dbg             ,a2d(i,n,1:4),left_bd(i,n),right_bd(i,n) &
!dbg             ,aleft_in,aright_in,aboth_over

        endif
      enddo

    if(iaer_src.eq.1) then
      ! case of mass transfering below the min bin boundary,
      !   put them into aerosol variables
!CDIR NODEP
      do n=1,L
        if(error_number(i,n)/=10.and.&
           left_bd(i,n)<binb8(1) ) then
          ibx=1

          bd1=left_bd(i,n)
          bd2=min(binb8(1),right_bd(i,n))
          em=0

          xA=bd2+bd1
          xB=bd2-bd1
          xC=bd2*bd1
          xD=bd2*bd2+bd1*bd1
          trans_dN=a2d(i,n,1)*xB &
                  +0.5d+0*a2d(i,n,2)*xA*xB &
                  +a2d(i,n,3)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,4)*xD*xA*xB


          trans_dM=0.5d+0*a2d(i,n,1)*xA*xB &
                  +a2d(i,n,2)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,3)*xD*xA*xB &
                  +0.2d+0*a2d(i,n,4)*xB*(xA*xA*(xD-xC)+xC*xC)

          trans_dN=max(0.0e+0_DS,trans_dN)
          trans_dM=max(0.0e+0_DS,trans_dM)

          eps_map=ratio_Mp(i,n,rmas_m)/ratio_Mp(i,n,rmat_m)
          i_epslt=0.5*(1.0-sign(1.0_PS,eps_map-lmt_frac))
          icat=i_epslt*2+(1-i_epslt)*1
          ap_dN(n,icat)=ap_dN(n,icat)+trans_dN
          ap_dM(n,icat)=ap_dM(n,icat)+trans_dM*ratio_Mp(i,n,rmat_m)
          ap_dMS(n,icat)=ap_dMS(n,icat)+trans_dM*ratio_Mp(i,n,rmas_m)
          ap_dNI(n,icat)=ap_dNI(n,icat)+actINF_p(i,n)*trans_dN
          ap_dMV(n,icat)=ap_dMV(n,icat)+trans_dM
          inevp(i,n)=icat

!!!          write(fid_alog,*) "i,n,ap_dn",i,n,ap_dn(n),ap_dm(n),trans_dn,trans_dm,left_bd(i,n),right_bd(i,n) &
!!!                  ,a2d(i,n,1:3),error_number(i,n)

        endif
      enddo
    endif

!    do in=1,nbin*L
!       n=(in-1)/nbin+1
!       ibx=in-(n-1)*nbin
    do k = 1,nmass
       do n = 1, L
       do ibx = 1, nbin
          if(icond1(ibx,n)==1) then
             new_M(ibx,n,k+1) = new_M(ibx,n,k+1) &
                              + ratio_Mp(i,n,k) * trans_dM_all(ibx,n)
          endif
       enddo
       enddo
    enddo

    enddo shift_bin_loop1

  else if(iphase==2) then
    !     --- case where the colletor is a solid hydromteor ---

    shift_bin_loop2: do i=1,nbin_s

       icond1(:,:) = 0
       trans_dm_all(:,:) = 0.0_RP

!      do in=1,nbin*L
!        n=(in-1)/nbin+1
!        ibx=in-(n-1)*nbin
       do n = 1, L
       do ibx = 1, nbin


        aleft_in=(left_bd(i,n)-binb8(ibx))*(binb8(ibx+1)-left_bd(i,n))
        aright_in=(right_bd(i,n)-binb8(ibx))*(binb8(ibx+1)-right_bd(i,n))
        aboth_over=(binb8(ibx)-left_bd(i,n))*(right_bd(i,n)-binb8(ibx+1))

        if(error_number(i,n)/=10.and.&
           (aleft_in>0.0_RP.or.aright_in>0.0_RP.or.aboth_over>=0.0_RP)) then
!!!           (aleft_in>-1.0e-30.or.aright_in>-1.0e-30.or.aboth_over>-1.0e-30)) then

          icond1(ibx,n)=1

          bd1=max(binb8(ibx), left_bd(i,n))
          bd2=min(binb8(ibx+1), right_bd(i,n))
          em=0

          xA=bd2+bd1
          xB=bd2-bd1
          xC=bd2*bd1
          xD=bd2*bd2+bd1*bd1
          trans_dN=a2d(i,n,1)*xB &
                  +0.5d+0*a2d(i,n,2)*xA*xB &
                  +a2d(i,n,3)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,4)*xD*xA*xB


          trans_dM=0.5d+0*a2d(i,n,1)*xA*xB &
                  +a2d(i,n,2)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,3)*xD*xA*xB &
                  +0.2d+0*a2d(i,n,4)*xB*(xA*xA*(xD-xC)+xC*xC)


          trans_dN=max(0.0e+0_DS,trans_dN)
          trans_dM=max(0.0e+0_DS,trans_dM)

          new_N(ibx,n) = new_N(ibx,n) + trans_dN
          new_M(ibx,n,imt) = new_M(ibx,n,imt) + trans_dM
          trans_dM_all(ibx,n)=trans_dM_all(ibx,n)+trans_dM
          new_mtend(ibx,n) = new_mtend(ibx,n) + mtend(i,n)*trans_dN
          trans_dM_dum=trans_dM*(1.0_RP-ratio_Mp(i,n,imw_m))
          trans_dvcs = trans_dM_dum/den_ip_p(i,n)



!org          call cal_transvolcs( dstype, a2d(i,n,1), den_ip_p(i,n), &
!org                  bd1,bd2, &
!org                  trans_dM_dum,trans_dvcs)
          new_Q(ibx,n,ivcs) = new_Q(ibx,n,ivcs) + trans_dvcs

!dbg        write(fid_alog,*) "ck transbin",i,n,ibx,trans_dn,trans_dm,new_n(ibx,n),new_m(ibx,n,imt),&
!dbg             new_Q(ibx,n,ivcs)

!!!          if(den_ic_p(i,n)<1.0e-30) then
!!!            write(fid_alog,*) "den_ic_p check",i,n,den_ic_p(i,n),den_ip_p(i,n),trans_dM,ratio_Mp(i,n,imw_m)
!!!          endif

          trans_dM_dum=trans_dM*ratio_Mp(i,n,imc_m)
!org          call cal_translen3( dstype, a2d(i,n,1), den_ic_p(i,n), &
!org                  axr_p(i,n,1),spx_p,habit_p(i,n), &
!org                  bd1,bd2, &
!org                  trans_dM_dum,n_exice_p(i,n),trans_dL)

          C1 = 1.0_PS/(4.18879020478639_PS*den_ic_p(i,n))
          axr1=axr_p(i,n,1)
          axr2=axr_p(i,n,2)
          trans_dL(1) = C1*trans_dM_dum/sqrt(1.0d+0+axr1*axr1)**3
          trans_dL(2) = C1*trans_dM_dum/sqrt(1.0d+0+1.0d+0/axr1/axr1)**3
          trans_dL(3) = trans_dL(1)*axr2*axr2*axr2

          new_Q(ibx,n,iacr) = new_Q(ibx,n,iacr)+trans_dL(1)
          new_Q(ibx,n,iacr+1) = new_Q(ibx,n,iacr+1)+trans_dL(2)
          new_Q(ibx,n,iacr+2) = new_Q(ibx,n,iacr+2)+trans_dL(3)
          new_Q(ibx,n,iag)=new_Q(ibx,n,iag)+trans_dL(1)*rag_p(i,n)
          new_Q(ibx,n,icg)=new_Q(ibx,n,icg)+trans_dL(2)*rcg_p(i,n)
          new_Q(ibx,n,inex)=new_Q(ibx,n,inex)+n_exice_p(i,n)*trans_dN

        endif
      enddo
      enddo
      ! case of mass transfering above the max bin boundary,
      !   put them into the largest bin
!CDIR NODEP
      do n=1,L
        if(error_number(i,n)/=10.and.&
           right_bd(i,n)>binb8(nbin+1) ) then
          ibx=nbin

          icond1(ibx,n)=1

          bd1=max(binb8(ibx+1), left_bd(i,n))
          bd2=right_bd(i,n)
          em=0

          xA=bd2+bd1
          xB=bd2-bd1
          xC=bd2*bd1
          xD=bd2*bd2+bd1*bd1
          trans_dN=a2d(i,n,1)*xB &
                  +0.5d+0*a2d(i,n,2)*xA*xB &
                  +a2d(i,n,3)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,4)*xD*xA*xB


          trans_dM=0.5d+0*a2d(i,n,1)*xA*xB &
                  +a2d(i,n,2)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,3)*xD*xA*xB &
                  +0.2d+0*a2d(i,n,4)*xB*(xA*xA*(xD-xC)+xC*xC)


          trans_dN=max(0.0e+0_DS,trans_dN)
          trans_dM=max(0.0e+0_DS,trans_dM)

          new_N(ibx,n) = new_N(ibx,n) + trans_dN
          new_M(ibx,n,imt) = new_M(ibx,n,imt) + trans_dM
          trans_dM_all(ibx,n)=trans_dM_all(ibx,n)+trans_dM
          new_mtend(ibx,n) = new_mtend(ibx,n) + mtend(i,n)*trans_dN
          trans_dM_dum=trans_dM*(1.0_RP-ratio_Mp(i,n,imw_m))

!org          call cal_transvolcs( dstype, a2d(i,n,1), den_ip_p(i,n), &
!org                  bd1,bd2, &
!org                  trans_dM_dum,trans_dvcs)
          trans_dvcs = trans_dM_dum/den_ip_p(i,n)
          new_Q(ibx,n,ivcs) = new_Q(ibx,n,ivcs) + trans_dvcs

          trans_dM_dum=trans_dM*ratio_Mp(i,n,imc_m)
!org          call cal_translen3( dstype, a2d(i,n,1), den_ic_p(i,n),&
!org                  axr_p(i,n,1),spx_p,habit_p(i,n), &
!org                  bd1,bd2, &
!org                  trans_dM_dum,n_exice_p(i,n),trans_dL)

          C1 = 1.0/(4.18879020478639*den_ic_p(i,n))
          axr1=axr_p(i,n,1)
          axr2=axr_p(i,n,2)
          trans_dL(1) = C1*trans_dM_dum/(1.0d+0+axr1*axr1)**1.5
          trans_dL(2) = C1*trans_dM_dum/(1.0d+0+1.0d+0/axr1/axr1)**1.5
          trans_dL(3) = trans_dL(1)*axr2*axr2*axr2

          new_Q(ibx,n,iacr) = new_Q(ibx,n,iacr)+trans_dL(1)
          new_Q(ibx,n,iacr+1) = new_Q(ibx,n,iacr+1)+trans_dL(2)
          new_Q(ibx,n,iacr+2) = new_Q(ibx,n,iacr+2)+trans_dL(3)
          new_Q(ibx,n,iag)=new_Q(ibx,n,iag)+trans_dL(1)*rag_p(i,n)
          new_Q(ibx,n,icg)=new_Q(ibx,n,icg)+trans_dL(2)*rcg_p(i,n)
          new_Q(ibx,n,inex)=new_Q(ibx,n,inex)+n_exice_p(i,n)*trans_dN

        endif
      enddo

    if(iaer_src.eq.1) then
      ! case of mass transfering below the min bin boundary,
      !   put them into aerosol variables
!CDIR NODEP
      do n=1,L
        if(error_number(i,n)/=10.and.&
           left_bd(i,n)<binb8(1) ) then
          ibx=1

          bd1=left_bd(i,n)
          bd2=min(binb8(1),right_bd(i,n))
          em=0

          xA=bd2+bd1
          xB=bd2-bd1
          xC=bd2*bd1
          xD=bd2*bd2+bd1*bd1
          trans_dN=a2d(i,n,1)*xB &
                  +0.5d+0*a2d(i,n,2)*xA*xB &
                  +a2d(i,n,3)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,4)*xD*xA*xB


          trans_dM=0.5d+0*a2d(i,n,1)*xA*xB &
                  +a2d(i,n,2)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,3)*xD*xA*xB &
                  +0.2d+0*a2d(i,n,4)*xB*(xA*xA*(xD-xC)+xC*xC)


          trans_dN=max(0.0e+0_DS,trans_dN)
          trans_dM=max(0.0e+0_DS,trans_dM)

          eps_map=ratio_Mp(i,n,imas_m)/ratio_Mp(i,n,imat_m)
          i_epslt=0.5*(1.0-sign(1.0_PS,eps_map-lmt_frac))
          icat=i_epslt*2+(1-i_epslt)*1
          ap_dN(n,icat)=ap_dN(n,icat)+trans_dN
          ap_dM(n,icat)=ap_dM(n,icat)+trans_dM*ratio_Mp(i,n,imat_m)
          ap_dMS(n,icat)=ap_dMS(n,icat)+trans_dM*ratio_Mp(i,n,imas_m)
!!!          ap_dNI(n,icat)=ap_dNI(n,icat)+actINF_p(i,n)*trans_dN
          ap_dNI(n,icat)=ap_dNI(n,icat)+trans_dN
          ap_dMV(n,icat)=ap_dMV(n,icat)+trans_dM
          inevp(i,n)=icat
        endif
      enddo
    endif

!    do in=1,nbin*L
!       n=(in-1)/nbin+1
!       ibx=in-(n-1)*nbin
    do k = 1,nmass
       do n = 1, L
       do ibx = 1, nbin
          if(icond1(ibx,n)==1) then
             new_M(ibx,n,k+1) = new_M(ibx,n,k+1) &
                              + ratio_Mp(i,n,k) * trans_dM_all(ibx,n)
          endif
       end do
       end do
    enddo

    enddo shift_bin_loop2
  end if
end subroutine cal_transbin_vec

!!$subroutine ck_errtrans(col, ng, ibx, trans_dN, trans_dM, em)
!!$  implicit none
!!$  real(PS) :: trans_dN, trans_dM
!!$  integer :: col,ng, ibx
!!$  integer :: em
!!$  if( trans_dN < 0.0 ) then
!!$     trans_dN=0.0
!!$     em = 1
!!$     !     stop
!!$  end if
!!$  if( trans_dM < 0.0 ) then
!!$     trans_dM=0.0
!!$     em = 1
!!$     !     stop
!!$  end if
!!$end subroutine ck_errtrans

!!$subroutine cal_transmom( dis_type, a, llmt, ulmt, tdn, tdm )
!!$  implicit none
!!$  !     **********************************************************************
!!$  !     Calculate the area under n(x) and xn(x) bounded by llmt and ulmt
!!$  !     **********************************************************************
!!$  !     type of the distribution for a bin
!!$  integer :: dis_type
!!$  !     lower limit and upper limit of integration
!!$  real(8) :: llmt, ulmt
!!$  !     parameters of the distribution
!!$  !     For linear distribution : n(x) = a(1) + a(3) ( x - a(2) )
!!$  !     For parabolic approx    : n(x) = a(1) x^2 + a(2) x + a(3)
!!$  real(8) :: a(*)
!!$  real(PS) :: tdn, tdm
!!$
!!$  if( dis_type == 1 ) then
!!$     !     --- case of linear distribution ---
!!$!tmp     tdn = (ulmt-llmt)*(a(1)-a(3)*(a(2)-(llmt+ulmt)/2.0))
!!$!tmp     tdm = (a(1)-a(3)*a(2))*0.5*(ulmt**2-llmt**2)&
!!$!tmp          + a(3)*(1.0/3.0)*(ulmt**3-llmt**3)
!!$
!!$     tdn = ulmt*a(1)-ulmt*a(3)*a(2)+ulmt*a(3)*0.5*(llmt+ulmt)&
!!$            -llmt*a(1)+llmt*a(3)*a(2)-llmt*a(3)*0.5*(llmt+ulmt)
!!$
!!$     tdm = a(1)*0.5*ulmt**2.0-a(1)*0.5*llmt**2.0 &
!!$           + a(3)/3.0*ulmt**3.0-a(3)/3.0*llmt**3.0 &
!!$           - a(3)*a(2)*0.5*ulmt**2.0+a(3)*a(2)*0.5*llmt**2.0
!!$
!!$
!!$     !       if( tdn < 0.0_RP) then
!!$     !          write(fid_alog,'("cal_transmom",8E15.8)') a(1),a(2),a(3),tdn,tdm,ulmt,llmt
!!$     !       end if
!!$     !      write(fid_alog,*) "tdn and tdm",tdn,tdm
!!$
!!$  else if( dis_type == 2 ) then
!!$     !     --- case of parabolic distribution ---
!!$     tdn = (ulmt**3.0 - llmt**3.0)*a(1)/3.0 &
!!$          + (ulmt**2.0 - llmt**2.0)*a(2)/2.0 &
!!$          + (ulmt - llmt)*a(3)
!!$     tdm =  (ulmt**4.0 - llmt**4.0)*a(1)/4.0 &
!!$          + (ulmt**3.0 - llmt**3.0)*a(2)/3.0 &
!!$          + (ulmt**2.0 - llmt**2.0)*a(3)/2.0
!!$  else if( dis_type == 3 ) then
!!$     !     --- case of cubic distribution ---
!!$     tdn = (ulmt**4.0 - llmt**4.0)*a(1)/4.0 &
!!$          + (ulmt**3.0 - llmt**3.0)*a(2)/3.0 &
!!$          + (ulmt**2.0 - llmt**2.0)*a(3)/2.0 &
!!$          + (ulmt - llmt)*a(4)
!!$     tdm = (ulmt**5.0 - llmt**5.0)*a(1)/5.0 &
!!$          + (ulmt**4.0 - llmt**4.0)*a(2)/4.0 &
!!$          + (ulmt**3.0 - llmt**3.0)*a(3)/3.0 &
!!$          + (ulmt**2.0 - llmt**2.0)*a(4)/2.0
!!$!org     if( tdn < 0.0 .OR. tdm < 0.0 ) then
!!$!org        write(fid_alog,*) "cal_transmom > tdn and tdm are negative."
!!$!org        stop
!!$!org     end if
!!$  end if
!!$end subroutine cal_transmom

!!$subroutine cal_transmom_lin(a, llmt, ulmt, tdn, tdm )
!!$  implicit none
!!$  !     **********************************************************************
!!$  !     Calculate the area under n(x) and xn(x) bounded by llmt and ulmt
!!$  !     **********************************************************************
!!$  !     lower limit and upper limit of integration
!!$  real(8) :: llmt, ulmt
!!$  !     parameters of the distribution
!!$  !     For linear distribution : n(x) = a(1) + a(3) ( x - a(2) )
!!$  real(8) :: a(3)
!!$  real(PS) :: tdn, tdm
!!$
!!$   !     --- case of linear distribution ---
!!$!tmp     tdn = (ulmt-llmt)*(a(1)-a(3)*(a(2)-(llmt+ulmt)/2.0))
!!$!tmp     tdm = (a(1)-a(3)*a(2))*0.5*(ulmt**2-llmt**2)&
!!$!tmp          + a(3)*(1.0/3.0)*(ulmt**3-llmt**3)
!!$
!!$     tdn = ulmt*a(1)-ulmt*a(3)*a(2)+ulmt*a(3)*0.5*(llmt+ulmt)&
!!$            -llmt*a(1)+llmt*a(3)*a(2)-llmt*a(3)*0.5*(llmt+ulmt)
!!$
!!$     tdm = a(1)*0.5*ulmt**2.0-a(1)*0.5*llmt**2.0 &
!!$           + a(3)/3.0*ulmt**3.0-a(3)/3.0*llmt**3.0 &
!!$           - a(3)*a(2)*0.5*ulmt**2.0+a(3)*a(2)*0.5*llmt**2.0
!!$
!!$end subroutine cal_transmom_lin

!!$subroutine cal_translen3(dis_type, a, den_ic_p,axr_p,spx_p,habit_p &
!!$     ,llmt, ulmt, tdM, n_exice_p,tdL)
!!$  implicit none
!!$  !     **********************************************************************
!!$  !     Calculate transferred lengths into a bin bounded by llmt and ulmt
!!$  !
!!$  !     This subroutine treats a and c lengths in terms of volume dimension.
!!$  !     The volume dimension of length is weighted by concentration.
!!$  !     **********************************************************************
!!$  !     type of the distribution for a bin
!!$  integer :: dis_type
!!$  !     lower limit and upper limit of integration
!!$  real(8) :: llmt, ulmt
!!$  !     transferred mass
!!$  real(PS) :: tdM
!!$  !     parameters of the distribution
!!$  !     For linear distribution : n(x) = a(1) + a(3) ( x - a(2) )
!!$  real(8) :: a(*)
!!$  real(RP) :: tdL(*)
!!$  !     apparent density in the shfited bin of the ice crystal in a shifted bin
!!$  real(PS) :: den_ic_p
!!$  !     axis ratio for ice crystals  in a shifted bin
!!$  !     1: c/a
!!$  !     2: d/a
!!$  !     3: r/a
!!$  !     4: e/a
!!$  real(PS) :: axr_p(*)
!!$  !     coefficients
!!$  real(PS) :: C1
!!$  !     habit in shifted bin
!!$  integer :: habit_p
!!$  !     aspect ratio of circumscribing cylinder
!!$  real(PS) :: spx_p
!!$  !  number of extra monocrystals
!!$  real(PS) ::  n_exice_p
!!$
!!$!!!  if( dis_type == 1 ) then
!!$     !     --- case of linear distribution ---
!!$
!!$     !     --- case of spheroidal ice crystals
!!$     !     C1 = 3.0/(4.0*PI*den_i)
!!$     !     C1 = 1.0/(den_i*sq_three*3.0)
!!$     !     --- case of plate-column-dendrite ice crystal model
!!$     !     C1 = 1.0/(den_i*sq_three*(1.0-psi_p)*(3.0-psi_p))
!!$     !     --- use the apparent density calculated with circumscribing sphere.
!!$     C1 = 1.0/(4.18879020478639*den_ic_p)
!!$
!!$     ! currently always habit_p==1
!!$!!!     if(habit_p==1) then
!!$        tdL(1) = C1*tdM/(1.0+axr_p(1)**2.0)**1.5
!!$        tdL(2) = C1*tdM/(1.0+axr_p(1)**(-2.0))**1.5
!!$        tdL(3) = tdL(1)*(axr_p(2)**3.0)
!!$        tdL(4) = tdL(1)*(axr_p(3)**3.0)
!!$        tdL(5) = tdL(1)*(axr_p(4)**3.0)
!!$!!!     elseif(habit_p==4) then
!!$!!!        tdL(1) = C1*tdM/(1.0+axr_p(1)**2.0)**1.5
!!$!!!        tdL(2) = C1*tdM/(1.0+axr_p(1)**(-2.0))**1.5
!!$!!!        tdL(3) = tdL(1)*(axr_p(2)**3.0)
!!$!!!        tdL(4) = tdL(1)*(axr_p(3)**3.0)
!!$!!!        tdL(5) = tdL(1)*(axr_p(4)**3.0)
!!$!!c        tdL(4) = C1*tdM/(1.0+spx_p**2.0)**1.5
!!$!!c        tdL(1) = tdL(4)/(axr_p(3)**3.0)
!!$!!c        tdL(2) = tDL(1)*(axr_p(1)**3.0)
!!$!!c        tdL(3) = tdL(1)*(axr_p(2)**3.0)
!!$!!c        tdL(5) = tdL(1)*(axr_p(4)**3.0)
!!$!!!     elseif(habit_p==5) then
!!$!!!        tdL(1) = C1*tdM/(1.0+axr_p(1)**2.0)**1.5
!!$!!!        tdL(2) = C1*tdM/(1.0+axr_p(1)**(-2.0))**1.5
!!$!!!        tdL(3) = tdL(1)*(axr_p(2)**3.0)
!!$!!!        tdL(4) = tdL(1)*(axr_p(3)**3.0)
!!$!!!        tdL(5) = tdL(1)*(axr_p(4)**3.0)
!!$!!c        tdL(5) = C1*tdM/(1.0+spx_p**2.0)**1.5
!!$!!c        tdL(1) = tdL(5)/(axr_p(4)**3.0)
!!$!!c        tdL(2) = tDL(1)*(axr_p(1)**3.0)
!!$!!c        tdL(3) = tdL(1)*(axr_p(2)**3.0)
!!$!!c        tdL(4) = tdL(1)*(axr_p(3)**3.0)
!!$!!!     end if
!!$!!!  end if
!!$end subroutine cal_translen3

!!$subroutine cal_transvolcs( dis_type, a, den_ip_p, &
!!$                           llmt, ulmt, &
!!$                           tdM,tdvcs)
!!$  implicit none
!!$  !     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$  !     Calculate transferred volume of circumscribing spheroid
!!$  !     into a bin bounded by llmt and ulmt
!!$  !
!!$  !     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$  !     type of the distribution for a bin
!!$  integer :: dis_type
!!$  !     lower limit and upper limit of integration
!!$  real(8) :: llmt, ulmt
!!$  !     transferred mass
!!$  real(PS) :: tdM
!!$  !     parameters of the distribution
!!$  !     For linear distribution : n(x) = a(1) + a(3) ( x - a(2) )
!!$  real(8) :: a(*)
!!$  real(PS) :: tdvcs
!!$  !     density in the shfited bin
!!$  real(PS) :: den_ip_p
!!$
!!$!!!  if( dis_type == 1 ) then
!!$     !     --- case of linear distribution ---
!!$     !     this is weighting by concentration
!!$     tdvcs = tdM/den_ip_p
!!$!!!  end if
!!$end subroutine cal_transvolcs

subroutine cal_mmass(mmassr,k1r,k2r,npr,nbr,ncr,qrv,n1,isnow)
  use par_amps
  implicit none
  integer,intent(in) :: k1r,k2r,n1,npr,nbr,ncr
  integer,intent(in) :: isnow
  real(MP_KIND),intent(in) :: qrv(npr,nbr,ncr,*)
  real(MP_KIND),intent(inout) :: mmassr(nbr,ncr,*)

  integer :: k,ibr,icr
  integer :: jmt_q, jcon_q

  select case(isnow)
  case(1)
    jmt_q=imt_q
    jcon_q=icon_q
  case default
    jmt_q=rmt_q
    jcon_q=rcon_q
  end select

  do k=k1r,k2r
    do icr=1,ncr
      do ibr=1,nbr
        mmassr(ibr,icr,k)=qrv(jmt_q,ibr,icr,k)/max(qrv(jcon_q,ibr,icr,k),1.0e-30_RP)
!!c        write(fid_alog,*) isnow,k,icr,ibr,mmassr(ibr,icr,k) &
!!c                ,qrv(jmt_q,ibr,icr,k),qrv(jcon_q,ibr,icr,k)
      end do
    end do
  end do
end subroutine cal_mmass
subroutine cal_mmass_scale(mmassr,llbin,k1r,k2r,npr,nbr,ncr,qrv,isnow)
  use par_amps
  implicit none
  integer,intent(in) :: llbin,k1r,k2r,npr,nbr,ncr
  integer,intent(in) :: isnow
  real(MP_KIND),intent(in) :: qrv(npr,nbr,ncr,llbin)
  real(MP_KIND),intent(inout) :: mmassr(nbr,ncr,llbin)

  integer :: k,ibr,icr
  integer :: jmt_q, jcon_q

  select case(isnow)
  case(1)
    jmt_q=imt_q
    jcon_q=icon_q
  case default
    jmt_q=rmt_q
    jcon_q=rcon_q
  end select

  do k=k1r,k2r
    do icr=1,ncr
      do ibr=1,nbr
        mmassr(ibr,icr,k)=qrv(jmt_q,ibr,icr,k)/max(qrv(jcon_q,ibr,icr,k),1.0e-30_RP)
!!c        write(fid_alog,*) isnow,k,icr,ibr,mmassr(ibr,icr,k) &
!!c                ,qrv(jmt_q,ibr,icr,k),qrv(jcon_q,ibr,icr,k)
      end do
    end do
  end do
end subroutine cal_mmass_scale
!!$subroutine cal_linprms( N, M, bd1, bd2, a, error_number,from)
!!$  !     **********************************************************************
!!$  !     Calculate the parameters of the linear distribution:
!!$  !     n(x) = n_0 + k ( x - x_0 )
!!$  !     = a(1)   + a(3) ( x - a(2))
!!$  !     **********************************************************************
!!$  use com_amps, only: fid_alog
!!$  real(8) :: N, M, bd1, bd2
!!$  real(8) :: a(*)
!!$  real(8) :: n_0, x_0, k
!!$  integer :: error_number
!!$  character (len=*) :: from
!!$
!!$  error_number = 0
!!$
!!$  if(N<=0.0_RP) then
!!$     write(fid_alog,*) "cal_linear_prms > ERROR! The concentration is negative."
!!$     write(fid_alog,*) "    bd1,bd2 :",bd1,bd2
!!$     write(fid_alog,*) "       M,N,from :",M,N,from
!!$     error_number=5
!!$     return
!!$     !         stop
!!$  elseif(M<=0.0_RP) then
!!$     write(fid_alog,*) "cal_linear_prms > ERROR! The mass is negative."
!!$     write(fid_alog,*) "    bd1,bd2 :",bd1,bd2
!!$     write(fid_alog,*) "       M,N,from :",M,N,from
!!$     error_number=6
!!$     return
!!$     !         stop
!!$  elseif((bd1<1.0e-30_RP.and.bd2<1.0e-30_RP).or.bd1==bd2) then
!!$     error_number=4
!!$     return
!!$     !         stop
!!$  end if
!!$
!!$  x_0 = (bd1+bd2)*0.5_RP
!!$  n_0 = N/(bd2-bd1)
!!$  k = 12.0_RP*(M-x_0*N)/(bd2-bd1)**3.0_RP
!!$
!!$
!!$  if( n_0 +k*( bd1 - x_0) < 0.0 ) then
!!$     x_0 = 3.0_RP*M/N - 2.0_RP*bd2
!!$     k = 2.0_RP*N/max((bd2 - x_0)*(bd2 - x_0), 1e-300_DS)
!!$     n_0 = -1.0_RP
!!$     if(x_0 > bd2) then
!!$        !!c            write(fid_alog,*) "cal_linear_prms > ERROR! x* > x_2.",from
!!$        k = 0.0_RP
!!$        x_0 = 0.0_RP
!!$        error_number = 1
!!$     end if
!!$  else if( n_0 +k*( bd2 - x_0) < 0.0_RP ) then
!!$     x_0 = 3.0_RP*M/N - 2.0*bd1
!!$     k = - 2.0_RP*N/max((bd1 - x_0)*(bd1 - x_0), 1e-300_DS)
!!$     n_0 = -2.0_RP
!!$     if(x_0 < 0.0_RP) then
!!$        !!c            write(fid_alog,*) "cal_linear_prms > ERROR! The x* is negative."
!!$        !!c     *           ,from
!!$        k = 0.0_RP
!!$        x_0 = 0.0_RP
!!$        error_number = 2
!!$     else if(x_0 < bd1 ) then
!!$        !!c            write(fid_alog,*) "cal_linear_prms > ERROR! x* < x_1."
!!$        !!c     *           ,from
!!$        k = 0.0_RP
!!$        x_0 = 0.0_RP
!!$        error_number = 3
!!$     end if
!!$  end if
!!$
!!$!!c  if(abs(k)>1.0e+20) then
!!$!!c     write(fid_alog,*) "cal_linear_prms > ERROR! The gradient is too big."
!!$!!c     write(fid_alog,'("    n_0,k,x_0,bd1,bd2 :",5ES15.6)') n_0,k,x_0,bd1,bd2
!!$!!c     write(fid_alog,*) "       M,N,from :",M,N,from
!!$!!c     error_number=4
!!$!!c     return
!!$!!c  end if
!!$
!!$  a(1) = n_0
!!$  a(2) = x_0
!!$  a(3) = k
!!$  return
!!$end subroutine cal_linprms

!!$subroutine cal_linprms_vec(nbin_mx,nbin,L, con, amass, binb3d, a, error_number,from)
!!$  !     **********************************************************************
!!$  !     Calculate the parameters of the linear distribution:
!!$  !     n(x) = n_0 + k ( x - x_0 )
!!$  !     = a(1)   + a(3) ( x - a(2))
!!$  !     **********************************************************************
!!$  integer,intent(in) :: nbin_mx,nbin,L
!!$  real(8),dimension(nbin_mx,*),intent(in) :: con,amass
!!$  real(PS),dimension(2,nbin_mx,*),intent(in) :: binb3d
!!$  real(8),dimension(3,nbin_mx,*) :: a
!!$  real(8) :: n_0, x_0, k,bd1,bd2
!!$  integer,dimension(nbin_mx,*) :: error_number
!!$  character (len=*) :: from
!!$  integer :: in,i,n
!!$
!!$  error_number(1:nbin,1:L) = 0
!!$
!!$!CDIR NODEP
!!$  do in=1,nbin*L
!!$    n=(in-1)/NBIN+1
!!$    i=in-(n-1)*NBIN
!!$
!!$    if( (binb3d(1,i,n)<1.0e-30.and.binb3d(2,i,n)<1.0e-30).or.&
!!$        (binb3d(1,i,n)==binb3d(2,i,n))) then
!!$      a(1,i,n)=0.0d+0
!!$      a(2,i,n)=0.0d+0
!!$      a(3,i,n)=0.0d+0
!!$      error_number(i,n)=4
!!$    elseif(con(i,n)>0.0d+0.and.amass(i,n)>0.0d+0) then
!!$
!!$      bd1=binb3d(1,i,n)
!!$      bd2=binb3d(2,i,n)
!!$      x_0=0.5d+0*(bd1+bd2)
!!$      a(2,i,n)=x_0
!!$      a(1,i,n)=con(i,n)/(bd2-bd1)
!!$      a(3,i,n)=12.0d+0*(amass(i,n)-x_0*con(i,n))/(bd2-bd1)**3
!!$    else
!!$      a(1,i,n)=0.0d+0
!!$      a(2,i,n)=0.0d+0
!!$      a(3,i,n)=0.0d+0
!!$      error_number(i,n)=10
!!$    endif
!!$  enddo
!!$
!!$!CDIR NODEP
!!$  do in=1,nbin*L
!!$    n=(in-1)/NBIN+1
!!$    i=in-(n-1)*NBIN
!!$
!!$    if( a(1,i,n) +a(3,i,n)*( real(binb3d(1,i,n),8) - a(2,i,n)) < 0.0d+0 ) then
!!$      x_0=3.0d+0*amass(i,n)/con(i,n) - 2.0d+0*real(binb3d(2,i,n),8)
!!$      a(2,i,n) =x_0
!!$      a(3,i,n) = 2.0d+0*con(i,n)/max((real(binb3d(2,i,n),8) - x_0)**2,1.0d-100)
!!$      a(1,i,n) = -1.0d+0
!!$
!!$      error_number(i,n)=floor(0.5*(1.0-sign(1.0d+0,real(binb3d(2,i,n),8)-x_0))+1.0e-5)  ! 1
!!$
!!$    elseif( a(1,i,n) +a(3,i,n)*( real(binb3d(2,i,n),8) - a(2,i,n)) < 0.0d+0 ) then
!!$      x_0=3.0d+0*amass(i,n)/con(i,n) - 2.0d+0*real(binb3d(1,i,n),8)
!!$      a(2,i,n) = x_0
!!$      a(3,i,n) = - 2.0d+0*con(i,n)/max((real(binb3d(1,i,n),8) - x_0)**2,1.0d-100)
!!$      a(1,i,n) = -2.0d+0
!!$
!!$      error_number(i,n)=floor(1.5*(1.0-sign(1.0d+0,x_0-real(binb3d(1,i,n),8)))+1.0e-5)   ! 3
!!$      error_number(i,n)=floor(1.0*(1.0-sign(1.0d+0,x_0))+1.0e-5)   ! 2
!!$
!!$    endif
!!$  enddo
!!$
!!$  return
!!$end subroutine cal_linprms_vec

subroutine cal_linprms_vec_s(con, amass, binb1,binb2, a, error_number)
  !     **********************************************************************
  !     Calculate the parameters of the linear distribution:
  !     n(x) = n_0 + k ( x - x_0 )
  !     = a(1)   + a(3) ( x - a(2))
  !     **********************************************************************
  real(8),intent(in) :: con,amass
  real(PS),intent(in) :: binb1,binb2
  real(8),dimension(*),intent(inout) :: a
  integer,intent(inout) :: error_number
  real(8) :: n_0, x_0, k,bd1,bd2
  integer :: i_neg_1,i_neg_2,i1,i2,i3

  error_number=0

  bd1=binb1
  bd2=binb2

  x_0=0.5d+0*(bd1+bd2)
  a(2)=x_0
  a(1)=con/(bd2-bd1)
  a(3)=12.0d+0*(amass-x_0*con)/(bd2-bd1)**3


!!!  if( a(1) +a(3)*( bd1 - a(2)) < 0.0d+0 ) then
  i_neg_1=int(0.5d+0*(1.0d+0-sign(1.0d+0, a(1) +a(3)*( bd1 - a(2))))) ! CHIARUI
!!!  elseif( a(1) +a(3)*( bd2 - a(2)) < 0.0d+0 ) then
  i_neg_2=int(0.5d+0*(1.0d+0-sign(1.0d+0, a(1) +a(3)*( bd2 - a(2))))) ! CHIARUI

  i1=(1-i_neg_1)*(1-i_neg_2)
  i2=(1-i1)*i_neg_1*(1-i_neg_2)
  i3=(1-i1)*(1-i_neg_1)*i_neg_2

  a(1)= real(i1,8)*a(1) +&
        real(i2,8)*(-1.0d+0) +&
        real(i3,8)*(-2.0d+0)

  x_0= real(i1,8)*( 0.5d+0*(bd1+bd2)) +&
       real(i2,8)*( 3.0d+0*amass/con - 2.0d+0*bd2 ) +&
       real(i3,8)*( 3.0d+0*amass/con - 2.0d+0*bd1 )

  a(2)= x_0

  a(3)= real(i1,8)*a(3) +&
        real(i2,8)*(  2.0d+0*con/max((bd2 - x_0)**2,1.0d-100) ) + &
        real(i3,8)*( -2.0d+0*con/max((bd1 - x_0)**2,1.0d-100) )

  error_number=i1*0 + &
               i2*floor(0.5*(1.0-sign(1.0d+0,bd2-x_0))+1.0e-5) +&
               i3*floor(1.5*(1.0-sign(1.0d+0,x_0-bd1))+1.0e-5)

  ! indication of linear distribution
  a(4)=-9.99d+100

  return
end subroutine cal_linprms_vec_s

subroutine cal_lincubprms_vec(nbin_mx,nbin,L, con, amass, binb3d, a, error_number,from)
  !     **********************************************************************
  !     Calculate the parameters of the linear distribution and cubic
  !     distributions:
  !     n(x) = n_0 + k ( x - x_0 )
  !          = max(0,a(1))   + a(3) ( x - a(2))
  !          = max(0,a(1))-a(2)*a(3) + a(3)*x
  !
  !    Reference
  !      Dinh and Durran (2012), Atm. Chem. Phys.
  !     **********************************************************************
  integer,intent(in) :: nbin_mx,nbin,L
  real(8),dimension(nbin_mx,*),intent(in) :: con,amass
  real(PS),dimension(nbin_mx,L,2),intent(in) :: binb3d
  real(8),dimension(nbin_mx,L,4) :: a
  integer,dimension(nbin_mx,*) :: error_number
  character (len=*) :: from
  !
  ! new space
  real(8),dimension(nbin_mx,L,4) :: a_cub
  real(8) :: n_0, x_0, k,bd1,bd2
  real(8) :: delta_m,am0
  real(8) :: b0,b1,b2,b3,a0_c,a1_c,a2_c,a3_c
  real(8) :: amb_L,anb_L,xib_L,A_L,B_L,C_L
  real(8) :: amb_R,anb_R,xib_R,A_R,B_R,C_R
  real(8) :: q,delta,am_e,an_me,an_ml,an_mr
  integer :: in,i,n
  integer :: i_con_gt0_left,i_con_gt0_right,i_delta_ge0 &
            ,i_ame_lt0,i_cubic,i_anl_lt0,i_anr_lt0 &
            ,i_anml_lt0,i_anmr_lt0,i_large


  error_number(1:nbin,1:L) = 0

!
!  first, fit a piece-wise linear
!

!CDIR NODEP
  do in=1,nbin*L
    n=(in-1)/NBIN+1
    i=in-(n-1)*NBIN

    if( (binb3d(i,n,1)<1.0e-30_PS.and.binb3d(i,n,2)<1.0e-30_PS).or.&
        (binb3d(i,n,1)==binb3d(i,n,2))) then
      a(i,n,1)=0.0d+0
      a(i,n,2)=0.0d+0
      a(i,n,3)=0.0d+0
      a(i,n,4)=0.0d+0
      error_number(i,n)=4
!!!      write(fid_alog,*) i,n,a(i,n,1:4)
    elseif(con(i,n)>0.0d+0.and.amass(i,n)>0.0d+0) then

      bd1=binb3d(i,n,1)
      bd2=binb3d(i,n,2)
      x_0=0.5d+0*(bd1+bd2)
      a(i,n,2)=x_0
      a(i,n,1)=con(i,n)/(bd2-bd1)
      a(i,n,3)=12.0d+0*(amass(i,n)-x_0*con(i,n))/(bd2-bd1)**3

      a(i,n,4)=-9.99d+100 ! indication of linear

    else
      a(i,n,1)=0.0d+0
      a(i,n,2)=0.0d+0
      a(i,n,3)=0.0d+0
      a(i,n,4)=0.0d+0
      error_number(i,n)=10
    endif
  enddo

!CDIR NODEP
  do in=1,nbin*L
    n=(in-1)/NBIN+1
    i=in-(n-1)*NBIN

    if( a(i,n,1) +a(i,n,3)*( binb3d(i,n,1) - a(i,n,2)) < 0.0d+0 ) then
      x_0=3.0d+0*amass(i,n)/con(i,n) - 2.0d+0*binb3d(i,n,2)
      a(i,n,2) =x_0
      a(i,n,3) = 2.0d+0*con(i,n)/max((binb3d(i,n,2) - x_0)**2,1.0d-100)
      a(i,n,1) = -1.0d+0

      error_number(i,n)=floor(0.5*(1.0-sign(1.0d+0,binb3d(i,n,2)-x_0))+1.0e-5)  ! 1

    elseif( a(i,n,1) +a(i,n,3)*( binb3d(i,n,2) - a(i,n,2)) < 0.0d+0 ) then
      x_0=3.0d+0*amass(i,n)/con(i,n) - 2.0d+0*binb3d(i,n,1)
      a(i,n,2) = x_0
      a(i,n,3) = - 2.0d+0*con(i,n)/max((binb3d(i,n,1) - x_0)**2,1.0d-100)
      a(i,n,1) = -2.0d+0

      error_number(i,n)=floor(1.5*(1.0-sign(1.0d+0,x_0-binb3d(i,n,1)))+1.0e-5)   ! 3
      error_number(i,n)=floor(1.0*(1.0-sign(1.0d+0,x_0))+1.0e-5)   ! 2

    endif
!    write(fid_alog,*) i,n,a(i,n,1:4)

!    if(error_number(i,n)==0) then
!      write(fid_alog,*) i,n,a(i,n,1),a(i,n,2),a(i,n,3)
    if( dabs(a(i,n,1))>1.0d+100 .or. dabs(a(i,n,2))>1.0d+100 .or. dabs(a(i,n,3))>1.0d+100) then
      error_number(i,n)=4
    endif
!    endif
  enddo

!!!  return

!
! second, calculate the cubic
!

!CDIR NODEP
  do in=1,nbin*L
    n=(in-1)/NBIN+1
    i=in-(n-1)*NBIN

!    if(( error_number(i,n)/=4.or.error_number(i,n)/=10) .and. &
!       ( error_number(max(1,i-1),n)/=4.or.error_number(max(1,i-1),n)/=10) .and. &
!       ( error_number(min(NBIN,i+1),n)/=4.or.error_number(min(NBIN,i+1),n)/=10) .and. &
    if( ( error_number(i,n)==0.and. &
          error_number(max(1,i-1),n)==0.and. &
          error_number(min(NBIN,i+1),n)==0  ).and. &
         (i > 1 .and. i < NBIN) ) then

      bd1=binb3d(i,n,1)
      bd2=binb3d(i,n,2)
      delta_m=bd2-bd1
      am0=0.5*(bd1+bd2)

      b0=con(i,n)/delta_m
      b1=6.0*(amass(i,n)-am0*con(i,n))/(delta_m*delta_m)

      ! negative flags
      i_con_gt0_left =int(0.5d+0*(1.0d+0-sign(1.0d+0,1.0d-30-con(i-1,n)))) ! CHIARUI
      i_con_gt0_right=int(0.5d+0*(1.0d+0-sign(1.0d+0,1.0d-30-con(i+1,n)))) ! CHIARUI

      ! left mean mass
      amb_L=amass(i-1,n)/con(i-1,n)
      ! number concentration at the mean mass
      anb_L=max(0.0d+0,a(i-1,n,1))+a(i-1,n,3)*(amb_L-a(i-1,n,2))
      ! mean Kai value
      xib_L=2.0d0*(amb_L-am0)/delta_m

      ! right mean mass
      amb_R=amass(i+1,n)/con(i+1,n)
      ! number concentration at the mean mass
      anb_R=max(0.0d+0,a(i+1,n,1))+a(i+1,n,3)*(amb_R-a(i+1,n,2))
      ! mean Kai value
      xib_R=2.0d0*(amb_R-am0)/delta_m

      A_L=0.5d0*(3.0d0*xib_L*xib_L-1.0d0)
      B_L=0.5d0*(5.0d0*xib_L*xib_L*xib_L-3.0d0*xib_L)
      C_L=anb_L-b0-b1*xib_L

      A_R=0.5d0*(3.0d0*xib_R*xib_R-1.0d0)
      B_R=0.5d0*(5.0d0*xib_R*xib_R*xib_R-3.0d0*xib_R)
      C_R=anb_R-b0-b1*xib_R

      b3=(C_L*A_R-A_L*C_R)/B_L*A_R-A_L*B_R
      b2=(C_L-b3*B_L)/A_L
!      b3=(C_L*A_R-A_L*C_R)/max(1.0d-150,B_L*A_R-A_L*B_R)
!      b2=(C_L-b3*B_L)/max(1.0d-150,A_L)

!      write(fid_alog,'("ck b:",2I5,20ES17.6e3)') i,n,b0,b1,b2,b3,amb_L,anb_L,amb_R,anb_R
!      write(fid_alog,'("ck b2:",2I5,20ES17.6e3)') &
!          i,n,C_L*A_R-A_L*C_R,B_L*A_R-A_L*B_R,C_L,A_R,A_L,C_R
!      write(fid_alog,*) "ck b3:",i,n,a(1:3,i-1,n)
!      write(fid_alog,*) "ck b4:",i,n,error_number(i-1,n)


      q=am0/delta_m
      a0_c=b0-2.0d0*q*b1+(6.0d0*q*q-0.5d0)*b2-(20.0d0*q*q*q-3.0d0*q)*b3
      a1_c=(2.0d0*b1-12.0d0*q*b2+(60.0d0*q*q-3.0d0)*b3)/delta_m
      a2_c=(6.0d0*b2-60.0d0*q*b3)/(delta_m*delta_m)
      a3_c=20.0d0*b3/(delta_m*delta_m*delta_m)

      ! criteria for negative valune in the bin
      delta=a2_c*a2_c-3.0*a1_c*a3_c
      am_e=(-a2_c+sqrt(max(0.0d+0,delta)))/(3.0d0*a3_c)
      an_me=a0_c+a1_c*am_e+a2_c*am_e*am_e+a3_c*am_e*am_e*am_e
      i_delta_ge0=int(0.5d+0*(1.0d+0+sign(1.0d+0,delta))) ! CHIARUI
      i_ame_lt0=int(0.5d+0*(1.0d+0-sign(1.0d+0,an_me))) ! CHIARUI

      an_ml=a0_c+a1_c*bd1+a2_c*bd1*bd1+a3_c*bd1*bd1*bd1
      an_mr=a0_c+a1_c*bd2+a2_c*bd2*bd2+a3_c*bd2*bd2*bd2
      i_anml_lt0=int(0.5d+0*(1.0d+0-sign(1.0d+0,an_ml))) ! CHIARUI
      i_anmr_lt0=int(0.5d+0*(1.0d+0-sign(1.0d+0,an_mr))) ! CHIARUI

      ! too large coefficients
      !   this seems to cause infinity in shift bin
      i_large=int(0.5d+0*(1.0d+0+sign(1.0d+0,dabs(a0_c)-1.0d+30))) ! CHIARUI

      i_cubic=i_con_gt0_left*i_con_gt0_right* &
           (1-i_delta_ge0*i_ame_lt0)*(1-i_anml_lt0)*(1-i_anmr_lt0)* &
           (1-i_large)


      a_cub(i,n,1)=dble(i_cubic)*a0_c +&
               dble(1-i_cubic)*a(i,n,1)

      a_cub(i,n,2)=dble(i_cubic)*a1_c +&
               dble(1-i_cubic)*a(i,n,2)

      a_cub(i,n,3)=dble(i_cubic)*a2_c +&
               dble(1-i_cubic)*a(i,n,3)

      a_cub(i,n,4)=dble(i_cubic)*a3_c +&
               dble(1-i_cubic)*(-9.99d+100)

    endif
  enddo

!
! Finally, update a parameters
!
!CDIR NODEP
  do in=1,nbin*L
    n=(in-1)/NBIN+1
    i=in-(n-1)*NBIN

    if( ( error_number(i,n)==0.and. &
          error_number(max(1,i-1),n)==0.and. &
          error_number(min(NBIN,i+1),n)==0  ).and. &
         (i > 1 .and. i < NBIN) ) then

      a(i,n,1)=a_cub(i,n,1)
      a(i,n,2)=a_cub(i,n,2)
      a(i,n,3)=a_cub(i,n,3)
      a(i,n,4)=a_cub(i,n,4)

    endif
  enddo

  return
end subroutine cal_lincubprms_vec

subroutine cdfnor ( which, p, q, x, mean, sd, status, bound )

!*******************************************************************************
!
!! CDFNOR evaluates the CDF of the Normal distribution.
!
!  Discussion:
!
!    A slightly modified version of ANORM from SPECFUN
!    is used to calculate the cumulative standard normal distribution.
!
!    The rational functions from pages 90-95 of Kennedy and Gentle
!    are used as starting values to a Newton iteration which
!    compute the inverse standard normal.  Therefore no searches are
!    necessary for any parameter.
!
!    For X < -15, the asymptotic expansion for the normal is used  as
!    the starting value in finding the inverse standard normal.
!
!    The normal density is proportional to
!    exp ( - 0.5D+00 * (( X - MEAN)/SD)**2)
!
!  Reference:
!
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions
!    1966, Formula 26.2.12.
!
!    W J Cody,
!    Algorithm 715: SPECFUN - A Portable FORTRAN Package of
!      Special Function Routines and Test Drivers,
!    ACM Transactions on Mathematical Software,
!    Volume 19, Number 1, pages 22-32, 1993.
!
!    Kennedy and Gentle,
!    Statistical Computing,
!    Marcel Dekker, NY, 1980,
!    QA276.4 K46
!
!  Parameters:
!
!    Input, integer WHICH, indicates which argument is to be calculated
!    from the others.
!    1: Calculate P and Q from X, MEAN and SD;
!    2: Calculate X from P, Q, MEAN and SD;
!    3: Calculate MEAN from P, Q, X and SD;
!    4: Calculate SD from P, Q, X and MEAN.
!
!    Input/output, real ( kind = 8 ) P, the integral from -infinity to X
!    of the Normal density.  If this is an input or output value, it will
!    lie in the range [0,1].
!
!    Input/output, real ( kind = 8 ) Q, equal to 1-P.  If Q is an input
!    value, it should lie in the range [0,1].  If Q is an output value,
!    it will lie in the range [0,1].
!
!    Input/output, real ( kind = 8 ) X, the upper limit of integration of
!    the Normal density.
!
!    Input/output, real ( kind = 8 ) MEAN, the mean of the Normal density.
!
!    Input/output, real ( kind = 8 ) SD, the standard deviation of the
!    Normal density.  If this is an input value, it should lie in the
!    range (0,+infinity).
!
!    Output, integer STATUS, the status of the calculation.
!    0, if calculation completed correctly;
!    -I, if input parameter number I is out of range;
!    1, if answer appears to be lower than lowest search bound;
!    2, if answer appears to be higher than greatest search bound;
!    3, if P + Q /= 1.
!
!    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
!    If STATUS is negative, then this is the value exceeded by parameter I.
!    if STATUS is 1 or 2, this is the search bound that was exceeded.
!
  implicit none

  real ( kind = 8 ) bound
!tp  real ( kind = 8 ) dinvnr
  real ( kind = 8 ) mean
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) sd
  integer status
  integer which
  real ( kind = 8 ) x
  real ( kind = 8 ) z
!
!  Check the arguments.
!
  status = 0

  if ( which < 1 ) then
    status = -1
    bound = 1.0D+00
    LOG_WARN("cdfnor",*) 'CDFNOR - Fatal error!'
    LOG_WARN_CONT(*) 'Input parameter WHICH is out of range.'
    return
  else if ( 4 < which ) then
    status = -1
    bound = 4.0D+00
    if ( debug ) then
       LOG_WARN("cdfnor",*) 'CDFNOR - Fatal error!'
       LOG_WARN_CONT(*) 'Input parameter WHICH is out of range.'
    end if
    return
  end if
!
!  Unless P is to be computed, make sure it is legal.
!
  if ( which /= 1 ) then
    if ( p < 0.0D+00 ) then
      status = -2
      bound = 0.0D+00
      LOG_WARN("cdfnor",*) 'CDFNOR - Fatal error!'
      LOG_WARN_CONT(*) 'Input parameter P is out of range.'
      return
    else if ( 1.0D+00 < p ) then
      status = -2
      bound = 1.0D+00
      LOG_WARN("cdfnor",*) 'CDFNOR - Fatal error!'
      LOG_WARN_CONT(*) 'Input parameter P is out of range.'
      return
    end if
  end if
!
!  Unless Q is to be computed, make sure it is legal.
!
  if ( which /= 1 ) then
    if ( q < 0.0D+00 ) then
      status = -3
      bound = 0.0D+00
      LOG_WARN("cdfnor",*) 'CDFNOR - Fatal error!'
      LOG_WARN_CONT(*) 'Input parameter Q is out of range.'
      return
    else if ( 1.0D+00 < q ) then
      status = -3
      bound = 1.0D+00
      LOG_WARN("cdfnor",*) 'CDFNOR - Fatal error!'
      LOG_WARN_CONT(*) 'Input parameter Q is out of range.'
      return
    end if
  end if
!
!  Check that P + Q = 1.
!
  if ( which /= 1 ) then
    if ( 3.0D+00 * epsilon ( p ) < abs ( ( p + q ) - 1.0D+00 ) ) then
      status = 3
      LOG_WARN("cdfnor",*) 'CDFNOR - Fatal error!'
      LOG_WARN_CONT(*) 'P + Q /= 1.'
      return
    end if
  end if

  if ( which /= 4 ) then
    if ( sd <= 0.0D+00 ) then
      bound = 0.0D+00
      status = -6
      LOG_WARN("cdfnor",*) 'CDFNOR - Fatal error!'
      LOG_WARN_CONT(*) 'Input parameter SD is out of range.'
      return
    end if
  end if
!
!  Calculate P and Q.
!
  if ( which == 1 ) then

    z = ( x - mean ) / sd
    call cumnor ( z, p, q )
!
!  Calculate X.
!
  else if ( which == 2 ) then

    z = dinvnr ( p, q )
    x = sd * z + mean
!
!  Calculate MEAN.
!
  else if ( which == 3 ) then

    z = dinvnr ( p, q )
    mean = x - sd * z
!
!  Calculate SD.
!
  else if ( which == 4 ) then

    z = dinvnr ( p, q )
    sd = ( x - mean ) / z

  end if

  return
end subroutine cdfnor
subroutine cumnor ( arg, cum, ccum )

!*******************************************************************************
!
!! CUMNOR computes the cumulative normal distribution.
!
!  Discussion:
!
!    This function evaluates the normal distribution function:
!
!                              / x
!                     1       |       -t*t/2
!          P(x) = ----------- |      e       dt
!                 sqrt(2 pi)  |
!                             /-oo
!
!    This transportable program uses rational functions that
!    theoretically approximate the normal distribution function to
!    at least 18 significant decimal digits.  The accuracy achieved
!    depends on the arithmetic system, the compiler, the intrinsic
!    functions, and proper selection of the machine dependent
!    constants.
!
!  Author:
!
!    W J Cody
!    Mathematics and Computer Science Division
!    Argonne National Laboratory
!    Argonne, IL 60439
!
!  Reference:
!
!    W J Cody,
!    Rational Chebyshev approximations for the error function,
!    Mathematics of Computation,
!    1969, pages 631-637.
!
!    W J Cody,
!    Algorithm 715:
!    SPECFUN - A Portable FORTRAN Package of Special Function Routines
!      and Test Drivers,
!    ACM Transactions on Mathematical Software,
!    Volume 19, Number 1, 1993, pages 22-32.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the upper limit of integration.
!
!    Output, real ( kind = 8 ) CUM, CCUM, the Normal density CDF and
!    complementary CDF.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) EPS, the argument below which anorm(x)
!    may be represented by 0.5 and above which  x*x  will not underflow.
!    A conservative value is the largest machine number X
!    such that   1.0D+00 + X = 1.0D+00   to machine precision.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 5 ) :: a = (/ &
    2.2352520354606839287D+00, &
    1.6102823106855587881D+02, &
    1.0676894854603709582D+03, &
    1.8154981253343561249D+04, &
    6.5682337918207449113D-02 /)
  real ( kind = 8 ) arg
  real ( kind = 8 ), parameter, dimension ( 4 ) :: b = (/ &
    4.7202581904688241870D+01, &
    9.7609855173777669322D+02, &
    1.0260932208618978205D+04, &
    4.5507789335026729956D+04 /)
  real ( kind = 8 ), parameter, dimension ( 9 ) :: c = (/ &
    3.9894151208813466764D-01, &
    8.8831497943883759412D+00, &
    9.3506656132177855979D+01, &
    5.9727027639480026226D+02, &
    2.4945375852903726711D+03, &
    6.8481904505362823326D+03, &
    1.1602651437647350124D+04, &
    9.8427148383839780218D+03, &
    1.0765576773720192317D-08 /)
  real ( kind = 8 ) ccum
  real ( kind = 8 ) cum
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    2.2266688044328115691D+01, &
    2.3538790178262499861D+02, &
    1.5193775994075548050D+03, &
    6.4855582982667607550D+03, &
    1.8615571640885098091D+04, &
    3.4900952721145977266D+04, &
    3.8912003286093271411D+04, &
    1.9685429676859990727D+04 /)
  real ( kind = 8 ) del
  real ( kind = 8 ) eps
  integer i
  real ( kind = 8 ), parameter, dimension ( 6 ) :: p = (/ &
    2.1589853405795699D-01, &
    1.274011611602473639D-01, &
    2.2235277870649807D-02, &
    1.421619193227893466D-03, &
    2.9112874951168792D-05, &
    2.307344176494017303D-02 /)
  real ( kind = 8 ), parameter, dimension ( 5 ) :: q = (/ &
    1.28426009614491121D+00, &
    4.68238212480865118D-01, &
    6.59881378689285515D-02, &
    3.78239633202758244D-03, &
    7.29751555083966205D-05 /)
  real ( kind = 8 ), parameter :: root32 = 5.656854248D+00
  real ( kind = 8 ), parameter :: sixten = 16.0D+00
  real ( kind = 8 ) temp
  real ( kind = 8 ), parameter :: sqrpi = 3.9894228040143267794D-01
  real ( kind = 8 ), parameter :: thrsh = 0.66291D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xden
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) xsq
!
!  Machine dependent constants
!
  eps = epsilon ( 1.0D+00 ) * 0.5D+00

  x = arg
  y = abs ( x )

  if ( y <= thrsh ) then
!
!  Evaluate  anorm  for  |X| <= 0.66291
!
    if ( eps < y ) then
      xsq = x * x
    else
      xsq = 0.0D+00
    end if

    xnum = a(5) * xsq
    xden = xsq
    do i = 1, 3
      xnum = ( xnum + a(i) ) * xsq
      xden = ( xden + b(i) ) * xsq
    end do
    cum = x * ( xnum + a(4) ) / ( xden + b(4) )
    temp = cum
    cum = 0.5D+00 + temp
    ccum = 0.5D+00 - temp
!
!  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
!
  else if ( y <= root32 ) then

    xnum = c(9) * y
    xden = y
    do i = 1, 7
      xnum = ( xnum + c(i) ) * y
      xden = ( xden + d(i) ) * y
    end do
    cum = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = aint ( y * sixten ) / sixten
    del = ( y - xsq ) * ( y + xsq )
    cum = exp ( - xsq * xsq * 0.5D+00 ) * exp ( -del * 0.5D+00 ) * cum
    ccum = 1.0D+00 - cum

    if ( 0.0D+00 < x ) then
      call d_swap ( cum, ccum )
    end if
!
!  Evaluate ANORM for sqrt(32) < |X|.
!
  else

    cum = 0.0D+00
    xsq = 1.0D+00 / ( x * x )
    xnum = p(6) * xsq
    xden = xsq
    do i = 1, 4
      xnum = ( xnum + p(i) ) * xsq
      xden = ( xden + q(i) ) * xsq
    end do

    cum = xsq * ( xnum + p(5) ) / ( xden + q(5) )
    cum = ( sqrpi - cum ) / y
    xsq = aint ( x * sixten ) / sixten
    del = ( x - xsq ) * ( x + xsq )
    cum = exp ( - xsq * xsq * 0.5D+00 ) &
      * exp ( - del * 0.5D+00 ) * cum
    ccum = 1.0D+00 - cum

    if ( 0.0D+00 < x ) then
      call d_swap ( cum, ccum )
    end if

  end if

  if ( cum < tiny ( cum ) ) then
    cum = 0.0D+00
  end if

  if ( ccum < tiny ( ccum ) ) then
    ccum = 0.0D+00
  end if

  return
end subroutine cumnor
subroutine d_swap ( x, y )

!*******************************************************************************
!
!! D_SWAP swaps two real ( kind = 8 ) values.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end subroutine d_swap
function dinvnr ( p, q )

!*******************************************************************************
!
!! DINVNR computes the inverse of the normal distribution.
!
!  Discussion:
!
!    This routine returns X such that
!
!      CUMNOR(X) = P,
!
!    that is, so that
!
!      P = integral ( -Infinity <= T <= X ) exp(-U*U/2)/sqrt(2*PI) dU
!
!    The rational function on page 95 of Kennedy and Gentle is used as a
!    starting value for the Newton method of finding roots.
!
!  Reference:
!
!    Kennedy and Gentle,
!    Statistical Computing,
!    Marcel Dekker, NY, 1980,
!    QA276.4 K46
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, Q, the probability, and the complementary
!    probability.
!
!    Output, real ( kind = 8 ) DINVNR, the argument X for which the
!    Normal CDF has the value P.
!
  implicit none

  real ( kind = 8 ) ccum
  real ( kind = 8 ) cum
  real ( kind = 8 ) dinvnr
  real ( kind = 8 ) dx
  real ( kind = 8 ), parameter :: eps = 1.0D-13
  integer i
  integer, parameter :: maxit = 100
  real ( kind = 8 ) p
  real ( kind = 8 ) pp
  real ( kind = 8 ) q
  real ( kind = 8 ), parameter :: r2pi = 0.3989422804014326D+00
  real ( kind = 8 ) strtx
!tmp  real ( kind = 8 ) stvaln
  real ( kind = 8 ) xcur

  pp = min ( p, q )
  strtx = stvaln ( pp )
  xcur = strtx
!
!  Newton iterations.
!
  do i = 1, maxit

    call cumnor ( xcur, cum, ccum )
    dx = ( cum - pp ) / ( r2pi * exp ( -0.5D+00 * xcur * xcur ) )
    xcur = xcur - dx

    if ( abs ( dx / xcur ) < eps ) then
      if ( p <= q ) then
        dinvnr = xcur
      else
        dinvnr = -xcur
      end if
      return
    end if

  end do

  if ( p <= q ) then
    dinvnr = strtx
  else
    dinvnr = -strtx
  end if

  return
end function dinvnr

!!$function gam1 ( a )
!!$
!!$!*******************************************************************************
!!$!
!!$!! GAM1 computes 1 / GAMMA(A+1) - 1 for -0.5 <= A <= 1.5
!!$!
!!$!  Reference:
!!$!
!!$!    A R DiDinato and A H Morris,
!!$!    Algorithm 708:
!!$!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!!$!    ACM Transactions on Mathematical Software,
!!$!    Volume 18, 1993, pages 360-373.
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real ( kind = 8 ) A, forms the argument of the Gamma function.
!!$!
!!$!    Output, real ( kind = 8 ) GAM1, the value of 1 / GAMMA ( A + 1 ) - 1.
!!$!
!!$  implicit none
!!$
!!$  real ( kind = 8 ) a
!!$  real ( kind = 8 ) bot
!!$  real ( kind = 8 ) d
!!$  real ( kind = 8 ) gam1
!!$  real ( kind = 8 ), parameter, dimension ( 7 ) :: p = (/ &
!!$     0.577215664901533D+00, -0.409078193005776D+00, &
!!$    -0.230975380857675D+00,  0.597275330452234D-01, &
!!$     0.766968181649490D-02, -0.514889771323592D-02, &
!!$     0.589597428611429D-03 /)
!!$  real ( kind = 8 ), dimension ( 5 ) :: q = (/ &
!!$    0.100000000000000D+01, 0.427569613095214D+00, &
!!$    0.158451672430138D+00, 0.261132021441447D-01, &
!!$    0.423244297896961D-02 /)
!!$  real ( kind = 8 ), dimension ( 9 ) :: r = (/ &
!!$    -0.422784335098468D+00, -0.771330383816272D+00, &
!!$    -0.244757765222226D+00,  0.118378989872749D+00, &
!!$     0.930357293360349D-03, -0.118290993445146D-01, &
!!$     0.223047661158249D-02,  0.266505979058923D-03, &
!!$    -0.132674909766242D-03 /)
!!$  real ( kind = 8 ), parameter :: s1 = 0.273076135303957D+00
!!$  real ( kind = 8 ), parameter :: s2 = 0.559398236957378D-01
!!$  real ( kind = 8 ) t
!!$  real ( kind = 8 ) top
!!$  real ( kind = 8 ) w
!!$
!!$  d = a - 0.5D+00
!!$
!!$  if ( 0.0D+00 < d ) then
!!$    t = d - 0.5D+00
!!$  else
!!$    t = a
!!$  end if
!!$
!!$  if ( t == 0.0D+00 ) then
!!$
!!$    gam1 = 0.0D+00
!!$
!!$  else if ( 0.0D+00 < t ) then
!!$
!!$    top = (((((    &
!!$            p(7)   &
!!$      * t + p(6) ) &
!!$      * t + p(5) ) &
!!$      * t + p(4) ) &
!!$      * t + p(3) ) &
!!$      * t + p(2) ) &
!!$      * t + p(1)
!!$
!!$    bot = ((( q(5) * t + q(4) ) * t + q(3) ) * t + q(2) ) * t &
!!$      + 1.0D+00
!!$
!!$    w = top / bot
!!$
!!$    if ( d <= 0.0D+00 ) then
!!$      gam1 = a * w
!!$    else
!!$      gam1 = ( t / a ) * ( ( w - 0.5D+00 ) &
!!$        - 0.5D+00 )
!!$    end if
!!$
!!$  else if ( t < 0.0D+00 ) then
!!$
!!$    top = (((((((  &
!!$            r(9)   &
!!$      * t + r(8) ) &
!!$      * t + r(7) ) &
!!$      * t + r(6) ) &
!!$      * t + r(5) ) &
!!$      * t + r(4) ) &
!!$      * t + r(3) ) &
!!$      * t + r(2) ) &
!!$      * t + r(1)
!!$
!!$    bot = ( s2 * t + s1 ) * t + 1.0D+00
!!$    w = top / bot
!!$
!!$    if ( d <= 0.0D+00 ) then
!!$      gam1 = a * ( ( w + 0.5D+00 ) + 0.5D+00 )
!!$    else
!!$      gam1 = t * w / a
!!$    end if
!!$
!!$  end if
!!$
!!$  return
!!$end function gam1

!!$function gamma_1 ( a )
!!$
!!$!*******************************************************************************
!!$!
!!$!! GAMMA evaluates the gamma function.
!!$!
!!$!  Author:
!!$!
!!$!    Alfred H Morris, Jr,
!!$!    Naval Surface Weapons Center,
!!$!    Dahlgren, Virginia.
!!$!
!!$!  Reference:
!!$!
!!$!    A R DiDinato and A H Morris,
!!$!    Algorithm 708:
!!$!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!!$!    ACM Transactions on Mathematical Software,
!!$!    Volume 18, 1993, pages 360-373.
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real ( kind = 8 ) A, the argument of the Gamma function.
!!$!
!!$!    Output, real ( kind = 8 ) GAMMA, the value of the Gamma function.
!!$!
!!$  implicit none
!!$
!!$  real ( kind = 8 ) a
!!$  real ( kind = 8 ) bot
!!$  real ( kind = 8 ), parameter :: d = 0.41893853320467274178D+00
!!$!tmp  real ( kind = 8 ) exparg
!!$  real ( kind = 8 ) g
!!$  real ( kind = 8 ) gamma_1
!!$  integer i
!!$  integer j
!!$  real ( kind = 8 ) lnx
!!$  integer m
!!$  integer n
!!$  real ( kind = 8 ), dimension ( 7 ) :: p = (/ &
!!$    0.539637273585445D-03, 0.261939260042690D-02, &
!!$    0.204493667594920D-01, 0.730981088720487D-01, &
!!$    0.279648642639792D+00, 0.553413866010467D+00, &
!!$    1.0D+00 /)
!!$  real ( kind = 8 ), parameter :: pi = 3.1415926535898D+00
!!$  real ( kind = 8 ), dimension ( 7 ) :: q = (/ &
!!$    -0.832979206704073D-03,  0.470059485860584D-02, &
!!$     0.225211131035340D-01, -0.170458969313360D+00, &
!!$    -0.567902761974940D-01,  0.113062953091122D+01, &
!!$     1.0D+00 /)
!!$  real ( kind = 8 ), parameter :: r1 =  0.820756370353826D-03
!!$  real ( kind = 8 ), parameter :: r2 = -0.595156336428591D-03
!!$  real ( kind = 8 ), parameter :: r3 =  0.793650663183693D-03
!!$  real ( kind = 8 ), parameter :: r4 = -0.277777777770481D-02
!!$  real ( kind = 8 ), parameter :: r5 =  0.833333333333333D-01
!!$  real ( kind = 8 ) s
!!$  real ( kind = 8 ) t
!!$  real ( kind = 8 ) top
!!$  real ( kind = 8 ) w
!!$  real ( kind = 8 ) x
!!$  real ( kind = 8 ) z
!!$
!!$  gamma_1 = 0.0D+00
!!$  x = a
!!$
!!$  if ( abs ( a ) < 15.0D+00 ) then
!!$!
!!$!  Evaluation of GAMMA(A) for |A| < 15
!!$!
!!$    t = 1.0D+00
!!$    m = int ( a ) - 1
!!$!
!!$!  Let T be the product of A-J when 2 <= A.
!!$!
!!$    if ( 0 <= m ) then
!!$
!!$      do j = 1, m
!!$        x = x - 1.0D+00
!!$        t = x * t
!!$      end do
!!$
!!$      x = x - 1.0D+00
!!$!
!!$!  Let T be the product of A+J WHEN A < 1
!!$!
!!$    else
!!$
!!$      t = a
!!$
!!$      if ( a <= 0.0D+00 ) then
!!$
!!$        m = - m - 1
!!$
!!$        do j = 1, m
!!$          x = x + 1.0D+00
!!$          t = x * t
!!$        end do
!!$
!!$        x = ( x + 0.5D+00 ) + 0.5D+00
!!$        t = x * t
!!$        if ( t == 0.0D+00 ) then
!!$          return
!!$        end if
!!$
!!$      end if
!!$!
!!$!  Check if 1/T can overflow.
!!$!
!!$      if ( abs ( t ) < 1.0D-30 ) then
!!$        if ( 1.0001D+00 < abs ( t ) * huge ( t ) ) then
!!$          gamma_1 = 1.0D+00 / t
!!$        end if
!!$        return
!!$      end if
!!$
!!$    end if
!!$!
!!$!  Compute Gamma(1 + X) for 0 <= X < 1.
!!$!
!!$    top = p(1)
!!$    bot = q(1)
!!$    do i = 2, 7
!!$      top = top * x + p(i)
!!$      bot = bot * x + q(i)
!!$    end do
!!$
!!$    gamma_1 = top / bot
!!$!
!!$!  Termination.
!!$!
!!$    if ( 1.0D+00 <= a ) then
!!$      gamma_1 = gamma_1 * t
!!$    else
!!$      gamma_1 = gamma_1 / t
!!$    end if
!!$!
!!$!  Evaluation of Gamma(A) FOR 15 <= ABS ( A ).
!!$!
!!$  else
!!$
!!$    if ( 1000.0D+00 <= abs ( a ) ) then
!!$      return
!!$    end if
!!$
!!$    if ( a <= 0.0D+00 ) then
!!$
!!$      x = -a
!!$      n = x
!!$      t = x - n
!!$
!!$      if ( 0.9D+00 < t ) then
!!$        t = 1.0D+00 - t
!!$      end if
!!$
!!$      s = sin ( pi * t ) / pi
!!$
!!$      if ( mod ( n, 2 ) == 0 ) then
!!$        s = -s
!!$      end if
!!$
!!$      if ( s == 0.0D+00 ) then
!!$        return
!!$      end if
!!$
!!$    end if
!!$!
!!$!  Compute the modified asymptotic sum.
!!$!
!!$    t = 1.0D+00 / ( x * x )
!!$
!!$    g = (((( r1 * t + r2 ) * t + r3 ) * t + r4 ) * t + r5 ) / x
!!$
!!$    lnx = log ( x )
!!$!
!!$!  Final assembly.
!!$!
!!$    z = x
!!$    g = ( d + g ) + ( z - 0.5D+00 ) &
!!$      * ( lnx - 1.0D+00 )
!!$    w = g
!!$    t = g - real ( w, kind = 8 )
!!$
!!$    if ( 0.99999D+00 * exparg ( 0 ) < w ) then
!!$      return
!!$    end if
!!$
!!$    gamma_1 = exp ( w )* ( 1.0D+00 + t )
!!$
!!$    if ( a < 0.0D+00 ) then
!!$      gamma_1 = ( 1.0D+00 / ( gamma_1 * s ) ) / x
!!$    end if
!!$
!!$  end if
!!$
!!$  return
!!$end function gamma_1

!!$subroutine cumgam ( x, a, cum, ccum )
!!$
!!$!*******************************************************************************
!!$!
!!$!! CUMGAM evaluates the cumulative incomplete gamma distribution.
!!$!
!!$!  Discussion:
!!$!
!!$!    This routine computes the cumulative distribution function of the
!!$!    incomplete gamma distribution, i.e., the integral from 0 to X of
!!$!
!!$!      (1/GAM(A))*EXP(-T)*T**(A-1) DT
!!$!
!!$!    where GAM(A) is the complete gamma function of A, i.e.,
!!$!
!!$!      GAM(A) = integral from 0 to infinity of EXP(-T)*T**(A-1) DT
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real ( kind = 8 ) X, the upper limit of integration.
!!$!
!!$!    Input, real ( kind = 8 ) A, the shape parameter of the incomplete
!!$!    Gamma distribution.
!!$!
!!$!    Output, real ( kind = 8 ) CUM, CCUM, the incomplete Gamma CDF and
!!$!    complementary CDF.
!!$!
!!$  implicit none
!!$
!!$  real ( kind = 8 ) a
!!$  real ( kind = 8 ) ccum
!!$  real ( kind = 8 ) cum
!!$  real ( kind = 8 ) x
!!$
!!$  if ( x <= 0.0D+00 ) then
!!$
!!$    cum = 0.0D+00
!!$    ccum = 1.0D+00
!!$
!!$  else
!!$
!!$    call gamma_inc ( a, x, cum, ccum, 0 )
!!$
!!$  end if
!!$
!!$  return
!!$end subroutine cumgam

!!$subroutine gamma_inc ( a, x, ans, qans, ind )
!!$
!!$!*******************************************************************************
!!$!
!!$!! GAMMA_INC evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
!!$!
!!$!  Discussion:
!!$!
!!$!    This is certified spaghetti code.
!!$!
!!$!  Modified:
!!$!
!!$!    16 April 2005
!!$!
!!$!  Author:
!!$!
!!$!    Alfred H Morris, Jr,
!!$!    Naval Surface Weapons Center,
!!$!    Dahlgren, Virginia.
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real ( kind = 8 ) A, X, the arguments of the incomplete
!!$!    gamma ratio.  A and X must be nonnegative.  A and X cannot
!!$!    both be zero.
!!$!
!!$!    Output, real ( kind = 8 ) ANS, QANS.  On normal output,
!!$!    ANS = P(A,X) and QANS = Q(A,X).  However, ANS is set to 2 if
!!$!    A or X is negative, or both are 0, or when the answer is
!!$!    computationally indeterminate because A is extremely large
!!$!    and X is very close to A.
!!$!
!!$!    Input, integer IND, indicates the accuracy request:
!!$!    0, as much accuracy as possible.
!!$!    1, to within 1 unit of the 6-th significant digit,
!!$!    otherwise, to within 1 unit of the 3rd significant digit.
!!$!
!!$!  Local Parameters:
!!$!
!!$!     ALOG10 = LN(10)
!!$!     RT2PIN = 1/SQRT(2*PI)
!!$!     RTPI   = SQRT(PI)
!!$!
!!$  implicit none
!!$
!!$  real ( kind = 8 ) a
!!$  real ( kind = 8 ) a2n
!!$  real ( kind = 8 ) a2nm1
!!$  real ( kind = 8 ) acc
!!$  real ( kind = 8 ), dimension ( 3 ) :: acc0 = (/ &
!!$    5.0D-15, 5.0D-07, 5.0D-04 /)
!!$  real ( kind = 8 ), parameter :: alog10 = 2.30258509299405D+00
!!$  real ( kind = 8 ) am0
!!$  real ( kind = 8 ) amn
!!$  real ( kind = 8 ) an
!!$  real ( kind = 8 ) an0
!!$  real ( kind = 8 ) ans
!!$  real ( kind = 8 ) apn
!!$  real ( kind = 8 ) b2n
!!$  real ( kind = 8 ) b2nm1
!!$  real ( kind = 8 ) big(3)
!!$  real ( kind = 8 ) c
!!$  real ( kind = 8 ) c0
!!$  real ( kind = 8 ) c1
!!$  real ( kind = 8 ) c2
!!$  real ( kind = 8 ) c3
!!$  real ( kind = 8 ) c4
!!$  real ( kind = 8 ) c5
!!$  real ( kind = 8 ) c6
!!$  real ( kind = 8 ) cma
!!$  real ( kind = 8 ) d0(13)
!!$  real ( kind = 8 ) d1(12)
!!$  real ( kind = 8 ) d2(10)
!!$  real ( kind = 8 ) d3(8)
!!$  real ( kind = 8 ) d4(6)
!!$  real ( kind = 8 ) d5(4)
!!$  real ( kind = 8 ) d6(2)
!!$  real ( kind = 8 ) d10
!!$  real ( kind = 8 ) d20
!!$  real ( kind = 8 ) d30
!!$  real ( kind = 8 ) d40
!!$  real ( kind = 8 ) d50
!!$  real ( kind = 8 ) d60
!!$  real ( kind = 8 ) d70
!!$  real ( kind = 8 ) e
!!$  real ( kind = 8 ) e0
!!$  real ( kind = 8 ) e00(3)
!!$  real ( kind = 8 ) erf
!!$!tmp  real ( kind = 8 ), external :: erfc1
!!$  real ( kind = 8 ) g
!!$!tmp  real ( kind = 8 ), external :: gam1
!!$!tmp  real ( kind = 8 ), external :: gamma_1
!!$  real ( kind = 8 ) h
!!$  integer i
!!$  integer ind
!!$  integer iop
!!$  real ( kind = 8 ) j
!!$  real ( kind = 8 ) l
!!$  integer m
!!$  integer n
!!$  integer n_max
!!$  real ( kind = 8 ) qans
!!$  real ( kind = 8 ) r
!!$!tmp  real ( kind = 8 ) rexp
!!$!tmp  real ( kind = 8 ) rlog
!!$  real ( kind = 8 ), parameter :: rt2pin = 0.398942280401433D+00
!!$  real ( kind = 8 ) rta
!!$  real ( kind = 8 ), parameter :: rtpi = 1.77245385090552D+00
!!$  real ( kind = 8 ) rtx
!!$  real ( kind = 8 ) s
!!$  real ( kind = 8 ) sum1
!!$  real ( kind = 8 ) t
!!$  real ( kind = 8 ) t1
!!$  real ( kind = 8 ) tol
!!$  real ( kind = 8 ) twoa
!!$  real ( kind = 8 ) u
!!$  real ( kind = 8 ) w
!!$  real ( kind = 8 ) wk(20)
!!$  real ( kind = 8 ) x
!!$  real ( kind = 8 ) x0
!!$  real ( kind = 8 ) x00(3)
!!$  real ( kind = 8 ) y
!!$  real ( kind = 8 ) z
!!$
!!$  data big(1)/20.0D+00/,big(2)/14.0D+00/,big(3)/10.0D+00/
!!$  data e00(1)/0.25D-03/,e00(2)/0.25D-01/,e00(3)/0.14D+00/
!!$  data x00(1)/31.0D+00/,x00(2)/17.0D+00/,x00(3)/9.7D+00/
!!$  data d0(1)/0.833333333333333D-01/
!!$  data d0(2)/-0.148148148148148D-01/
!!$  data d0(3)/0.115740740740741D-02/,d0(4)/0.352733686067019D-03/
!!$  data d0(5)/-0.178755144032922D-03/,d0(6)/0.391926317852244D-04/
!!$  data d0(7)/-0.218544851067999D-05/,d0(8)/-0.185406221071516D-05/
!!$  data d0(9)/0.829671134095309D-06/,d0(10)/-0.176659527368261D-06/
!!$  data d0(11)/0.670785354340150D-08/,d0(12)/0.102618097842403D-07/
!!$  data d0(13)/-0.438203601845335D-08/
!!$  data d10/-0.185185185185185D-02/,d1(1)/-0.347222222222222D-02/
!!$  data d1(2)/0.264550264550265D-02/,d1(3)/-0.990226337448560D-03/
!!$  data d1(4)/0.205761316872428D-03/,d1(5)/-0.401877572016461D-06/
!!$  data d1(6)/-0.180985503344900D-04/,d1(7)/0.764916091608111D-05/
!!$  data d1(8)/-0.161209008945634D-05/,d1(9)/0.464712780280743D-08/
!!$  data d1(10)/0.137863344691572D-06/,d1(11)/-0.575254560351770D-07/
!!$  data d1(12)/0.119516285997781D-07/
!!$  data d20/0.413359788359788D-02/,d2(1)/-0.268132716049383D-02/
!!$  data d2(2)/0.771604938271605D-03/,d2(3)/0.200938786008230D-05/
!!$  data d2(4)/-0.107366532263652D-03/,d2(5)/0.529234488291201D-04/
!!$  data d2(6)/-0.127606351886187D-04/,d2(7)/0.342357873409614D-07/
!!$  data d2(8)/0.137219573090629D-05/,d2(9)/-0.629899213838006D-06/
!!$  data d2(10)/0.142806142060642D-06/
!!$  data d30/0.649434156378601D-03/,d3(1)/0.229472093621399D-03/
!!$  data d3(2)/-0.469189494395256D-03/,d3(3)/0.267720632062839D-03/
!!$  data d3(4)/-0.756180167188398D-04/,d3(5)/-0.239650511386730D-06/
!!$  data d3(6)/0.110826541153473D-04/,d3(7)/-0.567495282699160D-05/
!!$  data d3(8)/0.142309007324359D-05/
!!$  data d40/-0.861888290916712D-03/,d4(1)/0.784039221720067D-03/
!!$  data d4(2)/-0.299072480303190D-03/,d4(3)/-0.146384525788434D-05/
!!$  data d4(4)/0.664149821546512D-04/,d4(5)/-0.396836504717943D-04/
!!$  data d4(6)/0.113757269706784D-04/
!!$  data d50/-0.336798553366358D-03/,d5(1)/-0.697281375836586D-04/
!!$  data d5(2)/0.277275324495939D-03/,d5(3)/-0.199325705161888D-03/
!!$  data d5(4)/0.679778047793721D-04/
!!$  data d60/0.531307936463992D-03/,d6(1)/-0.592166437353694D-03/
!!$  data d6(2)/0.270878209671804D-03/
!!$  data d70 / 0.344367606892378D-03/
!!$
!!$  e = epsilon ( 1.0D+00 )
!!$
!!$  if ( a < 0.0D+00 .or. x < 0.0D+00 ) then
!!$    ans = 2.0D+00
!!$    return
!!$  end if
!!$
!!$  if ( a == 0.0D+00 .and. x == 0.0D+00 ) then
!!$    ans = 2.0D+00
!!$    return
!!$  end if
!!$
!!$  if ( a * x == 0.0D+00 ) then
!!$    if ( x <= a ) then
!!$      ans = 0.0D+00
!!$      qans = 1.0D+00
!!$    else
!!$      ans = 1.0D+00
!!$      qans = 0.0D+00
!!$    end if
!!$    return
!!$  end if
!!$
!!$  iop = ind + 1
!!$  if ( iop /= 1 .and. iop /= 2 ) iop = 3
!!$  acc = max ( acc0(iop), e )
!!$  e0 = e00(iop)
!!$  x0 = x00(iop)
!!$!
!!$!  Select the appropriate algorithm.
!!$!
!!$  if ( 1.0D+00 <= a ) then
!!$    go to 10
!!$  end if
!!$
!!$  if ( a == 0.5D+00 ) then
!!$    go to 390
!!$  end if
!!$
!!$  if ( x < 1.1D+00 ) then
!!$    go to 160
!!$  end if
!!$
!!$  t1 = a * log ( x ) - x
!!$  u = a * exp ( t1 )
!!$
!!$  if ( u == 0.0D+00 ) then
!!$    ans = 1.0D+00
!!$    qans = 0.0D+00
!!$    return
!!$  end if
!!$
!!$  r = u * ( 1.0D+00 + gam1 ( a ) )
!!$  go to 250
!!$
!!$   10 continue
!!$
!!$  if ( big(iop) <= a ) then
!!$    go to 30
!!$  end if
!!$
!!$  if ( x < a .or. x0 <= x ) then
!!$    go to 20
!!$  end if
!!$
!!$  twoa = a + a
!!$  m = int ( twoa )
!!$
!!$  if ( twoa == real ( m, kind = 8 ) ) then
!!$    i = m / 2
!!$    if ( a == real ( i, kind = 8 ) ) then
!!$      go to 210
!!$    end if
!!$    go to 220
!!$  end if
!!$
!!$   20 continue
!!$
!!$  t1 = a * log ( x ) - x
!!$  r = exp ( t1 ) / gamma_1 ( a )
!!$  go to 40
!!$
!!$   30 continue
!!$
!!$  l = x / a
!!$
!!$  if ( l == 0.0D+00 ) then
!!$    ans = 0.0D+00
!!$    qans = 1.0D+00
!!$    return
!!$  end if
!!$
!!$  s = 0.5D+00 + ( 0.5D+00 - l )
!!$  z = rlog ( l )
!!$  if ( 700.0D+00 / a <= z ) then
!!$    go to 410
!!$  end if
!!$
!!$  y = a * z
!!$  rta = sqrt ( a )
!!$
!!$  if ( abs ( s ) <= e0 / rta ) then
!!$    go to 330
!!$  end if
!!$
!!$  if ( abs ( s ) <= 0.4D+00 ) then
!!$    go to 270
!!$  end if
!!$
!!$  t = ( 1.0D+00 / a )**2
!!$  t1 = ((( 0.75D+00 * t - 1.0D+00 ) * t + 3.5D+00 ) &
!!$    * t - 105.0D+00 ) / ( a * 1260.0D+00 )
!!$  t1 = t1 - y
!!$  r = rt2pin * rta * exp ( t1 )
!!$
!!$40    continue
!!$
!!$  if ( r == 0.0D+00 ) then
!!$    if ( x <= a ) then
!!$      ans = 0.0D+00
!!$      qans = 1.0D+00
!!$    else
!!$      ans = 1.0D+00
!!$      qans = 0.0D+00
!!$    end if
!!$    return
!!$  end if
!!$
!!$  if ( x <= max ( a, alog10 ) ) then
!!$    go to 50
!!$  end if
!!$
!!$  if ( x < x0 ) then
!!$    go to 250
!!$  end if
!!$
!!$  go to 100
!!$!
!!$!  Taylor series for P/R.
!!$!
!!$50    continue
!!$
!!$  apn = a + 1.0D+00
!!$  t = x / apn
!!$  wk(1) = t
!!$
!!$  n = 20
!!$
!!$  do i = 2, 20
!!$    apn = apn + 1.0D+00
!!$    t = t * ( x / apn )
!!$    if ( t <= 1.0D-03 ) then
!!$      n = i
!!$      exit
!!$    end if
!!$    wk(i) = t
!!$  end do
!!$
!!$  sum1 = t
!!$
!!$  tol = 0.5D+00 * acc
!!$
!!$  do
!!$
!!$    apn = apn + 1.0D+00
!!$    t = t * ( x / apn )
!!$    sum1 = sum1 + t
!!$
!!$    if ( t <= tol ) then
!!$      exit
!!$    end if
!!$
!!$  end do
!!$
!!$  n_max = n - 1
!!$  do m = 1, n_max
!!$    n = n - 1
!!$    sum1 = sum1 + wk(n)
!!$  end do
!!$
!!$  ans = ( r / a ) * ( 1.0D+00 + sum1 )
!!$  qans = 0.5D+00 + ( 0.5D+00 - ans )
!!$  return
!!$!
!!$!  Asymptotic expansion.
!!$!
!!$  100 continue
!!$
!!$  amn = a - 1.0D+00
!!$  t = amn / x
!!$  wk(1) = t
!!$
!!$  n = 20
!!$
!!$  do i = 2, 20
!!$    amn = amn - 1.0D+00
!!$    t = t * ( amn / x )
!!$    if ( abs ( t ) <= 1.0D-03 ) then
!!$      n = i
!!$      exit
!!$    end if
!!$    wk(i) = t
!!$  end do
!!$
!!$  sum1 = t
!!$
!!$  do
!!$
!!$    if ( abs ( t ) <= acc ) then
!!$      exit
!!$    end if
!!$
!!$    amn = amn - 1.0D+00
!!$    t = t * ( amn / x )
!!$    sum1 = sum1 + t
!!$
!!$  end do
!!$
!!$  n_max = n - 1
!!$  do m = 1, n_max
!!$    n = n - 1
!!$    sum1 = sum1 + wk(n)
!!$  end do
!!$  qans = ( r / x ) * ( 1.0D+00 + sum1 )
!!$  ans = 0.5D+00 + ( 0.5D+00 - qans )
!!$  return
!!$!
!!$!  Taylor series for P(A,X)/X**A
!!$!
!!$  160 continue
!!$
!!$  an = 3.0D+00
!!$  c = x
!!$  sum1 = x / ( a + 3.0D+00 )
!!$  tol = 3.0D+00 * acc / ( a + 1.0D+00 )
!!$
!!$  do
!!$
!!$    an = an + 1.0D+00
!!$    c = -c * ( x / an )
!!$    t = c / ( a + an )
!!$    sum1 = sum1 + t
!!$
!!$    if ( abs ( t ) <= tol ) then
!!$      exit
!!$    end if
!!$
!!$  end do
!!$
!!$  j = a * x * ( ( sum1 / 6.0D+00 - 0.5D+00 / &
!!$    ( a +  2.0D+00  ) ) * x + 1.0D+00 &
!!$    / ( a + 1.0D+00 ) )
!!$
!!$  z = a * log ( x )
!!$  h = gam1 ( a )
!!$  g = 1.0D+00 + h
!!$
!!$  if ( x < 0.25D+00 ) then
!!$    go to 180
!!$  end if
!!$
!!$  if ( a < x / 2.59D+00 ) then
!!$    go to 200
!!$  end if
!!$
!!$  go to 190
!!$
!!$  180 continue
!!$
!!$  if ( -0.13394D+00 < z ) then
!!$    go to 200
!!$  end if
!!$
!!$  190 continue
!!$
!!$  w = exp ( z )
!!$  ans = w * g * ( 0.5D+00 + ( 0.5D+00 - j ))
!!$  qans = 0.5D+00 + ( 0.5D+00 - ans )
!!$  return
!!$
!!$200   continue
!!$
!!$  l = rexp ( z )
!!$  w = 0.5D+00 + ( 0.5D+00 + l )
!!$  qans = ( w * j - l ) * g - h
!!$
!!$  if ( qans < 0.0D+00 ) then
!!$    ans = 1.0D+00
!!$    qans = 0.0D+00
!!$    return
!!$  end if
!!$
!!$  ans = 0.5D+00 + ( 0.5D+00 - qans )
!!$  return
!!$!
!!$!  Finite sums for Q when 1 <= A and 2*A is an integer.
!!$!
!!$210   continue
!!$
!!$  sum1 = exp ( - x )
!!$  t = sum1
!!$  n = 1
!!$  c = 0.0D+00
!!$  go to 230
!!$
!!$220   continue
!!$
!!$  rtx = sqrt ( x )
!!$  sum1 = erfc1 ( 0, rtx )
!!$  t = exp ( -x ) / ( rtpi * rtx )
!!$  n = 0
!!$  c = -0.5D+00
!!$
!!$  230 continue
!!$
!!$  do while ( n /= i )
!!$    n = n + 1
!!$    c = c + 1.0D+00
!!$    t = ( x * t ) / c
!!$    sum1 = sum1 + t
!!$  end do
!!$
!!$  240 continue
!!$
!!$  qans = sum1
!!$  ans = 0.5D+00 + ( 0.5D+00 - qans )
!!$  return
!!$!
!!$!  Continued fraction expansion.
!!$!
!!$250   continue
!!$
!!$  tol = max ( 5.0D+00 * e, acc )
!!$  a2nm1 = 1.0D+00
!!$  a2n = 1.0D+00
!!$  b2nm1 = x
!!$  b2n = x + ( 1.0D+00 - a )
!!$  c = 1.0D+00
!!$
!!$  do
!!$
!!$    a2nm1 = x * a2n + c * a2nm1
!!$    b2nm1 = x * b2n + c * b2nm1
!!$    am0 = a2nm1 / b2nm1
!!$    c = c + 1.0D+00
!!$    cma = c - a
!!$    a2n = a2nm1 + cma * a2n
!!$    b2n = b2nm1 + cma * b2n
!!$    an0 = a2n / b2n
!!$
!!$    if ( abs ( an0 - am0 ) < tol * an0 ) then
!!$      exit
!!$    end if
!!$
!!$  end do
!!$
!!$  qans = r * an0
!!$  ans = 0.5D+00 + ( 0.5D+00 - qans )
!!$  return
!!$!
!!$!  General Temme expansion.
!!$!
!!$270   continue
!!$
!!$  if ( abs ( s ) <= 2.0D+00 * e .and. 3.28D-03 < a * e * e ) then
!!$    ans =  2.0D+00
!!$    return
!!$  end if
!!$
!!$  c = exp ( - y )
!!$  w = 0.5D+00 * erfc1 ( 1, sqrt ( y ) )
!!$  u = 1.0D+00 / a
!!$  z = sqrt ( z + z )
!!$
!!$  if ( l < 1.0D+00 ) then
!!$    z = -z
!!$  end if
!!$
!!$  if ( iop < 2 ) then
!!$
!!$    if ( abs ( s ) <= 1.0D-03 ) then
!!$
!!$      c0 = ((((((     &
!!$              d0(7)   &
!!$        * z + d0(6) ) &
!!$        * z + d0(5) ) &
!!$        * z + d0(4) ) &
!!$        * z + d0(3) ) &
!!$        * z + d0(2) ) &
!!$        * z + d0(1) ) &
!!$        * z - 1.0D+00 / 3.0D+00
!!$
!!$      c1 = (((((      &
!!$              d1(6)   &
!!$        * z + d1(5) ) &
!!$        * z + d1(4) ) &
!!$        * z + d1(3) ) &
!!$        * z + d1(2) ) &
!!$        * z + d1(1) ) &
!!$        * z + d10
!!$
!!$      c2 = ((((d2(5)*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20
!!$
!!$      c3 = (((d3(4)*z+d3(3))*z+d3(2))*z+d3(1))*z + d30
!!$
!!$      c4 = ( d4(2) * z + d4(1) ) * z + d40
!!$      c5 = ( d5(2) * z + d5(1) ) * z + d50
!!$      c6 = d6(1) * z + d60
!!$
!!$      t = (((((( d70 &
!!$        * u + c6 ) &
!!$        * u + c5 ) &
!!$        * u + c4 ) &
!!$        * u + c3 ) &
!!$        * u + c2 ) &
!!$        * u + c1 ) &
!!$        * u + c0
!!$
!!$    else
!!$
!!$      c0 = (((((((((((( &
!!$              d0(13)   &
!!$        * z + d0(12) ) &
!!$        * z + d0(11) ) &
!!$        * z + d0(10) ) &
!!$        * z + d0(9)  ) &
!!$        * z + d0(8)  ) &
!!$        * z + d0(7)  ) &
!!$        * z + d0(6)  ) &
!!$        * z + d0(5)  ) &
!!$        * z + d0(4)  ) &
!!$        * z + d0(3)  ) &
!!$        * z + d0(2)  ) &
!!$        * z + d0(1)  ) &
!!$        * z - 1.0D+00 / 3.0D+00
!!$
!!$      c1 = ((((((((((( &
!!$                d1(12) &
!!$          * z + d1(11) &
!!$        ) * z + d1(10) &
!!$        ) * z + d1(9)  &
!!$        ) * z + d1(8)  &
!!$        ) * z + d1(7)  &
!!$        ) * z + d1(6)  &
!!$        ) * z + d1(5)  &
!!$        ) * z + d1(4)  &
!!$        ) * z + d1(3)  &
!!$        ) * z + d1(2)  &
!!$        ) * z + d1(1)  &
!!$        ) * z + d10
!!$
!!$      c2 = ((((((((( &
!!$                d2(10) &
!!$          * z + d2(9) &
!!$        ) * z + d2(8) &
!!$        ) * z + d2(7) &
!!$        ) * z + d2(6) &
!!$        ) * z + d2(5) &
!!$        ) * z + d2(4) &
!!$        ) * z + d2(3) &
!!$        ) * z + d2(2) &
!!$        ) * z + d2(1) &
!!$        ) * z + d20
!!$
!!$      c3 = ((((((( &
!!$                d3(8) &
!!$          * z + d3(7) &
!!$        ) * z + d3(6) &
!!$        ) * z + d3(5) &
!!$        ) * z + d3(4) &
!!$        ) * z + d3(3) &
!!$        ) * z + d3(2) &
!!$        ) * z + d3(1) &
!!$        ) * z + d30
!!$
!!$      c4 = ((((( d4(6)*z+d4(5))*z+d4(4))*z+d4(3))*z+d4(2))*z+d4(1))*z + d40
!!$
!!$      c5 = (((d5(4)*z+d5(3))*z+d5(2))*z+d5(1))*z + d50
!!$
!!$      c6 = ( d6(2) * z + d6(1) ) * z + d60
!!$
!!$      t = ((((((   &
!!$            d70    &
!!$        * u + c6 ) &
!!$        * u + c5 ) &
!!$        * u + c4 ) &
!!$        * u + c3 ) &
!!$        * u + c2 ) &
!!$        * u + c1 ) &
!!$        * u + c0
!!$
!!$    end if
!!$
!!$  else if ( iop == 2 ) then
!!$
!!$    c0 = (((((      &
!!$            d0(6)   &
!!$      * z + d0(5) ) &
!!$      * z + d0(4) ) &
!!$      * z + d0(3) ) &
!!$      * z + d0(2) ) &
!!$      * z + d0(1) ) &
!!$      * z - 1.0D+00 / 3.0D+00
!!$
!!$    c1 = ((( d1(4) * z + d1(3) ) * z + d1(2) ) * z + d1(1) ) * z + d10
!!$    c2 = d2(1) * z + d20
!!$    t = ( c2 * u + c1 ) * u + c0
!!$
!!$  else if ( 2 < iop ) then
!!$
!!$    t = (( d0(3) * z + d0(2) ) * z + d0(1) ) * z - 1.0D+00 / 3.0D+00
!!$
!!$  end if
!!$
!!$310   continue
!!$
!!$  if ( 1.0D+00 <= l ) then
!!$    qans = c * ( w + rt2pin * t / rta )
!!$    ans = 0.5D+00 + ( 0.5D+00 - qans )
!!$  else
!!$    ans = c * ( w - rt2pin * t / rta )
!!$    qans = 0.5D+00 + ( 0.5D+00 - ans )
!!$  end if
!!$
!!$  return
!!$!
!!$!  Temme expansion for L = 1
!!$!
!!$  330 continue
!!$
!!$  if ( 3.28D-03 < a * e * e ) then
!!$    ans =  2.0D+00
!!$    return
!!$  end if
!!$
!!$  c = 0.5D+00 + ( 0.5D+00 - y )
!!$  w = ( 0.5D+00 - sqrt ( y ) &
!!$    * ( 0.5D+00 &
!!$    + ( 0.5D+00 - y / 3.0D+00 ) ) / rtpi ) / c
!!$  u = 1.0D+00 / a
!!$  z = sqrt ( z + z )
!!$
!!$  if ( l < 1.0D+00 ) then
!!$    z = -z
!!$  end if
!!$
!!$  if ( iop < 2 ) then
!!$
!!$    c0 = ((((((     &
!!$            d0(7)   &
!!$      * z + d0(6) ) &
!!$      * z + d0(5) ) &
!!$      * z + d0(4) ) &
!!$      * z + d0(3) ) &
!!$      * z + d0(2) ) &
!!$      * z + d0(1) ) &
!!$      * z - 1.0D+00 / 3.0D+00
!!$
!!$    c1 = (((((      &
!!$            d1(6)   &
!!$      * z + d1(5) ) &
!!$      * z + d1(4) ) &
!!$      * z + d1(3) ) &
!!$      * z + d1(2) ) &
!!$      * z + d1(1) ) &
!!$      * z + d10
!!$
!!$    c2 = ((((d2(5)*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20
!!$
!!$    c3 = (((d3(4)*z+d3(3))*z+d3(2))*z+d3(1))*z + d30
!!$
!!$    c4 = ( d4(2) * z + d4(1) ) * z + d40
!!$    c5 = ( d5(2) * z + d5(1) ) * z + d50
!!$    c6 = d6(1) * z + d60
!!$
!!$    t = (((((( d70 &
!!$      * u + c6 ) &
!!$      * u + c5 ) &
!!$      * u + c4 ) &
!!$      * u + c3 ) &
!!$      * u + c2 ) &
!!$      * u + c1 ) &
!!$      * u + c0
!!$
!!$  else if ( iop == 2 ) then
!!$
!!$    c0 = ( d0(2) * z + d0(1) ) * z - 1.0D+00 / 3.0D+00
!!$    c1 = d1(1) * z + d10
!!$    t = ( d20 * u + c1 ) * u + c0
!!$
!!$  else if ( 2 < iop ) then
!!$
!!$    t = d0(1) * z - 1.0D+00 / 3.0D+00
!!$
!!$  end if
!!$
!!$  go to 310
!!$!
!!$!  Special cases
!!$!
!!$  390 continue
!!$
!!$  if ( x < 0.25D+00 ) then
!!$    ans = erf ( sqrt ( x ) )
!!$    qans = 0.5D+00 + ( 0.5D+00 - ans )
!!$  else
!!$    qans = erfc1 ( 0, sqrt ( x ) )
!!$    ans = 0.5D+00 + ( 0.5D+00 - qans )
!!$  end if
!!$
!!$  return
!!$
!!$  410 continue
!!$
!!$  if ( abs ( s ) <= 2.0D+00 * e ) then
!!$    ans =  2.0D+00
!!$    return
!!$  end if
!!$
!!$  if ( x <= a ) then
!!$    ans = 0.0D+00
!!$    qans = 1.0D+00
!!$  else
!!$    ans = 1.0D+00
!!$    qans = 0.0D+00
!!$  end if
!!$
!!$  return
!!$end subroutine gamma_inc

function stvaln ( p )

!*******************************************************************************
!
!! STVALN provides starting values for the inverse of the normal distribution.
!
!  Discussion:
!
!    The routine returns an X for which it is approximately true that
!      P = CUMNOR(X),
!    that is,
!      P = Integral ( -infinity < U <= X ) exp(-U*U/2)/sqrt(2*PI) dU.
!
!  Reference:
!
!    Kennedy and Gentle,
!    Statistical Computing,
!    Marcel Dekker, NY, 1980, page 95,
!    QA276.4 K46
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the probability whose normal deviate
!    is sought.
!
!    Output, real ( kind = 8 ) STVALN, the normal deviate whose probability
!    is approximately P.
!
  implicit none

!tmp  real ( kind = 8 ) eval_pol
  real ( kind = 8 ) p
  real ( kind = 8 ) sgn
  real ( kind = 8 ) stvaln
  real ( kind = 8 ), parameter, dimension(0:4) :: xden = (/ &
    0.993484626060D-01, &
    0.588581570495D+00, &
    0.531103462366D+00, &
    0.103537752850D+00, &
    0.38560700634D-02 /)
  real ( kind = 8 ), parameter, dimension(0:4) :: xnum = (/ &
    -0.322232431088D+00, &
    -1.000000000000D+00, &
    -0.342242088547D+00, &
    -0.204231210245D-01, &
    -0.453642210148D-04 /)
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  if ( p <= 0.5D+00 ) then

    sgn = -1.0D+00
    z = p

  else

    sgn = 1.0D+00
    z = 1.0D+00 - p

  end if

  y = sqrt ( -2.0D+00 * log ( z ) )
  stvaln = y + eval_pol ( xnum, 4, y ) / eval_pol ( xden, 4, y )
  stvaln = sgn * stvaln

  return
end function stvaln
function eval_pol ( a, n, x )

!*******************************************************************************
!
!! EVAL_POL evaluates a polynomial at X.
!
!  Discussion:
!
!    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
!
!  Modified:
!
!    15 December 1999
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(0:N), coefficients of the polynomial.
!
!    Input, integer N, length of A.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) EVAL_POL, the value of the polynomial at X.
!
  implicit none

  integer n

  real ( kind = 8 ) a(0:n)
  real ( kind = 8 ) eval_pol
  integer i
  real ( kind = 8 ) term
  real ( kind = 8 ) x

  term = a(n)
  do i = n - 1, 0, -1
    term = term * x + a(i)
  end do

  eval_pol = term

  return
end function eval_pol

!!c  function erf(x) result(p)
!!c    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!c    ! calculate the error function
!!c    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!c    real(PS),intent(in) :: x
!!c    real            :: xx,p,q,bound
!!c    integer :: status
!!c    status=0
!!c    xx=1.414213562*x
!!c    call cdfnor ( 1, p, q, xx, 0.0_RP, 1.0_RP, status, bound )
!!c    if(status/=0) then
!!c       write(*,*) "error in cdfnor"
!!c       stop
!!c    end if
!!c    p=2.0_PS*p-1.0_PS
!!c  end function erf

!=======================================================================================


!!$function erfc1 ( ind, x )
!!$
!!$!*******************************************************************************
!!$!
!!$!! ERFC1 evaluates the complementary error function.
!!$!
!!$!  Modified:
!!$!
!!$!    09 December 1999
!!$!
!!$!  Reference:
!!$!
!!$!    A R DiDinato and A H Morris,
!!$!    Algorithm 708:
!!$!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!!$!    ACM Transactions on Mathematical Software,
!!$!    Volume 18, 1993, pages 360-373.
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, integer IND, chooses the scaling.
!!$!    If IND is nonzero, then the value returned has been multiplied by
!!$!    EXP(X*X).
!!$!
!!$!    Input, real ( kind = 8 ) X, the argument of the function.
!!$!
!!$!    Output, real ( kind = 8 ) ERFC1, the value of the complementary
!!$!    error function.
!!$!
!!$  implicit none
!!$
!!$  real ( kind = 8 ), dimension ( 5 ) :: a = (/ &
!!$     0.771058495001320D-04,  -0.133733772997339D-02, &
!!$     0.323076579225834D-01,   0.479137145607681D-01, &
!!$     0.128379167095513D+00 /)
!!$  real ( kind = 8 ) ax
!!$  real ( kind = 8 ), dimension(3) :: b = (/ &
!!$    0.301048631703895D-02, &
!!$    0.538971687740286D-01, &
!!$    0.375795757275549D+00 /)
!!$  real ( kind = 8 ) bot
!!$  real ( kind = 8 ), parameter :: c = 0.564189583547756D+00
!!$  real ( kind = 8 ) e
!!$  real ( kind = 8 ) erfc1
!!$!tmp  real ( kind = 8 ) exparg
!!$  integer ind
!!$  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
!!$    -1.36864857382717D-07, 5.64195517478974D-01, &
!!$     7.21175825088309D+00, 4.31622272220567D+01, &
!!$     1.52989285046940D+02, 3.39320816734344D+02, &
!!$     4.51918953711873D+02, 3.00459261020162D+02 /)
!!$  real ( kind = 8 ), dimension ( 8 ) :: q = (/  &
!!$    1.00000000000000D+00, 1.27827273196294D+01, &
!!$    7.70001529352295D+01, 2.77585444743988D+02, &
!!$    6.38980264465631D+02, 9.31354094850610D+02, &
!!$    7.90950925327898D+02, 3.00459260956983D+02 /)
!!$  real ( kind = 8 ), dimension ( 5 ) :: r = (/ &
!!$    2.10144126479064D+00, 2.62370141675169D+01, &
!!$    2.13688200555087D+01, 4.65807828718470D+00, &
!!$    2.82094791773523D-01 /)
!!$  real ( kind = 8 ), dimension ( 4 ) :: s = (/ &
!!$    9.41537750555460D+01, 1.87114811799590D+02, &
!!$    9.90191814623914D+01, 1.80124575948747D+02 /)
!!$  real ( kind = 8 ) t
!!$  real ( kind = 8 ) top
!!$  real ( kind = 8 ) w
!!$  real ( kind = 8 ) x
!!$!
!!$!  ABS ( X ) <= 0.5
!!$!
!!$  ax = abs ( x )
!!$
!!$  if ( ax <= 0.5D+00 ) then
!!$
!!$    t = x * x
!!$
!!$    top = (((( a(1) * t + a(2) ) * t + a(3) ) * t + a(4) ) * t + a(5) ) &
!!$      + 1.0D+00
!!$
!!$    bot = (( b(1) * t + b(2) ) * t + b(3) ) * t + 1.0D+00
!!$
!!$    erfc1 = 0.5D+00 + ( 0.5D+00 &
!!$      - x * ( top / bot ) )
!!$
!!$    if ( ind /= 0 ) then
!!$      erfc1 = exp ( t ) * erfc1
!!$    end if
!!$
!!$    return
!!$
!!$  end if
!!$!
!!$!  0.5 < abs ( X ) <= 4
!!$!
!!$  if ( ax <= 4.0D+00 ) then
!!$
!!$    top = (((((( p(1) * ax + p(2)) * ax + p(3)) * ax + p(4)) * ax &
!!$      + p(5)) * ax + p(6)) * ax + p(7)) * ax + p(8)
!!$
!!$    bot = (((((( q(1) * ax + q(2)) * ax + q(3)) * ax + q(4)) * ax &
!!$      + q(5)) * ax + q(6)) * ax + q(7)) * ax + q(8)
!!$
!!$    erfc1 = top / bot
!!$!
!!$!  4 < ABS ( X )
!!$!
!!$  else
!!$
!!$    if ( x <= -5.6D+00 ) then
!!$
!!$      if ( ind == 0 ) then
!!$        erfc1 =  2.0D+00
!!$      else
!!$        erfc1 =  2.0D+00  * exp ( x * x )
!!$      end if
!!$
!!$      return
!!$
!!$    end if
!!$
!!$    if ( ind == 0 ) then
!!$
!!$      if ( 100.0D+00 < x ) then
!!$        erfc1 = 0.0D+00
!!$        return
!!$      end if
!!$
!!$      if ( -exparg ( 1 ) < x * x ) then
!!$        erfc1 = 0.0D+00
!!$        return
!!$      end if
!!$
!!$    end if
!!$
!!$    t = ( 1.0D+00 / x )**2
!!$
!!$    top = ((( r(1) * t + r(2) ) * t + r(3) ) * t + r(4) ) * t + r(5)
!!$
!!$    bot = ((( s(1) * t + s(2) ) * t + s(3) ) * t + s(4) ) * t &
!!$      + 1.0D+00
!!$
!!$    erfc1 = ( c - t * top / bot ) / ax
!!$
!!$  end if
!!$!
!!$!  Final assembly.
!!$!
!!$  if ( ind /= 0 ) then
!!$
!!$    if ( x < 0.0D+00 ) then
!!$      erfc1 =  2.0D+00  * exp ( x * x ) - erfc1
!!$    end if
!!$
!!$  else
!!$
!!$    w = x * x
!!$    t = w
!!$    e = w - t
!!$    erfc1 = (( 0.5D+00 &
!!$      + ( 0.5D+00 - e )) * exp ( - t ) ) * erfc1
!!$
!!$    if ( x < 0.0D+00 ) then
!!$      erfc1 =  2.0D+00  - erfc1
!!$    end if
!!$
!!$  end if
!!$
!!$  return
!!$end function erfc1

!!$function erf ( x )
!!$
!!$!*******************************************************************************
!!$!
!!$!! ERF evaluates the error function.
!!$!
!!$!  Reference:
!!$!
!!$!    A R DiDinato and A H Morris,
!!$!    Algorithm 708:
!!$!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!!$!    ACM Transactions on Mathematical Software,
!!$!    Volume 18, 1993, pages 360-373.
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real ( kind = 8 ) X, the argument.
!!$!
!!$!    Output, real ( kind = 8 ) ERF, the value of the error function at X.
!!$!
!!$  implicit none
!!$
!!$  real ( kind = 8 ), parameter, dimension ( 5 ) :: a = (/ &
!!$    0.771058495001320D-04, &
!!$    -0.133733772997339D-02, &
!!$    0.323076579225834D-01, &
!!$    0.479137145607681D-01, &
!!$    0.128379167095513D+00 /)
!!$  real ( kind = 8 ) ax
!!$  real ( kind = 8 ), parameter, dimension ( 3 ) :: b = (/ &
!!$    0.301048631703895D-02, &
!!$    0.538971687740286D-01, &
!!$    0.375795757275549D+00 /)
!!$  real ( kind = 8 ) bot
!!$  real ( kind = 8 ), parameter :: c = 0.564189583547756D+00
!!$  real ( kind = 8 ) erf
!!$  real ( kind = 8 ), dimension ( 8 ) :: p = (/   &
!!$   -1.36864857382717D-07, 5.64195517478974D-01, &
!!$    7.21175825088309D+00, 4.31622272220567D+01, &
!!$    1.52989285046940D+02, 3.39320816734344D+02, &
!!$    4.51918953711873D+02, 3.00459261020162D+02 /)
!!$  real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
!!$    1.00000000000000D+00, 1.27827273196294D+01, &
!!$    7.70001529352295D+01, 2.77585444743988D+02, &
!!$    6.38980264465631D+02, 9.31354094850610D+02, &
!!$    7.90950925327898D+02, 3.00459260956983D+02 /)
!!$  real ( kind = 8 ), dimension ( 5 ) :: r = (/ &
!!$    2.10144126479064D+00, 2.62370141675169D+01, &
!!$    2.13688200555087D+01, 4.65807828718470D+00, &
!!$    2.82094791773523D-01 /)
!!$  real ( kind = 8 ), parameter, dimension ( 4 ) :: s = (/ &
!!$    9.41537750555460D+01, 1.87114811799590D+02, &
!!$    9.90191814623914D+01, 1.80124575948747D+02 /)
!!$  real ( kind = 8 ) t
!!$  real ( kind = 8 ) top
!!$  real ( kind = 8 ) x
!!$  real ( kind = 8 ) x2
!!$
!!$  ax = abs ( x )
!!$
!!$  if ( ax <= 0.5D+00 ) then
!!$
!!$    t = x * x
!!$
!!$    top = (((( a(1)   * t &
!!$             + a(2) ) * t &
!!$             + a(3) ) * t &
!!$             + a(4) ) * t &
!!$             + a(5) ) + 1.0D+00
!!$
!!$    bot = (( b(1) * t + b(2) ) * t + b(3) ) * t + 1.0D+00
!!$    erf = ax * ( top / bot )
!!$
!!$  else if ( ax <= 4.0D+00 ) then
!!$
!!$    top = (((((( p(1)   * ax &
!!$               + p(2) ) * ax &
!!$               + p(3) ) * ax &
!!$               + p(4) ) * ax &
!!$               + p(5) ) * ax &
!!$               + p(6) ) * ax &
!!$               + p(7) ) * ax &
!!$               + p(8)
!!$
!!$    bot = (((((( q(1) * ax + q(2) ) * ax + q(3) ) * ax + q(4) ) * ax &
!!$      + q(5) ) * ax + q(6) ) * ax + q(7) ) * ax + q(8)
!!$
!!$    erf = 0.5D+00 &
!!$      + ( 0.5D+00 - exp ( - x * x ) * top / bot )
!!$
!!$  else if ( ax < 5.8D+00 ) then
!!$
!!$    x2 = x * x
!!$    t = 1.0D+00 / x2
!!$
!!$    top = ((( r(1) * t + r(2) ) * t + r(3) ) * t + r(4) ) * t + r(5)
!!$
!!$    bot = ((( s(1) * t + s(2) ) * t + s(3) ) * t + s(4) ) * t &
!!$      + 1.0D+00
!!$
!!$    erf = ( c - top / ( x2 * bot )) / ax
!!$    erf = 0.5D+00 &
!!$      + ( 0.5D+00 - exp ( - x2 ) * erf )
!!$
!!$  else
!!$
!!$    erf = 1.0D+00
!!$
!!$  end if
!!$
!!$  if ( x < 0.0D+00 ) then
!!$    erf = -erf
!!$  end if
!!$
!!$  return
!!$end function erf

!!$function exparg ( l )
!!$
!!$!*******************************************************************************
!!$!
!!$!! EXPARG returns the largest or smallest legal argument for EXP.
!!$!
!!$!  Discussion:
!!$!
!!$!    Only an approximate limit for the argument of EXP is desired.
!!$!
!!$!  Modified:
!!$!
!!$!    09 December 1999
!!$!
!!$!  Reference:
!!$!
!!$!    A R DiDinato and A H Morris,
!!$!    Algorithm 708:
!!$!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!!$!    ACM Transactions on Mathematical Software,
!!$!    Volume 18, 1993, pages 360-373.
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, integer L, indicates which limit is desired.
!!$!    If L = 0, then the largest positive argument for EXP is desired.
!!$!    Otherwise, the largest negative argument for EXP for which the
!!$!    result is nonzero is desired.
!!$!
!!$!    Output, real ( kind = 8 ) EXPARG, the desired value.
!!$!
!!$  implicit none
!!$
!!$  integer b
!!$  real ( kind = 8 ) exparg
!!$!tmp  integer ipmpar
!!$  integer l
!!$  real ( kind = 8 ) lnb
!!$  integer m
!!$!
!!$!  Get the arithmetic base.
!!$!
!!$  b = ipmpar(4)
!!$!
!!$!  Compute the logarithm of the arithmetic base.
!!$!
!!$  if ( b == 2 ) then
!!$    lnb = 0.69314718055995D+00
!!$  else if ( b == 8 ) then
!!$    lnb = 2.0794415416798D+00
!!$  else if ( b == 16 ) then
!!$    lnb = 2.7725887222398D+00
!!$  else
!!$    lnb = log ( real ( b, kind = 8 ) )
!!$  end if
!!$
!!$  if ( l /= 0 ) then
!!$    m = ipmpar(9) - 1
!!$    exparg = 0.99999D+00 * ( m * lnb )
!!$  else
!!$    m = ipmpar(10)
!!$    exparg = 0.99999D+00 * ( m * lnb )
!!$  end if
!!$
!!$  return
!!$end function exparg

!!$function rexp ( x )
!!$
!!$!*******************************************************************************
!!$!
!!$!! REXP evaluates the function EXP(X) - 1.
!!$!
!!$!  Modified:
!!$!
!!$!    09 December 1999
!!$!
!!$!  Reference:
!!$!
!!$!    A R DiDinato and A H Morris,
!!$!    Algorithm 708:
!!$!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!!$!    ACM Transactions on Mathematical Software,
!!$!    Volume 18, 1993, pages 360-373.
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real ( kind = 8 ) X, the argument of the function.
!!$!
!!$!    Output, real ( kind = 8 ) REXP, the value of EXP(X)-1.
!!$!
!!$  implicit none
!!$
!!$  real ( kind = 8 ), parameter :: p1 =  0.914041914819518D-09
!!$  real ( kind = 8 ), parameter :: p2 =  0.238082361044469D-01
!!$  real ( kind = 8 ), parameter :: q1 = -0.499999999085958D+00
!!$  real ( kind = 8 ), parameter :: q2 =  0.107141568980644D+00
!!$  real ( kind = 8 ), parameter :: q3 = -0.119041179760821D-01
!!$  real ( kind = 8 ), parameter :: q4 =  0.595130811860248D-03
!!$  real ( kind = 8 ) rexp
!!$  real ( kind = 8 ) w
!!$  real ( kind = 8 ) x
!!$
!!$  if ( abs ( x ) <= 0.15D+00 ) then
!!$
!!$    rexp = x * ( ( ( p2 * x + p1 ) * x + 1.0D+00 ) &
!!$      / ( ( ( ( q4 * x + q3 ) * x + q2 ) * x + q1 ) * x + 1.0D+00 ) )
!!$
!!$  else
!!$
!!$    w = exp ( x )
!!$
!!$    if ( x <= 0.0D+00 ) then
!!$      rexp = ( w - 0.5D+00 ) - 0.5D+00
!!$    else
!!$      rexp = w * ( 0.5D+00 + ( 0.5D+00 - 1.0D+00 / w ) )
!!$    end if
!!$
!!$  end if
!!$
!!$  return
!!$end function rexp

!!$function rlog ( x )
!!$
!!$!*******************************************************************************
!!$!
!!$!! RLOG computes X - 1 - LN(X).
!!$!
!!$!  Modified:
!!$!
!!$!    06 August 2004
!!$!
!!$!  Reference:
!!$!
!!$!    A R DiDinato and A H Morris,
!!$!    Algorithm 708:
!!$!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!!$!    ACM Transactions on Mathematical Software,
!!$!    Volume 18, 1993, pages 360-373.
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, real ( kind = 8 ) X, the argument of the function.
!!$!
!!$!    Output, real ( kind = 8 ) RLOG, the value of the function.
!!$!
!!$  implicit none
!!$
!!$  real ( kind = 8 ), parameter :: a  =  0.566749439387324D-01
!!$  real ( kind = 8 ), parameter :: b  =  0.456512608815524D-01
!!$  real ( kind = 8 ), parameter :: half = 0.5D+00
!!$  real ( kind = 8 ), parameter :: p0 =  0.333333333333333D+00
!!$  real ( kind = 8 ), parameter :: p1 = -0.224696413112536D+00
!!$  real ( kind = 8 ), parameter :: p2 =  0.620886815375787D-02
!!$  real ( kind = 8 ), parameter :: q1 = -0.127408923933623D+01
!!$  real ( kind = 8 ), parameter :: q2 =  0.354508718369557D+00
!!$  real ( kind = 8 ) r
!!$  real ( kind = 8 ) rlog
!!$  real ( kind = 8 ) t
!!$  real ( kind = 8 ), parameter :: two =  2.0D+00
!!$  real ( kind = 8 ) u
!!$  real ( kind = 8 ) w
!!$  real ( kind = 8 ) w1
!!$  real ( kind = 8 ) x
!!$
!!$  if ( x < 0.61D+00 ) then
!!$
!!$    r = ( x - 0.5D+00 ) - 0.5D+00
!!$    rlog = r - log ( x )
!!$
!!$  else if ( x < 1.57D+00 ) then
!!$
!!$    if ( x < 0.82D+00 ) then
!!$
!!$      u = x - 0.7D+00
!!$      u = u / 0.7D+00
!!$      w1 = a - u * 0.3D+00
!!$
!!$    else if ( x < 1.18D+00 ) then
!!$
!!$      u = ( x - half ) - half
!!$      w1 = 0.0D+00
!!$
!!$    else if ( x < 1.57D+00 ) then
!!$
!!$      u = 0.75D+00 * x - 1.0D+00
!!$      w1 = b + u / 3.0D+00
!!$
!!$    end if
!!$
!!$    r = u / ( u + two )
!!$    t = r * r
!!$    w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
!!$    rlog = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1
!!$
!!$  else if ( 1.57D+00 <= x ) then
!!$
!!$    r = ( x - half ) - half
!!$    rlog = r - log ( x )
!!$
!!$  end if
!!$
!!$  return
!!$end function rlog

!!$function ipmpar ( i )
!!$
!!$!*******************************************************************************
!!$!
!!$!! IPMPAR returns integer machine constants.
!!$!
!!$!  Discussion:
!!$!
!!$!    Input arguments 1 through 3 are queries about integer arithmetic.
!!$!    We assume integers are represented in the N-digit, base A form
!!$!
!!$!      sign * ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
!!$!
!!$!    where 0 <= X(0:N-1) < A.
!!$!
!!$!    Then:
!!$!
!!$!      IPMPAR(1) = A, the base of integer arithmetic;
!!$!      IPMPAR(2) = N, the number of base A digits;
!!$!      IPMPAR(3) = A**N - 1, the largest magnitude.
!!$!
!!$!    It is assumed that the single and real ( kind = 8 ) floating
!!$!    point arithmetics have the same base, say B, and that the
!!$!    nonzero numbers are represented in the form
!!$!
!!$!      sign * (B**E) * (X(1)/B + ... + X(M)/B**M)
!!$!
!!$!    where X(1:M) is one of { 0, 1,..., B-1 }, and 1 <= X(1) and
!!$!    EMIN <= E <= EMAX.
!!$!
!!$!    Input argument 4 is a query about the base of real arithmetic:
!!$!
!!$!      IPMPAR(4) = B, the base of single and real ( kind = 8 ) arithmetic.
!!$!
!!$!    Input arguments 5 through 7 are queries about single precision
!!$!    floating point arithmetic:
!!$!
!!$!     IPMPAR(5) = M, the number of base B digits for single precision.
!!$!     IPMPAR(6) = EMIN, the smallest exponent E for single precision.
!!$!     IPMPAR(7) = EMAX, the largest exponent E for single precision.
!!$!
!!$!    Input arguments 8 through 10 are queries about real ( kind = 8 )
!!$!    floating point arithmetic:
!!$!
!!$!     IPMPAR(8) = M, the number of base B digits for real ( kind = 8 ).
!!$!     IPMPAR(9) = EMIN, the smallest exponent E for real ( kind = 8 ).
!!$!     IPMPAR(10) = EMAX, the largest exponent E for real ( kind = 8 ).
!!$!
!!$!  Reference:
!!$!
!!$!    Fox, Hall, and Schryer,
!!$!    Algorithm 528,
!!$!    Framework for a Portable FORTRAN Subroutine Library,
!!$!    ACM Transactions on Mathematical Software,
!!$!    Volume 4, 1978, pages 176-188.
!!$!
!!$!  Parameters:
!!$!
!!$!    Input, integer I, the index of the desired constant.
!!$!
!!$!    Output, integer IPMPAR, the value of the desired constant.
!!$!
!!$  implicit none
!!$
!!$  integer i
!!$  integer imach(10)
!!$  integer ipmpar
!!$!
!!$!     MACHINE CONSTANTS FOR AMDAHL MACHINES.
!!$!
!!$!     data imach( 1) /   2 /
!!$!     data imach( 2) /  31 /
!!$!     data imach( 3) / 2147483647 /
!!$!     data imach( 4) /  16 /
!!$!     data imach( 5) /   6 /
!!$!     data imach( 6) / -64 /
!!$!     data imach( 7) /  63 /
!!$!     data imach( 8) /  14 /
!!$!     data imach( 9) / -64 /
!!$!     data imach(10) /  63 /
!!$!
!!$!     Machine constants for the AT&T 3B SERIES, AT&T
!!$!     PC 7300, AND AT&T 6300.
!!$!
!!$!     data imach( 1) /     2 /
!!$!     data imach( 2) /    31 /
!!$!     data imach( 3) / 2147483647 /
!!$!     data imach( 4) /     2 /
!!$!     data imach( 5) /    24 /
!!$!     data imach( 6) /  -125 /
!!$!     data imach( 7) /   128 /
!!$!     data imach( 8) /    53 /
!!$!     data imach( 9) / -1021 /
!!$!     data imach(10) /  1024 /
!!$!
!!$!     Machine constants for the BURROUGHS 1700 SYSTEM.
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   33 /
!!$!     data imach( 3) / 8589934591 /
!!$!     data imach( 4) /    2 /
!!$!     data imach( 5) /   24 /
!!$!     data imach( 6) / -256 /
!!$!     data imach( 7) /  255 /
!!$!     data imach( 8) /   60 /
!!$!     data imach( 9) / -256 /
!!$!     data imach(10) /  255 /
!!$!
!!$!     Machine constants for the BURROUGHS 5700 SYSTEM.
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   39 /
!!$!     data imach( 3) / 549755813887 /
!!$!     data imach( 4) /    8 /
!!$!     data imach( 5) /   13 /
!!$!     data imach( 6) /  -50 /
!!$!     data imach( 7) /   76 /
!!$!     data imach( 8) /   26 /
!!$!     data imach( 9) /  -50 /
!!$!     data imach(10) /   76 /
!!$!
!!$!     Machine constants for the BURROUGHS 6700/7700 SYSTEMS.
!!$!
!!$!     data imach( 1) /      2 /
!!$!     data imach( 2) /     39 /
!!$!     data imach( 3) / 549755813887 /
!!$!     data imach( 4) /      8 /
!!$!     data imach( 5) /     13 /
!!$!     data imach( 6) /    -50 /
!!$!     data imach( 7) /     76 /
!!$!     data imach( 8) /     26 /
!!$!     data imach( 9) / -32754 /
!!$!     data imach(10) /  32780 /
!!$!
!!$!     Machine constants for the CDC 6000/7000 SERIES
!!$!     60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
!!$!     ARITHMETIC (NOS OPERATING SYSTEM).
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   48 /
!!$!     data imach( 3) / 281474976710655 /
!!$!     data imach( 4) /    2 /
!!$!     data imach( 5) /   48 /
!!$!     data imach( 6) / -974 /
!!$!     data imach( 7) / 1070 /
!!$!     data imach( 8) /   95 /
!!$!     data imach( 9) / -926 /
!!$!     data imach(10) / 1070 /
!!$!
!!$!     Machine constants for the CDC CYBER 995 64 BIT
!!$!     ARITHMETIC (NOS/VE OPERATING SYSTEM).
!!$!
!!$!     data imach( 1) /     2 /
!!$!     data imach( 2) /    63 /
!!$!     data imach( 3) / 9223372036854775807 /
!!$!     data imach( 4) /     2 /
!!$!     data imach( 5) /    48 /
!!$!     data imach( 6) / -4096 /
!!$!     data imach( 7) /  4095 /
!!$!     data imach( 8) /    96 /
!!$!     data imach( 9) / -4096 /
!!$!     data imach(10) /  4095 /
!!$!
!!$!     Machine constants for the CRAY 1, XMP, 2, AND 3.
!!$!
!!$!     data imach( 1) /     2 /
!!$!     data imach( 2) /    63 /
!!$!     data imach( 3) / 9223372036854775807 /
!!$!     data imach( 4) /     2 /
!!$!     data imach( 5) /    47 /
!!$!     data imach( 6) / -8189 /
!!$!     data imach( 7) /  8190 /
!!$!     data imach( 8) /    94 /
!!$!     data imach( 9) / -8099 /
!!$!     data imach(10) /  8190 /
!!$!
!!$!     Machine constants for the data GENERAL ECLIPSE S/200.
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   15 /
!!$!     data imach( 3) / 32767 /
!!$!     data imach( 4) /   16 /
!!$!     data imach( 5) /    6 /
!!$!     data imach( 6) /  -64 /
!!$!     data imach( 7) /   63 /
!!$!     data imach( 8) /   14 /
!!$!     data imach( 9) /  -64 /
!!$!     data imach(10) /   63 /
!!$!
!!$!     Machine constants for the HARRIS 220.
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   23 /
!!$!     data imach( 3) / 8388607 /
!!$!     data imach( 4) /    2 /
!!$!     data imach( 5) /   23 /
!!$!     data imach( 6) / -127 /
!!$!     data imach( 7) /  127 /
!!$!     data imach( 8) /   38 /
!!$!     data imach( 9) / -127 /
!!$!     data imach(10) /  127 /
!!$!
!!$!     Machine constants for the HONEYWELL 600/6000
!!$!     AND DPS 8/70 SERIES.
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   35 /
!!$!     data imach( 3) / 34359738367 /
!!$!     data imach( 4) /    2 /
!!$!     data imach( 5) /   27 /
!!$!     data imach( 6) / -127 /
!!$!     data imach( 7) /  127 /
!!$!     data imach( 8) /   63 /
!!$!     data imach( 9) / -127 /
!!$!     data imach(10) /  127 /
!!$!
!!$!     Machine constants for the HP 2100
!!$!     3 WORD real ( kind = 8 ) OPTION WITH FTN4
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   15 /
!!$!     data imach( 3) / 32767 /
!!$!     data imach( 4) /    2 /
!!$!     data imach( 5) /   23 /
!!$!     data imach( 6) / -128 /
!!$!     data imach( 7) /  127 /
!!$!     data imach( 8) /   39 /
!!$!     data imach( 9) / -128 /
!!$!     data imach(10) /  127 /
!!$!
!!$!     Machine constants for the HP 2100
!!$!     4 WORD real ( kind = 8 ) OPTION WITH FTN4
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   15 /
!!$!     data imach( 3) / 32767 /
!!$!     data imach( 4) /    2 /
!!$!     data imach( 5) /   23 /
!!$!     data imach( 6) / -128 /
!!$!     data imach( 7) /  127 /
!!$!     data imach( 8) /   55 /
!!$!     data imach( 9) / -128 /
!!$!     data imach(10) /  127 /
!!$!
!!$!     Machine constants for the HP 9000.
!!$!
!!$!     data imach( 1) /     2 /
!!$!     data imach( 2) /    31 /
!!$!     data imach( 3) / 2147483647 /
!!$!     data imach( 4) /     2 /
!!$!     data imach( 5) /    24 /
!!$!     data imach( 6) /  -126 /
!!$!     data imach( 7) /   128 /
!!$!     data imach( 8) /    53 /
!!$!     data imach( 9) / -1021 /
!!$!     data imach(10) /  1024 /
!!$!
!!$!     Machine constants for the IBM 360/370 SERIES,
!!$!     THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
!!$!     5/7/9 AND THE SEL SYSTEMS 85/86.
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   31 /
!!$!     data imach( 3) / 2147483647 /
!!$!     data imach( 4) /   16 /
!!$!     data imach( 5) /    6 /
!!$!     data imach( 6) /  -64 /
!!$!     data imach( 7) /   63 /
!!$!     data imach( 8) /   14 /
!!$!     data imach( 9) /  -64 /
!!$!     data imach(10) /   63 /
!!$!
!!$!     Machine constants for the IBM PC.
!!$!
!!$!      data imach(1)/2/
!!$!      data imach(2)/31/
!!$!      data imach(3)/2147483647/
!!$!      data imach(4)/2/
!!$!      data imach(5)/24/
!!$!      data imach(6)/-125/
!!$!      data imach(7)/128/
!!$!      data imach(8)/53/
!!$!      data imach(9)/-1021/
!!$!      data imach(10)/1024/
!!$!
!!$!     Machine constants for the MACINTOSH II - ABSOFT
!!$!     MACFORTRAN II.
!!$!
!!$!     data imach( 1) /     2 /
!!$!     data imach( 2) /    31 /
!!$!     data imach( 3) / 2147483647 /
!!$!     data imach( 4) /     2 /
!!$!     data imach( 5) /    24 /
!!$!     data imach( 6) /  -125 /
!!$!     data imach( 7) /   128 /
!!$!     data imach( 8) /    53 /
!!$!     data imach( 9) / -1021 /
!!$!     data imach(10) /  1024 /
!!$!
!!$!     Machine constants for the MICROVAX - VMS FORTRAN.
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   31 /
!!$!     data imach( 3) / 2147483647 /
!!$!     data imach( 4) /    2 /
!!$!     data imach( 5) /   24 /
!!$!     data imach( 6) / -127 /
!!$!     data imach( 7) /  127 /
!!$!     data imach( 8) /   56 /
!!$!     data imach( 9) / -127 /
!!$!     data imach(10) /  127 /
!!$!
!!$!     Machine constants for the PDP-11 FORTRAN SUPPORTING
!!$!     32-BIT integer ARITHMETIC.
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   31 /
!!$!     data imach( 3) / 2147483647 /
!!$!     data imach( 4) /    2 /
!!$!     data imach( 5) /   24 /
!!$!     data imach( 6) / -127 /
!!$!     data imach( 7) /  127 /
!!$!     data imach( 8) /   56 /
!!$!     data imach( 9) / -127 /
!!$!     data imach(10) /  127 /
!!$!
!!$!     Machine constants for the SEQUENT BALANCE 8000.
!!$!
!!$!     data imach( 1) /     2 /
!!$!     data imach( 2) /    31 /
!!$!     data imach( 3) / 2147483647 /
!!$!     data imach( 4) /     2 /
!!$!     data imach( 5) /    24 /
!!$!     data imach( 6) /  -125 /
!!$!     data imach( 7) /   128 /
!!$!     data imach( 8) /    53 /
!!$!     data imach( 9) / -1021 /
!!$!     data imach(10) /  1024 /
!!$!
!!$!     Machine constants for the SILICON GRAPHICS IRIS-4D
!!$!     SERIES (MIPS R3000 PROCESSOR).
!!$!
!!$!     data imach( 1) /     2 /
!!$!     data imach( 2) /    31 /
!!$!     data imach( 3) / 2147483647 /
!!$!     data imach( 4) /     2 /
!!$!     data imach( 5) /    24 /
!!$!     data imach( 6) /  -125 /
!!$!     data imach( 7) /   128 /
!!$!     data imach( 8) /    53 /
!!$!     data imach( 9) / -1021 /
!!$!     data imach(10) /  1024 /
!!$!
!!$!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
!!$!     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
!!$!     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
!!$!
!!$  data imach( 1) /     2 /
!!$  data imach( 2) /    31 /
!!$  data imach( 3) / 2147483647 /
!!$  data imach( 4) /     2 /
!!$  data imach( 5) /    24 /
!!$  data imach( 6) /  -125 /
!!$  data imach( 7) /   128 /
!!$  data imach( 8) /    53 /
!!$  data imach( 9) / -1021 /
!!$  data imach(10) /  1024 /
!!$!
!!$!     Machine constants for the UNIVAC 1100 SERIES.
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   35 /
!!$!     data imach( 3) / 34359738367 /
!!$!     data imach( 4) /    2 /
!!$!     data imach( 5) /   27 /
!!$!     data imach( 6) / -128 /
!!$!     data imach( 7) /  127 /
!!$!     data imach( 8) /   60 /
!!$!     data imach( 9) /-1024 /
!!$!     data imach(10) / 1023 /
!!$!
!!$!     Machine constants for the VAX 11/780.
!!$!
!!$!     data imach( 1) /    2 /
!!$!     data imach( 2) /   31 /
!!$!     data imach( 3) / 2147483647 /
!!$!     data imach( 4) /    2 /
!!$!     data imach( 5) /   24 /
!!$!     data imach( 6) / -127 /
!!$!     data imach( 7) /  127 /
!!$!     data imach( 8) /   56 /
!!$!     data imach( 9) / -127 /
!!$!     data imach(10) /  127 /
!!$!
!!$  ipmpar = imach(i)
!!$
!!$  return
!!$end function ipmpar

  function acos_m(x) result(out)
    real(PS),intent(in) :: x
    real(PS) :: out
    if(x>0.999999_PS) then
       out=0.0_PS
    elseif(x<-0.999999_PS) then
       out=PI
    else
       out=acos(x)
    end if
  end function acos_m

  function fGM( p, x_min, x_max) result(out)
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate the gamma distribution
    !
    ! fGM(p, x_min, x_max) = int_x_min^x_max x^(p-1) exp(-x) dx
    !
    ! NOTE:
    !      This version only can calculate it for cases
    !           p : integer
    !           x_min = 0 and x_max = + infinity
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real(PS)    :: p
    real(PS), optional :: x_min, x_max
    real(PS)    :: out

    if( present(x_min)) then
    else; x_min = 0.0_PS
    end if
    if( present(x_max)) then
    else; x_max = 0.0_PS
    end if

    if( x_min == 0.0_PS .and. x_max == 0.0_PS ) then
       if( p - int(p) == 0.0_PS ) then
          out = factorial(int( p) - 1)
       end if
    else
       if( p  == 4.0_PS ) then
          out = (x_min**3.0+3.0_PS*x_min**2.0+6.0_PS*x_min+6.0_PS)*exp(-x_min)-&
               (x_max**3.0+3.0_PS*x_max**2.0+6.0_PS*x_max+6.0_PS)*exp(-x_max)
       end if
    end if
  end function fGM

  RECURSIVE FUNCTION factorial( n ) RESULT(res)
    INTEGER, INTENT(IN) :: n
    INTEGER :: res
    IF( n==1 ) THEN
       res = 1
    ELSE
       res = n*factorial( n-1 )
    END IF
  END FUNCTION factorial

  function zbrent(iwhich,p_mode,nx,h,d_0,mu,sig,m_l,d_coal,x1,x2,tol) result(out)
    implicit none
    integer :: iwhich
    integer :: itmax
    real(ds) :: out,tol,x1,x2,eps,p_mode,nx,h,d_0,mu,sig,m_l,d_coal

    parameter (itmax=100,eps=3.e-8)
    !  using brent's method, find the root of a function func known to lie between x1 and x2.
    !  the root, returned as zbrent, will be refined until its accuracy is tol.

    !  parameters: maximum allowed number of iterations, and machine floating-point precision.
    integer :: iter
    real(ds) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    a=x1
    b=x2
    fa=func(a)
    fb=func(b)
    if((fa.gt.0.0_ds.and.fb.gt.0.0_ds).or.(fa.lt.0.0_ds.and.fb.lt.0.0_ds)) then
!!c     write(*,*) "root must be bracketed for zbrent"
       out=-999.9_ds
       return
    end if

    c=b
    fc=fb
    do iter=1,itmax
       if((fb.gt.0.0_ds.and.fc.gt.0.0_ds).or.(fb.lt.0.0_ds.and.fc.lt.0.0_ds))then
          c=a ! rename a, b, c and adjust bounding interval d.
          fc=fa
          d=b-a
          e=d
       endif
       if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       endif
       tol1=2.0_ds*eps*abs(b)+0.5_ds*tol ! convergence check.
       xm=0.5_ds*(c-b)
       if(abs(xm).le.tol1 .or. fb.eq.0.)then
          out=b
          return
       endif
       if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa ! attempt inverse quadratic interpolation.
          if(a.eq.c) then
             p=2.0_ds*xm*s
             q=1.0_ds-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.0_ds))
             q=(q-1.0_ds)*(r-1.0_ds)*(s-1.0_ds)
          endif
          if(p.gt.0.) q=-q ! check whether in bounds.
          p=abs(p)
          if(2.0_ds*p .lt. min(3.0_ds*xm*q-abs(tol1*q),abs(e*q))) then
             e=d ! accept interpolation.
             d=p/q
          else
             d=xm ! interpolation failed, use bisection.
             e=d
          endif
       else ! bounds decreasing too slowly, use bisection.
          d=xm
          e=d
       endif
       a=b ! move last best guess to a.
       fa=fb
       if(abs(d) .gt. tol1) then ! evaluate new trial root.
          b=b+d
       else
          b=b+sign(tol1,xm)
       endif
       fb=func(b)
    enddo
    if(debug) write(*,*) "brent exceeding maximum iterations"

    out=b
    return
  contains
    function func(x) result(out)
      real(ds) :: x,out
      real(ds) :: c1,c2,a1,mass,z_coal,z_0
      if(iwhich==1) then
!!c      out=1.0_ds-(getznorm2((d_coal-mu_s)/x )-getznorm2((d_0-mu_s)/x))
!!c         out=1.0_ds-getznorm2((d_coal-mu_s)/x )

         out=x-nx/h/coedsq2p/(&
              1.0_ds-&
              min(0.99999_DS,getznorm2((d_0-mu)/max(x,1.0e-20_DS))))
!org              min(0.99999999_DS,getznorm2((d_0-mu)/max(x,1.0e-20_DS))))
!tempei bug              getznorm2((d_0-mu)/x))
      elseif(iwhich==2) then
         h=p_mode*mu*exp(0.5_ds*x**2)
         c1=log(mu)+x**2
         out=x-nx/h/coedsq2p/(&
              1.0_ds-&
!tempei bug              getznorm2((dlog(d_0)-c1)/x))
              min(0.99999999_DS,getznorm2((dlog(d_0)-c1)/max(x,1.0e-20_DS))))
      elseif(iwhich==3) then
         z_coal=(d_coal-x)/sig
         z_0=-x/sig

         c1=sig*(sig**2*(z_0**2+2.0_ds)+3.0_ds*x*(sig*z_0+x))/coedsq2p
         c2=sig*(sig**2*(z_coal**2+2.0_ds)+3.0_ds*x*(sig*z_coal+x))/coedsq2p
         a1=x*(3.0_ds*sig**2.0+x**2)

         if(abs(z_coal)>1.0e+30_DS) then
            c2=0.0_ds
         else
            c2=c2*dexp(-z_coal**2/2.0_ds)
         end if

         if(abs(z_0)>1.0e+30_DS) then
            c1=0.0_ds
         else
            c1=c1*dexp(-z_0**2/2.0_ds)
         end if

         mass=-c2+c1+a1*(getznorm2(z_coal)-getznorm2(z_0))
         out=(pi/6.0_ds)*mass-m_l*(getznorm2(z_coal)-getznorm2(z_0))

      end if
    end function func

  end function zbrent

  function get_cmod_inh(x) result(y)
    !
    ! get a coefficient to modifiy gamma to get better length growth.
    !
    real(PS),intent(in) :: x
    real(PS) :: y
    real(PS),parameter :: a1=-8.727618e-6,a2=4.806144e-4,a3=-9.762830e-3,&
                          a4=9.723944e-2,a5=5.013389e-1

    y=max(0.5_PS,min(1.0_PS,&
      a1*x**4+a2*x**3+a3*x**2+a4*x+a5))

  end function get_cmod_inh

  subroutine init_osmo_par
  use com_amps
  implicit none

    real(PS) :: x,dx
    real(PS) :: xmin,xmax
    real(PS),dimension(100) :: y
    integer :: i,n

    ! (NH4)2 SO4
    xmin=0.0
    xmax=5.5
    dx=0.1

    n=(xmax-xmin)/dx+1
    i=1
    x=xmin
    do
      if(x>xmax) exit
      y(i)=osm_ammsul(x)
      x=x+dx
      i=i+1
    enddo
    n_osm_nh42so4=n
    xs_osm_nh42so4=xmin
    do i=1,n
      y_osm_nh42so4(i)=y(i)
    enddo
    dx_osm_nh42so4=dx
    if ( IsMaster .and. debug ) then
      write(fid_alog,*) "osm NH42SO4"
      write(fid_alog,*) "N",n
      write(fid_alog,*) "xs",xmin
      write(fid_alog,*) "y",y(1:n)
    end if

    ! NaCl2 (chrolide sodium)
    xmin=0.0
    xmax=6.0
    dx=0.1

    n=(xmax-xmin)/dx+1
    i=1
    x=xmin
    do
      if(x>xmax) exit
      y(i)=osm_sodchl(x)
      x=x+dx
      i=i+1
    enddo
    n_osm_sodchl=n
    xs_osm_sodchl=xmin
    do i=1,n
      y_osm_sodchl(i)=y(i)
    enddo
    dx_osm_sodchl=dx
    if ( IsMaster .and. debug ) then
      write(fid_alog,*) "osm sodchl"
      write(fid_alog,*) "N",n
      write(fid_alog,*) "xs",xmin
      write(fid_alog,*) "y",y(1:n)
    end if

  end subroutine init_osmo_par

  function osm_ammsul(molality) result(out)
    use com_amps, only: fid_alog
    integer, parameter :: nt=23
!!$      integer, parameter :: nt=6
    real(PS) :: molality,out
    real(PS), dimension(nt) :: x,y
    integer :: i

    x = (/ 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5/)
    y = (/ 1.0,0.767,0.731,0.707,0.690,0.677,0.667,0.658,0.652,0.646,0.640,0.632,0.628,0.624,0.623,0.623 &
          ,0.626,0.635,0.647,0.660,0.673,0.686,0.699/)
!!$  x = (/ 0.0,0.1,0.5,1.0,2.0,5.0 /)
!!$  y = (/ 1.0,0.767,0.677,0.640,0.623,0.672 /)


    if(molality>=x(nt)) then
!!$       out=(y(nt)-y(nt-1))/(x(nt)-x(nt-1))*(molality-x(nt-1))+y(nt-1)
       out=y(nt)
       return
    elseif(molality<x(1)) then
       out=y(1)
       return
    else
       do i=1,nt-1
          if(x(i)<=molality.and.x(i+1)>molality) then
             out=(y(i+1)-y(i))/(x(i+1)-x(i))*(molality-x(i))+y(i)
             return
          end if
       end do
    end if

    if(debug) write(fid_alog,*) "ERROR! osm_ammsul: Molality is out of bound",molality



!!$  if(molality<x(2)) then
!!$     out=(y(2)-y(1))/(x(2)-x(1))*(molality-x(1))+y(1)
!!$  elseif(molality>=x(nt)) then
!!$     out=(y(nt)-y(nt-1))/(x(nt)-x(nt-1))*(molality-x(nt-1))+y(nt-1)
!!$  else
!!$     do i=1,nt-2
!!$        if(x(i+1)<=molality.and.x(i+2)>molality) then
!!$           out=y(i)*(molality-x(i+1))*(molality-x(i+2))/(x(i)-x(i+1))/(x(i)-x(i+2))+&
!!$                y(i+1)*(molality-x(i))*(molality-x(i+2))/(x(i+1)-x(i))/(x(i+1)-x(i+2))+&
!!$                y(i+2)*(molality-x(i))*(molality-x(i+1))/(x(i+2)-x(i))/(x(i+2)-x(i+1))
!!$           return
!!$        end if
!!$     end do
!!$  end if
!!$  if(molality>=x(nt)) then
!!$     out=(y(nt)-y(nt-1))/(x(nt)-x(nt-1))*(molality-x(nt-1))+y(nt-1)
!!$  else
!!$     do i=1,nt-2
!!$        if(x(i)<=molality.and.x(i+2)>molality) then
!!$        if(x(i)<=molality.and.x(i+1)>molality) then
!!$           out=y(i)*(molality-x(i+1))*(molality-x(i+2))/(x(i)-x(i+1))/(x(i)-x(i+2))+&
!!$                y(i+1)*(molality-x(i))*(molality-x(i+2))/(x(i+1)-x(i))/(x(i+1)-x(i+2))+&
!!$                y(i+2)*(molality-x(i))*(molality-x(i+1))/(x(i+2)-x(i))/(x(i+2)-x(i+1))
!!$           return
!!$        end if
!!$     end do
!!$  end if

  end function osm_ammsul


  function osm_sodchl(molality) result(out)
    ! for chrolide sodium (NaCl2)
    integer, parameter :: nt=24
!!$  integer, parameter :: nt=6
    real(PS) :: molality,out
    real(PS), dimension(nt) :: x,y
    integer :: i

    x = (/ 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0/)
    y = (/ 1.0,0.932,0.925,0.922,0.920,0.921,0.923,0.926,0.929,0.932,0.936,0.943,0.951,0.962,0.972,0.983,1.013 &
          ,1.045,1.080,1.116,1.153,1.192,1.231,1.271/)

!!$  x = (/ 0.0,0.1,0.5,1.0,2.0,5.0 /)
!!$  y = (/ 1.0,0.767,0.677,0.640,0.623,0.672 /)


    if(molality>=x(nt)) then
!!$       out=(y(nt)-y(nt-1))/(x(nt)-x(nt-1))*(molality-x(nt-1))+y(nt-1)
       out=y(nt)
       return
    elseif(molality<x(1)) then
       out=y(1)
       return
    else
       do i=1,nt-1
          if(x(i)<=molality.and.x(i+1)>molality) then
             out=(y(i+1)-y(i))/(x(i+1)-x(i))*(molality-x(i))+y(i)
             return
          end if
       end do
    end if
  end function osm_sodchl

  subroutine init_normal_lut
  use com_amps
  implicit none
    real(8) :: x,dx
    real(8) :: xmin,xmax
    real(PS),dimension(501) :: y
    integer :: i,n

    real(8) :: mean,sd,p,q,bound
    integer :: status

    ! standard normal distribution
    mean=0.0d+0
    sd=1.0d+0

    xmin=0.0d+0
    xmax=5.0d+0
    dx=0.01d+0

    n=(xmax-xmin)/dx+1
    i=1
    x=xmin
    do
      if(x>xmax) exit
      call cdfnor(1,p,q,x,mean,sd,status,bound)
      y(i)=1.0-p
      x=x+dx
      i=i+1
    enddo
    if ( IsMaster .and. debug ) then
      write(fid_alog,*) "normal_lut"
      write(fid_alog,*) "N",n
      write(fid_alog,*) "xs",xmin
      write(fid_alog,*) "y",y(1:n)
    end if

    n_snrml=n
    xs_snrml=xmin
    dx_snrml=dx
    y_snrml(1:n)=y(1:n)


  end subroutine init_normal_lut

  subroutine init_inv_normal_lut
!
! lookup tables of inverse normal distribution
!
  use com_amps
  implicit none
    real(8) :: y,dy
    real(8) :: ymin,ymax
    real(PS),dimension(501) :: x
    real(PS),dimension(501) :: ye
    integer :: i,n

    real(8) :: mean,sd,p,q
    integer :: status

    ymin=1.0e-30
    ymax=0.5d+0

    n=501
    dy=(log10(ymax)-log10(ymin))/real(n-1,PS_KIND)
    if ( IsMaster .and. debug ) then
      write(fid_alog,*) "dy is ",dy
    end if
    i=1
    y=ymin
    do
      if(y>ymax) exit

      x(i)= dinvnr( y, 1.0d+0-y )
      y=10.0**(log10(y)+dy)
      ye(i)=y
      i=i+1
    enddo
    n=i-1
!    x(1)=-11.46402469_RP   ! corresponding to 1.0e-30 CDF
!    x(n)=7.941444487_RP    ! very close to 1.0
!    do i=2,n
!      if(x(i)==0.0_RP) then
!        n=n-1
!        x(n)=7.941444487_RP    ! very close to 1.0
!        exit
!      endif
!    enddo
    if ( IsMaster .and. debug ) then
      write(fid_alog,*) "inv_normal_lut"
      write(fid_alog,*) "N",n
      write(fid_alog,*) "ys",ymin
      write(fid_alog,*) "ye",ye(1:n)
      write(fid_alog,*) "x",x(1:n)
    end if

    n_isnrml=n
    xs_isnrml=log10(ymin)
    dx_isnrml=dy
    y_isnrml(1:n)=x(1:n)

  end subroutine init_inv_normal_lut

!!$  SUBROUTINE LINSYS( MAT, B, N, IER)
!!$    !     -----------------
!!$    !
!!$    !     THIS ROUTINE SOLVES A SYSTEM OF LINEAR EQUATIONS
!!$    !             A*X = B
!!$    !     THE METHOD USED IS GAUSSIAN ELIMINATION WITH
!!$    !     PARTIAL PIVOTING.
!!$    !
!!$    !     INPUT:
!!$    !     THE COEFFICIENT MATRIX  A  IS STORED IN THE ARRAY MAT.
!!$    !     THE RIGHT SIDE CONSTANTS ARE IN THE ARRAY  B.
!!$    !     THE ORDER OF THE LINEAR SYSTEM IS  N.
!!$    !     THE VARIABLE  MD  IS THE NUMBER OF ROWS THAT  MAT
!!$    !     IS DIMENSIONED AS HAVING IN THE CALLING PROGRAM.
!!$    !     THE SIZE OF  MAXPIV, GIVEN BELOW, MUST BE GREATER
!!$    !     THAN  N.  IF NOT, IT IS A FATAL ERROR.
!!$    !
!!$    !     OUTPUT:
!!$    !     THE ARRAY  B  CONTAINS THE SOLUTION  X.
!!$    !     MAT CONTAINS THE UPPER TRIANGULAR MATRIX  U
!!$    !     OBTAINED BY ELIMINATION.  THE ROW MULTIPLIERS
!!$    !     USED IN THE ELIMINATION ARE STORED IN THE
!!$    !     LOWER TRIANGULAR PART OF  MAT.
!!$    !     IER=0 MEANS THE MATRIX  A  WAS COMPUTATIONALLY
!!$    !     NONSINGULAR, AND THE GAUSSIAN ELIMINATION
!!$    !     WAS COMPLETED SATISFACTORILY.
!!$    !     IER=1 MEANS THAT THE MATRIX  A  WAS
!!$    !     COMPUTATIONALLY SINGULAR.
!!$    !
!!$    integer, intent(in)             :: N
!!$    real(PS), dimension(:,:)        :: MAT
!!$    real(PS), dimension(:)          :: B
!!$    real(PS)                        :: MULT
!!$    integer, dimension(N)           :: PIVOT
!!$    integer, intent(out)            :: IER
!!$    !
!!$    integer                         :: K, I, J, var_Status
!!$    real(PS)                        :: AMAX, ABSA, TEMP, SUM
!!$    !
!!$
!!$    !     BEGIN ELIMINATION STEPS.
!!$    DO K=1,N-1
!!$       !       CHOOSE PIVOT ROW.
!!$       PIVOT(K) = K
!!$       AMAX = ABS(MAT(K,K))
!!$       DO I=K+1,N
!!$          ABSA = ABS(MAT(I,K))
!!$          IF(ABSA .GT. AMAX) THEN
!!$             PIVOT(K) = I
!!$             AMAX = ABSA
!!$          END IF
!!$       end do
!!$       !
!!$       IF(AMAX .EQ. 0.0) THEN
!!$          !           COEFFICIENT MATRIX IS SINGULAR.
!!$          IER = 1
!!$          RETURN
!!$       END IF
!!$       !
!!$       IF(PIVOT(K) .NE. K) THEN
!!$          !           SWITCH ROWS K AND PIVOT(K).
!!$          I = PIVOT(K)
!!$          TEMP = B(K)
!!$          B(K) = B(I)
!!$          B(I) = TEMP
!!$          DO J=K,N
!!$             TEMP = MAT(K,J)
!!$             MAT(K,J) = MAT(I,J)
!!$             MAT(I,J) = TEMP
!!$          end do
!!$       END IF
!!$       !
!!$       !       PERFORM STEP #K OF ELIMINATION.
!!$       DO I=K+1,N
!!$          MULT = MAT(I,K)/MAT(K,K)
!!$          MAT(I,K) = MULT
!!$          B(I) = B(I) - MULT*B(K)
!!$          DO J=K+1,N
!!$             MAT(I,J) = MAT(I,J) - MULT*MAT(K,J)
!!$          end do
!!$       end do
!!$    end do
!!$
!!$
!!$    !
!!$    IF(MAT(N,N) .EQ. 0.0) THEN
!!$       !         COEFFICIENT MATRIX IS SINGULAR.
!!$       IER = 1
!!$       RETURN
!!$    END IF
!!$    !
!!$    !     SOLVE FOR SOLUTION X USING BACK SUBSTITUTION.
!!$    DO I=N,1,-1
!!$       SUM = 0.0
!!$       DO J=I+1,N
!!$          SUM = SUM + MAT(I,J)*B(J)
!!$       end do
!!$       B(I) = (B(I) - SUM)/MAT(I,I)
!!$    end do
!!$    IER = 0
!!$    RETURN
!!$  end SUBROUTINE LINSYS

!!$  SUBROUTINE DIANT(DIA,ENT,RX,DN0,FCTMS,PWMAS,RDIST)
!!$    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! calculate mean diameter for the given total concentration and distribution.
!!$    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! mean diameter
!!$    real(PS), intent(inout)        :: DIA
!!$    ! total concentration
!!$    real(PS), intent(in)        :: ENT
!!$    ! mixing ratio
!!$    real(PS), intent(in)        :: RX
!!$    ! density of air
!!$    real(PS), intent(in)        :: DN0
!!$    ! coefficient of mass-dimension relation: m=FCTMS*(D**PWMAS)
!!$    real(PS), intent(in)        :: FCTMS, PWMAS
!!$    ! type of distribution
!!$    integer, intent(in)         :: RDIST
!!$    real(PS)                    :: GAMFCN, POW
!!$    IF( RDIST .EQ. 2) THEN
!!$       GAMFCN = fGM( PWMAS+1.0)
!!$    ELSE
!!$       GAMFCN = 1.0_PS
!!$    ENDIF
!!$    POW=1.0_PS/PWMAS
!!$
!!$    if( RX > 0.0_PS) then
!!$          DIA = ( DN0 * RX/(FCTMS*MAX(1.0E-20_PS,ENT)*GAMFCN))**POW
!!$    else
!!$          DIA =0.
!!$    endif
!!$
!!$  END SUBROUTINE DIANT

end module mod_amps_utility

!!$module tic_toc
!!$! abandoned for scale
!!$! replaced by scale_prof
!!$! Define global constants for timing increments
!!$use acc_amps
!!$implicit none
!!$private
!!$!  integer, private :: start   ! current value of system clock
!!$!  integer, private :: rate    ! system clock counts/sec
!!$!  integer, private :: finish  ! ending value of system clock
!!$! Useage: use tic_toc      ! for access to data
!!$!         call tic         ! start clock
!!$!         ...              ! use some cpu time
!!$!         cputime = toc () ! for increment in seconds
!!$
!!$contains ! access to start, rate, finish

!!$subroutine  tic(s1,r1)
!!$! -------------------------------------------------
!!$! Model the matlab tic function, for use with toc
!!$! -------------------------------------------------
!!$  real(PS),intent(inout) :: s1,r1
!!$  real(PS),dimension(2) :: tarray
!!$  real(PS) :: result
!!$! CHIARUI >>>>>>>>>> 8/4/2019
!!$  !real(4),external :: etime
!!$! CHIARUI <<<<<<<<<<
!!$!!!  call system_clock ( COUNT = start, COUNT_RATE = rate )
!!$!!!  call system_clock ( COUNT = s1, COUNT_RATE = r1)
!!$  !result=etime(tarray)
!!$!uzushio  call etime(tarray,result)
!!$  s1=result
!!$  r1=0.0
!!$end subroutine  tic

!!$real(PS) function  toc (s1,r1)
!!$! -------------------------------------------------
!!$! Model the matlab toc function, for use with tic
!!$! -------------------------------------------------
!!$! Useage:   call tic         ! start clock
!!$!           cputime = toc () ! for increment
!!$  real(PS),intent(in) :: s1,r1
!!$  real(PS),dimension(2) :: tarray
!!$  real(PS) :: result
!!$! CHIARUI >>>>>>>>>> 8/4/2019
!!$  !real(4),external :: etime
!!$! CHIARUI <<<<<<<<<<
!!$!!!  call system_clock ( COUNT = finish ) ! stop execution timer
!!$!!!  toc = float( finish - start ) / float( rate )
!!$!!!  call system_clock ( COUNT = f1 ) ! stop execution timer
!!$!!!  toc = float( f1 - s1 ) / float( r1 )
!!$  !result=etime(tarray)
!!$!uzushio  call etime(tarray,result)
!!$
!!$  toc = result-s1
!!$end function  toc

!!$end module tic_toc
