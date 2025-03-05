#include "scalelib.h"
!OCL SERIAL
module mod_amps_lib
  use scale_io
  use acc_amps
!  use mod_amps_const
  implicit none

  private

  public :: ifc_cloud_micro
  public :: ini_aerosols_mpace
  public :: ini_aerosols_const
  public :: ini_aerosols_sheba
  public :: ini_aerosols_profile_mpace
  public :: ini_aerosols_profile_const
  public :: ini_aerosols_profile_sheba
  public :: binmicrosetup_scale
  public :: get_sat_vapor_pres_lk

  real(RP), parameter, private :: p0 = 100000.0_RP
!  real(RP), parameter, private :: r  = R_d / 10000_RP
  real(RP), parameter, private :: r  = 287.04_PS
!  real(RP), parameter, private :: rm = R_v / 10000_RP
  real(RP), parameter, private :: rm = 461.50_RP
!  real(RP), parameter, private :: cp = C_pa / 10000_RP
  real(RP), parameter, private :: cp = 1004.64_PS
contains

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ifc_cloud_micro.f90
!   Interface for the cloud microphyscs program "class_cloud_micro.f90"
!
!   Version 2.0
!                                        2004/05/28      Tempei Hashino
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!tmpSUBROUTINE ifc_cloud_micro( &
!tmp     XC, vtcloud  &
SUBROUTINE ifc_cloud_micro( &
     CM &
    ,XC, vtcloud  &
    ,XR, NRTYPE,NRBIN, NRCAT &
    ,XS, NSTYPE,NSBIN, NSCAT &
    ,XA, NATYPE,NABIN, NACAT &
    ,RV, DEN, PT, T, W &
    ,L,ID,JD,KD,iproc_t,iproc,istrt &
    ,qtp,thil &
    ,jseed,ifrst,isect_seed &
    ,nextn &
    ,dmtendl,dcontendl,dbintendl &
!    ,mpirank_w,mpirank_e,mpirank_s,mpirank_n &
! added for sheba
!    ,dmtend,tcon,tmass,tvar,dctend_in &
!    ,zstv,zz,xx,yy,nzp,nxp,nyp,nfpt,nfpth &
!    ,dmtend_z,tcon_z,tmass_z,n1mx &
! end added for sheba
! <<< 2014/10 T. Hashino added for KiD
!    ,dM_auto_liq,dM_accr_liq &
!    ,dM_auto_ice,dM_accr_ice &
!    ,dM_auto_rim,dM_accr_rim &
!    ,cptacc_ifc,cptacc_mtd &
! >>> 2014/10 T. Hashino added for KiD
     )
  !use tic_toc
  use maxdims
  use class_Cloud_Micro, only: &
     Cloud_Micro, &
     ini_cloud_micro, &
     reality_check, &
     check_water_apmass, &
     cal_micro_tendency, &
     return_output, &
     update_terminal_vel
!  use class_group, only: denchk1
!!!  use mod_amps_utility, only: random_genvar
  implicit none

  type (Cloud_Micro),intent(inout)       :: CM
!tmp  type (Cloud_Micro)      :: CM
!tmp  COMMON /ICMPRIV/CM
!!!!$omp threadprivate(/ICMPRIV/)

  ! C: cloud_drop
  ! R: rain
  ! S: solid hydrometeor
  ! NXBIN       : number of bins
  ! NXTYPE      : number of variables that each bin contain
  ! X(IBIN,ITYPE)    : quantity for ITYPE of IBIN-th bin
  ! DXDT(ITYPE) : tendency of the quantity
  ! ITYPE : 1. mixing ratio (g/g)
  !         2. concentration (#/g)
  !         3. volume of circumscribing sphere (cm^3/g)
  !         4. con. weighted a-axis length^3 (cm^3/g)
  !         5. con. weighted c-axis length^3 (cm^3/g)
  !         6. con. weighted d-axis length^3 (cm^3/g)
  !         7. con. weighted r-axis length^3 (cm^3/g)
  !         8. mixing ratio of riming production (g/g)
  !         9. mixing ratio of ice crystals (g/g)
  !        10. mixing ratio of melt water (g/g)
  !        11. heat (erg/g)
  !
!tmp  integer, intent(in)    :: NRBIN, NSBIN,NABIN
!tmp  integer, intent(in)    :: NRTYPE, NSTYPE,NATYPE
!tmp  integer, intent(in)    :: NRCAT, NSCAT,NACAT
  integer    :: NRBIN, NSBIN,NABIN
  integer    :: NRTYPE, NSTYPE,NATYPE
  integer    :: NRCAT, NSCAT,NACAT
  ! length of microphysics grid vector, thread number, flag for initial do loop
  !  iproc_t includes mpi process and openMP threads
  !  iproc   includes only openMP threads
  integer    :: L,iproc_t,iproc,istrt
  real(MP_KIND)  :: XC(L),vtcloud(L)
  real(MP_KIND)  :: XR(NRTYPE,NRBIN,NRCAT,L), XS(NSTYPE,NSBIN,NSCAT,L),XA(NATYPE,NABIN,NACAT,L)
  real(MP_KIND)  :: RV(L), DEN(L), PT(L), T(L), W(L),qtp(L),thil(L)
  integer   :: ID(L),JD(L),KD(L)

  ! random generator vars
!!!  type (random_genvar),intent(in) :: rdsd0
  integer,intent(in) ::  jseed,ifrst,isect_seed,nextn

  ! tendencies output
  real(MP_KIND), intent(inout) :: dmtendl(10,2,L),dcontendl(10,2,L),dbintendl(3,2,mxnbin,L)

  ! process ranks for MPI parallerization
  !integer,intent(in) :: mpirank_w,mpirank_e,mpirank_s,mpirank_n

  ! current model time in second
!tmp  real,intent(in) :: c_time
! added for sheba
!  real :: dmtend(11,9,*),tcon(9,*),tmass(9,*),zstv(*),zz(*),xx(*),yy(*)
!  real :: tvar(4,*),dctend_in(7,*)
!  integer :: nfpt,nfpth,nzp,nxp,nyp
!  integer :: n1mx
!  real :: dmtend_z(15,n1mx,*),tcon_z(2,n1mx,*),tmass_z(2,n1mx,*)
! end added for sheba
! <<< 2014/10 T. Hashino added for KiD
!  real, dimension(*) :: &
!                  dM_auto_liq,dM_accr_liq &
!                 ,dM_auto_ice,dM_accr_ice &
!                 ,dM_auto_rim,dM_accr_rim
! >>> 2014/10 T. Hashino added for KiD
  !real,dimension(20,*) :: cptacc_ifc,cptacc_mtd


  !real,dimension(20) :: cptime
  !real :: cpsum
!tmp  character(len=32),dimension(20)  :: proname
  integer :: i,n

  real(PS) :: s1,r1

!tmp  proname(1)="ini_cloud_micro"
!tmp  proname(2)="reality_check"
!tmp  proname(3)="check_water_apmass 1"
!tmp  proname(4)="cal_micro_tendency"
!tmp  proname(5)="print_cloud_tendency"
!tmp  proname(6)="return_output"
!tmp  proname(7)="check_water_apmass 2"
!tmp  proname(8)="update_terminal_vel"
!tmp  proname(9)="cal_dmtend"


!!c  write(fid_alog,*) "I am inside! nbin",L,KD(1),ID(1),JD(1)
!!c  write(fid_alog,'(4ES15.8)') RV(1),RV(20),T(1),T(20)

  !cptime=0.0

  ! set random number vars
  CM%rdsd%jseed=jseed
  CM%rdsd%ifrst=ifrst
  CM%rdsd%isect_seed=isect_seed
  CM%rdsd%nextn=nextn

  ! 1. Initialize the Cloud_Micro object


  !call tic(s1,r1)
!tmp  call ini_cloud_micro( &
!tmp        XC &
  call ini_cloud_micro( &
        CM &
       ,XC &
       ,XR, NRTYPE,NRBIN, NRCAT &
       ,XS, NSTYPE,NSBIN, NSCAT &
       ,XA, NATYPE,NABIN, NACAT &
       ,RV, DEN, PT, T, W &
       ,L,qtp,ID,JD,KD)
  !cptime(1)=toc(s1,r1)

!!c       n=16
!!c       write(fid_alog,*) "ck af inicloud:",kd(n),id(n),jd(n),n,CM%air%TV(n)%rv,CM%air%TV(n)%t

!!c  write(fid_alog,*) "af make cloud"
!!  n=178
!!  i=7
!!  call denchk1(CM%solid_hydro,i,n,id,jd,kd,'af_inicloud')

  ! 2. check whether the advected prognostic variables are
  !    in realistic, physical range.
  !call tic(s1,r1)

  call reality_check( CM,qtp,ID,JD,KD)
  !cptime(2)=toc(s1,r1)

!!  n=178
!!  i=7
!!  call denchk1(CM%solid_hydro,i,n,id,jd,kd,'af_real')

  !call tic(s1,r1)

  call check_water_apmass(CM,NATYPE,NABIN,NACAT,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
                          XA,XR,XS,qtp,RV,"make_cloud",0,ID,JD,KD)
  !cptime(3)=toc(s1,r1)

!!c  write(fid_alog,*) "af water apmass"
!!c       n=16
!!c       write(fid_alog,*) "ck af waer apmass:",kd(n),id(n),jd(n),n,CM%air%TV(n)%rv,CM%air%TV(n)%t

  ! 3. calculate the tendency of cloud microphysical variables
  !call tic(s1,r1)
!  call PROF_rapstart("amps_tendency",3)
!!  n=178
!!  i=7
!!  call denchk1(CM%solid_hydro,i,n,id,jd,kd,'bf_microtend')

  call cal_micro_tendency(CM,ID,JD,KD,qtp,thil,iproc_t,istrt &
                         ,L,dmtendl,dcontendl,dbintendl &
! <<< 2014/10 T. Hashino added for KiD
!                 ,dM_auto_liq,dM_accr_liq &
!                 ,dM_auto_ice,dM_accr_ice &
!                 ,dM_auto_rim,dM_accr_rim)
!                 ,cptacc_mtd(1,iproc+1))
                         )
! >>> 2014/10 T. Hashino added for KiD

  !cptime(4)=toc(s1,r1)

!!c  write(fid_alog,*) "af cal_mic_tend"

  ! 4. print tendencies
  !call tic(s1,r1)

  !call print_cloud_tendency(CM,ID,JD,KD,iproc_t,istrt)
  !cptime(5)=toc(s1,r1)
!  call PROF_rapend("amps_tendency",3)

!!c  write(fid_alog,*) "bf ro"

  ! 5. return the calculated tendencies to model prognostic variables.
  !call tic(s1,r1)

  call return_output(CM,L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       NATYPE,NABIN,NACAT,&
       RV,XC,XR,XS,XA,&
       qtp)
  !cptime(6)=toc(s1,r1)

!!c  write(fid_alog,*) "bf cktp1"
  !call tic(s1,r1)

  call check_water_apmass(CM,NATYPE,NABIN,NACAT,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       XA,XR,XS,qtp,RV,"return_model",1,ID,JD,KD)
  !cptime(7)=toc(s1,r1)


!!c  write(fid_alog,*) "bf update"

  ! 6. calculate terminal velocity
  !call tic(s1,r1)

  call update_terminal_vel(CM, L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       vtcloud, XR, XS)
  !cptime(8)=toc(s1,r1)

  ! added for sheba
  ! 7. calculate domain integrated variables for analysis
  !call tic(s1,r1)

! the following is for 3D dynamical simulation
!!c  call cal_dmtend(dmtend,tcon,tmass,tvar,dctend_in,iproc, &
!!c         CM,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT, &
!!c         XR,XS,KD,ID,JD,zstv,zz,xx,yy,nzp,nxp,nyp,nfpt,nfpth, &
!!c         mpirank_w,mpirank_e,mpirank_s,mpirank_n)

! the following is for KiD
  !call cal_dmtend_kid(dmtend,tcon,tmass,iproc,&
  !     CM,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
  !     XR,XS,KD,ID,JD,zstv,zz,xx,yy,nzp,nxp,nyp,nfpt,nfpth)

  !if( mod(int(CM%cur_time), 60) == 0 ) then
  !  call cal_dmtend_z(dmtend_z,tcon_z,tmass_z,iproc,&
  !       CM,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
  !       dM_auto_liq,dM_accr_liq, &
  !       XR,XS,KD,ID,JD,zstv,zz,xx,yy,nzp,nxp,nyp,nfpt,nfpth,&
  !       n1mx,&
  !       mpirank_w,mpirank_e,mpirank_s,mpirank_n)
  !endif
  !cptime(9)=toc(s1,r1)

  !cptacc_ifc(1:20,iproc+1)=cptacc_ifc(1:20,iproc+1)+cptime(1:20)

!tmp  write(fid_alog,*) "iproc+1:",iproc+1
!tmp  write(fid_alog,*) "ifc, cptime",cptime(1:20)
!tmp  write(fid_alog,*) "ifc, cptacc",cptacc_ifc(1:20,iproc+1)

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!tmp  if(mod(CM%cur_time,600.0)==0) then
!!!  if(mod(CM%cur_time,1800.0)==0) then
!tmp    write(fid_alog,*) "CPU time (s)"
!tmp    cpsum=sum(CM%cptime_acc_ifc)
!tmp    if(cpsum>0.0) then
!tmp      do i=1,9
!tmp        write(fid_alog,100) proname(i),CM%cptime_acc_ifc(i),CM%cptime_acc_ifc(i)*100.0/cpsum
!tmp      end do
!tmp      write(fid_alog,*) " "
!tmp100     format(A20,":",ES13.4,"s,",F6.2,"%")
!tmp      CM%cptime_acc_ifc(1:20)=0.0
!tmp    end if
!tmp  end if

  ! end added for sheba

!!! stop


END SUBROUTINE ifc_cloud_micro

subroutine binmicrosetup_scale(npr,nbr,ncr,npi,nbi,nci, &
                               npa,nba,nca,model_k, model_i, model_j,nbhzcl, &
                               estbar,esitbar)
  use scale_prc, only: &
     PRC_abort
  use maxdims
  use com_amps
  use mod_amps_utility, only: &
       read_AMPSTASK, &
       RDCETB, &
       RDAPTB, &
       RDSTTB, &
       init_inherent_growth_par, &
       init_osmo_par, &
       init_normal_lut, &
       init_inv_normal_lut, &
       read_seed
  use par_amps, only: &
       isplit_bin_liq, &
       isplit_bin_ice
  implicit none
!  arguments
      integer, intent(in)  :: npr, nbr, ncr, npi, nbi, nci, npa, nba, nca, model_k, model_i, model_j
      integer, intent(out) :: nbhzcl
      real(DS), DIMENSION(*) ::  ESTBAR
      real(DS), DIMENSION(*) ::  ESITBAR

      integer :: i, j

      integer :: ierr

!
!     IXCNFL ...... level of complexity for X hydrometeor
!     default total number of bins
      integer :: idbnum(10)
      DATA idbnum/0, 3, 5, 10, 15, 20, 30, 40, 60, 80/
!
!     setting for terminal velocity
      real(PS) :: VRR(11)
      DATA VRR/8.545e-1,8.59895E-01,8.65272E-01,8.75851E-01 &
           ,8.92279E-01  &
           ,9.11267E-01,9.28026E-01,9.40583E-01,9.49129E-01 &
           ,9.55522E-01,9.58816E-01/
      real(PS) :: VSR(5)
      DATA VSR/6.12185E-01,6.30215E-01,7.35575E-01,8.51963E-01, &
           8.95877E-01/

!     pi/6.0, 4.0*pi/3.0, 2.0*pi
      real(PS), parameter :: coef2=0.523598776, coef3=4.18879020478639, coef4=6.28318530717959
      real(PS) :: m0, m1

!     the minimum cloud droplet mass
      real(PS) :: c_mnm, c_max
!     1 um
!      parameter(c_mnm=4.188790205e-12)
!     0.1 um
      parameter(c_mnm=4.188790205e-15,c_max=6.54498e-8)
!     0.01 um
!     c_max -> rad=25 micron
!      parameter(c_mnm=4.188790205e-18,c_max=6.54498e-8)

      real(PS) :: dsrat, maxmass_r, dsrat_h
!     corresponds to 1 cm diameter of a rain drop
      parameter(maxmass_r=5.2359870E-01)

!     the mass boundaries to define pristince crystal, aggregates, hail
      real(PS) :: mg1, mg2, mg3, mg4
!ccc      parameter(mg1=4.18879020478639E-12,mg2=1.0e-7,mg3=1.0e-1
!ccc     *     ,mg4=1.0e+1)
      parameter(mg1=4.18879020478639E-12,mg2=1.0e-6,mg3=1.0e-2 &
           ,mg4=1.0e+1)
!ccc      parameter(mg1=4.18879020478639E-12,mg2=1.0e-6,mg3=1.0e+0
!ccc     *     ,mg4=1.0e+1)


      integer :: NBIN20(3), NBIN10(3), NBIN40(3)
!ccc      NBIN20=(/4,14,2/)
!ccc      NBIN10=(/2,7,1/)
!ccc      NBIN40=(/8,28,4/)
      NBIN20=(/4,10,6/)
      NBIN10=(/2,5,3/)
      NBIN40=(/8,20,12/)

!
!   Create the file id for log
!
      call set_ampslog

      !open(fid_alog,file=fname_ampslog)
      !write(*,*) "fid_alog",fid_alog

      ! set grid maximum limits
      n1mx = model_k
      n2mx = model_i
      n3mx = model_j
      nxpmax = n2mx
      nypmax = n3mx
      nzpmax = n1mx
      nzgmax = n4mx
      mxln = max(n2mx,n1mx,n3mx)
      mxlh = max(n2mx,n3mx)
      mxlv = n1mx

      LMAX = n1mx + 1

      ! set 6D microphysical variable maximum limits
      if (nbr /= 40 .and. nbr /= 80) then
         LOG_ERROR("binmicrosetup_scale",*) "# of bins for liquid is not 40 or 80", nbr
         call PRC_abort
      endif
      if (nbi /= 10 .and. nbi /= 20 .and. nbi /= 40) then
         LOG_ERROR("binmicrosetup_scale",*) "# of bins for ice is not 10 or 20 or 40", nbi
         call PRC_abort
      endif
      if (nbr == 40) then
         allocate(bu_fd(2,17835))
         allocate(bu_tmass(435))
      else if (nbr == 80) then
         allocate(bu_fd(2,62400))
         allocate(bu_tmass(780))
      endif
      bu_fd = 0.0_PS
      bu_tmass = 0.0_PS

      nprmx = npr
      npimx = npi
      npamx = npa
      nbrmx = nbr
      nbimx = nbi
      nbamx = nba
      ncrmx = ncr
      ncimx = nci
      ncamx = nca

      mxnbinr = nbrmx
      mxnbini = nbimx
      mxnbina = nbamx

      mxnbin = max(nbrmx,nbimx)
      mxnbinb = mxnbin + 1
!
!  Read the AMPSTASK.F
!
         call read_AMPSTASK

!     add haze and cloud droplets category
!tmp tempei 2012/10         nbr = nbr+nbin_h
         nbhzcl=nbin_h
         if ( IsMaster .and. debug ) then
           write(fid_alog,*) "in binmic",nbr,nbin_h,nbhzcl
         end if
!
!     calculation of bin boundaries
         if(size(binbr).lt.nbr) then
            LOG_ERROR("binmicrosetup_scale",*) "# of bins for liquid is not enough in binbr", size(binbr), nbr
            call PRC_abort
         endif
         if(size(binbi).lt.nbi) then
            LOG_ERROR("binmicrosetup_scale",*) "# of bins for ice is not enough in binbi", size(binbi), nbi
            call PRC_abort
         endif

!     for rain category, including cloud cat.
         binbr(1) = c_mnm
!sheba         dsrat_h=(minmass_r/c_mnm)**(1.0/(nbin_h))
         dsrat_h=(c_max/c_mnm)**(1.0/(nbin_h))
         do i=2,nbin_h
            binbr(i)=binbr(i-1)*dsrat_h
         enddo
!sheba         binbr(nbin_h+1) = minmass_r
         binbr(nbin_h+1) = c_max

!sheba         dsrat=(maxmass_r/ &
!sheba              minmass_r*sth_r**(0.5* &
!sheba              real((nbr-nbin_h-1)*(nbr-nbin_h-0)))) &
!sheba              **(1.0/real(nbr-nbin_h-0))

         dsrat=(maxmass_r/c_max)**(1.0/(nbr-nbin_h))

         binbr(nbr+1)=maxmass_r

         if ( debug ) then
            write(fid_alog,*) "bin boundaries of rain"
            write(fid_alog,*) "c_mnm: ",c_mnm, c_max
!sheba         write(fid_alog,*) "minmass_r: ",minmass_r
            write(fid_alog,*) "maxmass_r: ",maxmass_r
            write(fid_alog,*) "nbin_h, dsrat_h: ",nbin_h,dsrat_h
            write(fid_alog,*) "theta: ",sth_r
            write(fid_alog,*) "bin#,binbr,srat,rad"
            do i=1,nbin_h+1
               write(fid_alog,'(I5,3ES15.6)') i,binbr(i),dsrat_h,(binbr(i)/coef3)**(1.0_RP/3.0_RP)
            end do
         end if
         do i=nbin_h+2,nbr+1
            binbr(i)=binbr(i-1)*dsrat
            if(debug) write(fid_alog,'(I5,3ES15.6)') i,binbr(i),dsrat,(binbr(i)/coef3)**(1.0_RP/3.0_RP)
!sheba            dsrat=dsrat/sth_r
         end do


!cc         do i=2,nbr+1
!cc            write(fid_alog,'(I5,2ES15.6)') i,binbr(i),binbr(i)/binbr(i-1)
!cc         end do
!ccc         do i=3,nbr+1
!ccc            binbr(i)=srat_r*binbr(i-1)+sadd_r
!ccc         end do


!     for ice category
         if(nbi==40) then
            srat_s=(mg2/mg1)**(1.0/(real(NBIN40(1),PS_KIND)))
            binbi(1) = minmass_s
            do i=2,NBIN40(1)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
            srat_s=(mg3/mg2)**(1.0/(real(NBIN40(2),PS_KIND)))
            do i=NBIN40(1)+2,NBIN40(1)+NBIN40(2)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
            srat_s=(mg4/mg3)**(1.0/(real(NBIN40(3),PS_KIND)))
            do i=NBIN40(1)+NBIN40(2)+2,NBIN40(1)+NBIN40(2)+NBIN40(3)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
         elseif(nbi==20) then
            srat_s=(mg2/mg1)**(1.0/(real(NBIN20(1),PS_KIND)))
            binbi(1) = minmass_s
            do i=2,NBIN20(1)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
            srat_s=(mg3/mg2)**(1.0/(real(NBIN20(2),PS_KIND)))
            do i=NBIN20(1)+2,NBIN20(1)+NBIN20(2)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
            srat_s=(mg4/mg3)**(1.0/(real(NBIN20(3),PS_KIND)))
            do i=NBIN20(1)+NBIN20(2)+2,NBIN20(1)+NBIN20(2)+NBIN20(3)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
         elseif(nbi==10) then
            srat_s=(mg2/mg1)**(1.0/(real(NBIN10(1),PS_KIND)))
            binbi(1) = minmass_s
            do i=2,NBIN10(1)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
            srat_s=(mg3/mg2)**(1.0/(real(NBIN10(2),PS_KIND)))
            do i=NBIN10(1)+2,NBIN10(1)+NBIN10(2)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
            srat_s=(mg4/mg3)**(1.0/(real(NBIN10(3),PS_KIND)))
            do i=NBIN10(1)+NBIN10(2)+2,NBIN10(1)+NBIN10(2)+NBIN10(3)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
         else
            LOG_ERROR("binmicrosetup_scale",*) "# of bins for ice is not 40, 20, 10", nbi
            call PRC_abort
            binbi(1) = minmass_s
            do i=2,nbi+1
               binbi(i)=srat_s*binbi(i-1)+sadd_s
            end do
         end if
!     write out the bin boundaries
         !if ( IsMaster ) then
         if ( debug ) then
           write(fid_alog,*) "bin boundaries of ice"
           write(fid_alog,*) "minmass_s: ",minmass_s
           write(fid_alog,*) "bin#,binbi,srat"
           write(fid_alog,'(I5,2ES15.6)') 1,binbi(1)
           do i=2,nbi+1
              write(fid_alog,'(I5,2ES15.6)') i,binbi(i),binbi(i)/binbi(i-1)
           end do
        end if

!
!     define split bin
!
         do i=1,nbr+1
!!!           if((binbr(i)/coef3)**(1.0/3.0).gt.25.0e-4) then
           if((binbr(i)/coef3)**(1.0/3.0).gt.20.0e-4) then
!           if((binbr(i)/coef3)**(1.0/3.0).gt.14.0e-4) then
             isplit_bin_liq=i
             exit
           endif
         end do
         do i=1,nbi+1
            if( binbi(i) > binbr(isplit_bin_liq) ) then
               isplit_bin_ice = i-1
               exit
            end if
         end do
         if ( IsMaster .and. debug ) then
           write(fid_alog,*) "Split bin for liquid and ice: ",isplit_bin_liq, isplit_bin_ice
         end if


!     Read lookup tables
         if ( IsMaster .and. debug ) then
           write(fid_alog,*) "Directory for collision coefficiency files:" &
                , DRCETB
         end if
         call RDCETB(nbr,ncr,nbi,nci)

         if ( IsMaster .and. debug ) then
           write(fid_alog,*) "Directory for aerosol activation file:" &
                , DRAPTB
         end if
         call RDAPTB()

         if ( IsMaster .and. debug ) then
           write(fid_alog,*) "Directory for statistical lookup file:" &
                , DRSTTB
         end if
         call RDSTTB()

!     set up the factor for terminal velocity
!     (concentration-weighted Vt)/(mass-weighted Vt)
!     for rain
!         do i=1,nnbr(ngrid)
!            do j=1,nnpr(ngrid)
!               if(j.eq.1.or.j.eq.9.or.j.eq.10) then
!                  VTRFCT(j,i)=1.0
!               else
!                  VTRFCT(j,i)=VRR(i)
!               end if
!            end do
!         end do

!     for solid hydrometeors
!         do i=1,nnbi(ngrid)
!            do j=1,nnpi(ngrid)
!               if(j.eq.1.or.j.eq.9.or.j.eq.10) then
!                  VTSFCT(j,i)=1.0
!               else
!                  VTSFCT(j,i)=VSR(i)
!               end if
!            end do
!         end do

!     calculate collisional-breakup fragment distribution
         call cal_breakfragment(nbr,estbar,esitbar)



!     find smallest and largest bin numbers necessary for breakup calculation
         r_frg0=50.0e-4
         r_frg1=0.45
         ibinr_frg0=0
         ibinr_frg1=0

         m0=coef3*r_frg0**3.0
         do i=1, nbr+1
            if(m0<binbr(i)) then
               ibinr_frg0=i-1
               if ( IsMaster .and. debug ) then
                 write(fid_alog,*) "ibinr_frg0",i-1,m0,binbr(i)
               end if
               exit
            end if
         end do
         m1=coef3*r_frg1**3.0
         do i = 1, nbr+1
            if(m1<binbr(i)) then
               ibinr_frg1=i-1
               if ( IsMaster .and. debug ) then
                 write(fid_alog,*) "ibinr_frg1",i-1,m1,binbr(i)
               end if
               exit
            end if
         end do
         if(ibinr_frg0==0.or.ibinr_frg1==0) then
            LOG_ERROR("binmicrosetup_scale",*) "ibinr_frg0.or.ibinr_frg1 is 0. check the bin boundary"
            call PRC_abort
         end if
         if ( IsMaster .and. debug ) then
           write(fid_alog,*) "ibinr_frg0,ibinr_frg1",ibinr_frg0,ibinr_frg1
         end if

!         sizeseed=2
         call random_seed(size=sizeseed)
         if ( IsMaster .and. debug ) then
           write(fid_alog,*) "size of seeds for random generator",sizeseed
         end if
!rand         call read_seed(0.0)

         ! initialize inherent growth parameterization
         call init_inherent_growth_par

         ! initialize osmotic coefficient LUT
         call init_osmo_par

         ! initialize normal distribution LUT
         call init_normal_lut

         ! initialize inverse normal distribution LUT
         call init_inv_normal_lut

      return
      end subroutine binmicrosetup_scale

      subroutine ini_aerosols_ctr(qapv,nz,npa,nba,nca &
           ,z,zz,denv,thp,piv,estbar,esitbar,HDL_switch)
!-------------------------------------------------------------------
! This is for KiD
!-------------------------------------------------------------------
        use scale_prc, only: &
           PRC_abort
        use maxdims
        use par_amps
        use com_amps
        use mod_amps_utility, only: &
             get_inact_tropic2,cdfnor
        implicit none

      integer nz,npa,nba,nca,ns,HDL_switch
      integer nin,i,j,k,iba,ica
      real(PS) :: qapv(npamx,nbamx,ncamx,*)
      real(MP_KIND) ::  z(*),zz(*) &
          ,denv(*),thp(*) &
          ,piv(*)
      real(DS) :: estbar(*),esitbar(*)
      real(PS) ::  AMCCN,AMIN,w,pi_chiarui,ptot,til,wt,es,esi,sw,si
!tmp,get_inact_tropic2
      real(PS) ::  c1,c2,c3,c4,z1,z2,z3,z4,a,b,V,r1,r2,r3

      real(RP) :: cv, ep, gama, gamai, cpr, rcp
      !
      real(8) :: zn
      real(8) :: p,q,mean,sd,bound
      integer :: status
      real(PS),parameter :: factIN=0.0e+0_PS

!
!     category 1: accumulation mode (CN)
!     category 2: accumulation mode (IN)
!     category 3: coarse mode (CCN)
!     category 4: coarse mode (CCN)
!
!     NOTE: Evaporation of hydrometeors transfers aerosol masses into
!           either category 1 or 2. category 3 and 4 are used for
!           activation only.
!

      cv=cp-r
      ep=r/rm
      gama=cp/cv
      gamai=1./gama
      cpr=cp/r
      rcp=r/cp

      do ica=1,nca
!     set up the mean mass for a AP particle
         V=4.188790205_PS*(ap_mean(ica)**3)
         den_apt(ica)=den_api(ica)/(1.0_PS-eps_ap(ica)*(1.0_PS-den_api(ica) &
              /den_aps(ica)))
         ap_mass(ica)=den_apt(ica)*V

!    calculate coef for getting total mass from total concentration
         if(dtype_a(ica)==3) then
!           lognormal distribution
            if(debug) write(fid_alog,*) "AP lognormal",ica
            coef_ap(ica)=4.188790205_PS*den_apt(ica) &
                 *exp(4.5_PS*ap_lnsig(ica)**2+3.0_PS*log(ap_mean(ica)))
         elseif(dtype_a(ica)==4) then
!           uniform distribution
            if(debug) write(fid_alog,*) "AP uniform",ica
            coef_ap(ica)=4.188790205_PS*den_apt(ica) &
                 *ap_mean(ica)**3
         end if
      end do
      !
      !    calculate the base and normalization factors for theta-PDF scheme
      !
      do ica=1,nca
        mean=ap_mean_cp(ica)
        sd=ap_sig_cp(ica)
        zn=0.0d+0
        if(debug) write(fid_alog,*) "AP contact parameter mean and sd:",ica,mean,sd
        if(mean<=0.0.or.sd<=0.0) then
           LOG_ERROR("ini_aerosols_ctr",*) "Need to specify contact parameter variables"
           call PRC_abort
        endif
        call cdfnor(1,p,q,zn,mean,sd,status,bound)
        cdf_cp_0(ica)=p
        zn=180.0d+0
        call cdfnor(1,p,q,zn,mean,sd,status,bound)
        cdf_cp_180m0(ica)=p-cdf_cp_0(ica)
        if(debug) write(fid_alog,*) "AP cdf of 0 and 180 deg",ica,cdf_cp_0(ica),cdf_cp_180m0(ica)
      enddo

      if (HDL_switch == 0)  return

!    CCN concentration is given in cm^(-3)
!    calculate the coefficients of exponential distribution on vertical
      c1=100.0_PS
      z1=0.0_PS
      c2=500.0_PS
      z2=10000.0_PS
      c3=500.0_PS
      z3=14000.0_PS
      c4=2.0_PS
      z4=20000.0_PS

!tmp      b=log(c1/c2)/(z2-z1)
!tmp      a=c1*(c1/c2)**(z1/(z2-z1))

!   assume the ratio of each distributions.
      r1=66.6_PS/(66.6_PS+133.0_PS+3.1_PS)
      r2=133.0_PS/(66.6_PS+133.0_PS+3.1_PS)
      r3=3.1_PS/(66.6_PS+133.0_PS+3.1_PS)

      if(debug) write(fid_alog,*) "CCN setup"
      do i=1,nz
         if(z(i)<=z2) then
            qapv(acon_q,1,1,i)=log10(c1)+(log10(c2)-log10(c1)) &
                         /(z2-z1)*(z(i)-z1)
            qapv(acon_q,1,1,i)=10.0_PS**qapv(acon_q,1,1,i)


         elseif(z(i)<=z3) then
            qapv(acon_q,1,1,i)=c2
         else
            qapv(acon_q,1,1,i)=log10(c3)+(log10(c4)-log10(c3)) &
                         /(z4-z3)*(z(i)-z3)
            qapv(acon_q,1,1,i)=10.0_PS**qapv(acon_q,1,1,i)
         endif

         qapv(acon_q,1,3,i)=qapv(acon_q,1,1,i)*r2
         qapv(acon_q,1,4,i)=qapv(acon_q,1,1,i)*r3
         qapv(acon_q,1,1,i)=qapv(acon_q,1,1,i)*r1

         if(debug) write(fid_alog,254) z(i),qapv(acon_q,1,1,i),qapv(acon_q,1,3,i) &
              ,qapv(acon_q,1,4,i) &
              ,qapv(acon_q,1,1,i)+qapv(acon_q,1,3,i)+qapv(acon_q,1,4,i)

      end do
254   format("z,CCN 1 (#/cm^3), CCN 2, CCN 3, CCN T",10ES15.6)
!
!    IN concentration is given in cm^(-3)
!    calculate the coefficients of exponential distribution on vertical
      c1=2.5e-3_PS
      z1=1000.0_PS
!     set background concentration and height
      c2=2.5e-3_PS
      z2=2000.0_PS

!tmp      b=log(c1/c2)/(z2-z1)
!tmp      a=c1*(c1/c2)**(z1/(z2-z1))
      if(debug) write(fid_alog,*) "IN setup 1"
      do k=1,nz
         pi_chiarui=piv(k)/cp
         ptot=pi_chiarui**cpr*1.e5_RP
         til=min(273.16_RP,thp(k)*pi_chiarui)
         i=max(1,min(int(til)-163,149))
         wt=max(min(til-(i+163),1.0_RP),0.0_RP)
         es=estbar(i)*(1.0_RP-wt)+estbar(i+1)*wt

!         qs=0.622_RP*es/max(ptot-es,es)
!         sw=ptot*qv(k)/(0.622_RP+qv(k))/es-1.0_RP

         J=MAX(1,MIN(int(til)-163,110))
         wt=MAX(MIN(til-FLOAT(J+163),1.0_RP),0.0_RP)
         esi=esitbar(J)*(1.-wt)+esitbar(J+1)*wt

!         si=es/esi*(sw+1.0)-1.0_RP
         si=es/esi*(0.0_RP+1.0_RP)-1.0_RP

         qapv(acon_q,1,2,k)=get_inact_tropic2(si,til)
         if(debug) write(fid_alog,255) z(k),qapv(acon_q,1,2,k),til-273.16_RP,si
      end do

255   format("z,IN (#/cm^3), T(C),si",4ES15.6)
!      nin=131

!     impose the dust layer in lower troposphere
!     peak elevation
!tmp      c1=1.0_RP
!tmp      z1=2000.0_RP
!     set background concentration and height
!tmp      c2=1.0_RP
!tmp      z2=3000.0_RP
!tmp      z3=5000.0_RP
!tmp      a=(c2-c1)/(z2-z1)**2
!tmp      write(fid_alog,*) "IN setup 2"
!tmp      do i=1,nzp
!tmp         if(Z(i)<z1) then
!tmp            QAP01D(2,1,2,i)=QAP01D(2,1,2,i)
!tmp     *          +(c1-0.0_RP)/(z1-0.0_RP)*(Z(i)-0.0_RP)+0.0_RP
!tmp         elseif(Z(i)<=z2) then
!tmp            QAP01D(2,1,2,i)=QAP01D(2,1,2,i)+c1
!tmp         elseif(Z(i)<=z3) then
!tmp            QAP01D(2,1,2,i)=QAP01D(2,1,2,i)
!tmp     *          +(0.0_RP-c2)/(z3-z2)*(Z(i)-z2)+c2
!tmp         end if
!tmp         write(fid_alog,255) z(i),QAP01D(2,1,2,i)
!tmp      end do

      ! IN concentration is given in cm^(-3)
!      open(1,file='avgIN.txt')
!      do i=1,nin
!         read(1,*) PIN(i),CIN(i)
!         PIN(i)=PIN(i)*1.0e+2
!         CIN(i)=CIN(i)
!      end do
!      close(1)

!      do i=1,ns
!         CINP(i)=1.0_RP
!         CINP(i)=0.0_RP
!         do j=1,nin-1
!c            if(pin(j)<ps(i).and.ps(i)<=pin(j+1)) then
!            if(pin(j+1)<ps(i).and.ps(i)<=pin(j)) then
!              w=(ps(i)-pin(j))/(pin(j+1)-pin(j))
!              CINP(i)=(1.0_RP-w)*CIN(j)+w*CIN(j+1)
!               exit
!            end if
!         end do
!      end do

!      do i=1,ns
!         CCCNP(i)=25.0_RP
!      end do

!     get interpolated values
!      CALL HTINT(NS,CCCNP,ZINP,NZP,CAP01D(1,1),Z)
!     get interpolated values
!      CALL HTINT(NS,CINP,ZINP,NZP,CAP01D(1,2),Z)


      iba=1
      do ica=1,nca
         do k=1,nz
            qapv(acon_q,1,ica,k)=qapv(acon_q,1,ica,k)/(denv(k)*1.0e-3_RP)
!     realation between total mass and concentration for log normal distribution
!     is given by M=N (4*PI/3)*den * exp(9*lnsig**2/2+3*ln(am))
            qapv(amt_q,1,ica,k)=qapv(acon_q,1,ica,k)*coef_ap(ica)
            if(level_comp>=4) then
              qapv(ams_q,iba,ica,k)=qapv(amt_q,iba,ica,k)* &
                                    eps_ap(ica)
            end if
         end do
      end do
      return
      end subroutine ini_aerosols_ctr

      subroutine ini_aerosols_w2(qapv,nz,npa,nba,nca &
           ,z,zz,denv,thp,piv,estbar,esitbar,HDL_switch)
!-------------------------------------------------------------------
! This is for KiD
!-------------------------------------------------------------------
        use scale_prc, only: &
           PRC_abort
        use maxdims
        use par_amps
        use com_amps
        use mod_amps_utility, only: &
           get_inact_tropic2, &
           cdfnor
        implicit none
      integer, intent(in)    ::  nz,npa,nba,nca,HDL_switch
      real(MP_KIND), intent(inout)    ::  qapv(npamx,nbamx,ncamx,*)
      real(MP_KIND), intent(in)       ::  z(*),zz(*) &
                                ,denv(*),thp(*) &
                                ,piv(*)
      real(DS) :: estbar(*),esitbar(*)
      integer ::  nin,i,j,k,iba,ica,ns
      real(PS) :: AMCCN,AMIN,w,pi_chiarui,ptot,til,wt,es,esi,sw,si
!tmp,get_inact_tropic2
      real(PS) :: c1,c2,c3,c4,z1,z2,z3,z4,a,b,V,r1,r2,r3

      real(RP) :: cv, ep, gama, gamai, cpr, rcp
      real,parameter :: factIN=0.0e+0
      !
      real(8) :: zn
      real(8) :: p,q,mean,sd,bound
      integer :: status

!
!     category 1: accumulation mode (CN)
!     category 2: accumulation mode (IN)
!     category 3: coarse mode (CCN)
!     category 4: coarse mode (CCN)
!
!     NOTE: Evaporation of hydrometeors transfers aerosol masses into
!           either category 1 or 2. category 3 and 4 are used for
!           activation only.
!

      cv=cp-r
      ep=r/rm
      gama=cp/cv
      gamai=1./gama
      cpr=cp/r
      rcp=r/cp


    if (HDL_switch == 0)  then
      do ica=1,nca
!     set up the mean mass for a AP particle
         V=4.188790205*(ap_mean(ica)**3.0)
         den_apt(ica)=den_api(ica)/(1.0-eps_ap(ica)*(1.0-den_api(ica) &
              /den_aps(ica)))
         ap_mass(ica)=den_apt(ica)*V

!    calculate coef for getting total mass from total concentration
         if(dtype_a(ica)==3) then
!           lognormal distribution
            coef_ap(ica)=4.188790205*den_apt(ica) &
                 *exp(4.5*ap_lnsig(ica)**2.0+3.0*log(ap_mean(ica)))
            if(debug) write(fid_alog,*) "AP lognormal",ica,coef_ap(ica)
         elseif(dtype_a(ica)==4) then
!           uniform distribution
            coef_ap(ica)=4.188790205*den_apt(ica) &
                 *ap_mean(ica)**3.0
            if(debug) write(fid_alog,*) "AP uniform",ica,coef_ap(ica)
         end if
      end do

      !
      !    calculate the base and normalization factors for theta-PDF scheme
      !
      do ica=1,nca
        mean=ap_mean_cp(ica)
        sd=ap_sig_cp(ica)
        zn=0.0d+0
        if(debug) write(fid_alog,*) "AP contact parameter mean and sd:",ica,mean,sd
        if(mean<=0.0.or.sd<=0.0) then
           LOG_ERROR("ini_aerosols_w2",*) "Need to specify contact parameter variables"
           call PRC_abort
        endif
        call cdfnor(1,p,q,zn,mean,sd,status,bound)
        cdf_cp_0(ica)=p
        zn=180.0d+0
        call cdfnor(1,p,q,zn,mean,sd,status,bound)
        cdf_cp_180m0(ica)=p-cdf_cp_0(ica)
        if(debug) write(fid_alog,*) "AP cdf of 0 and 180 deg",ica,cdf_cp_0(ica),cdf_cp_180m0(ica)
      enddo

    elseif (HDL_switch == 1)  then

!    CCN concentration is given in cm^(-3)
!    constant over the vertical
!    only one category here

      if(debug) write(fid_alog,*) "CCN setup"
      do i=1,nz
         if(N_ap_ini(1)/=0.0) then
           qapv(acon_q,1,1,i)=N_ap_ini(1)
         else
           qapv(acon_q,1,1,i)=50.0
         endif

         if(debug) write(fid_alog,254) z(i),qapv(acon_q,1,1,i)

      end do
254   format("z,CCN 1 (#/cm^3)",10ES15.6)
!
!    IN concentration is given in cm^(-3)
!    calculate the coefficients of exponential distribution on vertical
      c1=2.5e-3
      z1=1000.0
!     set background concentration and height
      c2=2.5e-3
      z2=2000.0

!tmp      b=log(c1/c2)/(z2-z1)
!tmp      a=c1*(c1/c2)**(z1/(z2-z1))
      if(debug) write(fid_alog,*) "IN setup 1"
      do k=1,nz
         pi_chiarui=piv(k)/cp
         ptot=pi_chiarui**cpr*1.e5_RP
         til=min(273.16_RP,thp(k)*pi_chiarui)
         i=max(1,min(int(til)-163,149))
         wt=max(min(til-(i+163),1.0_RP),0.0_RP)
         es=estbar(i)*(1.-wt)+estbar(i+1)*wt

!         qs=0.622_RP*es/max(ptot-es,es)
!         sw=ptot*qv(k)/(0.622+qv(k))/es-1.0_RP

         J=MAX(1,MIN(int(til)-163,110))
         wt=MAX(MIN(til-(J+163),1.0_RP),0.0_RP)
         esi=esitbar(J)*(1.-wt)+esitbar(J+1)*wt

!         si=es/esi*(sw+1.0)-1.0
         si=es/esi*(0.0+1.0)-1.0

         qapv(acon_q,1,2,k)=get_inact_tropic2(si,til)
         if(debug) write(fid_alog,255) z(k),qapv(acon_q,1,2,k),til-273.16,si,es,esi
      end do

255   format("z,IN (#/cm^3), T(C),si,es,esi",6ES15.6)
!      nin=131

!     impose the dust layer in lower troposphere
!     peak elevation
!tmp      c1=1.0_RP
!tmp      z1=2000.0_RP
!     set background concentration and height
!tmp      c2=1.0_RP
!tmp      z2=3000.0_RP
!tmp      z3=5000.0_RP
!tmp      a=(c2-c1)/(z2-z1)**2
!tmp      write(fid_alog,*) "IN setup 2"
!tmp      do i=1,nzp
!tmp         if(Z(i)<z1) then
!tmp            QAP01D(2,1,2,i)=QAP01D(2,1,2,i)
!tmp     *          +(c1-0.0_RP)/(z1-0.0_RP)*(Z(i)-0.0_RP)+0.0_RP
!tmp         elseif(Z(i)<=z2) then
!tmp            QAP01D(2,1,2,i)=QAP01D(2,1,2,i)+c1
!tmp         elseif(Z(i)<=z3) then
!tmp            QAP01D(2,1,2,i)=QAP01D(2,1,2,i)
!tmp     *          +(0.0_RP-c2)/(z3-z2)*(Z(i)-z2)+c2
!tmp         end if
!tmp         write(fid_alog,255) z(i),QAP01D(2,1,2,i)
!tmp      end do

      ! IN concentration is given in cm^(-3)
!      open(1,file='avgIN.txt')
!      do i=1,nin
!         read(1,*) PIN(i),CIN(i)
!         PIN(i)=PIN(i)*1.0e+2_RP
!         CIN(i)=CIN(i)
!      end do
!      close(1)

!      do i=1,ns
!         CINP(i)=1.0_RP
!         CINP(i)=0.0_RP
!         do j=1,nin-1
!c            if(pin(j)<ps(i).and.ps(i)<=pin(j+1)) then
!            if(pin(j+1)<ps(i).and.ps(i)<=pin(j)) then
!              w=(ps(i)-pin(j))/(pin(j+1)-pin(j))
!              CINP(i)=(1.0-w)*CIN(j)+w*CIN(j+1)
!               exit
!            end if
!         end do
!      end do

!      do i=1,ns
!         CCCNP(i)=25.0_RP
!      end do

!     get interpolated values
!      CALL HTINT(NS,CCCNP,ZINP,NZP,CAP01D(1,1),Z)
!     get interpolated values
!      CALL HTINT(NS,CINP,ZINP,NZP,CAP01D(1,2),Z)


      iba=1
      do ica=1,nca
         do k=1,nz
            qapv(acon_q,1,ica,k)=qapv(acon_q,1,ica,k)/(denv(k)*1.0e-3_RP)
!     realation between total mass and concentration for log normal distribution
!     is given by M=N (4*PI/3)*den * exp(9*lnsig**2/2+3*ln(am))
            qapv(amt_q,1,ica,k)=qapv(acon_q,1,ica,k)*coef_ap(ica)
            if(level_comp>=4) then
              qapv(ams_q,iba,ica,k)=qapv(amt_q,iba,ica,k)* &
                                    eps_ap(ica)
            end if
         end do
      end do
    endif
    return
  end subroutine ini_aerosols_w2


  subroutine ini_aerosols_const(nca)
    use scale_prc, only: &
       PRC_abort
    use maxdims
    use par_amps
    use com_amps
    use mod_amps_utility, only: &
!       get_inact_pamarcmip, &
       cdfnor
    implicit none

      integer, intent(in)    :: nca
      integer ::  ica
      real(PS) :: pi_chiarui,til
!tmp,get_inact_tropic2
      real(PS) :: c1,c2,c3,c4,z1,z2,z3,z4,a,b,V,r1,r2,r3

      real(RP) :: cv, ep, gama, gamai, cpr, rcp
      real(PS),parameter :: factIN=0.0e+0
      !
      real(8) :: zn
      real(8) :: p,q,mean,sd,bound
      integer :: status

      integer :: ierr
!
!     category 1: accumulation mode (CN)
!     category 2: accumulation mode (IN)
!     category 3: coarse mode (CCN)
!     category 4: coarse mode (CCN)
!
!     NOTE: Evaporation of hydrometeors transfers aerosol masses into
!           either category 1 or 2. category 3 and 4 are used for
!           activation only.
!

      cv=cp-r
      ep=r/rm
      gama=cp/cv
      gamai=1./gama
      cpr=cp/r
      rcp=r/cp


      do ica=1,nca
!     set up the mean mass for a AP particle
         V=4.188790205_PS*(ap_mean(ica)**3)
         den_apt(ica)=den_api(ica)/(1.0_PS-eps_ap(ica)*(1.0_PS-den_api(ica) &
              /den_aps(ica)))
         ap_mass(ica)=den_apt(ica)*V

!    calculate coef for getting total mass from total concentration
         if(dtype_a(ica)==3) then
!           lognormal distribution
            coef_ap(ica)=4.188790205_PS*den_apt(ica) &
                 *exp(4.5_PS*ap_lnsig(ica)**2+3.0_PS*log(ap_mean(ica)))
            if(debug) write(fid_alog,*) "AP lognormal",ica,coef_ap(ica)
         elseif(dtype_a(ica)==4) then
!           uniform distribution
            coef_ap(ica)=4.188790205_PS*den_apt(ica) &
                 *ap_mean(ica)**3
            if(debug) write(fid_alog,*) "AP uniform",ica,coef_ap(ica)
         end if
      end do

      !
      !    calculate the base and normalization factors for theta-PDF scheme
      !
      do ica=1,nca
        mean=ap_mean_cp(ica)
        sd=ap_sig_cp(ica)
        zn=0.0d+0
        if(debug) write(fid_alog,*) "AP contact parameter mean and sd:",ica,mean,sd
        if(mean<=0.0_PS.or.sd<=0.0_PS) then
          LOG_ERROR("ini_aerosols_const",*) "Need to specify contact parameter variables"
          call PRC_abort
        endif
        call cdfnor(1,p,q,zn,mean,sd,status,bound)
        cdf_cp_0(ica)=p
        zn=180.0d+0
        call cdfnor(1,p,q,zn,mean,sd,status,bound)
        cdf_cp_180m0(ica)=p-cdf_cp_0(ica)
        if(debug) write(fid_alog,*) "AP cdf of 0 and 180 deg",ica,cdf_cp_0(ica),cdf_cp_180m0(ica)
      enddo

      return
    end subroutine ini_aerosols_const

    subroutine ini_aerosols_profile_const(qapv,nz,nca,z,denv)
      use maxdims
      use par_amps
      use com_amps
      use mod_amps_utility, only: &
!         get_inact_pamarcmip, &
         cdfnor
      implicit none
      integer,  intent(in)    :: nz, nca
      real(MP_KIND), intent(inout) :: qapv(npamx,nbamx,ncamx,*)
      real(MP_KIND), intent(in)    :: z(*), denv(*)
      integer ::  i,k,ica,iba
      real(PS) :: pi_chiarui!,til
!tmp,get_inact_tropic2
      real(PS) :: c1,c2,c3,c4,z1,z2,z3,z4,a,b,V,r1,r2,r3

      real(RP) :: cv, ep, gama, gamai, cpr, rcp
      real(PS),parameter :: factIN=0.0e+0


!    CCN concentration is given in cm^(-3)
!    constant over the vertical
!    only one category here

      if(debug) write(fid_alog,*) "CCN setup"
      do i=1,nz
         qapv(acon_q,1,1,i)=N_ap_ini(1)
         qapv(acon_q,1,3,i)=N_ap_ini(3)
         qapv(acon_q,1,4,i)=N_ap_ini(4)

         if(debug) write(fid_alog,254) z(i),qapv(acon_q,1,1,i), qapv(acon_q,1,3,i), qapv(acon_q,1,4,i)

      end do
254   format("z, CCN 1, CCN 3, CCN 4(#/cm^3)",10ES15.6)

!
!    IN concentration is given in cm^(-3)
!    calculate the coefficients of exponential distribution on vertical
      c1=0.15e-3
      z1=1000.0
!     set background concentration and height
      c2=0.15e-3
      z2=2000.0

      if(debug) write(fid_alog,*) "IN setup 1"
      do k=1,nz
         qapv(acon_q,1,2,k)=N_ap_ini(2)
         !pi_chiarui=piv(k)/cp
         !til=min(273.16,thp(k)*pi_chiarui)
         !qapv(acon_q,1,2,k)=get_inact_pamarcmip(til)
         if(debug) write(fid_alog,255) z(k),qapv(acon_q,1,2,k)!, til
      end do

255   format("z,IN (#/cm^3), T(C)",3ES15.6)

      iba=1
      do ica=1,nca
         do k=1,nz
            qapv(acon_q,1,ica,k)=qapv(acon_q,1,ica,k)/(denv(k)*1.0e-3_PS)
!     realation between total mass and concentration for log normal distribution
!     is given by M=N (4*PI/3)*den * exp(9*lnsig**2/2+3*ln(am))
            qapv(amt_q,1,ica,k)=qapv(acon_q,1,ica,k)*coef_ap(ica)
            if(level_comp>=4) then
              qapv(ams_q,iba,ica,k)=qapv(amt_q,iba,ica,k)* &
                                    eps_ap(ica)
            end if
         end do
      end do

      return
    end subroutine ini_aerosols_profile_const


    subroutine ini_aerosols_mpace(nca)
      use scale_prc, only: &
         PRC_abort
      use maxdims
      use par_amps
      use com_amps
      use mod_amps_utility, only: &
         get_inact_tropic2, &
         cdfnor
      implicit none

      integer, intent(in) :: nca
      integer :: ica
      real(PS) :: AMCCN,AMIN,w,pi_chiarui,ptot,til,wt,es,esi,sw,si
!tmp,get_inact_tropic2
      real(PS) :: c1,c2,c3,c4,z1,z2,z3,z4,a,b,V,r1,r2,r3

      real(RP) :: cv, ep, gama, gamai, cpr, rcp
      real(PS),parameter :: factIN=0.0e+0
      !
      real(8) :: zn
      real(8) :: p,q,mean,sd,bound
      integer :: status

      integer :: ierr
!
!     category 1: accumulation mode (CN)
!     category 2: accumulation mode (IN)
!     category 3: coarse mode (CCN)
!     category 4: coarse mode (CCN)
!
!     NOTE: Evaporation of hydrometeors transfers aerosol masses into
!           either category 1 or 2. category 3 and 4 are used for
!           activation only.
!

      cv=cp-r
      ep=r/rm
      gama=cp/cv
      gamai=1./gama
      cpr=cp/r
      rcp=r/cp

      do ica=1,nca
!     set up the mean mass for a AP particle
         V=4.188790205_PS*(ap_mean(ica)**3)
         den_apt(ica)=den_api(ica)/(1.0_PS-eps_ap(ica)*(1.0_PS-den_api(ica) &
              /den_aps(ica)))
         ap_mass(ica)=den_apt(ica)*V

!    calculate coef for getting total mass from total concentration
         if(dtype_a(ica)==3) then
!           lognormal distribution
            if(debug) write(fid_alog,*) "AP lognormal",ica
            coef_ap(ica)=4.188790205_PS*den_apt(ica) &
                 *exp(4.5_PS*ap_lnsig(ica)**2.0+3.0*log(ap_mean(ica)))
         elseif(dtype_a(ica)==4) then
!           uniform distribution
            if(debug) write(fid_alog,*) "AP uniform",ica
            coef_ap(ica)=4.188790205_PS*den_apt(ica) &
                 *ap_mean(ica)**3
         end if
      end do
      !
      !    calculate the base and normalization factors for theta-PDF scheme
      !
      do ica=1,nca
        mean=ap_mean_cp(ica)
        sd=ap_sig_cp(ica)
        zn=0.0d+0
        if(debug) write(fid_alog,*) "AP contact parameter mean and sd:",ica,mean,sd
        if(mean<=0.0_PS.or.sd<=0.0_PS) then
            LOG_ERROR("ini_aerosols_mpace",*) "Need to specify contact parameter variables"
          call PRC_abort
        endif
        call cdfnor(1,p,q,zn,mean,sd,status,bound)
        cdf_cp_0(ica)=p
        zn=180.0d+0
        call cdfnor(1,p,q,zn,mean,sd,status,bound)
        cdf_cp_180m0(ica)=p-cdf_cp_0(ica)
        if(debug) write(fid_alog,*) "AP cdf of 0 and 180 deg",ica,cdf_cp_0(ica),cdf_cp_180m0(ica)
      enddo

      return
    end subroutine ini_aerosols_mpace


    subroutine ini_aerosols_profile_mpace(qapv,nz,nca,z,denv)
      use scale_prc, only: &
         PRC_abort
      use maxdims
      use par_amps
      use com_amps
      use mod_amps_utility, only: &
         get_inact_tropic2, &
         cdfnor
      implicit none

      integer, intent(in) ::  nz, nca
      real(MP_KIND), intent(inout) :: qapv(npamx,nbamx,ncamx,*)
      real(MP_KIND), intent(in) :: z(*), denv(*)
      integer :: i,k,ica,iba
      real(PS) :: AMCCN,AMIN,w,pi_chiarui,ptot,til,wt,es,esi,sw,si
!tmp,get_inact_tropic2
      real(PS) :: c1,c2,c3,c4,z1,z2,z3,z4,a,b,V,r1,r2,r3

      real(RP) :: cv, ep, gama, gamai, cpr, rcp
      real(PS),parameter :: factIN=0.0e+0


!    CCN concentration is given in cm^(-3)
!    calculate the coefficients of exponential distribution on vertical
! changed for sheba
      c1=1.8_PS
!      c2=72.2
      c2=72.2_PS
!      c2=650.0

!      z1=5000.0
!     set background concentration and height

!      z2=10500.0
!      b=log(c1/c2)/(z2-z1)
!      a=c1*(c1/c2)**(z1/(z2-z1))

      if(debug) write(fid_alog,*) "CCN setup"
      do i=1,nz
         qapv(acon_q,1,1,i)=c2
         qapv(acon_q,1,3,i)=c1
         qapv(acon_q,1,4,i)=0.0_PS

         if(debug) write(fid_alog,254) z(i),qapv(acon_q,1,1,i),qapv(acon_q,1,3,i) &
              ,qapv(acon_q,1,4,i) &
              ,qapv(acon_q,1,1,i)+qapv(acon_q,1,3,i)+qapv(acon_q,1,4,i)
      end do
! end changed for sheba
254   format("z,CCN 1 (#/cm^3), CCN 2, CCN 3, CCN T",10ES15.6)
!
!    IN concentration is given in cm^(-3)
!    calculate the coefficients of exponential distribution on vertical
!GCSS      c1=1.7e-3
!!!      c1=0.29e-3  ! Fridlind et al (2011)
      c1=0.16e-3_PS  ! Fridlind et al (2011)
      z1=1000.0_PS
!     set background concentration and height
!GCSS      c2=1.7e-3
!!!      c2=0.29e-3 ! Fridlind et al (2011)
      c2=0.16e-3_PS ! Fridlind et al (2011)
      z2=2000.0_PS
!sheba      b=log(c1/c2)/(z2-z1)
!sheba      a=c1*(c1/c2)**(z1/(z2-z1))
      if(debug) write(fid_alog,*) "IN setup 1"
      do i=1,nz
         qapv(acon_q,1,2,i)=c1
!sheba         QAP01D(acon_q,1,2,i)=a*exp(-b*Z(i))
         if(z(i)>=z2) then
            qapv(acon_q,1,2,i)=c2
         end if
         if(debug) write(fid_alog,255) z(i),qapv(acon_q,1,2,i)
      end do

255   format("z,IN (#/cm^3)",2ES15.6)
!      nin=131

!     impose the dust layer in lower troposphere
!     peak elevation
!sheba      c1=1.0
!sheba      z1=2000.0
!shebac     set background concentration and height
!sheba      c2=1.0
!sheba      z2=3000.0
!sheba      z3=5000.0
!sheba      a=(c2-c1)/(z2-z1)**2.0
!sheba      write(fid_alog,*) "IN setup 2"
!sheba      do i=1,nzp
!sheba         if(Z(i)<z1) then
!sheba            QAP01D(acon_q,1,2,i)=QAP01D(acon_q,1,2,i)
!sheba     *          +(c1-0.0)/(z1-0.0)*(Z(i)-0.0)+0.0
!sheba         elseif(Z(i)<=z2) then
!sheba            QAP01D(acon_q,1,2,i)=QAP01D(acon_q,1,2,i)+c1
!sheba         elseif(Z(i)<=z3) then
!sheba            QAP01D(acon_q,1,2,i)=QAP01D(acon_q,1,2,i)
!sheba     *          +(0.0-c2)/(z3-z2)*(Z(i)-z2)+c2
!sheba         end if
!sheba         write(fid_alog,255) z(i),QAP01D(acon_q,1,2,i)
!sheba      end do

      ! IN concentration is given in cm^(-3)
!      open(1,file='avgIN.txt')
!      do i=1,nin
!         read(1,*) PIN(i),CIN(i)
!         PIN(i)=PIN(i)*1.0e+2
!         CIN(i)=CIN(i)
!      end do
!      close(1)

!      do i=1,ns
!         CINP(i)=1.0
!         CINP(i)=0.0
!         do j=1,nin-1
!c            if(pin(j)<ps(i).and.ps(i)<=pin(j+1)) then
!            if(pin(j+1)<ps(i).and.ps(i)<=pin(j)) then
!              w=(ps(i)-pin(j))/(pin(j+1)-pin(j))
!              CINP(i)=(1.0-w)*CIN(j)+w*CIN(j+1)
!               exit
!            end if
!         end do
!      end do

!      do i=1,ns
!         CCCNP(i)=25.0
!      end do

!     get interpolated values
!      CALL HTINT(NS,CCCNP,ZINP,NZP,CAP01D(1,1),Z)
!     get interpolated values
!      CALL HTINT(NS,CINP,ZINP,NZP,CAP01D(1,2),Z)


      do ica=1,nca
         do k=1,nz
            qapv(acon_q,1,ica,k)=qapv(acon_q,1,ica,k)/(denv(k)*1.0e-3_PS)
!     realation between total mass and concentration for log normal distribution
!     is given by M=N (4*PI/3)*den * exp(9*lnsig**2/2+3*ln(am))
            qapv(amt_q,1,ica,k)=qapv(acon_q,1,ica,k)*coef_ap(ica)

            if(level_comp>=4) then
              qapv(ams_q,1,ica,k)=qapv(amt_q,1,ica,k)*eps_ap(ica)

!     activated IN cocnetration
              !qapv(acnc_q,1,ica,k)=qapv(acon_q,1,ica,k)*factIN
              if(debug) write(fid_alog,*) "qapv check",ica,k,qapv(acon_q,1,ica,k)&
                        ,qapv(amt_q,1,ica,k),qapv(ams_q,1,ica,k)!&
                        !,qap01d(acnc_q,1,ica,k)
            endif
         end do
      end do

      return
    end subroutine ini_aerosols_profile_mpace


    subroutine ini_aerosols_sheba(nca)
      use scale_prc, only: &
         PRC_abort
      use com_amps
      use par_amps
      use maxdims
      use mod_amps_utility, only: &
         cdfnor
      implicit none

      integer, intent(in) :: nca
      integer :: ica
      real(PS) :: AMCCN,AMIN,w
      real(PS) :: c1,c2,c3,c4,z1,z2,z3,z4,a,b,V,r1,r2,r3
      real(PS) ,parameter :: factIN=1.6e-4
      !
      real(8) :: zn
      real(8) :: p,q,mean,sd,bound
      integer :: status

      integer :: ierr
!
!     category 1: accumulation mode (CN)
!     category 2: accumulation mode (IN)
!     category 3: nucleation mode (CCN)
!     category 4: coarse mode (CCN)
!
!     NOTE: Evaporation of hydrometeors transfers aerosol masses into
!           either category 1 or 2. category 3 and 4 are used for
!           activation only.
!
      do ica=1,nca
!     set up the mean mass for a AP particle
         V=4.188790205_PS*(ap_mean(ica)**3)
         den_apt(ica)=den_api(ica)/(1.0_PS-eps_ap(ica)*(1.0_PS-den_api(ica) &
              /den_aps(ica)))
         ap_mass(ica)=den_apt(ica)*V

         if(debug) write(fid_alog,*) "ap_mass ini",ica,ap_mass(ica)

!    calculate coef for getting total mass from total concentration
         if(dtype_a(ica)==3) then
!           lognormal distribution
            if(debug) write(fid_alog,*) "AP lognormal",ica
            coef_ap(ica)=4.188790205_PS*den_apt(ica) &
                 *exp(4.5_PS*ap_lnsig(ica)**2+3.0_PS*log(ap_mean(ica)))
         elseif(dtype_a(ica)==4) then
!           uniform distribution
            if(debug) write(fid_alog,*) "AP uniform",ica
            coef_ap(ica)=4.188790205_PS*den_apt(ica) &
                 *ap_mean(ica)**3
         end if
      end do
      !
      !    calculate the base and normalization factors for theta-PDF scheme
      !
      do ica=1,nca
        mean=ap_mean_cp(ica)
        sd=ap_sig_cp(ica)
        zn=0.0d+0
        if(debug) write(fid_alog,*) "AP contact parameter mean and sd:",ica,mean,sd
        if(mean<=0.0.or.sd<=0.0) then
          LOG_ERROR("ini_aerosols_sheba",*) "Need to specify contact parameter variables"
          call PRC_abort
        endif
        call cdfnor(1,p,q,zn,mean,sd,status,bound)
        cdf_cp_0(ica)=p
        zn=180.0d+0
        call cdfnor(1,p,q,zn,mean,sd,status,bound)
        cdf_cp_180m0(ica)=p-cdf_cp_0(ica)
        if(debug) write(fid_alog,*) "AP cdf of 0 and 180 deg",ica,cdf_cp_0(ica),cdf_cp_180m0(ica)
      enddo

      return
      end subroutine ini_aerosols_sheba


      subroutine ini_aerosols_profile_sheba(qapv,nz,nca,z,denv)
      use com_amps
      use par_amps
      use maxdims
      implicit none

      integer, intent(in) :: nz, nca
      real(MP_KIND), intent(in) :: z(*), denv(*)
      real(MP_KIND), intent(inout) ::  qapv(npamx,nbamx,ncamx,*)
      integer :: i,k,ica,iba
      real(PS) :: AMCCN,AMIN,w
      real(PS) :: c1,c2,c3,c4,z1,z2,z3,z4,a,b,V,r1,r2,r3
      real(PS) ,parameter :: factIN=1.6e-4


!    CCN concentration is given in cm^(-3)
!    calculate the coefficients of exponential distribution on vertical
! changed for sheba
      c1=1.8
!      c2=72.2
      c2=350.0
!      c2=650.0

!      z1=5000.0
!     set background concentration and height

!      z2=10500.0
!      b=log(c1/c2)/(z2-z1)
!      a=c1*(c1/c2)**(z1/(z2-z1))

      if(debug) write(fid_alog,*) "CCN setup"
      do i=1,nz
         qapv(acon_q,1,1,i)=c2
         qapv(acon_q,1,3,i)=c1
         if (ncamx .ge. 4) then
            qapv(acon_q,1,4,i)=0.0_MP_KIND
            if(debug) write(fid_alog,254) z(i),qapv(acon_q,1,1,i),qapv(acon_q,1,3,i), &
                 qapv(acon_q,1,4,i), &
                 qapv(acon_q,1,1,i)+qapv(acon_q,1,3,i)+qapv(acon_q,1,4,i)
         else
            if(debug) write(fid_alog,254) z(i),qapv(acon_q,1,1,i),qapv(acon_q,1,3,i), &
                 qapv(acon_q,1,1,i)+qapv(acon_q,1,3,i)
         end if
      end do
! end changed for sheba
254   format("z,CCN 1 (#/cm^3), CCN 2, CCN 3, CCN T",10ES15.6)
!
!    IN concentration is given in cm^(-3)
!    calculate the coefficients of exponential distribution on vertical
!GCSS      c1=1.7e-3
!!!      c1=0.29e-3  ! Fridlind et al (2011)
      c1=0.15e-3_PS  ! Fridlind et al (2011)
      z1=1000.0_PS
!     set background concentration and height
!GCSS      c2=1.7e-3
!!!      c2=0.29e-3 ! Fridlind et al (2011)
      c2=0.15e-3_PS ! Fridlind et al (2011)
      z2=2000.0_PS
!sheba      b=log(c1/c2)/(z2-z1)
!sheba      a=c1*(c1/c2)**(z1/(z2-z1))
      if(debug) write(fid_alog,*) "IN setup 1"
      do i=1,nz
         qapv(acon_q,1,2,i)=c1
!sheba         QAP01D(acon_q,1,2,i)=a*exp(-b*Z(i))
         if(z(i)>=z2) then
            qapv(acon_q,1,2,i)=c2
         end if
         if(debug) write(fid_alog,255) z(i),qapv(acon_q,1,2,i)
      end do

255   format("z,IN (#/cm^3)",2ES15.6)
!      nin=131

!     impose the dust layer in lower troposphere
!     peak elevation
!sheba      c1=1.0
!sheba      z1=2000.0
!shebac     set background concentration and height
!sheba      c2=1.0
!sheba      z2=3000.0
!sheba      z3=5000.0
!sheba      a=(c2-c1)/(z2-z1)**2.0
!sheba      write(fid_alog,*) "IN setup 2"
!sheba      do i=1,nzp
!sheba         if(Z(i)<z1) then
!sheba            QAP01D(acon_q,1,2,i)=QAP01D(acon_q,1,2,i)
!sheba     *          +(c1-0.0)/(z1-0.0)*(Z(i)-0.0)+0.0
!sheba         elseif(Z(i)<=z2) then
!sheba            QAP01D(acon_q,1,2,i)=QAP01D(acon_q,1,2,i)+c1
!sheba         elseif(Z(i)<=z3) then
!sheba            QAP01D(acon_q,1,2,i)=QAP01D(acon_q,1,2,i)
!sheba     *          +(0.0-c2)/(z3-z2)*(Z(i)-z2)+c2
!sheba         end if
!sheba         write(fid_alog,255) z(i),QAP01D(acon_q,1,2,i)
!sheba      end do

      ! IN concentration is given in cm^(-3)
!      open(1,file='avgIN.txt')
!      do i=1,nin
!         read(1,*) PIN(i),CIN(i)
!         PIN(i)=PIN(i)*1.0e+2
!         CIN(i)=CIN(i)
!      end do
!      close(1)

!      do i=1,ns
!         CINP(i)=1.0
!         CINP(i)=0.0
!         do j=1,nin-1
!c            if(pin(j)<ps(i).and.ps(i)<=pin(j+1)) then
!            if(pin(j+1)<ps(i).and.ps(i)<=pin(j)) then
!              w=(ps(i)-pin(j))/(pin(j+1)-pin(j))
!              CINP(i)=(1.0-w)*CIN(j)+w*CIN(j+1)
!               exit
!            end if
!         end do
!      end do

!      do i=1,ns
!         CCCNP(i)=25.0
!      end do

!     get interpolated values
!      CALL HTINT(NS,CCCNP,ZINP,NZP,CAP01D(1,1),Z)
!     get interpolated values
!      CALL HTINT(NS,CINP,ZINP,NZP,CAP01D(1,2),Z)


      do ica=1,nca
         do k=1,nz
            qapv(acon_q,1,ica,k)=qapv(acon_q,1,ica,k)/(denv(k)*1.0e-3)
!     realation between total mass and concentration for log normal distribution
!     is given by M=N (4*PI/3)*den * exp(9*lnsig**2/2+3*ln(am))
            qapv(amt_q,1,ica,k)=qapv(acon_q,1,ica,k)*coef_ap(ica)

            if(level_comp>=4) then
              qapv(ams_q,1,ica,k)=qapv(amt_q,1,ica,k)*eps_ap(ica)

!     activated IN cocnetration
              !qapv(acnc_q,1,ica,k)=qapv(acon_q,1,ica,k)*factIN
              if(debug) write(fid_alog,*) "qapv check",ica,k,qapv(acon_q,1,ica,k)&
                        ,qapv(amt_q,1,ica,k),qapv(ams_q,1,ica,k)!&
                        !,qap01d(acnc_q,1,ica,k)
            endif
         end do
      end do

      return
    end subroutine ini_aerosols_profile_sheba

    function get_sat_vapor_pres_lk(phase, T, estbar, esitbar ) result(e_sat)
      use class_Thermo_Var, only: &
         get_sat_vapor_pres_lk_upper => get_sat_vapor_pres_lk
      implicit none
      integer,intent(in)  :: phase
      real(PS),intent(in) :: T
      real(DS),intent(in) :: estbar(150),esitbar(111)
      real(PS) :: e_sat
      e_sat = get_sat_vapor_pres_lk_upper(phase, T, estbar, esitbar)

      return
    end function get_sat_vapor_pres_lk

    subroutine cal_breakfragment(NRBIN,estbar,esitbar)
      use scale_prc, only: &
         PRC_abort
      use com_amps
      use class_AirGroup, only: &
         AirGroup, &
         make_AirGroup_2
      use class_Group, only: &
         Group, &
         make_group, &
         ini_group_MP
      use mod_amps_utility, only: &
         random_genvar
      use mod_amps_core, only: &
         cal_Coalescence_Efficiency, &
         cal_breakup_dis_LL, &
         diag_pq
      implicit none

  integer,intent(in) :: NRBIN

  ! saturation vapor lookup table
!tmp  real  :: estbar(150),esitbar(111)
  real(DS),intent(in)  :: estbar(150),esitbar(111)

  ! steady growth problem: phase over which saturation is calculated.
  ! phase 1: water, 2: ice
  integer    :: phase

  ! random generator vars (not used)
  type (random_genvar) :: rdsd

  ! steady growth problem: relative humidity
  real(PS)    :: RH

  ! define the vapor pressure
  real(PS)              :: e

  ! definition for reading command line
  character(len=20)      :: chbuf
  real(PS)               :: temp ! ambient temperature in Celcius


  real(PS),dimension(1) :: rv,den,pt,T,W,z


  type (Group)               :: liquid
  type (AirGroup)            :: steady


  ! coalescence efficiency
  real (PS)                 :: E_coal
  real (PS) :: D,D_L,D_S,S_T,S_C,DS_S,CKE
  ! low-diameter cut off related to the resolution of the experiments (cm)
  real(DS),parameter :: D_0=0.01

  integer                                :: i, j,n



  ! message from reality_check
!tmp  integer,pointer,dimension(:)   :: mes_rc
  integer,dimension(1)   :: mes_rc
  integer,dimension(1)   :: ID,JD,KD

  ! optional status check
  integer  :: var_Status,em

  ! coefficient: sqrt(2*pi), pi/6.0
  real(DS),parameter  :: coef1=2.506628275,coef2=0.523598776

  ! maximum number of elements for lookup table, bu_tmass and bu_fd
  integer :: i1d_pair_max,kk_max

  integer :: ierr

!tp  allocate(mes_rc(1))


  !if ( IsMaster ) then
  !  write(fid_alog,*) "Now calculating collisional breakup fragments."
  !end if

  bu_fd=0.0_PS
  bu_tmass=0.0_PS

  T=278.6795
  PT=850.0e+2_PS
  W=0.0
  phase=1
  RH=100.0_PS

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! construct liquid hydrometeor
!tmp  liquid = make_Group ( 1, 0.0, 1, NRBIN, &
!tmp       srat_r, sadd_r, minmass_r, 1,1.0,1.0,1.0,1.0,binbr)
  call make_Group (liquid,1, 0.0_PS, 1, NRBIN, &
       srat_r, sadd_r, minmass_r, 1,1.0_PS,den_aps(1),den_api(1),eps_ap(1),binbr)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! construt thermo_var object
  if ( IsMaster .and. debug ) then
     write(fid_alog,*) "ck estbar",estbar(1:10)
  end if

  steady=make_AirGroup_2(1,estbar,esitbar,RV,DEN,PT,T,W,phase,RH)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! initialize the mass and con of liquid
  call ini_group_MP(liquid)


  ! diagnose the property
  call diag_pq(liquid,steady,1,em, mes_rc,ID,JD,KD,rdsd,ihabit_gm_random &
              ,eps_ap(1),nu_aps,phi_aps,m_aps)


  ! find bin that has the minimum size for possible breakup.
  do i=1,NRBIN
     if(liquid%MS(i,1)%len>=D_0) then
        jmin_bk=i
        exit
     end if
  end do

  imin_bk=jmin_bk+1
  imax_bk=NRBIN
  jmax_bk=NRBIN-1

  i1d_pair_max=(imax_bk-1)-jmin_bk+1+(imax_bk-imin_bk)*(1+imax_bk-imin_bk)/2
  kk_max=i1d_pair_max*liquid%N_BIN

  if ( IsMaster .and. debug ) then
    write(fid_alog,*) "i1d_pair_max,kk_max",i1d_pair_max,kk_max

    write(fid_alog,*) "size of bu_tmass, bu_fd",size(bu_tmass),size(bu_fd)/2
  end if

  if(size(bu_tmass)<i1d_pair_max) then
    LOG_ERROR("cal_breakfragment",*) "Error.  Increase the bu_tmass array"
    call PRC_abort
  endif
  if(size(bu_fd)/2<kk_max) then
    LOG_ERROR("cal_breakfragment",*) "Error.  Increase the bu_fd array"
    call PRC_abort
  endif



  do i=imin_bk,imax_bk
     do j=jmin_bk,i-1
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! determine the coalescence efficiency
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        call cal_Coalescence_Efficiency(liquid,i,liquid,j,1,steady%TV(1),E_coal,&
             D_L,D_S,S_T,S_C,DS_S,CKE)

        if(CKE<=1.0e-20) cycle
!!c             call cal_breakup_dis(liquid,i,j,imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,bu_fd,&
!!c                                  1.0,D_L,D_S,S_T,S_C,dS_S,CKE)
        call cal_breakup_dis_LL(liquid,i,j,imin_bk,jmin_bk,bu_tmass,bu_fd,&
                                D_L,D_S,S_T,S_C,CKE)


     end do

     if(D_L*1.0e+2>0.4.and.D_L*1.0e+2<0.5) then
        call print_bufd(i,liquid,"check_bkup0p4.dat")
     elseif(D_L*1.0e+2>0.3.and.D_L*1.0e+2<0.4) then
        call print_bufd(i,liquid,"check_bkup0p3.dat")
     end if

  end do

!!c  write(fid_alog,*) "cal_break bu_fd"
!!c  write(fid_alog,*) bu_fd
!!c  write(fid_alog,*) "cal_break bu_tmass"
!!c  write(fid_alog,*) bu_tmass


  ! +++ deallocate memory used for the Group Object +++
!tmp  call delete_Group(liquid)
!tmp  call delete_AirGroup (steady)

!tmp  deallocate(mes_rc)
end subroutine cal_breakfragment

subroutine print_bufd(ii,g,infile)
  use com_amps
  use class_Group, only: &
     Group
  type (Group)        :: g
  character (len=*) :: infile
  integer           :: ii,jj,kk,k,i1d_pair, io_print
  real(PS) :: d2


  if ( IsMaster ) then
    io_print = IO_get_available_fid()
    open(io_print,file=infile)
    do jj=jmin_bk,ii-1

       i1d_pair=jj-jmin_bk+1+(ii-imin_bk)*(1+ii-imin_bk)/2
       do k=1,g%n_bin
          d2=(0.5*(binbr(k)+binbr(k+1))/coefpi6)**(1.0/3.0)
          kk=(i1d_pair-1)*g%n_bin+k
          if(debug) write(io_print,'(5es15.6)') g%ms(ii,1)%len,g%ms(jj,1)%len,d2,bu_fd(1,kk),bu_fd(2,kk)
       end do
    end do
    close(io_print)
  end if

end subroutine print_bufd


! AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS
! AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS
! AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS
! AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS
! AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS
! AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS


subroutine wave_imitator(W,NGRIDS,c_time,dZ,w0)
  implicit none
  real(8),intent(in) :: c_time
  integer,intent(in) :: NGRIDS
  real(MP_KIND), dimension(NGRIDS)               :: w
  real(PS),intent(in) :: dZ,w0
  !
  integer :: i
  real(PS) :: tau_cycle

  tau_cycle=3.141592653*dZ/w0
  ! reduce the total water based on a time scale constant
  do i=1,NGRIDS

    w(i)=w0*sin(2.0*PI*c_time/tau_cycle)
  end do
end subroutine wave_imitator
! AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS
! AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS AMPS

subroutine cal_molality(mo_h,mo_c,mo_r,npr,nbr,ncr,ngrids,XR)
  use par_amps
  implicit none
  integer,intent(in) :: NGRIDS,npr,nbr,ncr
  real(MP_KIND), dimension(NGRIDS)               :: mo_h,mo_c,mo_r
  real(MP_KIND), dimension(npr,nbr,ncr,ngrids)  :: XR
  real(MP_KIND) :: cc,cr,ch
  real(PS),parameter :: M_s=115.11
  !

  integer :: i,ibr,ibi,j
  real(PS),parameter :: cdmn=1.0e-4,cdmx=50.0e-4

  do i=1,NGRIDS
     mo_c(i)=0.0_RP
     mo_h(i)=0.0_RP
     mo_r(i)=0.0_RP
     cc=0.0_RP
     cr=0.0_RP
     ch=0.0_RP
     do ibr=1,nbr
        if( (XR(rmt_q,ibr,1,i)/XR(rcon_q,ibr,1,i)*6.0_RP/3.14_RP)**(1.0_RP/3.0_RP)<cdmx.and.&
            (XR(rmt_q,ibr,1,i)/XR(rcon_q,ibr,1,i)*6.0_RP/3.14_RP)**(1.0_RP/3.0_RP)>=cdmn) then
           cc=cc+XR(rcon_q,ibr,1,i)
           mo_c(i)=mo_c(i)+XR(rcon_q,ibr,1,i)*max(0.0_RP,1.0e+3_RP*XR(rmas_q,ibr,1,i)/M_s/(XR(rmt_q,ibr,1,i)-XR(rmat_q,ibr,1,i)))
        elseif( (XR(rmt_q,ibr,1,i)/XR(rcon_q,ibr,1,i)*6.0_RP/3.14_RP)**(1.0_RP/3.0_RP)<cdmn) then
           ch=ch+XR(rcon_q,ibr,1,i)
           mo_h(i)=mo_h(i)+XR(rcon_q,ibr,1,i)*max(0.0_RP,1.0e+3_RP*XR(rmas_q,ibr,1,i)/M_s/(XR(rmt_q,ibr,1,i)-XR(rmat_q,ibr,1,i)))
        else
           cr=cr+XR(rcon_q,ibr,1,i)
           mo_r(i)=mo_r(i)+XR(rcon_q,ibr,1,i)*max(0.0_RP,1.0e+3_RP*XR(rmas_q,ibr,1,i)/M_s/(XR(rmt_q,ibr,1,i)-XR(rmat_q,ibr,1,i)))

        endif
     enddo
     mo_c(i)=mo_c(i)/cc
     mo_h(i)=mo_h(i)/ch
     mo_r(i)=mo_r(i)/cr

  end do
end subroutine cal_molality
subroutine cal_initqvar(qt,qrb,qib,qh,qc,qr,qrt,ch,cc,cr,ci,radc,radi,ac1,ac2,ac3,&
                        am1,am2,am3,am_l,am_i,&
                        NGRIDS,npr,nbr,ncr,npi,nbi,nci,npa,nba,nca,rv,den,XR,XS,XA)
  use par_amps
  implicit none
  integer,intent(in) :: NGRIDS,npr,nbr,ncr,npi,nbi,nci,npa,nba,nca
  real(MP_KIND), dimension(NGRIDS)               :: RV,qt,qrb,qib,qh,qc,qr,ch,cc,cr,ci,den &
                              ,radc,radi,ac1,ac2,ac3,qrt &
                              ,am1,am2,am3,am_l,am_i
  real(MP_KIND), dimension(npr,nbr,ncr,ngrids)  :: XR
  real(MP_KIND), dimension(npi,nbi,nci,ngrids)  :: XS
  real(MP_KIND), dimension(npa,nba,nca,ngrids)  :: XA
  integer :: i,ibr,ibi,icr,ici,iba,ica
  real(PS),parameter :: cdmn=1.0e-4,cdmx=50.0e-4
  real(MP_KIND) :: dia

  do i=1,NGRIDS

    qr(i)=0.0_RP
    qc(i)=0.0_RP
    qh(i)=0.0_RP
    qrb(i)=0.0_RP
    qrt(i)=0.0_RP
    am_l(i)=0.0_RP
    do icr=1,ncr
      do ibr=1,nbr
        qrb(i)=qrb(i)+max(0.0_RP,XR(rmt_q,ibr,icr,i)-XR(rmat_q,ibr,icr,i))
        qrt(i)=qrt(i)+max(0.0_RP,XR(rmt_q,ibr,icr,i))
        am_l(i)=am_l(i)+max(0.0_RP,XR(rmat_q,ibr,icr,i)*den(i)*1.0e-3_RP)
        if(XR(rcon_q,ibr,icr,i)>1.0e-30_RP) then
          dia=(XR(rmt_q,ibr,icr,i)/XR(rcon_q,ibr,icr,i)*6.0_RP/3.14_RP)**(1.0_RP/3.0_RP)
          if( dia<cdmx.and.dia>=cdmn) then
             qc(i)=qc(i)+max(0.0_RP,XR(rmt_q,ibr,icr,i)-XR(rmat_q,ibr,icr,i))
          elseif( dia<cdmn) then
             qh(i)=qh(i)+max(0.0_RP,XR(rmt_q,ibr,icr,i)-XR(rmat_q,ibr,icr,i))
          else
             qr(i)=qr(i)+max(0.0_RP,XR(rmt_q,ibr,icr,i)-XR(rmat_q,ibr,icr,i))
          endif
        endif
      enddo
    enddo
    cr(i)=0.0_RP
    cc(i)=0.0_RP
    ch(i)=0.0_RP
    do icr=1,ncr
      do ibr=1,nbr
        if(XR(rcon_q,ibr,icr,i)>1.0e-30_RP) then
          dia=(XR(rmt_q,ibr,icr,i)/XR(rcon_q,ibr,icr,i)*6.0/3.14)**(1.0_RP/3.0_RP)
          if( dia<cdmx.and.dia>=cdmn) then
             cc(i)=cc(i)+max(0.0_RP,XR(rcon_q,ibr,icr,i)*DEN(i)*1.0e-3_RP)
          elseif( dia<cdmn) then
             ch(i)=ch(i)+max(0.0_RP,XR(rcon_q,ibr,icr,i)*DEN(i)*1.0e-3_RP)
          else
             cr(i)=cr(i)+max(0.0_RP,XR(rcon_q,ibr,icr,i)*DEN(i)*1.0e-3_RP)
          endif
        endif
      enddo
    enddo

    qib(i)=0.0_RP
    ci(i)=0.0_RP
    radi(i)=0.0_RP
    am_i(i)=0.0_RP
    do ici=1,nci
      do ibi=1,nbi
        qib(i)=qib(i)+max(0.0_RP,XS(imt_q,ibi,ici,i)-XS(imat_q,ibi,ici,i))
        am_i(i)=am_i(i)+max(0.0_RP,XS(imat_q,ibr,icr,i)*DEN(i)*1.0e-3_RP)
        if(XS(icon_q,ibi,ici,i)>1.0e-30_RP) then
          ci(i)=ci(i)+max(0.0_RP,XS(icon_q,ibi,ici,i)*DEN(i)*1.0e-3_RP)
          dia=(XS(ivcs_q,ibi,ici,i)/XS(icon_q,ibi,ici,i)*6.0_RP/3.14_RP)**(1.0_RP/3.0_RP)
          radi(i)=radi(i)+max(0.0_RP,XS(icon_q,ibi,ici,i)*DEN(i)*1.0e-3_RP)*dia
        endif
      enddo
    enddo
    radi(i)=radi(i)/max(ci(i),1.0e-30_RP)*0.5_RP

    ac1(i)=0.0_RP
    am1(i)=0.0_RP
    ica=1
    do iba=1,nba
      ac1(i)=ac1(i)+XA(acon_q,iba,ica,i)*DEN(i)*1.0e-3_RP
      am1(i)=am1(i)+XA(amt_q,iba,ica,i)*DEN(i)*1.0e-3_RP
    enddo
    ac2(i)=0.0_RP
    am2(i)=0.0_RP
    ica=2
    do iba=1,nba
      ac2(i)=ac2(i)+XA(acon_q,iba,ica,i)*DEN(i)*1.0e-3_RP
      am2(i)=am2(i)+XA(amt_q,iba,ica,i)*DEN(i)*1.0e-3_RP
    enddo
    ac3(i)=0.0_RP
    am3(i)=0.0_RP
    ica=3
    do iba=1,nba
      ac3(i)=ac3(i)+XA(acon_q,iba,ica,i)*DEN(i)*1.0e-3_RP
      am3(i)=am3(i)+XA(amt_q,iba,ica,i)*DEN(i)*1.0e-3_RP
    enddo

    qt(i)=rv(i)+qrb(i)+qib(i)
  end do

end subroutine cal_initqvar
subroutine cal_til(til,NGRIDS,thetail,pt)
  implicit none
  integer,intent(in) :: NGRIDS
  real(MP_KIND), dimension(NGRIDS)               :: TIL,PT,thetail
  !

  integer :: i,ibr,ibi

  real(MP_KIND) :: pit,theta

  ! this is the easiest parcel model
  do i=1,NGRIDS
     pit=cp*(pt(i)/p0)**(r/cp)
     til(i)=pit/cp*thetail(i)
  enddo
end subroutine cal_til
subroutine cal_thetail(thetail,NGRIDS,qr,qi,pt,T)
  implicit none
  integer,intent(in) :: NGRIDS
  real(MP_KIND), dimension(NGRIDS)               :: TIL,PT, T,qr,qi,thetail
  !

  integer :: i,ibr,ibi
  ! latent heat of condensation in J/kg
  real(MP_KIND), parameter             :: L_e = 2.5e+6
  ! latent heat of sublimation in J/kg
  real(MP_KIND), parameter             :: L_s = 2.8337e+6
  ! gas constant of air
  real(MP_KIND), parameter             :: Ra = 287.0
  real(MP_KIND),parameter :: g=9.8
  real(MP_KIND),parameter :: cp=1004.0

  real(MP_KIND) :: pit,theta

  ! this is the easiest parcel model
  do i=1,NGRIDS
     pit=cp*(pt(i)/p0)**(Ra/cp)
     theta=T(i)*cp/pit
     thetail(i)=theta/(1.0_RP+(L_e*qr(i)+L_s*qi(i))/(cp*max(T(i),253.0_RP)))
  enddo
end subroutine cal_thetail

subroutine set_ampslog
  use scale_prc, only: &
     PRC_abort, &
     PRC_myrank, &
     PRC_IsMaster
  use com_amps

  implicit none
  character(len=256) :: word
  integer :: ierr
!org   use module_mpinms, only: mpi_getfid &
!org                              ,mpi_create_fname

!org  call mpi_create_fname(fname_ampslog,'amps_','.log')
!org  fid_alog=mpi_getfid(111)

  !call MPI_COMM_rank(MPI_COMM_WORLD, myrank, ierr)
  write(word,'(I0)') PRC_myrank
  !IsMaster = myrank == 0
  IsMaster = PRC_IsMaster

  if ( debug ) then
     fname_ampslog='amps_000000_'//trim(word)//'.log'
     !fname_ampslog='amps_000000.log'
     string_fid_amps='000000_'//trim(word)
     fid_alog = IO_get_available_fid()
     open( unit   = fid_alog, &
           file   = fname_ampslog, &
           form   = 'formatted', &
           iostat = ierr )
     if ( ierr /= 0 ) then
        LOG_ERROR("set_ampslog",*) "AMPS log file cannot be open, check! ", trim(fname_ampslog)
        call PRC_abort
     end if
  end if

end subroutine set_ampslog

end module mod_amps_lib
