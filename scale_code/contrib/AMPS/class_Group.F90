#include "scalelib.h"
!OCL SERIAL
MODULE class_Group
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! class_group.f90
!
! version 4.0
!
! The object Group may consist of mass bins, shape bins, or
! aerosol bins.
!
! Also, the object Group may represent one category of
! Bulk Water Parameterization.
! The spectrum is devided into bins and dealt with the same way as bin model.
! This version passes the mass and volume components.
!
! written by Tempei Hashino
!
! All units are CGS.
! ergs = g cm^2/s^2
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use scale_io
  use maxdims
  use mod_amps_const
  use acc_amps
  use par_amps
  use com_amps, only: &
     debug
  use class_Mass_Bin, only: &
     Mass_Bin
  use class_AirGroup, only: &
     AirGroup
  use class_Thermo_Var, only: &
     Thermo_Var
  use class_Ice_Shape, only: &
     Ice_Shape
  implicit none
  private

  public :: make_Group
  public :: make_col_lut
  public :: ini_group_all
  public :: ini_tendency
  public :: ini_tendency_ag
  public :: ini_group_MP
  public :: update_group_all
  public :: update_modelvars_all
  public :: update_mesrc
  public :: get_totalmass
  public :: print_tendency
  public :: print_tendency_ap
  public :: cal_surface_temp2_vec
  public :: cal_capacitance_vec
  public :: cal_ventilation_coef_vec
  public :: cal_terminal_vel_vec
  public :: cal_wterm_vel_v3_vec
  public :: cal_growth_mode_vec
  public :: cal_coef_vapdep_ap_vec
  public :: cal_coef_vapdep2_vec
  public :: cal_needgive
  public :: cal_den_aclen_vec

  public :: Group
  public :: col_lut_aux
  public :: vap_igp_aux

  TYPE Group
     ! memory is aligned as declared below for use in common block
     sequence
!!c     private

     ! mass object which has property of each bin
     ! first argument : bin number
     ! second argument : grid number
     !type (Mass_Bin), dimension(mxnbin,LMAX)       :: MS
     type (Mass_Bin),  allocatable :: MS(:,:)

     ! ice shape object
     ! first argument : bin number
     ! second argument : grid number
     !type (Ice_Shape), dimension(mxnbin,LMAX)       :: IS
     type (Ice_Shape), allocatable :: IS(:,:)

     ! group mark for each grid
     !   indicates if the concentration, or/and mass variables are fixed by
     !   reality_check routine.
     !integer,dimension(LMAX)         :: mark_cm    ! 4*100
     integer, allocatable          :: mark_cm(:)    ! 4*100

     ! error mark for each grid
     !   indicates any errors occuring in the schemes
     !    0 : no error
     !    1 : activation scheme not converging
     !integer,dimension(LMAX)         :: mark_er    ! 4*100
     integer, allocatable          :: mark_er(:)    ! 4*100

     ! boundaries of bins in the ascending order
     !real(PS), dimension(mxnbinb)            :: binb   ! 4*max(30,20)
     real(PS), allocatable         :: binb(:)   ! 4*max(30,20)

     ! coefficients defining bin boundaries binb(i+1)=a*binb(i)+b
     real(PS)        :: a_b,b_b

     ! minimum bin boundary
     real(PS)            :: mbinb


     ! time step in second
     real(PS)            :: dt

     ! token indicating phase of water, ice, or aerosol in the grid cell
     !    Three types for bin model
     !    water        : token = 1
     !    ice          : token = 2
     !    aerosol      : token = 3
     !
     !    Two types for bulk water parameterization
     !    cloud droplets   : token = 11
     !    rain             : token = 12
     integer           :: token


     ! total number of mass component
     integer         :: N_masscom

     ! total number of mass to describe the vapor deposition
     integer         :: N_axis

     ! total number of mass to describe the riming process
     integer         :: N_rimemass

     ! total number of mass to describe the aggregation
     integer         :: N_aggmass

     ! total number of mass to describe an aerosol per hydrometeor
     integer         :: N_apmass

     ! total number of mass to describe the melt water
     integer         :: N_meltmass

     ! total number of mass to describe the mass by freezing supercooled liquids.
     integer         :: N_frozenmass

     ! total number of non-mass variables
     integer         :: N_nonmass

     ! total number of volume variables
     integer         :: N_vol

     ! total number of processes for tendency
     integer         :: N_tendpros

     ! total number of bins
     integer         :: N_BIN

     ! total number of boundaries
     integer         :: N_binb


     ! concentration on the bin boundary
     ! If you use 3 parameter distribution in a bin, this has to
     ! be given
!!c     real(PS), pointer, dimension(:)            :: con_onbinb

     ! +++ type of distribution for a bin +++
     !     0 : constant-size disbtribution
     !     1 : linear distribution (default)
     !     2 : Marshall-Palmer distribution
     integer                                    :: org_dtype

     ! number of grids that have microphyscis in it.
     integer     :: L

     integer :: align_Group ! CHIARUI

  endtype Group

  type col_lut_aux
    sequence
    real(PS) :: xs,dx,ys,dy
    integer :: nr,nc
  end type col_lut_aux

  ! Inherent Growth parameterization
  integer,parameter :: nok_max=23
  type vap_igp_aux
    sequence
    ! number of knots
    integer :: nok
    integer :: align_vap_igp_aux ! CHIARUI
    ! parameters for cubic approximation
    real(PS),dimension(nok_max)   :: x
    real(PS),dimension(nok_max,4) :: a,b
  end type vap_igp_aux


CONTAINS

  ! ************* Optional Constructor for a group type **********************
  subroutine make_Group(g, dL, ddt, tok, n_bin, dab, dbb, dmbinb, dtype,&
       dden,dden_as,dden_ai,deps_map,dbinb)
    use class_Mass_Bin, only: &
       make_mass_bin

    integer, intent(in) :: dL, tok, n_bin, dtype
!tmp    integer, optional,intent(in) :: dNtp
    real(PS),intent(in)    :: ddt, dab, dbb, dmbinb,dden,dden_as,dden_ai,deps_map
!tmp    real(PS),optional,dimension(n_bin+1) :: dbinb
    real(PS),dimension(*) :: dbinb
    type (Group)  :: g

    ! The boundary of the bin binb_i is given by
    ! binb_i = a_b * binb_{i-1} + b_b
    !real(PS)            :: a_b, b_b
    integer :: i,j


    allocate( g%MS(mxnbin,LMAX) )
    allocate( g%IS(mxnbin,LMAX) )
    allocate( g%mark_cm(LMAX) )
    allocate( g%mark_er(LMAX) )
    allocate( g%binb(mxnbinb) )

    ! definition of total number of grids point with microphysics in it.
    g%L = dL

    ! definition of time step
    g%dt = ddt

    ! definition of token
    g%token = tok

    ! definition of parameters defining bin property
    g%N_BIN = n_bin
    g%N_binb = n_bin + 1

    g%a_b = dab
    g%b_b = dbb
    g%mbinb = dmbinb


    ! +++ definition of concentration on the boundary of bins +++
    !     allocate memory to type of distribution +++
!!c    allocate( g%con_onbinb(g%N_binb), stat = var_Status)
!!c    if(var_Status /= 0 ) stop "Memory not available for con_onbinb &
!!c         in class_group"
!!c    if( present(conb)) then; g%con_onbinb = conb
!!c    else;                  g%con_onbinb   = 0.0_PS
!!c    end if


    ! +++ definition of distribution
    g%org_dtype = dtype

    ! +++ memory allocation +++
!!c    call allocate_group(g)


!!c    write(*,*) "make_group",dL,ddt,tok,n_bin,dab,dbb,dmbinb,dtype
!!c    write(*,*) "makegrou2",dden,dden_as,dden_ai,deps_map
!!c    write(*,*) "makegroup3",dbinb

    ! +++ calculation of bin boundaries +++
       do i = 1, g%N_binb
          g%binb(i) = dbinb(i)
!!c          write(*,*) "gbinb",i,g%binb(i)
       end do
!tmp    if( present(dbinb)) then
!tmp       do i = 1, g%N_binb
!tmp          g%binb(i) = dbinb(i)
!tmp!!c          write(*,*) "gbinb",i,g%binb(i)
!tmp       end do
!tmp    else
!tmp       g%binb(1) = g%mbinb
!tmp       do i = 2, g%N_binb
!tmp          g%binb(i) = g%a_b*g%binb(i-1) + g%b_b
!tmp       end do
!tmp    end if
    ! +++ definition of number of process for tendency +++
       g%N_tendpros = 12
!tmp    if( present(dNtp)) then
!tmp       g%N_tendpros = dNtp
!tmp    else
!tmp       g%N_tendpros = 12
!tmp    end if

    if( g%token == 1 ) then
      ! +++ construct the liquid objects +++
      g%N_axis = 0
      g%N_rimemass = 0
      g%N_aggmass = 0
      g%N_apmass=3
      g%N_meltmass=0
      g%N_masscom=g%N_apmass
      g%N_vol = 0
      g%N_nonmass = g%N_vol

      call initialize

    else if( g%token == 2 ) then
      ! +++ construct the ice objects +++
      g%N_axis = 5
      g%N_rimemass = 1
      g%N_aggmass = 1
      g%N_apmass=3
      g%N_meltmass=1
      g%N_frozenmass=1
      g%N_masscom = g%N_rimemass+g%N_aggmass+1+g%N_apmass+g%N_meltmass+g%N_frozenmass
      g%N_vol = 0
      g%N_nonmass = 1 + g%N_axis + 1 + g%N_vol

      call initialize

      call initialize_shape

    else if( g%token == 3 ) then
      ! +++ construct the aerosol objects +++
      g%N_axis = 0
      g%N_rimemass = 0
      g%N_aggmass = 0
      g%N_apmass=2
      g%N_meltmass=0
      g%N_masscom=g%N_apmass
      g%N_vol = 0
      g%N_nonmass = g%N_vol

      call initialize

    else if( g%token == 11 ) then
       ! +++ construct cloud droplets category +++
       g%N_axis = 0
       g%N_rimemass = 0
       g%N_aggmass = 0
       g%N_apmass=2
       g%N_meltmass=0
       g%N_masscom=g%N_apmass
       g%N_vol = 0
       g%N_nonmass = 0


       do j=1,g%L
          do i=1,g%N_BIN
             ! The default is 300 spherical cloud drops with 1 um radius.
!tmp             g%MS(i,j) = make_Mass_Bin( g%binb(i), g%binb(i+1), g%org_dtype, &
             call make_Mass_Bin( g%MS(i,j), g%binb(i), g%binb(i+1), g%org_dtype, &
                  g%N_axis, g%N_rimemass, g%N_aggmass, g%N_apmass, g%N_meltmass, &
                  g%N_vol, g%N_tendpros, &
                  300.0_PS, 4.18879e-12_PS, 0.0_PS, 0.0_PS, &
                  2.0e-4_PS, 1.0_PS, dden_as, dden_ai,&
                  1.0e-4_PS, 1.0e-4_PS, 0.0_PS, 0.0_PS,&
                  deps_map,&
                  0.0_PS, 0.0_PS, 0.0_PS, 0.0_PS,&
                  0.0_PS, 0.0_PS, 0.0_PS, 0.0_PS,&
                  0)
          end do
       end do

    else if( g%token == 12 ) then
       ! +++ construct rain category +++
       g%N_axis = 0
       g%N_rimemass = 0
       g%N_aggmass = 0
       g%N_apmass=2
       g%N_meltmass=0
       g%N_masscom=g%N_apmass
       g%N_vol = 0
       g%N_nonmass = 0
       g%N_apmass=2

       do j=1,g%L
          do i=1,g%N_BIN
             ! The default is a cloud drops with 1 um radius.
!tmp             g%MS(i,j) = make_Mass_Bin( g%binb(i), g%binb(i+1), g%org_dtype, &
             call make_Mass_Bin( g%MS(i,j), g%binb(i), g%binb(i+1), g%org_dtype, &
                  g%N_axis, g%N_rimemass, g%N_aggmass, g%N_apmass, g%N_meltmass,&
                  g%N_vol, g%N_tendpros, &
                  1.0_PS, 4.18879e-12_PS, 0.0_PS, 0.0_PS, &
                  2.0e-4_PS, 1.0_PS,dden_as, dden_ai,&
                  1.0e-4_PS, 1.0e-4_PS, 0.0_PS, 0.0_PS,&
                  deps_map,&
                  0.0_PS, 0.0_PS, 0.0_PS, 0.0_PS,&
                  0.0_PS, 0.0_PS, 0.0_PS, 0.0_PS,&
                  0)
          end do
       end do

    end if

    g%mark_cm = 0
    g%mark_er = 0




  contains
    subroutine initialize
      integer  :: i,j,k,ij,l

!      do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN

          g%MS(i,j)%lmt = g%binb(i)
          g%MS(i,j)%umt = g%binb(i+1)

          g%MS(i,j)%con = 0.0_PS
          g%MS(i,j)%len = 0.0_PS
          g%MS(i,j)%den = dden
          g%MS(i,j)%den_as = dden_as
          g%MS(i,j)%den_ai = dden_ai

          g%MS(i,j)%a_len = 0.0_PS
          g%MS(i,j)%c_len = 0.0_PS

          g%MS(i,j)%semi_a = 0.0_PS
          g%MS(i,j)%semi_c = 0.0_PS

          g%MS(i,j)%eps_map = deps_map
          g%MS(i,j)%vtm = 0.0_PS

          g%MS(i,j)%tmp = 0.0_PS
          g%MS(i,j)%fv = 0.0_PS

          g%MS(i,j)%fkn = 0.0_PS

          g%MS(i,j)%fac = 0.0_PS
          g%MS(i,j)%fh = 0.0_PS

          g%MS(i,j)%Nre = 0.0_PS
          g%MS(i,j)%CAP = 0.0_PS

          g%MS(i,j)%CAP_hex=0.0_PS

          g%MS(i,j)%e_sat=0.0_PS

          g%MS(i,j)%dis_type = g%org_dtype

          g%MS(i,j)%mean_mass=0.0_PS

          g%MS(i,j)%r_act=0.0_PS
          g%MS(i,j)%r_crt=0.0_PS

          g%MS(i,j)%th00_cp=0.0_PS

          g%MS(i,j)%mark = 0

          g%MS(i,j)%inmlt=0

          g%MS(i,j)%inevp=0

          g%MS(i,j)%Ldmassdt(1) = 0.0_PS
          g%MS(i,j)%Ldmassdt(2) = 0.0_PS
          g%MS(i,j)%coef(1)=0.0_PS
          g%MS(i,j)%coef(2)=0.0_PS
      enddo
      enddo

!      do k=1,1+g%N_masscom
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN
      do k=1,1+g%N_masscom
         g%MS(i,j)%mass(k)=0.0_PS
      enddo
      enddo
      enddo

!      do k=1,g%N_vol
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN
      do k=1,g%N_vol
         g%MS(i,j)%vol(k)=0.0_PS
      enddo
      enddo
      enddo

!      do k=1,g%N_tendpros
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN
      do k=1,g%N_tendpros
         g%MS(i,j)%dcondt(k) = 0.0_PS
      enddo
      enddo
      enddo
!      do l=1,1+g%N_masscom
!      do k=1,g%N_tendpros
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN
      do k=1,g%N_tendpros
      do l=1,1+g%N_masscom
         g%MS(i,j)%dmassdt(l,k) = 0.0_PS
      enddo
      enddo
      enddo
      enddo

!      do l=1,g%N_nonmass
!      do k=1,g%N_tendpros
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN
      do k=1,g%N_tendpros
      do l=1,g%N_nonmass
         g%MS(i,j)%dvoldt(l,k) = 0.0_PS
      enddo
      enddo
      enddo
      enddo

    end subroutine initialize

    subroutine initialize_shape
      integer  :: i,j

!      do ij=1,g%N_BIN*g%L
!        j=(ij-1)/g%N_BIN+1
!        i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN

        g%IS(i,j)%a = 0.0_PS
        g%IS(i,j)%d = 0.0_PS
        g%IS(i,j)%c = 0.0_PS
        g%IS(i,j)%r = 0.0_PS
        g%IS(i,j)%e = 0.0_PS
        g%IS(i,j)%ag = 0.0_PS
        g%IS(i,j)%cg = 0.0_PS
        g%IS(i,j)%n_exice = 0.0_PS

        ! +++ calculate the volume +++
        g%IS(i,j)%V_ic = 0.0_PS
        ! circumscribing volume is calculated later.
        g%IS(i,j)%V_cs = 0.0_PS

        g%IS(i,j)%semi_aip = 0.0_PS
        g%IS(i,j)%semi_cip = 0.0_PS
        g%IS(i,j)%den_ip = 0.0_PS
        g%IS(i,j)%den_ic = 0.0_PS
        g%IS(i,j)%V_csw = 0.0_PS
        g%IS(i,j)%phi_ic = 0.0_PS
        g%IS(i,j)%psi_ic = 0.0_PS
        g%IS(i,j)%gam_ic = 0.0_PS
        g%IS(i,j)%eta_ic = 0.0_PS
        g%IS(i,j)%phi_cs = 0.0_PS
        g%IS(i,j)%q_e = 0.0_PS


        ! +++ categorize the ice hydrometeor +++
        g%IS(i,j)%habit = 0
        g%IS(i,j)%sh_type = 0
        g%IS(i,j)%growth_mode=0
        g%IS(i,j)%init_growth=0

        g%IS(i,j)%is_mod(1) = 0
        g%IS(i,j)%is_mod(2) = 0
      enddo
      enddo

    end subroutine initialize_shape

  end subroutine make_Group

  function make_col_lut(nr_in,nc_in,xs_in,ys_in,dx_in,dy_in) &
       result (c)
    implicit none
    integer, intent(in) :: nr_in,nc_in
    real(PS),intent(in) :: xs_in,ys_in,dx_in,dy_in
    type (col_lut_aux)  :: c

    c%nr = nr_in
    c%nc = nc_in
    c%xs = xs_in
    c%dx = dx_in
    c%ys = ys_in
    c%dy = dy_in
  end function make_col_lut


  subroutine ini_group_all(level, ag,gr,gs,ga, con_c,&
       L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,NATYPE,NABIN,NACAT,&
       XC,XR,XS,XA,flagp_c,flagp_r,flagp_s,flagp_a,coef_ap,eps_ap0,den_apt0,den_aps0,den_api0, &
       ID,JD,KD)
    use scale_prc, only: &
       PRC_abort
    use mod_amps_utility, only: &
       get_len_s1, &
       get_len_c2a, &
       get_len_s3
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! initialize the group
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! level of complexity
    integer, intent(in)  :: level
    type (Group), intent(inout)  :: gr,gs
!!c    type (Group), pointer, dimension(:) :: ga
    type (AirGroup), intent(in)   :: ag
    integer,intent(in)  :: L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,NATYPE,NABIN,NACAT
    real(MP_KIND) :: XC(*),XR(NRTYPE,NRBIN,NRCAT,*),XS(NSTYPE,NSBIN,NSCAT,*),XA(NATYPE,NABIN,NACAT,*)
    real(PS) :: coef_ap(*),eps_ap0(*),den_apt0(*),den_aps0(*),den_api0(*)
    integer :: ID(*),JD(*),KD(*)
    type (Group), dimension(NACAT) :: ga
    real(PS) :: con_c
    integer,intent(in)  :: flagp_c,flagp_r,flagp_s,flagp_a
    integer       :: i,j,k,m
    integer :: var_status
!    integer,dimension(mxnbin*LMAX) :: ierror
!!!    integer,dimension(:),allocatable :: ierror

    ! minimum, and maximum volume of circumscribing sphere that corresponds to 1 um of
    ! hex ice crystal and sphere of 10cm radius
    real(PS),parameter :: V_csmin_hex=1.184768784d-11,V_csmax=4188.79020478639D0
!!c    real,parameter :: V_csmin_hex=1.184768784e-11,V_csmax=4188.79020478639e+10
    real(PS) :: V_csmin2,V_csmin3
    real(PS) :: tmass1,tmass2,tmass3

    ! minimum mass of ice crystal with 1 um radius
    real(PS),parameter :: m_icmin=4.763209003d-12
    ! density minimum with circumscribing cylinder to hexagonal ice
    real(PS),parameter :: den_max=0.758088258D0
    real(PS) :: old_con,phi,psi,pag,pcg,con1
    real(PS) :: a_min,c_min,d_min,ag_min,cg_min
    ! maximum and minimum aspect ratio
    real(PS),parameter  :: phi_max=2.0d+1,phi_min=5.0d-3

    ! minimum possible ratio of water mass on the cloud or rain drop.
    real(PS),parameter :: minr_watermass=1.0d-2

    ! mass fraction of aerosol particles to total mass of hydrometeors and
    ! mass fraction of soluble aerosol particles to total aerosol mass
    ! Those are used to maintain positive aerosol masses, and use fraction
    ! of smaller bin.
    real(PS) :: fapt_r,faps_r,fapt_s,faps_s

    real(PS) :: m_lmt
    !     The minimum possible background concentration for large particle
    !     is assumed to be 1.0e-5 cm^-3
    real(PS),parameter :: n_lmt_ap=1.0d-5
!parcel model    real(PS),parameter :: n_lmt_ap=1.0e-15
    !     The minimum radius possible for accumulation particles
    real(PS),parameter :: r3_lmt=1.0d-18
    !     The minimum radius possible for nucleation particles
    real(PS),parameter :: r3_lmt_nuc=1.0d-21

    ! minimum possible fractions
    real(PS),parameter :: min_fapt_r=1.0d-18,min_faps_r=1.0d-5&
!!c    real(PS),parameter :: min_fapt_r=1.0e-18_PS,min_faps_r=0.0_PS&
         ,min_fapt_s=1.0d-18,min_faps_s=0.0D0,sep_faps=1.0d-5


    ! minimum limit for prediction
!org    real(PS),parameter :: m_lmt_ap=1.0e-25
    real(PS),parameter :: m_lmt_ap=1.0d-30

    !
    real(PS),parameter :: mx_aprat=0.95D0

    real(PS) :: max_rad
    ! factor
    real(PS),parameter :: fct1=0.1,fct2=0.9D0

    ! realistic bounds for length predictions
    real(ps),parameter :: max_exice=1.0D0


    real(PS) :: rat
    real(PS),dimension(gs%N_BIN,gs%L) :: mass_left,mass_stot


!!c    real(PS) :: alen_n,clen_n,den_n

!!c    integer :: omp_in_parallel,omp_get_num_threads,omp_get_thread_num


!!c    write(*,'("omp in ini_group",3I5)') omp_in_parallel(),omp_get_num_threads(),omp_get_thread_num()
!!c    write(*,'("memory loc in ini_group",4I15)') loc(ag),loc(gr),loc(gs),loc(XR),loc(XS)
!!c    allocate(&
!!c      ierror(max(maxval(ga(1:nacat)%N_bin),gr%N_bin,gs%N_bin)* &
!!c             max(maxval(ga(1:nacat)%L),gr%L,gs%L)), &
!!c             stat=var_status)
!!c    if(var_status/=0) then
!!c      write(*,*) "ini_group_all: allocation failed 1"
!!c      stop
!!c    endif


    ! initialize tendencies
    if( flagp_r > 0 ) then
      call initialize_ig(gr,1,1)
    endif

    if( flagp_s > 0 ) then
      call initialize_ig(gs,2,2)
      call initialize_shape_ig(gs)
    endif

    if( flagp_a /= 0 ) then
      do k=1,NACAT
        call initialize_ig(ga(k),3,k)
      enddo
    endif


!do ij=1,gr%N_BIN*gr%L
!j=(ij-1)/gr%N_BIN+1
!i=ij-(j-1)*gr%N_BIN
!if (JD(j) == 1 .and. (XS(imt_q,ibi,ici,j)/qipv(icon_q,ibi,1,j) > 1.0D0) .and. (qipv(icon_q,ibi,ici,k+1) /= 0.0D0)) write(*,*) ID(j), JD(j), j, i, qipv(imt_q,ibi,ici,k+1), qipv(icon_q,ibi,ici,k+1), qipv(imt_q,ibi,ici,k+1)/qipv(icon_q,ibi,ici,k+1)
!if (JD(j) == 1 .and. XS(imt_q,i,1,j) > 0.0D0 .and. (qipv(icon_q,i,1,j) <= 1.0e-25)) write(*,*) "FK!", ID(j), JD(j), j, i, qipv(imt_q,ibi,ici,k+1), qipv(icon_q,ibi,ici,k+1), qipv(imt_q,ibi,ici,k+1)/qipv(icon_q,ibi,ici,k+1)
!enddo
!!c    call qiprint2(XS,NSTYPE,NSBIN,NSCAT,L,"bf ini")


    !
    ! liquid spectrum
    !
    if( flagp_r > 0 ) then

       if(debug) then
!        ierror(1:gr%N_BIN*gr%L)=0
!        do ij=1,gr%N_BIN*gr%L
 !         j=(ij-1)/gr%N_BIN+1
!          i=ij-(j-1)*gr%N_BIN
          do j = 1, gr%L
          do i = 1, gr%N_BIN
             if(XR(rcon_q,i,NRCAT,j)>0.0) then
                write(*,'("ini_group_all",4I5,50ES15.6)') i,ID(j),JD(j),KD(j), &
                     XR(rmt_q,i,NRCAT,j),XR(rcon_q,i,NRCAT,j)
             endif
          enddo
          enddo
       endif

!      ierror(1:gr%N_BIN*gr%L)=0
!      do ij=1,gr%N_BIN*gr%L
!        j=(ij-1)/gr%N_BIN+1
!        i=ij-(j-1)*gr%N_BIN
       do j = 1, gr%L
       do i = 1, gr%N_BIN

          gr%MS(i,j)%mass(rmt)=XR(rmt_q,i,NRCAT,j)*ag%TV(j)%den
          if(XR(rmt_q,i,NRCAT,j)*ag%TV(j)%den>0.0) then
             gr%MS(i,j)%con=XR(rcon_q,i,NRCAT,j)*ag%TV(j)%den
          else
             gr%MS(i,j)%con=0.0_PS
             if(debug .and. XR(rmt_q,i,NRCAT,j)*ag%TV(j)%den<0.0_PS) then
                write(*,*) "rain mass is negative at",i,j
                write(*,'(5ES15.6)') XR(rmt_q,i,NRCAT,j),XR(rcon_q,i,NRCAT,j),&
                     gr%MS(i,j)%mass(rmt),gr%MS(i,j)%con,ag%TV(j)%den
             endif
          end if
       enddo
       enddo

      if(level>=4) then
!        do ij=1,gr%N_BIN*gr%L
!          j=(ij-1)/gr%N_BIN+1
!          i=ij-(j-1)*gr%N_BIN
         do j = 1, gr%L
         do i = 1, gr%N_BIN
          if(gr%MS(i,j)%con>0.0_PS) then
            if(XR(rmat_q,i,NRCAT,j)>0.0) then
              ! total aerosol mass
              gr%MS(i,j)%mass(rmat)=min( &
                     max(m_lmt_ap,min_fapt_r*gr%MS(i,j)%mass(rmt) &
                    ,XR(rmat_q,i,NRCAT,j)*ag%TV(j)%den),mx_aprat*gr%MS(i,j)%mass(rmt) )
!org              gr%MS(i,j)%mass(rmat)=min(gr%MS(i,j)%mass(rmt),max(m_lmt_ap,&
!org                     max(&
!org                     min_fapt_r*gr%MS(i,j)%mass(rmt)&
!org                    ,min(XR(rmat_q,i,NRCAT,j)*ag%TV(j)%den,mx_aprat*gr%MS(i,j)%mass(rmt))) &
!org                      ))
!tmp              fapt_r=min(1.0_PS,max(min_fapt_r,gr%MS(i,j)%mass(rmat)/gr%MS(i,j)%mass(rmt)))
!!c              write(*,*) "ck rmatini:",i,j,&
!!c                   gr%MS(i,j)%mass(rmt),gr%MS(i,j)%mass(rmat),&
!!c                   min_fapt_r*gr%MS(i,j)%mass(rmt),&
!!c                   XR(rmat_q,i,NRCAT,j)*ag%TV(j)%den,&
!!c                   mx_aprat*gr%MS(i,j)%mass(rmt)
            else
              gr%MS(i,j)%mass(rmat)=min(max(m_lmt_ap,min_fapt_r*gr%MS(i,j)%mass(rmt)),&
                     gr%MS(i,j)%mass(rmt))
            end if
          endif
        enddo
        enddo
!        do ij=1,gr%N_BIN*gr%L
!          j=(ij-1)/gr%N_BIN+1
!          i=ij-(j-1)*gr%N_BIN
        do j = 1, gr%L
        do i = 1, gr%N_BIN
          if(gr%MS(i,j)%con>0.0_PS) then
            if(XR(rmas_q,i,NRCAT,j)>0.0) then
              ! soluble mass of aerosols
              gr%MS(i,j)%mass(rmas)=min(gr%MS(i,j)%mass(rmat),max(m_lmt_ap,&
                     max(&
                     min_faps_r*gr%MS(i,j)%mass(rmat)&
                    ,min(XR(rmas_q,i,NRCAT,j)*ag%TV(j)%den,gr%MS(i,j)%mass(rmat))) &
                       ))
!tmp              faps_r=min(1.0_PS,max(min_faps_r,gr%MS(i,j)%mass(rmas)/gr%MS(i,j)%mass(rmat)))
            else
              gr%MS(i,j)%mass(rmas)=min(max(m_lmt_ap,min_faps_r*gr%MS(i,j)%mass(rmat)),&
                     gr%MS(i,j)%mass(rmat))
            end if

            ! insoluble mass is diagnosed in diag_pq
            gr%MS(i,j)%mass(rmai)=0.0_PS
          else
            gr%MS(i,j)%mass(rmat)=0.0_PS
            gr%MS(i,j)%mass(rmas)=0.0_PS
            gr%MS(i,j)%mass(rmai)=0.0_PS
          end if
        enddo
        enddo

      end if

!debug      do ij=1,gr%N_BIN*gs%L
!debug        j=(ij-1)/gr%N_BIN+1
!debug        i=ij-(j-1)*gr%N_BIN
!debug        if (gr%MS(i,j)%con /= 0.0_PS .or. XR(rcon_q,i,NRCAT,j) /= 0.0 .or. XR(rmt_q,i,NRCAT,j) /= 0.0) then
!debug          write(fid_alog,'(a,4I3,3ES17.9)') "R1", JD(j), ID(j), i, j, &
!debug               gr%MS(i,j)%mass(rmt), &
!debug               gr%MS(i,j)%con, ag%TV(j)%den
!debug          write(fid_alog,'(a,4I3,2ES17.9)') "R2", JD(j), ID(j), i, j, &
!debug               XR(rmt_q,i,NRCAT,j), &
!debug               XR(rcon_q,i,NRCAT,j)
!debug          write(fid_alog,'(a,4I3,7ES17.9)') "R3", JD(j), ID(j), i, j, &
!debug               XR(rmt_q,i,NRCAT,j)*ag%TV(j)%den, &
!debug               XR(rcon_q,i,NRCAT,j)*ag%TV(j)%den
!debug        end if
!debug      end do

      if(debug) then
!        ierror(1:gr%L)=0
         do j=1,gr%L
            if(sum(gr%MS(1:gr%N_BIN,j)%con)>200.0) then
               do i=1,gr%N_BIN
                  write(*,'("ini_group_all",4I5,50ES15.6)') i,ID(j),JD(j),KD(j), &
                       XR(rmt_q,i,NRCAT,j),XR(rcon_q,i,NRCAT,j), &
                       gr%MS(i,j)%con
               enddo
            endif
         end do
      end if

   end if

    !
    ! ice spectrum
    !
    if( flagp_s > 0 ) then

      if(debug) then
!        ierror(1:gs%N_BIN*gs%L)=0
!        do ij=1,gs%N_BIN*gs%L
!          j=(ij-1)/gs%N_BIN+1
!          i=ij-(j-1)*gs%N_BIN
         do j = 1, gs%L
         do i = 1, gs%N_BIN
            if(XS(icon_q,i,NSCAT,j)>0.0) then
               write(*,'("ini_group_all,ice:",4I5,50ES15.6)') i,ID(j),JD(j),KD(j) &
                    ,XS(imt_q,i,NSCAT,j),XS(icon_q,i,NSCAT,j)
            endif
         end do
         end do
      endif

!      do ij=1,gs%N_BIN*gs%L
!        j=(ij-1)/gs%N_BIN+1
!        i=ij-(j-1)*gs%N_BIN
      do j = 1, gs%L
      do i = 1, gs%N_BIN

        gs%MS(i,j)%mass(imt) = XS(imt_q,i,NSCAT,j)*ag%TV(j)%den
        if(XS(imt_q,i,NSCAT,j)*ag%TV(j)%den>0.0) then
          gs%MS(i,j)%con=XS(icon_q,i,NSCAT,j)*ag%TV(j)%den
        else
          gs%MS(i,j)%con=0.0_PS
        end if

      enddo
      enddo

!dbg      do j=1,gs%L
!dbg        write(*,'("ck ini_group con_s",I5,30ES15.6)') j,(gs%MS(i,j)%con,i=1,gs%N_BIN)
!dbg      enddo

      if(level>=4) then
        !
        ! aerosol prediction
        !
!        do ij=1,gs%N_BIN*gs%L
!          j=(ij-1)/gs%N_BIN+1
!          i=ij-(j-1)*gs%N_BIN
         do j = 1, gs%L
         do i = 1, gs%N_BIN
          if( gs%MS(i,j)%con > 0.0_PS ) then
            if(XS(imat_q,i,NSCAT,j)>0.0_PS) then
              ! total aerosol mass

              gs%MS(i,j)%mass(imat)=min( &
                    max(m_lmt_ap,min_fapt_s*gs%MS(i,j)%mass(imt) &
                   ,XS(imat_q,i,NSCAT,j)*ag%TV(j)%den),mx_aprat*gs%MS(i,j)%mass(imt) )

!org              gs%MS(i,j)%mass(imat)=min(gs%MS(i,j)%mass(imt),max(m_lmt_ap,&
!org                    max(&
!org                    min_fapt_s*gs%MS(i,j)%mass(imt)&
!org                   ,min(XS(imat_q,i,NSCAT,j)*ag%TV(j)%den,mx_aprat*gs%MS(i,j)%mass(imt))) &
!org                         ))
!tmp              fapt_s=min(1.0_PS,max(min_fapt_s,gs%MS(i,j)%mass(imat)/gs%MS(i,j)%mass(imt)))
            else
              gs%MS(i,j)%mass(imat)=min(max(m_lmt_ap,min_fapt_s*gs%MS(i,j)%mass(imt)),&
                    gs%MS(i,j)%mass(imt))
            end if
          endif
        enddo
        enddo
!        do ij=1,gs%N_BIN*gs%L
!          j=(ij-1)/gs%N_BIN+1
!          i=ij-(j-1)*gs%N_BIN
        do j = 1, gs%L
        do i = 1, gs%N_BIN
          if( gs%MS(i,j)%con > 0.0_PS ) then
            if(XS(imas_q,i,NSCAT,j)>=0.0_PS) then
              ! soluble mass of aerosols
              gs%MS(i,j)%mass(imas)=min(gs%MS(i,j)%mass(imat),max(m_lmt_ap, &
                    max(&
                    min_faps_s*gs%MS(i,j)%mass(imat)&
                   ,min(XS(imas_q,i,NSCAT,j)*ag%TV(j)%den,gs%MS(i,j)%mass(imat))) &
                        ))
!tmp              faps_s=min(1.0_PS,max(min_faps_s,gs%MS(i,j)%mass(imas)/gs%MS(i,j)%mass(imat)))
            else
              gs%MS(i,j)%mass(imas)=min(max(m_lmt_ap,min_faps_s*gs%MS(i,j)%mass(imat)),&
                    gs%MS(i,j)%mass(imat))
            end if

            ! insoluble mass is diagnosed in diag_pq
            gs%MS(i,j)%mass(imai)=0.0_PS
          else

            gs%MS(i,j)%mass(imat)=0.0_PS
            gs%MS(i,j)%mass(imas)=0.0_PS
            gs%MS(i,j)%mass(imai)=0.0_PS
          end if
        enddo
        enddo

      endif

!      do ij=1,gs%N_BIN*gs%L
!        j=(ij-1)/gs%N_BIN+1
!        i=ij-(j-1)*gs%N_BIN
      do j = 1, gs%L
      do i = 1, gs%N_BIN
        if( gs%MS(i,j)%con > 0.0_PS ) then
          !
          ! crystal mass component
          !
          if( XS(imt_q,i,NSCAT,j) <= XS(imc_q,i,NSCAT,j)) then
            gs%MS(i,j)%mass(imc) = gs%MS(i,j)%mass(imt)
          else
            if(gs%MS(i,j)%mass(imt)>=gs%MS(i,j)%con*m_icmin) then
              gs%MS(i,j)%mass(imc)=max( &
                   min(XS(imc_q,i,NSCAT,j)*ag%TV(j)%den,gs%MS(i,j)%mass(imt)) &
                  ,gs%MS(i,j)%con*m_icmin)
            else
              gs%MS(i,j)%mass(imc) = min(XS(imc_q,i,NSCAT,j)*ag%TV(j)%den,gs%MS(i,j)%mass(imt))
            end if
          end if
        else
          gs%MS(i,j)%mass(imc) = 0.0_PS

        endif
      enddo
      enddo
!      do ij=1,gs%N_BIN*gs%L
!        j=(ij-1)/gs%N_BIN+1
!        i=ij-(j-1)*gs%N_BIN
      do j = 1, gs%L
      do i = 1, gs%N_BIN
        if( gs%MS(i,j)%con > 0.0_PS ) then
          !
          ! rime mass component
          !
          if( gs%MS(i,j)%mass(imt)-gs%MS(i,j)%mass(imc) < XS(imr_q,i,NSCAT,j)*ag%TV(j)%den ) then
            gs%MS(i,j)%mass(imr) = max(gs%MS(i,j)%mass(imt)-gs%MS(i,j)%mass(imc),0.0_PS)
          else
            gs%MS(i,j)%mass(imr) = max(XS(imr_q,i,NSCAT,j)*ag%TV(j)%den,0.0_PS)
          end if
          !
          ! aggregation mass component
          !
          if( gs%MS(i,j)%mass(imt)-gs%MS(i,j)%mass(imc)<XS(ima_q,i,NSCAT,j)*ag%TV(j)%den ) then
            gs%MS(i,j)%mass(ima) = max(gs%MS(i,j)%mass(imt)-gs%MS(i,j)%mass(imc),0.0_PS)
          else
            gs%MS(i,j)%mass(ima) = max(XS(ima_q,i,NSCAT,j)*ag%TV(j)%den,0.0_PS)
          end if
        else
          gs%MS(i,j)%mass(imr) = 0.0_PS
          gs%MS(i,j)%mass(ima) = 0.0_PS
        endif
      enddo
      enddo

      if(level>=6) then
!        do ij=1,gs%N_BIN*gs%L
!          j=(ij-1)/gs%N_BIN+1
!          i=ij-(j-1)*gs%N_BIN
         do j = 1, gs%L
         do i = 1, gs%N_BIN
          if( gs%MS(i,j)%con > 0.0_PS ) then
            !
            ! melt water mass
            !
            if( gs%MS(i,j)%mass(imt)-gs%MS(i,j)%mass(imc)<XS(imw_q,i,NSCAT,j)*ag%TV(j)%den ) then
              gs%MS(i,j)%mass(imw) = max(gs%MS(i,j)%mass(imt)-gs%MS(i,j)%mass(imc),0.0_PS)
            else
              gs%MS(i,j)%mass(imw) = max(XS(imw_q,i,NSCAT,j)*ag%TV(j)%den,0.0_PS)
            end if
            !
            ! freezing nucleation mass
            ! Contact freezing, immersion freezing, and homogeneous freezing produce this mass.
            !
            if( gs%MS(i,j)%mass(imc)<XS(imf_q,i,NSCAT,j)*ag%TV(j)%den ) then
              gs%MS(i,j)%mass(imf) = gs%MS(i,j)%mass(imc)
            else
              gs%MS(i,j)%mass(imf) = max(XS(imf_q,i,NSCAT,j)*ag%TV(j)%den,0.0_PS)
            end if
          else
            gs%MS(i,j)%mass(imf) = 0.0_PS
          end if
        enddo
        enddo
      endif

!debug      do ij=1,gs%N_BIN*gs%L
!debug        j=(ij-1)/gs%N_BIN+1
!debug        i=ij-(j-1)*gs%N_BIN
!debug        if (gs%MS(i,j)%con /= 0.0_PS .or. XS(icon_q,i,NSCAT,j) /= 0.0 .or. XS(imt_q,i,NSCAT,j) /= 0.0) then
!debug          write(fid_alog,'(a,4I3,8ES17.9)') "1", JD(j), ID(j), i, j, &
!debug               gs%MS(i,j)%mass(imt), &
!debug               gs%MS(i,j)%mass(imr), &
!debug               gs%MS(i,j)%mass(ima), &
!debug               gs%MS(i,j)%mass(imc), &
!debug               gs%MS(i,j)%mass(imw), &
!debug               gs%MS(i,j)%mass(imf), &
!debug               gs%MS(i,j)%con, ag%TV(j)%den
!debug          write(fid_alog,'(a,4I3,7ES17.9)') "2", JD(j), ID(j), i, j, &
!debug               XS(imt_q,i,NSCAT,j), &
!debug               XS(imr_q,i,NSCAT,j), &
!debug               XS(ima_q,i,NSCAT,j), &
!debug               XS(imc_q,i,NSCAT,j), &
!debug               XS(imw_q,i,NSCAT,j), &
!debug               XS(imf_q,i,NSCAT,j), &
!debug               XS(icon_q,i,NSCAT,j)
!debug          write(fid_alog,'(a,4I3,7ES17.9)') "3", JD(j), ID(j), i, j, &
!debug               XS(imt_q,i,NSCAT,j)*ag%TV(j)%den, &
!debug               XS(imr_q,i,NSCAT,j)*ag%TV(j)%den, &
!debug               XS(ima_q,i,NSCAT,j)*ag%TV(j)%den, &
!debug               XS(imc_q,i,NSCAT,j)*ag%TV(j)%den, &
!debug               XS(imw_q,i,NSCAT,j)*ag%TV(j)%den, &
!debug               XS(imf_q,i,NSCAT,j)*ag%TV(j)%den, &
!debug               XS(icon_q,i,NSCAT,j)*ag%TV(j)%den
!debug        end if
!debug      end do

!      do ij=1,gs%N_BIN*gs%L
!        j=(ij-1)/gs%N_BIN+1
!        i=ij-(j-1)*gs%N_BIN
      do j = 1, gs%L
      do i = 1, gs%N_BIN
        if( gs%MS(i,j)%con > 0.0_PS ) then
          mass_left(i,j)=max(0.0_PS,gs%MS(i,j)%mass(imt)-gs%MS(i,j)%mass(imc))
          mass_stot(i,j)=gs%MS(i,j)%mass(imr)+gs%MS(i,j)%mass(ima)+gs%MS(i,j)%mass(imw)
        endif
      enddo
      enddo

!      do ij=1,gs%N_BIN*gs%L
!        j=(ij-1)/gs%N_BIN+1
!        i=ij-(j-1)*gs%N_BIN
      do j = 1, gs%L
      do i = 1, gs%N_BIN
        if( gs%MS(i,j)%con > 0.0_PS ) then
          !
          ! maintain the relative relation of predicted mass component
          !
          if(mass_stot(i,j)>1.0e-25_PS) then
            rat=mass_left(i,j)/mass_stot(i,j)
            gs%MS(i,j)%mass(imr)=gs%MS(i,j)%mass(imr)*rat
            gs%MS(i,j)%mass(ima)=gs%MS(i,j)%mass(ima)*rat
            gs%MS(i,j)%mass(imw)=gs%MS(i,j)%mass(imw)*rat
          else
            gs%MS(i,j)%mass(imc)=gs%MS(i,j)%mass(imt)
            gs%MS(i,j)%mass(imr)=0.0_PS
            gs%MS(i,j)%mass(ima)=0.0_PS
            gs%MS(i,j)%mass(imw)=0.0_PS
          end if
        endif
      enddo
      enddo

!      do ij=1,gs%N_BIN*gs%L
!        j=(ij-1)/gs%N_BIN+1
!        i=ij-(j-1)*gs%N_BIN
      do j = 1, gs%L
      do i = 1, gs%N_BIN
        if( gs%MS(i,j)%con > 0.0_PS ) then
          !
          ! a and c are considered together, because imbalance causes
          ! rediculous aspect ratio.
          ! a and c axis lengths of ice crystals: max 5 cm.
          !
          if( XS(iacr_q,i,NSCAT,j) <= 0.0 .or. XS(iccr_q,i,NSCAT,j) <= 0.0) then

            if(gs%MS(i,j)%mass(imt)<=1.0e-30.or.gs%MS(i,j)%con<=1.0e-30) then
              gs%MS(i,j)%mass(imt)=0.0
              gs%MS(i,j)%con=0.0
            else
              call gen_length(gs%MS(i,j)%mass(imc)/gs%MS(i,j)%con,& !ag%TV(j)%T,
                      gs%MS(i,j)%a_len,gs%MS(i,j)%c_len)

              gs%IS(i,j)%d=0.0_PS
              gs%IS(i,j)%ag=0.0_PS
              gs%IS(i,j)%cg=0.0_PS
              gs%IS(i,j)%n_exice=0.0_PS

!           write(*,'("lengths generated",4I5,14ES15.6)') i,ID(j),JD(j),KD(j), &
!                    gs%MS(i,j)%mass(imc), gs%MS(i,j)%con, XS(iacr_q,i,NSCAT,j), XS(iccr_q,i,NSCAT,j), XS(ivcs_q,i,NSCAT,j), &
!                    ag%TV(j)%T,gs%MS(i,j)%mass(imc)/gs%MS(i,j)%con,&
!                    gs%MS(i,j)%a_len,gs%MS(i,j)%c_len,gs%IS(i,j)%d,&
!                    gs%IS(i,j)%ag,gs%IS(i,j)%cg,gs%IS(i,j)%n_exice
            endif

          else
!!c            write(*,'("ck_in_lengths:",4I5,10ES15.6)') i,ID(j),JD(j),KD(j),&
!!c                    ag%TV(j)%T,gs%MS(i,j)%mass(imc)/gs%MS(i,j)%con,&
!!c                    gs%MS(i,j)%a_len,gs%MS(i,j)%c_len,gs%IS(i,j)%d,&
!!c                    XS(iacr_q,i,nscat,j),XS(iccr_q,i,nscat,j),XS(idcr_q,i,nscat,j)
!!c            if(isnan(XS(iacr_q,i,nscat,j))) then
!!c               write(*,*) "NaN"
!!c               stop
!!c            endif

            gs%MS(i,j)%a_len=min(max(1.0e-4_PS,&
                   (XS(iacr_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j))**(1.0/3.0)),5.0_PS)
            gs%MS(i,j)%c_len=min(max(1.0e-4_PS,&
                   (XS(iccr_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j))**(1.0/3.0)),5.0_PS)
          endif
        endif
      enddo
      enddo

!      do ij=1,gs%N_BIN*gs%L
!        j=(ij-1)/gs%N_BIN+1
!        i=ij-(j-1)*gs%N_BIN
      do j = 1, gs%L
      do i = 1, gs%N_BIN
        if( gs%MS(i,j)%con > 0.0_PS ) then
          if( XS(iacr_q,i,NSCAT,j) > 0.0.and.XS(iccr_q,i,NSCAT,j) > 0.0) then
            !
            ! reality check of axis ratio
            !
            if(gs%MS(i,j)%c_len>gs%MS(i,j)%a_len*phi_max) then
!!c              write(*,'("axis ratio too big",4I5)') i,ID(j),JD(j),KD(j)
              gs%MS(i,j)%c_len=gs%MS(i,j)%a_len*phi_max
            elseif(gs%MS(i,j)%c_len<gs%MS(i,j)%a_len*phi_min) then
!!c              write(*,'("axis ratio too small",4I5)') i,ID(j),JD(j),KD(j)
              gs%MS(i,j)%a_len=gs%MS(i,j)%c_len/phi_min
            end if

            !
            ! d axis lengths of ice crystals should be bounded by a-axis length
            !
            if( level >= 2 ) then
              if( XS(idcr_q,i,NSCAT,j) > 0.0.and.gs%MS(i,j)%a_len>=10.0e-4_PS ) then
                gs%IS(i,j)%d = max(min( &
                        (XS(idcr_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j))**(1.0/3.0) &
                        ,0.9_PS*gs%MS(i,j)%a_len),0.0_PS)
              else
!!c               write(*,*) "ini_group > neg length for d, modified to 0.0e-4"
                gs%IS(i,j)%d = 0.0_PS
              end if
            end if
            !
            ! number of extra single ice crystals
            !
            if( level==7 ) then
              gs%IS(i,j)%n_exice=max(0.0_PS,min(max_exice,XS(inex_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j)))
            end if
!tmp            if(gs%IS(i,j)%n_exice>1.0) then
!tmp              write(*,*) "ini_group > nexice",gs%IS(i,j)%n_exice
!tmp            end if

            !
            ! position of center of gravity: a-axis coordinate
            !
            if( level==3.or.level==5.or. level==7 ) then
              if( XS(iag_q,i,NSCAT,j) > 0.0 ) then
                gs%IS(i,j)%ag=max(min( &
                     (XS(iag_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j))**(1.0/3.0) &
                     ,gs%MS(i,j)%a_len),0.0_PS)
              else
                gs%IS(i,j)%ag=0.0_PS
              end if
            end if
            !
            ! position of center of gravity: c-axis coordinate
            !
            if( level==3.or.level==5.or.level==7 ) then
              if( XS(icg_q,i,NSCAT,j) > 0.0_PS ) then
                gs%IS(i,j)%cg=max(min( &
                     (XS(icg_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j))**(1.0_PS/3.0_PS) &
                     ,gs%MS(i,j)%c_len),0.0_PS)
              else
                gs%IS(i,j)%cg=0.0_PS
              end if
            end if

            ! reality check of density
            ! calculate minimum lengths, assuming hexagonal plates
!!c                   if(gs%IS(i,j)%n_exice<0.5_PS) then
            phi=gs%MS(i,j)%c_len/gs%MS(i,j)%a_len
            psi=gs%IS(i,j)%d/gs%MS(i,j)%a_len
            pag=gs%IS(i,j)%ag/gs%MS(i,j)%a_len
            pcg=gs%IS(i,j)%cg/gs%MS(i,j)%c_len

            if(debug .and. phi>2.0_PS) then
              write(*,*) "ini_group_all:phi>2.",i,j, &
                gs%MS(i,j)%con,gs%MS(i,j)%mass,gs%MS(i,j)%a_len,gs%MS(i,j)%c_len
            endif

            a_min=(gs%MS(i,j)%mass(imc)/gs%MS(i,j)%con/&
                        (3.0_PS*sq_three*phi*(1.0-psi)*den_i))**0.33333333
!!c                   a_min = ((gs%MS(i,j)%mass(imc)/gs%MS(i,j)%con/(1.0_PS+gs%IS(i,j)%n_exice))/&
!!c                        (3.0_PS*sq_three*phi*den_i))**0.33333333
            c_min=a_min*phi
            d_min=a_min*psi
            ag_min=a_min*pag
            cg_min=c_min*pcg

            ! reality check of density
            if(a_min>gs%MS(i,j)%a_len) then
!!c                      write(*,*) "ini_group>a_len is too small, so modified from:",gs%MS(i,j)%a_len,a_min
              gs%MS(i,j)%c_len=c_min
              gs%MS(i,j)%a_len=a_min
              gs%IS(i,j)%d=d_min
              gs%IS(i,j)%ag=ag_min
              gs%IS(i,j)%cg=cg_min
            elseif(c_min>gs%MS(i,j)%c_len) then
!!c                      write(*,*) "ini_group>c_len is too small, so modified from:",gs%MS(i,j)%c_len,c_min
              gs%MS(i,j)%c_len=c_min
              gs%MS(i,j)%a_len=a_min
              gs%IS(i,j)%d=d_min
              gs%IS(i,j)%ag=ag_min
              gs%IS(i,j)%cg=cg_min
            endif
          end if
        endif
      enddo
      enddo

!      do ij=1,gs%N_BIN*gs%L
!        j=(ij-1)/gs%N_BIN+1
!        i=ij-(j-1)*gs%N_BIN
      do j = 1, gs%L
      do i = 1, gs%N_BIN
        if( gs%MS(i,j)%con > 0.0_PS ) then
          !
          ! circumscribing spheroid: max with 1.0e+1 (cm) radius
          !                          min with 1.0e-4 (cm) radius
          if(XS(ivcs_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j) <= 0.0_PS ) then
!!c                   write(*,*) "ini_group > Vcs is non-positive!"
!!c                   write(*,100) gs%MS(i,j)%con, gs%MS(i,j)%mass(imt), &
!!c                        gs%MS(i,j)%mass(imr),gs%MS(i,j)%mass(ima), gs%IS(i,j)%V_cs
            gs%MS(i,j)%con = 0.0_PS
            gs%MS(i,j)%mass(imt) = 0.0_PS
            gs%MS(i,j)%mass(imc)=0.0_PS
            gs%MS(i,j)%mass(imr)=0.0_PS
            gs%MS(i,j)%mass(ima)=0.0_PS
            gs%MS(i,j)%mass(imw)=0.0_PS
            gs%MS(i,j)%mass(imf)=0.0_PS
            gs%MS(i,j)%mass(imat)=0.0_PS
            gs%MS(i,j)%mass(imas)=0.0_PS
            gs%IS(i,j)%V_cs=0.0_PS
            gs%MS(i,j)%a_len=0.0_PS
            gs%MS(i,j)%c_len=0.0_PS
            gs%IS(i,j)%d=0.0_PS
            gs%IS(i,j)%ag=0.0_PS
            gs%IS(i,j)%cg=0.0_PS
            gs%IS(i,j)%n_exice=0.0_PS
          elseif(XS(ivcs_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j)<V_csmin_hex) then
!!c                   write(*,*) "Volume is smaller than V_csmin",i,j
            gs%IS(i,j)%V_cs=V_csmin_hex
            gs%MS(i,j)%a_len = 1.0e-4_PS
            gs%MS(i,j)%c_len = 1.0e-4_PS
            gs%IS(i,j)%d = 0.0_PS
            gs%IS(i,j)%ag=0.0_PS
            gs%IS(i,j)%cg=0.0_PS
            gs%IS(i,j)%n_exice=0.0_PS
            gs%MS(i,j)%mass(imr)=0.0_PS
            tmass1=m_icmin*gs%MS(i,j)%con
            gs%MS(i,j)%mass(imc)=tmass1
            gs%MS(i,j)%mass(imt)=tmass1
            tmass2=min(max(m_lmt_ap,min_fapt_s*tmass1),tmass1)
            gs%MS(i,j)%mass(imat)=tmass2
            gs%MS(i,j)%mass(imas)=min(max(m_lmt_ap,min_faps_s*tmass2),tmass2)
            gs%MS(i,j)%mass(imw)=0.0_PS
            gs%MS(i,j)%mass(imf)=0.0_PS

          elseif(XS(ivcs_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j)>V_csmax) then
            gs%IS(i,j)%V_cs=V_csmax
          else
            gs%IS(i,j)%V_cs=XS(ivcs_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j)
          end if
        else
          gs%MS(i,j)%mass(imat)=0.0_PS
          gs%MS(i,j)%mass(imas)=0.0_PS
          gs%MS(i,j)%mass(imc)=0.0_PS
          gs%MS(i,j)%mass(imr)=0.0_PS
          gs%MS(i,j)%mass(ima)=0.0_PS
          gs%MS(i,j)%mass(imw)=0.0_PS
          gs%MS(i,j)%mass(imf)=0.0_PS
          gs%MS(i,j)%a_len=0.0_PS
          gs%MS(i,j)%c_len=0.0_PS
          gs%IS(i,j)%d=0.0_PS
          gs%IS(i,j)%n_exice=0.0_PS
          gs%IS(i,j)%ag=0.0_PS
          gs%IS(i,j)%cg=0.0_PS
          gs%IS(i,j)%V_cs=0.0_PS
        end if
      enddo
      enddo

      if ( debug ) then
!      ierror(1:gs%N_BIN*gs%L)=0
!      do ij=1,gs%N_BIN*gs%L
!        j=(ij-1)/gs%N_BIN+1
!        i=ij-(j-1)*gs%N_BIN
         do j = 1, gs%L
         do i = 1, gs%N_BIN
            if(gs%MS(i,j)%mass(imt)<0.0_PS) then
               write(*,17) i,j,XS(imt_q,i,NSCAT,j),XS(icon_q,i,NSCAT,j),ag%TV(j)%den
17             format("solid mass is negative at bin,grid,mass,con,den",2I5,3ES15.6)
            end if
            if(gs%MS(i,j)%con<0.0_PS) then
               write(*,16) i,j,XS(imt_q,i,NSCAT,j),XS(icon_q,i,NSCAT,j),ag%TV(j)%den
16             format("solid con is negative at bin,grid,mass,con,den",2I5,3ES15.6)
            end if
            if(gs%MS(i,j)%mass(imt)<gs%MS(i,j)%mass(imc)) then
               write(*,*) "m1<m4",i,j,gs%MS(i,j)%mass(imt),gs%MS(i,j)%mass(imc)
            endif
!!!        if(gs%MS(i,j)%con>1.0e-30.and.gs%IS(i,j)%V_cs>1.0e-30.and.gs%MS(i,j)%mass(imt)>1.0e-30.and.&
!!!               gs%MS(i,j)%mass(imt)/max(1.0e-30,gs%IS(i,j)%V_cs)/max(1.0e-30,gs%MS(i,j)%con)<1.0e-2) then
!!!          ierror(ij)=7
!            write(*,*) "ini_group_all: den<1.0e-2",i,KD(j),ID(j),JD(j),gs%MS(i,j)%con,gs%MS(i,j)%mass(imt) &
!                    ,gs%IS(i,j)%V_cs &
!                    ,gs%MS(i,j)%mass(imt)/gs%IS(i,j)%V_cs/max(1.0e-30_RP,gs%MS(i,j)%con)
!!!        endif
         enddo
         enddo
      end if
    end if

    !
    ! aerosol
    !
    if( flagp_a==1 ) then
      do k=1,NACAT
!        do ij=1,ga(k)%N_BIN*ga(k)%L
!          j=(ij-1)/ga(k)%N_BIN+1
!          i=ij-(j-1)*ga(k)%N_BIN
         do j = 1, ga(k)%L
         do i = 1, ga(k)%N_BIN

          con1=XA(acon_q,i,k,j)*ag%TV(j)%den
          if(con1<0.0_PS) then
!!c                   write(*,*) "ap con is negative at",i,j
!!c                   write(*,'(5ES15.6)') XA(amt_q,i,k,j),XA(acon_q,i,k,j),&
!!c                        ga(k)%MS(i,j)%mass(amt),ga(k)%MS(i,j)%con,ag%TV(j)%den
            ga(k)%MS(i,j)%con=0.0_PS
            ga(k)%MS(i,j)%mass(amt)=0.0_PS
            ga(k)%MS(i,j)%mass(ams)=0.0_PS
          else
            ga(k)%MS(i,j)%con=con1
            tmass1=con1*coef_ap(k)
            ga(k)%MS(i,j)%mass(amt)=tmass1
            ga(k)%MS(i,j)%mass(ams)=tmass1*ga(k)%MS(i,j)%eps_map
          end if
          !write(fid_alog,'(a,3I3,7ES15.6)'), "aerosol", k, j, i, ga(k)%MS(i,j)%con, ga(k)%MS(i,j)%mass(amt), XA(acon_q,i,k,j), XA(acon_q,i,k,j)*ag%TV(j)%den, ag%TV(j)%den, tmass1, coef_ap(k)
        end do
        end do
      enddo
    elseif(flagp_a/=0) then
      do k=1,NACAT
!        ierror(1:ga(k)%N_BIN*ga(k)%L)=0

! <<< 2014/12 T. Hashino bug found
        m_lmt=coef4pi3*r3_lmt*den_apt0(k)
!test        if(k==1.or.k==2.or.k==4) then
!test          m_lmt=coef4pi3*r3_lmt*den_apt0(k)
!test        elseif(k==3) then
!test          m_lmt=coef4pi3*r3_lmt_nuc*den_apt0(k)
! >>> 2014/12 T. Hashino bug found
!test        else
!test          write(*,*) "this category not defined in ini_group_all"
!test          stop
!test        endif

!        do ij=1,ga(k)%N_BIN*ga(k)%L
!          j=(ij-1)/ga(k)%N_BIN+1
!          i=ij-(j-1)*ga(k)%N_BIN
        do j = 1, ga(k)%L
        do i = 1, ga(k)%N_BIN

          ga(k)%MS(i,j)%con=max(n_lmt_ap,XA(acon_q,i,k,j)*ag%TV(j)%den)

          tmass1=XA(amt_q,i,k,j)*ag%TV(j)%den
          ga(k)%MS(i,j)%mass(amt)=tmass1
!tmp                write(*,*) "ck ini_group0:",k,j,amt,amt_q,ga(k)%MS(i,j)%mass(amt),XA(amt_q,i,k,j),ag%TV(j)%den

          if(level>=4) then
            tmass2=max(0.0_PS,min(XA(ams_q,i,k,j)*ag%TV(j)%den,tmass1))

            ga(k)%MS(i,j)%mass(ams)=tmass2
            if(k==1) then
              ! for CCN accumulation mode (category 1)
              ga(k)%MS(i,j)%mass(ams)=min(max(tmass1*sep_faps,&
                          tmass2),tmass1)
            elseif(k==2) then
              ! for IN accumulation mode (category 2)
              ga(k)%MS(i,j)%mass(ams)=max(min(tmass1*0.99*sep_faps,&
                          tmass2),0.0_PS)
            elseif(k==3) then
              ! for CCN nucleation mode (category 3)
              ga(k)%MS(i,j)%mass(ams)=min(max(tmass1*sep_faps,&
                          tmass2),tmass1)
            elseif(k==4) then
              ! for CCN accumulation mode (category 4)
              ga(k)%MS(i,j)%mass(ams)=min(max(tmass1*sep_faps,&
                          ga(k)%MS(i,j)%mass(ams)),tmass1)
            end if

          end if
        enddo
        enddo

!        do ij=1,ga(k)%N_BIN*ga(k)%L
!          j=(ij-1)/ga(k)%N_BIN+1
!          i=ij-(j-1)*ga(k)%N_BIN
        do j = 1, ga(k)%L
        do i = 1, ga(k)%N_BIN

          tmass1=m_lmt*ga(k)%MS(i,j)%con
          if(ga(k)%MS(i,j)%mass(amt)<tmass1) then
            ga(k)%MS(i,j)%mass(amt)=tmass1
            ga(k)%MS(i,j)%mass(ams)=eps_ap0(k)*tmass1
          end if
        enddo
        enddo
!        ierror(1:ga(k)%N_BIN*ga(k)%L)=0
!        do ij=1,ga(k)%N_BIN*ga(k)%L
!          j=(ij-1)/ga(k)%N_BIN+1
!          i=ij-(j-1)*ga(k)%N_BIN
        do j = 1, ga(k)%L
        do i = 1, ga(k)%N_BIN
          if(ga(k)%MS(i,j)%mass(amt)<1.0e-30_PS) then
             LOG_ERROR("ini_group_all",*) "ap mass is negative at",i,j
             LOG_ERROR_CONT(*) XA(amt_q,i,k,j),XA(acon_q,i,k,j),&
                    ga(k)%MS(i,j)%mass(amt),ga(k)%MS(i,j)%con,ag%TV(j)%den
              call PRC_abort
          endif
          if(ga(k)%MS(i,j)%con<n_lmt_ap*0.1) then
             LOG_ERROR("ini_group_all",*) "ap con is negative at",i,j
             LOG_ERROR_CONT(*) XA(amt_q,i,k,j),XA(acon_q,i,k,j),&
                    ga(k)%MS(i,j)%mass(amt),ga(k)%MS(i,j)%con,ag%TV(j)%den
              call PRC_abort
          endif
        enddo
        enddo
      enddo
    endif

!!c    deallocate(ierror)

100 format("ini_group >", 5(2X,E15.8))
!!c    call qiprint3(gs,L,"af ini")

  CONTAINS
    subroutine initialize_ig(g,iphase,ica)
      type (Group), intent(inout)  :: g
      integer,intent(in) :: iphase,ica
      integer :: i,j,k,m

      ! Note: initialize the tendencies with 0, but others with
      ! undef since these are calculated later.

!      do k=1,g%N_tendpros
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN
      do k=1,g%N_tendpros
         g%MS(i,j)%dcondt(k) = 0.0_PS
      enddo
      enddo
      enddo
!      do m=1,1+g%N_masscom
!      do k=1,g%N_tendpros
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN
      do k=1,g%N_tendpros
      do m=1,1+g%N_masscom
         g%MS(i,j)%dmassdt(m,k) = 0.0_PS
      enddo
      enddo
      enddo
      enddo

!      do m=1,g%N_nonmass
!      do k=1,g%N_tendpros
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN
      do m=1,g%N_nonmass
      do k=1,g%N_tendpros
         g%MS(i,j)%dvoldt(m,k) = 0.0_PS
      enddo
      enddo
      enddo
      enddo

!      do k=1,1+g%N_masscom
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN
      do k=1,1+g%N_masscom
         g%MS(i,j)%mass(k)=undef
      enddo
      enddo
      enddo

!      do k=1,g%N_vol
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN
      do k=1,g%N_vol
         g%MS(i,j)%vol(k)=undef
      enddo
      enddo
      enddo

!      do ij=1,g%N_BIN*g%L
!        j=(ij-1)/g%N_BIN+1
!        i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN

        g%MS(i,j)%Ldmassdt(1) = 0.0_PS
        g%MS(i,j)%Ldmassdt(2) = 0.0_PS

        g%MS(i,j)%coef(1)=undef
        g%MS(i,j)%coef(2)=undef

        g%MS(i,j)%p(1)=undef
        g%MS(i,j)%p(2)=undef
        g%MS(i,j)%p(3)=undef

        g%MS(i,j)%con=undef
        g%MS(i,j)%mean_mass=undef
        g%MS(i,j)%len=undef

        g%MS(i,j)%den=undef
        g%MS(i,j)%den_as=den_aps0(ica)
        g%MS(i,j)%den_ai=den_api0(ica)

        g%MS(i,j)%a_len=undef
        g%MS(i,j)%c_len=undef

        g%MS(i,j)%eps_map=eps_ap0(ica)

        g%MS(i,j)%vtm=undef
        g%MS(i,j)%tmp=undef
        g%MS(i,j)%fv=undef
        g%MS(i,j)%fkn=undef
        g%MS(i,j)%fac=undef
        g%MS(i,j)%fh=undef
        g%MS(i,j)%Nre=undef
        g%MS(i,j)%CAP=undef
        g%MS(i,j)%CAP_hex=undef
        g%MS(i,j)%semi_a=undef
        g%MS(i,j)%semi_c=undef
        g%MS(i,j)%e_sat=undef
        g%MS(i,j)%r_act=undef
        g%MS(i,j)%r_crt=undef
        g%MS(i,j)%th00_cp=undef

        g%MS(i,j)%mark=0
        g%MS(i,j)%inevp=0
        g%MS(i,j)%inmlt=0
      enddo
      enddo


      do j=1,g%L
        g%mark_cm(j)=0
        g%mark_er(j) = 0
      enddo

    end subroutine initialize_ig

    subroutine initialize_shape_ig(g)
      type (Group), intent(inout)  :: g
      integer  :: i,j,k,l

!      do ij=1,g%N_BIN*g%L
!        j=(ij-1)/g%N_BIN+1
!        i=ij-(j-1)*g%N_BIN
      do j = 1, g%L
      do i = 1, g%N_BIN

        g%IS(i,j)%a = undef
        g%IS(i,j)%d = undef
        g%IS(i,j)%c = undef
        g%IS(i,j)%r = undef
        g%IS(i,j)%e = undef
        g%IS(i,j)%ag = undef
        g%IS(i,j)%cg = undef
        g%IS(i,j)%n_exice = undef

        g%IS(i,j)%V_ic = undef
        g%IS(i,j)%V_cs = undef
        g%IS(i,j)%semi_aip = undef
        g%IS(i,j)%semi_cip = undef
        g%IS(i,j)%den_ip = undef
        g%IS(i,j)%den_ic = undef
        g%IS(i,j)%V_csw = undef
        g%IS(i,j)%phi_ic = undef
        g%IS(i,j)%psi_ic = undef
        g%IS(i,j)%gam_ic = undef
        g%IS(i,j)%eta_ic = undef
        g%IS(i,j)%phi_cs = undef
        g%IS(i,j)%q_e = undef

        g%IS(i,j)%habit = iundef
        g%IS(i,j)%sh_type = iundef
        g%IS(i,j)%is_mod(1) = iundef
        g%IS(i,j)%is_mod(2) = iundef
        g%IS(i,j)%growth_mode=iundef
        g%IS(i,j)%init_growth=iundef

      enddo
      enddo


    end subroutine initialize_shape_ig

  end subroutine ini_group_all

  subroutine update_group_all(gr,gs,ga,ag,&
       level,nacat,update_vapor,mes_rc,&
       qtp,rv2,&
       iupdate_gr,iupdate_gs,iupdate_ga,&
       flagp_r,flagp_s,flagp_a,eps_ap0,ID,JD,KD)
    use scale_prc, only: &
       PRC_abort
    use class_Thermo_Var, only: &
       get_sat_vapor_pres_lk
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Note: time increment has to be specified for each group
    !       before this one is called.
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! level of complexity
    integer, intent(in)  :: level
    integer,intent(in)  :: update_vapor,nacat
    integer,intent(in) :: ID(*),JD(*),KD(*)
    integer,dimension(*),intent(in) :: iupdate_gr,iupdate_gs
    integer,dimension(mxntend,*),intent(in) :: iupdate_ga
    type (group), intent(inout)  :: gr,gs
    type (group), dimension(*),intent(inout) :: ga
    type (airgroup), intent(inout) :: ag
    ! message from reality-check
    integer,dimension(*)   :: mes_rc
    real(MP_KIND),dimension(*),intent(in) :: qtp
    real(PS),dimension(*),intent(inout) :: rv2
    real(PS),dimension(*),intent(in)  :: eps_ap0
    integer,intent(in)  :: flagp_r,flagp_s,flagp_a
    ! new space
    real(8),dimension(mxnbin,LMAX)     :: ndxdt,ndxdt2
    !real(PS) :: ndxdt3
    real (ps)     :: var_n1,var_n2
    ! original mass concentration
    real(8),dimension(mxnbin,LMAX)  :: om
    ! original number concentration
    real(8),dimension(mxnbin,LMAX)  :: oc
    ! original total volume variables
    real(8),dimension(mxnbin,LMAX)  :: oQ
    real(8),dimension(mxnbin,LMAX)  :: oq_alen!,phi_dbg,psi_dbg
    real(8),dimension(mxnbin,LMAX)  :: oq_clen,ndxdt_alen,ndxdt_clen

    ! mass fraction of aerosol particles to total mass of hydrometeors and
    ! mass fraction of soluble aerosol particles to total aerosol mass
    ! those are used to maintain positive aerosol masses, and use fraction
    ! of smaller bin.
    real(ps),dimension(max(gr%n_bin,gs%n_bin),LMAX) :: fap
    !real(ps) :: fapt_r,faps_r,fapt_s,faps_s, fap3
    ! minimum possible fractions
    real(ps),parameter :: min_fapt_r=1.0e-18_ps,min_faps_r=1.0e-5_ps&
!!c    real(ps),parameter :: min_fapt_r=1.0e-18_ps,min_faps_r=0.0_ps&
         ,min_fapt_s=1.0e-18_ps,min_faps_s=0.0_ps,sep_faps=1.0e-5_PS
    ! possible ratio of mean mass to bin boundaries
    real(ps),parameter :: brat1=1.001,brat2=0.999
    real(ps),parameter :: mlmt=1.0e-30,nlmt=1.0e-30,m_lmt_ap=1.0e-25
    !     the minimum possible background concentration for large particle
    !     is assumed to be 1.0e-5 cm^-3
    real(ps),parameter :: n_lmt_ap=1.0e-5,n_max_ap=1.0e+4
!parcel model    real(ps),parameter :: n_lmt_ap=1.0e-15,n_max_ap=1.0e+4
    !     The minimum radius possible for accumulation particles
    real(ps),parameter :: r3_lmt=1.0e-18
    !     The minimum radius possible for nucleation particles
    real(PS),parameter :: r3_lmt_nuc=1.0e-21


!org    real(ps),parameter :: mx_aprat=0.99999
    real(ps),parameter :: mx_aprat=0.95
    real(ps) :: m_lmt
    ! realistic bounds for length predictions
    real(ps),parameter :: min_hexlen=1.0e-4,max_hexlen=3.0,&
                          max_roslen=3.0,max_irrlen=3.0,max_exice=1.0
    ! realistic bounds for crystal mass
    real(ps),parameter :: m_icmin=4.763209003e-12

    real(PS) :: tmass1,tmass2!,tmass3
    real(PS) :: phi,psi,pag,pcg,den_ip1
    real(PS) :: a_min,c_min,d_min,ag_min,cg_min

    ! minimum, and maximum volume of circumscribing sphere that corresponds to 1 um of
    ! hex ice crystal and sphere of 10cm radius
    real(PS),parameter :: V_csmin_hex=1.184768784e-11,V_csmax=4188.79020478639
!!c    real,parameter :: V_csmin_hex=1.184768784e-11,V_csmax=4188.79020478639e+10


    ! maximum and minimum aspect ratio
    real(PS),parameter  :: phi_max=2.0e+1,phi_min=5.0e-3

    ! accuracy limit to prevent dubious crystal features by numerical errors
    real(PS),parameter :: accuracy_lmt=1.0e-6

    ! vapor specific humidity cacluated with activation scheme
    real(PS) :: rv_n
    ! total mixing ratio of hydrometeors
    real(PS),dimension(LMAX) :: totmixr,totmixs
    real(PS),dimension(LMAX) :: s_n
    real(PS) :: e_n

    real(8) :: dtd,con8

!    integer,dimension(mxnbin*LMAX*mxntend) :: ierror
    integer,dimension(mxnbin,LMAX) :: icond1,icond2,icond3 !,ierror2,ierror3

    integer  :: i,j,k,n,ic,L!,ij,var_status
    !integer :: ie
    !integer :: ic1,ic2
    integer :: iq1

    !integer :: imark

    integer  :: maxice_kidx, maxice_bidx, maxice_kidx2, maxice_bidx2
    real(PS) :: maxice_s, maxice_s2, maxice_d(gs%n_tendpros), maxice_m(16)

    L=max(gr%L,gs%L,ga(1)%L)


    do n=1,L
      totmixs(n)=0.0_PS
      totmixr(n)=0.0_PS
    enddo

    !
    ! +++ for rain +++
    !
    if( flagp_r > 0 ) then

      dtd=gr%dt

!      do in=1,gr%n_bin*L
!        n=(in-1)/gr%N_BIN+1
!        i=in-(n-1)*gr%N_BIN
      do n = 1, L
      do i = 1, gr%N_BIN
        ndxdt(i,n) = 0.0d+0
        ndxdt2(i,n) = 0.0d+0
      enddo
      enddo

      do j=1,gr%n_tendpros
        if(iupdate_gr(j)==0) cycle
!        do in=1,gr%n_bin*L
!          n=(in-1)/gr%N_BIN+1
!          i=in-(n-1)*gr%N_BIN
        do n = 1, L
        do i = 1, gr%N_BIN
          ndxdt(i,n) = ndxdt(i,n)+gr%MS(i,n)%dmassdt(rmt,j)
          ndxdt2(i,n) = ndxdt2(i,n)+gr%MS(i,n)%dcondt(j)
        enddo
        enddo
      enddo

      if ( debug ) then

!      ierror(1:gr%n_bin*gr%n_tendpros*L)=0
!      do ijn=1,gr%n_bin*gr%n_tendpros*L
!        in=(ijn-1)/gr%n_tendpros+1
!        j=ijn-(in-1)*gr%n_tendpros
!        n=(in-1)/gr%n_bin+1
!        i=in-(n-1)*gr%n_bin
         do n = 1, L
         do i = 1, gr%n_bin
         do j = 1, gr%n_tendpros

            if(iupdate_gr(j)==0) cycle

            if(gr%MS(i,n)%dmassdt(rmt,j)>1.0e+03.or.&
               gr%MS(i,n)%dmassdt(rmt,j)<-1.0e+03) then
               write(*,201) n,i,j,gr%MS(i,n)%dmassdt(rmt,j)
201            format("Warning mt_rain>grid,bin,process,dmassdt",3i5,es15.6)
            endif
         enddo
         enddo
         enddo

!      ierror(1:gr%n_bin*gr%n_tendpros*L)=0
!      do ijn=1,gr%n_bin*gr%n_tendpros*L
!        in=(ijn-1)/gr%n_tendpros+1
!        j=ijn-(in-1)*gr%n_tendpros
!        n=(in-1)/gr%n_bin+1
!        i=in-(n-1)*gr%n_bin
         do n = 1, L
         do i = 1, gr%n_bin
         do j = 1, gr%n_tendpros

            if(iupdate_gr(j)==0) cycle

            if(gr%MS(i,n)%dcondt(j)<-1.0e+8 .or. &
               gr%MS(i,n)%dcondt(j)>1.0e+8 ) then
               write(*,*) "Warning con_r:bin,grid,process",i,n,j,gr%MS(i,n)%dcondt(j)
            endif
         enddo
         enddo
         enddo

      end if

!      do in=1,gr%n_bin*L
!        n=(in-1)/gr%N_BIN+1
!        i=in-(n-1)*gr%N_BIN
      do n = 1, L
      do i = 1, gr%N_BIN
        om(i,n)=gr%MS(i,n)%mass(rmt)
        gr%MS(i,n)%mass(rmt)=gr%MS(i,n)%mass(rmt)+ndxdt(i,n)*dtd

        oc(i,n)=gr%MS(i,n)%con
        gr%MS(i,n)%con=gr%MS(i,n)%con+ndxdt2(i,n)*dtd
      enddo
      enddo
      !
      ! check
      !
      icond1(1:gr%n_bin,1:L)=0
      icond2(1:gr%n_bin,1:L)=0
!      do in=1,gr%n_bin*L
!        n=(in-1)/gr%N_BIN+1
!        i=in-(n-1)*gr%N_BIN
      do n = 1, L
      do i = 1, gr%N_BIN
        if(gr%MS(i,n)%mass(rmt)<1.0e-30.or.gr%MS(i,n)%mass(rmt)<om(i,n)*1.0e-6) then
          ! case of zero out.  Due to single-presicion of mass variables.
          icond1(i,n)=1
        elseif(gr%MS(i,n)%mass(rmt)>1.0e-1) then
           icond1(i,n)=2
        end if

        if(gr%MS(i,n)%con<1.0e-30.or.gr%MS(i,n)%con<oc(i,n)*1.0e-6) then
          ! case of zero out.  Due to single-presicion of con variable.
          icond2(i,n)=1
        else
          icond2(i,n)=3
        endif
      enddo
      enddo

!      do in=1,gr%n_bin*L
!        n=(in-1)/gr%N_BIN+1
!        i=in-(n-1)*gr%N_BIN
      do n = 1, L
      do i = 1, gr%N_BIN
        if(icond1(i,n)==1.or.icond2(i,n)==1) then
          gr%MS(i,n)%con=0.0_PS
          gr%MS(i,n)%mass(rmt)=0.0_PS
        endif
      end do
      end do
!      do in=1,gr%n_bin*L
!        n=(in-1)/gr%N_BIN+1
!        i=in-(n-1)*gr%N_BIN
      do n = 1, L
      do i = 1, gr%N_BIN
        if(icond1(i,n)/=1.and.icond2(i,n)==3) then
          gr%MS(i,n)%con= &
                  gr%MS(i,n)%mass(rmt)&
                  /max(brat1*gr%binb(i),min(gr%MS(i,n)%mass(rmt)/gr%MS(i,n)%con&
                 ,brat2*gr%binb(i+1)))
        endif
      enddo
      enddo

      if ( debug ) then
!        do in=1,gr%n_bin*L
!          n=(in-1)/gr%N_BIN+1
!          i=in-(n-1)*gr%N_BIN
         do n = 1, L
         do i = 1, gr%N_BIN
            if(icond1(i,n)==2) then
               write(*,251) i,n,gr%MS(i,n)%mass(rmt),gr%MS(i,n)%con
251            format("warning:qr is large at bin,grid:",2i5,2es15.6)
            endif
         end do
         end do
      end if

      if(level>=4) then
!        do in=1,gr%n_bin*L
!          n=(in-1)/gr%N_BIN+1
!          i=in-(n-1)*gr%N_BIN
         do n = 1, L
         do i = 1, gr%N_BIN
          ndxdt(i,n) = 0.0d+0
          ndxdt2(i,n) = 0.0d+0
        enddo
        enddo
        do j=1,gr%n_tendpros
          if(iupdate_gr(j)==0) cycle
!          do in=1,gr%n_bin*L
!            n=(in-1)/gr%N_BIN+1
!            i=in-(n-1)*gr%N_BIN
          do n = 1, L
          do i = 1, gr%N_BIN
            ! mass by total mass of aerosols (g/g)
            ndxdt(i,n)=ndxdt(i,n)+gr%MS(i,n)%dmassdt(rmat,j)

            ! mass by soluble mass of aerosols (g/g)
            ndxdt2(i,n)=ndxdt2(i,n)+gr%MS(i,n)%dmassdt(rmas,j)
          enddo
          enddo
        enddo

        if ( debug ) then
!        ierror(1:gr%n_bin*gr%n_tendpros*L)=0
!        do ijn=1,gr%n_bin*gr%n_tendpros*L
!          in=(ijn-1)/gr%n_tendpros+1
!          j=ijn-(in-1)*gr%n_tendpros
!          n=(in-1)/gr%n_bin+1
!          i=in-(n-1)*gr%n_bin
           do n = 1, L
           do i = 1, gr%n_bin
           do j = 1, gr%n_tendpros
              if(iupdate_gr(j)==0) cycle
              if(gr%MS(i,n)%dmassdt(rmat,j)>1.0e+08.or.&
                 gr%MS(i,n)%dmassdt(rmat,j)<-1.0e+08)then
                 write(*,236) i,n,j,gr%MS(i,n)%dmassdt(rmat,j)
236              format("Warning aptm in rain>bin,grid,process,dmassdt",3i5,es15.6)
              endif
           enddo
           enddo
           enddo

!        ierror(1:gr%n_bin*gr%n_tendpros*L)=0
!        do ijn=1,gr%n_bin*gr%n_tendpros*L
!          in=(ijn-1)/gr%n_tendpros+1
!          j=ijn-(in-1)*gr%n_tendpros
!          n=(in-1)/gr%n_bin+1
!          i=in-(n-1)*gr%n_bin
           do n = 1, L
           do i = 1, gr%n_bin
           do j = 1, gr%n_tendpros
              if(iupdate_gr(j)==0) cycle
              if(gr%MS(i,n)%dmassdt(rmas,j)>1.0e+08.or.&
                 gr%MS(i,n)%dmassdt(rmas,j)<-1.0e+08)then
                 write(*,237) i,n,j,gr%MS(i,n)%dmassdt(rmas,j)
237              format("Warning apsm in solid>bin,grid,process,dmassdt",3i5,es15.6)
              endif
           enddo
           enddo
           enddo
        end if

!        do in=1,gr%n_bin*L
!          n=(in-1)/gr%N_BIN+1
!          i=in-(n-1)*gr%N_BIN
        do n = 1, L
        do i = 1, gr%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            gr%MS(i,n)%mass(rmat)=gr%MS(i,n)%mass(rmat)+ndxdt(i,n)*dtd
            fap(i,n)=gr%MS(i,n)%mass(rmat)/gr%MS(i,n)%mass(rmt)
          endif
        enddo
        enddo
!        do in=1,gr%n_bin*L
!          n=(in-1)/gr%N_BIN+1
!          i=in-(n-1)*gr%N_BIN
        do n = 1, L
        do i = 1, gr%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            if(fap(i,n)>1.0_PS) then
              gr%MS(i,n)%mass(rmt)=gr%MS(i,n)%mass(rmat)/mx_aprat
              gr%MS(i,n)%con= &
                  gr%MS(i,n)%mass(rmt)&
                  /max(brat1*gr%binb(i),min(gr%MS(i,n)%mass(rmt)/gr%MS(i,n)%con&
                 ,brat2*gr%binb(i+1)))

            elseif(fap(i,n)<=min_fapt_r) then

              gr%MS(i,n)%mass(rmat)=min(max(m_lmt_ap, &
                                   min_fapt_r*gr%MS(i,n)%mass(rmt))&
                                  ,mx_aprat*gr%MS(i,n)%mass(rmt))
            else
              gr%MS(i,n)%mass(rmat)=min(mx_aprat*gr%MS(i,n)%mass(rmt), &
                                    max(m_lmt_ap, &
                                    gr%MS(i,n)%mass(rmat)))

            endif

          endif
        enddo
        enddo
!        do in=1,gr%n_bin*L
!          n=(in-1)/gr%N_BIN+1
!          i=in-(n-1)*gr%N_BIN
        do n = 1, L
        do i = 1, gr%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            gr%MS(i,n)%mass(rmas)=gr%MS(i,n)%mass(rmas)+ndxdt2(i,n)*dtd
            fap(i,n)=gr%MS(i,n)%mass(rmas)/gr%MS(i,n)%mass(rmat)
          endif
        enddo
        enddo
!        do in=1,gr%n_bin*L
!          n=(in-1)/gr%N_BIN+1
!          i=in-(n-1)*gr%N_BIN
        do n = 1, L
        do i = 1, gr%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            if(fap(i,n)<min_faps_r) then
              gr%MS(i,n)%mass(rmas)=min(max(m_lmt_ap,&
                                 min_faps_r*gr%MS(i,n)%mass(rmat)) &
                                ,gr%MS(i,n)%mass(rmat))
            endif
          endif
        enddo
        enddo

      endif

      do i=1,gr%n_bin
        do n=1,L
          totmixr(n)=totmixr(n)+max(0.0_PS,gr%MS(i,n)%mass(rmt)-gr%MS(i,n)%mass(rmat))
        enddo
      end do

      if(debug) then
        do n=1,L
!         if(sum(XR(rcon_q,1:gr%n_bin,nrcat,n))*ag%tv(n)%den>100.0) then
         if(sum(gr%MS(1:gr%n_bin,n)%con)>0.0) then
         write(*,'("update_group:bin,i,j,k,total rcon, others",4i5,100es15.6)') n,ID(n),JD(n),KD(n) &
                      ,sum(gr%MS(1:gr%n_bin,n)%con) &
                      ,gr%MS(1:gr%n_bin,n)%con

         write(*,'("update_group:bin,i,j,k,total act, total vap mass tendency",4i5,100es15.6)') n,ID(n),JD(n),KD(n) &
                      ,sum(gr%MS(1:gr%n_bin,n)%dmassdt(rmt,3)) &  ! activation
                      ,sum(gr%MS(1:gr%n_bin,n)%dmassdt(rmt,1))    ! vapor
         write(*,*) "cal_model_tendency:con-bin,i,j,k,total rim, total melt, total mossop, total coal"
         write(*,'(4I5,10ES15.6)') n,ID(n),JD(n),KD(n) &
                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(4)) &  ! riming
                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(10)) &  ! melting_shedding
                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(9)) &  ! mossop and hallet
                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(2))    ! collision-coal
         endif
       enddo
      endif

    endif
    !
    ! ice
    !
    if( flagp_s > 0 ) then

      dtd=gs%dt

!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        ndxdt(i,n) = 0.0d+0
        ndxdt2(i,n) = 0.0d+0
      enddo
      enddo
      maxice_d = 0.0D0
      do j=1,gs%n_tendpros
         if(iupdate_gs(j)==0) cycle
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
         do n = 1, L
         do i = 1, gs%N_BIN
            ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(imt,j)
            ndxdt2(i,n)=ndxdt2(i,n)+gs%MS(i,n)%dcondt(j)
            if (abs(gs%MS(i,n)%dmassdt(imt,j)) > abs(maxice_d(j))) then
               maxice_d(j) = gs%MS(i,n)%dmassdt(imt,j)
            endif
         enddo
         enddo
      enddo

      maxice_s = 0.0D0
      maxice_s2 = 0.0D0
      maxice_m = 0.0D0
      maxice_kidx = 0
      maxice_bidx = 0
      maxice_kidx2 = 0
      maxice_bidx2 = 0
!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        if (abs(ndxdt2(i,n)) > abs(maxice_s)) then
          maxice_s = ndxdt2(i,n)
          maxice_kidx = n
          maxice_bidx = i
        endif
        if (abs(gs%MS(i,n)%mass(imt)) > abs(maxice_s2)) then
          maxice_s2 = gs%MS(i,n)%mass(imt)
          maxice_kidx2 = n
          maxice_bidx2 = i
        endif
        if (abs(gs%MS(i,n)%mass(imt)) > abs(maxice_m(1))) then
          maxice_m(1) = gs%MS(i,n)%mass(imt)
        endif
        if (abs(gs%MS(i,n)%mass(imr)) > abs(maxice_m(2))) then
          maxice_m(2) = gs%MS(i,n)%mass(imr)
        endif
        if (abs(gs%MS(i,n)%mass(ima)) > abs(maxice_m(3))) then
          maxice_m(3) = gs%MS(i,n)%mass(ima)
        endif
        if (abs(gs%MS(i,n)%mass(imc)) > abs(maxice_m(4))) then
          maxice_m(4) = gs%MS(i,n)%mass(imc)
        endif
        if (abs(gs%MS(i,n)%mass(imw)) > abs(maxice_m(5))) then
          maxice_m(5) = gs%MS(i,n)%mass(imw)
        endif
        if (abs(gs%MS(i,n)%mass(imf)) > abs(maxice_m(6))) then
          maxice_m(6) = gs%MS(i,n)%mass(imf)
        endif
        if (abs(gs%MS(i,n)%mass(imat)) > abs(maxice_m(7))) then
          maxice_m(7) = gs%MS(i,n)%mass(imat)
        endif
        if (abs(gs%MS(i,n)%mass(imas)) > abs(maxice_m(8))) then
          maxice_m(8) = gs%MS(i,n)%mass(imas)
        endif
        if (abs(gs%MS(i,n)%con) > abs(maxice_m(9))) then
          maxice_m(9) = gs%MS(i,n)%con
        endif
        if (abs(gs%IS(i,n)%v_cs) > abs(maxice_m(10))) then
          maxice_m(10) = gs%IS(i,n)%v_cs
        endif
        if (abs(gs%MS(i,n)%a_len) > abs(maxice_m(11))) then
          maxice_m(11) = gs%MS(i,n)%a_len
        endif
        if (abs(gs%MS(i,n)%c_len) > abs(maxice_m(12))) then
          maxice_m(12) = gs%MS(i,n)%c_len
        endif
        if (abs(gs%IS(i,n)%d) > abs(maxice_m(13))) then
          maxice_m(13) = gs%IS(i,n)%d
        endif
        if (abs(gs%IS(i,n)%ag) > abs(maxice_m(14))) then
          maxice_m(14) = gs%IS(i,n)%ag
        endif
        if (abs(gs%IS(i,n)%cg) > abs(maxice_m(15))) then
          maxice_m(15) = gs%IS(i,n)%cg
        endif
        if (abs(gs%IS(i,n)%n_exice) > abs(maxice_m(16))) then
          maxice_m(16) = gs%IS(i,n)%n_exice
        endif
      enddo
      enddo

      !if (maxice_bidx2 /= 0 .and. maxice_kidx2 /= 0 ) then
      !  write(fid_alog,'(a,4I5,4ES22.12)') "CONT", JD(1), ID(1), maxice_bidx, maxice_kidx, maxice_s, maxice_s*dtd, maxice_s2, gs%MS(maxice_bidx2,maxice_kidx2)%con
        !write(fid_alog,'(a,12ES16.9)') "CONT1", maxice_d
      !  write(fid_alog,'(a,8ES16.9)') "CONT1", maxice_m(1:8)
      !  write(fid_alog,'(a,8ES16.9)') "CONT2", maxice_m(9:16)
      !endif

!      do n=1,L
!        write(*,*) "ck ice,c,m,dm1,dm11:",n, &
!           sum(gs%MS(1:gs%n_bin,n)%con),sum(gs%MS(1:gs%n_bin,n)%mass(1)), &
!           sum(gs%MS(1:gs%n_bin,n)%dmassdt(imt,1)),sum(gs%MS(1:gs%n_bin,n)%dmassdt(imt,11))
!      enddo

      if ( debug ) then
!      ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!      do ijn=1,gs%n_bin*gs%n_tendpros*L
!        in=(ijn-1)/gs%n_tendpros+1
!        j=ijn-(in-1)*gs%n_tendpros
!        n=(in-1)/gs%n_bin+1
!        i=in-(n-1)*gs%n_bin
         do n = 1, L
         do i = 1, gs%n_bin
         do j = 1, gs%n_tendpros
            if(iupdate_gs(j)==0) cycle
            if(gs%MS(i,n)%dmassdt(imt,j)>1.0e+08.or.&
               gs%MS(i,n)%dmassdt(imt,j)<-1.0e+08)then
               write(*,202) n,i,j,gs%MS(i,n)%dmassdt(imt,j)
202            format("Warning mt_s>grid,bin,process,dmassdt",3i5,es15.6)
            endif
         enddo
         enddo
         enddo

!      ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!      do ijn=1,gs%n_bin*gs%n_tendpros*L
!        in=(ijn-1)/gs%n_tendpros+1
!        j=ijn-(in-1)*gs%n_tendpros
!        n=(in-1)/gs%n_bin+1
!        i=in-(n-1)*gs%n_bin
         do n = 1, L
         do i = 1, gs%n_bin
         do j = 1, gs%n_tendpros
            if(iupdate_gs(j)==0) cycle
            if(gs%MS(i,n)%dcondt(j)<-1.0e+8 .or. &
               gs%MS(i,n)%dcondt(j)>1.0e+8 ) then
               write(*,*) "Warning con_s:bin,grid,process",i,n,j,gs%MS(i,n)%dcondt(j)
            endif
         enddo
         enddo
         enddo
      end if

!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        om(i,n)=gs%MS(i,n)%mass(imt)
        gs%MS(i,n)%mass(imt)=gs%MS(i,n)%mass(imt)+ndxdt(i,n)*dtd
        oc(i,n)=gs%MS(i,n)%con
        gs%MS(i,n)%con=gs%MS(i,n)%con+ndxdt2(i,n)*dtd
      enddo
      enddo


!dbg      do n=1,gs%L
!dbg        write(*,'("ck update_group con_s",I5,30ES15.6)') n,(gs%MS(i,n)%con,i=1,gs%N_BIN)
!dbg      enddo
!!!!CDIR NOVECTOR
!!!      do in=1,gs%n_bin*L
!!!        n=(in-1)/gs%N_BIN+1
!!!        i=in-(n-1)*gs%N_BIN
!!!        if(kd(n).eq.21.and.jd(n).eq.4.and.i.eq.3) then
!!!          write(*,*) "ck con1",oc(i,n),gs%MS(i,n)%con,ndxdt2(i,n)
!!!        endif
!!!      enddo
      !
      ! check
      !
      icond1(1:gs%n_bin,1:L)=0
      icond2(1:gs%n_bin,1:L)=0
!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        if(gs%MS(i,n)%mass(imt)<1.0e-30.or.gs%MS(i,n)%mass(imt)<om(i,n)*1.0e-6) then
          ! case of zero out.  Due to single-presicion of mass variables.
          icond1(i,n)=1
        elseif(gs%MS(i,n)%mass(imt)>1.0e-1) then
          icond1(i,n)=2
        endif
        if(gs%MS(i,n)%con<1.0e-30.or.gs%MS(i,n)%con<oc(i,n)*1.0e-6) then
          ! case of zero out.  Due to single-presicion of con variable.
          icond2(i,n)=1
        else
          icond2(i,n)=3
        endif
      enddo
      enddo
!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        if(icond1(i,n)==1.or.icond2(i,n)==1) then
          gs%MS(i,n)%con=0.0_PS
          gs%MS(i,n)%mass(imt)=0.0_PS
        endif
      end do
      end do
!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        if(icond1(i,n)/=1.and.icond2(i,n)==3) then
           gs%MS(i,n)%con=&
               gs%MS(i,n)%mass(imt)&
               /max(brat1*gs%binb(i),&
                min(gs%MS(i,n)%mass(imt)/gs%MS(i,n)%con&
                 ,brat2*gs%binb(i+1)))
        endif
      enddo
      enddo
      if ( debug ) then
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
         do n = 1, L
         do i = 1, gs%N_BIN
            if(icond1(i,n)==2) then
               write(*,252) i,n,gs%MS(i,n)%mass(imt)
252            format("warning:qs is large at bin,grid:",2i5,2es15.6)
            endif
         enddo
         enddo
      end if

      ! volume of circumscribing sphere (cm^3/g) and
      ! volume of axis length production (cm^3/g)
      iq1=ivcs_q
      nonmass_loop: do k=1,gs%N_nonmass
        !      1. volume of circumscribing sphere (cm^3/g)
        !      2. con. weighted a-axis length^3 (cm^3/g)
        !      3. con. weighted c-axis length^3 (cm^3/g)
        !      4. con. weighted d-axis length^3 (cm^3/g)
        !      5. con. weighted ag-axis length^3 (cm^3/g) (polycrystals)
        !      6. con. weighted cg-axis length^3 (cm^3/g) (polycrystals)
        !      7. con. weighted # of poly        (cm^3/g) (polycrystals)

        if(k==ivcs) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            oq(i,n)=oc(i,n)*gs%IS(i,n)%v_cs
          enddo
          enddo
        elseif(k==iacr) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            oq(i,n)=oc(i,n)*gs%MS(i,n)%a_len*gs%MS(i,n)%a_len*gs%MS(i,n)%a_len
            oq_alen(i,n)=oc(i,n)*gs%MS(i,n)%a_len*gs%MS(i,n)%a_len*gs%MS(i,n)%a_len
          enddo
          enddo
        elseif(k==iccr) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            oq(i,n)=oc(i,n)*gs%MS(i,n)%c_len*gs%MS(i,n)%c_len*gs%MS(i,n)%c_len
            oq_clen(i,n)=oc(i,n)*gs%MS(i,n)%c_len*gs%MS(i,n)%c_len*gs%MS(i,n)%c_len
          enddo
          enddo
        elseif(k==idcr) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            oq(i,n)=oc(i,n)*gs%IS(i,n)%d*gs%IS(i,n)%d*gs%IS(i,n)%d
          enddo
          enddo
        elseif(k==iag) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            oq(i,n)=oc(i,n)*gs%IS(i,n)%ag*gs%IS(i,n)%ag*gs%IS(i,n)%ag
          enddo
          enddo
        elseif(k==icg) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            oq(i,n)=oc(i,n)*gs%IS(i,n)%cg*gs%IS(i,n)%cg*gs%IS(i,n)%cg
          enddo
          enddo
        elseif(k==inex) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            oq(i,n)=oc(i,n)*gs%IS(i,n)%n_exice
          enddo
          enddo
        endif

        ndxdt(:,:) = 0.0_DP
        do j=1,gs%n_tendpros
          if(iupdate_gs(j)==0) cycle
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
          do n = 1, L
          do i = 1, gs%n_bin
            ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dvoldt(k,j)
          enddo
          enddo
        enddo

        if ( debug ) then
!        ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!        do ijn=1,gs%n_bin*gs%n_tendpros*L
!          in=(ijn-1)/gs%n_tendpros+1
!          j=ijn-(in-1)*gs%n_tendpros
!          n=(in-1)/gs%n_bin+1
!          i=in-(n-1)*gs%n_bin
           do n = 1, L
           do i = 1, gs%n_bin
           do j = 1, gs%n_tendpros
              if( ( iupdate_gs(j) .ne. 0 ) &
                   .and. (    gs%MS(i,n)%dvoldt(k,j)>1.0e+20_RP &
                         .or. gs%MS(i,n)%dvoldt(k,j)<-1.0e+20_RP ) ) then
                 write(*,231) n,k,j,gs%MS(i,n)%dvoldt(k,j)
231              format("axis>grid,axis#,process,dvoldt",3i5,es15.6)
                 write(*,'("phi_cs,den",2es15.6)') gs%IS(i,n)%phi_cs,gs%MS(i,n)%den
                 write(*,'("psi_ic,phi_ic,alen",3es15.6)') gs%IS(i,n)%psi_ic,gs%IS(i,n)%phi_ic,gs%MS(i,n)%a_len
                 write(*,'("i,ht,sh",3i4)') i,gs%IS(i,n)%habit,gs%IS(i,n)%sh_type

              endif
           enddo
           enddo
           enddo
        end if


        if(k==ivcs) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN

            icond3(i,n)=0
            if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
              con8=gs%MS(i,n)%con
              gs%IS(i,n)%v_cs=max(0.0_DS,(oq(i,n)+ndxdt(i,n)*dtd)/con8)
              icond3(i,n)=1
              if(gs%IS(i,n)%v_cs==0.0_PS) then
                icond3(i,n)=2
              elseif(gs%IS(i,n)%v_cs<V_csmin_hex) then
                icond3(i,n)=3
                gs%IS(i,n)%V_cs=V_csmin_hex
              elseif(gs%IS(i,n)%v_cs>V_csmax) then
                icond3(i,n)=4
                gs%IS(i,n)%V_cs=V_csmax
              endif
            else
              gs%IS(i,n)%v_cs=0.0_PS
            endif
          enddo
          enddo
        elseif(k==iacr) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
              ndxdt_alen(i,n)=gs%MS(i,n)%a_len
              con8=gs%MS(i,n)%con
              gs%MS(i,n)%a_len=(max(0.0_DS,oq(i,n)+ndxdt(i,n)*dtd)/con8)**(1.0/3.0)

            else
              gs%MS(i,n)%a_len=0.0
            endif
          enddo
          enddo
        elseif(k==iccr) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
              ndxdt_clen(i,n)=gs%MS(i,n)%c_len
              con8=gs%MS(i,n)%con
              gs%MS(i,n)%c_len=(max(0.0_DS,oq(i,n)+ndxdt(i,n)*dtd)/con8)**(1.0/3.0)
            else
              gs%MS(i,n)%c_len=0.0
            endif
          enddo
          enddo
        elseif(k==idcr) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
           if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
              con8=gs%MS(i,n)%con
              gs%IS(i,n)%d=(max(0.0_DS,oq(i,n)+ndxdt(i,n)*dtd)/con8)**(1.0/3.0)
            else
              gs%IS(i,n)%d=0.0_PS
            endif
          enddo
          enddo
        elseif(k==iag) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
              con8=gs%MS(i,n)%con
              gs%IS(i,n)%ag=(max(0.0_DS,oq(i,n)+ndxdt(i,n)*dtd)/con8)**(1.0/3.0)
            else
              gs%IS(i,n)%ag=0.0_PS
            endif
          enddo
          enddo
        elseif(k==icg) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
              con8=gs%MS(i,n)%con
              gs%IS(i,n)%cg=(max(0.0_DS,oq(i,n)+ndxdt(i,n)*dtd)/con8)**(1.0/3.0)
            else
              gs%IS(i,n)%cg=0.0_PS
            endif
          enddo
          enddo
        elseif(k==inex) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
              con8=gs%MS(i,n)%con
              gs%IS(i,n)%n_exice=max(0.0_DS,oq(i,n)+ndxdt(i,n)*dtd)/con8
            else
              gs%IS(i,n)%n_exice=0.0_PS
            endif
          enddo
          enddo
        endif

      enddo nonmass_loop

!dbg      call check_phi('af_nonmass_loop')

      ! bound with maximum and minium possible length.
!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
          if(gs%MS(i,n)%a_len<=min_hexlen.or.&
             gs%MS(i,n)%c_len<=min_hexlen) then
            gs%MS(i,n)%a_len=min_hexlen
            gs%MS(i,n)%c_len=min_hexlen
          elseif(gs%MS(i,n)%a_len>=max_hexlen.or.&
             gs%MS(i,n)%c_len>=max_hexlen) then
            gs%MS(i,n)%a_len=max_hexlen
            gs%MS(i,n)%c_len=max_hexlen
          endif
!org          gs%MS(i,n)%a_len=max(min_hexlen,&
!org                           min(max_hexlen,gs%MS(i,n)%a_len))
!org          gs%MS(i,n)%c_len=max(min_hexlen,&
!org                           min(max_hexlen,gs%MS(i,n)%c_len))

          !
          ! d axis lengths of ice crystals should be bounded by a-axis length
          !
          gs%IS(i,n)%d=max(0.0_RP,min(0.9_PS*gs%MS(i,n)%a_len,gs%IS(i,n)%d))

          gs%IS(i,n)%ag=max(0.0_RP,min(max_roslen,gs%IS(i,n)%ag))
          gs%IS(i,n)%cg=max(0.0_RP,min(max_irrlen,gs%IS(i,n)%cg))
          gs%IS(i,n)%n_exice=max(0.0_RP,min(max_exice,gs%IS(i,n)%n_exice))

        endif
      enddo
      enddo

!dbg      call check_phi('af_con1')

      ! another constraints on d-axis length and a and c coordinates according to
      ! the a and c axis lengths
!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        if(gs%MS(i,n)%a_len<=1.0e-3_PS.or.&
           gs%MS(i,n)%c_len<=1.0e-4_PS) then
          gs%IS(i,n)%d=0.0_PS
        end if
        if(gs%MS(i,n)%a_len<=min_hexlen.or.&
           gs%MS(i,n)%c_len<=min_hexlen) then
          gs%IS(i,n)%ag=0.0_PS
          gs%IS(i,n)%cg=0.0_PS
          gs%IS(i,n)%n_exice=0.0_PS
        end if
      enddo
      enddo

!dbg      call check_phi('af_massup1')
      !
      ! mass by ice crystals (g/g)
      !
!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        ndxdt(i,n) = 0.0d+0
      enddo
      enddo
      do j=1,gs%n_tendpros
        if(iupdate_gs(j)==0) cycle
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
        do n = 1, L
        do i = 1, gs%N_BIN
          ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(imc,j)
        enddo
        enddo
      enddo

      if ( debug ) then
!      ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!      do ijn=1,gs%n_bin*gs%n_tendpros*L
!        in=(ijn-1)/gs%n_tendpros+1
!        j=ijn-(in-1)*gs%n_tendpros
!        n=(in-1)/gs%n_bin+1
!        i=in-(n-1)*gs%n_bin
         do n = 1, L
         do i = 1, gs%n_bin
         do j = 1, gs%n_tendpros
            if(iupdate_gs(j)==0) cycle
            if(gs%MS(i,n)%dmassdt(imc,j)>1.0e+08.or.&
               gs%MS(i,n)%dmassdt(imc,j)<-1.0e+08)then
               write(*,233) n,j,gs%MS(i,n)%dmassdt(imc,j)
233            format("Warning icemass>grid,process,dmassdt",2i5,es15.6)
            endif
         enddo
         enddo
         enddo
      end if

!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
          gs%MS(i,n)%mass(imc)=min(gs%MS(i,n)%mass(imt),&
                     max(gs%MS(i,n)%mass(imc)+real(ndxdt(i,n)*dtd,PS),&
                     m_icmin*gs%MS(i,n)%con))
        endif
      enddo
      enddo

!dbg      call check_phi('af_massup2')

      !
      ! mass by riming process (g/g)
      !
!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        ndxdt(i,n) = 0.0d+0
      enddo
      enddo
      do j=1,gs%n_tendpros
        if(iupdate_gs(j)==0) cycle
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
        do n = 1, L
        do i = 1, gs%N_BIN
           ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(imr,j)
        enddo
        enddo
      enddo

      if ( debug ) then
!      ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!      do ijn=1,gs%n_bin*gs%n_tendpros*L
!        in=(ijn-1)/gs%n_tendpros+1
!        j=ijn-(in-1)*gs%n_tendpros
!        n=(in-1)/gs%n_bin+1
!        i=in-(n-1)*gs%n_bin
         do n = 1, L
         do i = 1, gs%n_bin
         do j = 1, gs%n_tendpros
            if(iupdate_gs(j)==0) cycle
            if(gs%MS(i,n)%dmassdt(imr,j)>1.0e+03.or.&
               gs%MS(i,n)%dmassdt(imr,j)<-1.0e+03)then
               write(*,232) n,j,gs%MS(i,n)%dmassdt(imr,j)
232            format("Warning rimemass>grid,process,dmassdt",2i5,es15.6)
            endif
         enddo
         enddo
         enddo
      end if

!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
 !       i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
          gs%MS(i,n)%mass(imr)=min(gs%MS(i,n)%mass(imt)-gs%MS(i,n)%mass(imc),&
                     max(gs%MS(i,n)%mass(imr)+real(ndxdt(i,n)*dtd,ps),0.0_ps))
        endif
      enddo
      enddo
!dbg      call check_phi('af_massup3')
      !
      ! mass by agg process (g/g)
      !
!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        ndxdt(i,n) = 0.0d+0
      enddo
      enddo
      do j=1,gs%n_tendpros
        if(iupdate_gs(j)==0) cycle
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
        do n = 1, L
        do i = 1, gs%N_BIN
          ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(ima,j)
        enddo
        enddo
      enddo

      if ( debug ) then
!      ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!      do ijn=1,gs%n_bin*gs%n_tendpros*L
!        in=(ijn-1)/gs%n_tendpros+1
!        j=ijn-(in-1)*gs%n_tendpros
!        n=(in-1)/gs%n_bin+1
!        i=in-(n-1)*gs%n_bin
         do n = 1, L
         do i = 1, gs%n_bin
         do j = 1, gs%n_tendpros
            if(iupdate_gs(j)==0) cycle
            if(gs%MS(i,n)%dmassdt(ima,j)>1.0e+03.or.&
               gs%MS(i,n)%dmassdt(ima,j)<-1.0e+03)then
               write(*,247) n,j,gs%MS(i,n)%dmassdt(ima,j)
247            format("Warning aggmass>grid,process,dmassdt",2i5,es15.6)
            endif
         enddo
         enddo
         enddo
      end if

!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
          gs%MS(i,n)%mass(ima)=min(gs%MS(i,n)%mass(imt)-gs%MS(i,n)%mass(imc),&
                   max(gs%MS(i,n)%mass(ima)+real(ndxdt(i,n)*dtd,ps),0.0_ps))
        endif
      enddo
      enddo

!dbg      call check_phi('af_massup4')
      if(level>=6) then
        !
        ! melt water mass  (g/g)
        !
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
         do n = 1, L
         do i = 1, gs%N_BIN
          ndxdt(i,n) = 0.0d+0
        enddo
        enddo

        do j=1,gs%n_tendpros
          if(iupdate_gs(j)==0) cycle
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
          do n = 1, L
          do i = 1, gs%N_BIN
            ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(imw,j)
          enddo
          enddo
        enddo

        if ( debug ) then
!        ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!        do ijn=1,gs%n_bin*gs%n_tendpros*L
!          in=(ijn-1)/gs%n_tendpros+1
!          j=ijn-(in-1)*gs%n_tendpros
!          n=(in-1)/gs%n_bin+1
!          i=in-(n-1)*gs%n_bin
           do n = 1, L
           do i = 1, gs%n_bin
           do j = 1, gs%n_tendpros
              if(iupdate_gs(j)==0) cycle
              if(gs%MS(i,n)%dmassdt(imw,j)>1.0e+08.or.&
                 gs%MS(i,n)%dmassdt(imw,j)<-1.0e+08)then
                 write(*,239) n,j,gs%MS(i,n)%dmassdt(imw,j)
239              format("Warning meltmass>grid,process,dmassdt",2i5,es15.6)
              endif
           enddo
           enddo
           enddo
        end if

!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
        do n = 1, L
        do i = 1, gs%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            gs%MS(i,n)%mass(imw)=min(max(gs%MS(i,n)%mass(imt)-&
                       gs%MS(i,n)%mass(imc)-gs%MS(i,n)%mass(imr),0.0_PS),&
                       max(gs%MS(i,n)%mass(imw)+real(ndxdt(i,n)*dtd,ps),0.0_ps))

          endif
        enddo
        enddo
!dbg      call check_phi('af_massup5')
        !
        ! freezing nucleation mass  (g/g)
        !
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
        do n = 1, L
        do i = 1, gs%N_BIN
          ndxdt(i,n) = 0.0d+0
        enddo
        enddo
        do j=1,gs%n_tendpros
          if(iupdate_gs(j)==0) cycle
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
          do n = 1, L
          do i = 1, gs%N_BIN
            ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(imf,j)
          enddo
          enddo
        enddo

        if ( debug ) then
!        ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!        do ijn=1,gs%n_bin*gs%n_tendpros*L
!          in=(ijn-1)/gs%n_tendpros+1
!          j=ijn-(in-1)*gs%n_tendpros
!          n=(in-1)/gs%n_bin+1
!          i=in-(n-1)*gs%n_bin
           do n = 1, L
           do i = 1, gs%n_bin
           do j = 1, gs%n_tendpros
              if(iupdate_gs(j)==0) cycle
              if(gs%MS(i,n)%dmassdt(imf,j)>1.0e+08.or.&
                 gs%MS(i,n)%dmassdt(imf,j)<-1.0e+08)then
                 write(*,249) n,j,gs%MS(i,n)%dmassdt(imf,j)
249              format("Warning freezmass>grid,process,dmassdt",2i5,es15.6)
              endif
           enddo
           enddo
           enddo
        end if

!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
        do n = 1, L
        do i = 1, gs%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            gs%MS(i,n)%mass(imf)=min(gs%MS(i,n)%mass(imc),&
                       max(gs%MS(i,n)%mass(imf)+real(ndxdt(i,n)*dtd,ps),0.0_ps))

          endif
        enddo
        enddo
      end if

!dbg      call check_phi('af_massup6')
      if(level>=4) then
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
         do n = 1, L
         do i = 1, gs%N_BIN
          ndxdt(i,n) = 0.0d+0
          ndxdt2(i,n) = 0.0d+0
        enddo
        enddo
        do j=1,gs%n_tendpros
          if(iupdate_gs(j)==0) cycle
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
          do n = 1, L
          do i = 1, gs%N_BIN
            ! mass by total mass of aerosols (g/g)
            ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(imat,j)

            ! mass by soluble mass of aerosols (g/g)
            ndxdt2(i,n)=ndxdt2(i,n)+gs%MS(i,n)%dmassdt(imas,j)
          enddo
          enddo
        enddo

        if ( debug ) then
!        ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!        do ijn=1,gs%n_bin*gs%n_tendpros*L
!          in=(ijn-1)/gs%n_tendpros+1
!          j=ijn-(in-1)*gs%n_tendpros
!          n=(in-1)/gs%n_bin+1
!          i=in-(n-1)*gs%n_bin
           do n = 1, L
           do i = 1, gs%n_bin
           do j = 1, gs%n_tendpros
              if(iupdate_gs(j)==0) cycle
              if(gs%MS(i,n)%dmassdt(imat,j)>1.0e+03.or.&
                 gs%MS(i,n)%dmassdt(imat,j)<-1.0e+03)then
                 write(*,234) KD(n),ID(n),JD(n),i,j,gs%MS(i,n)%dmassdt(imat,j)
234              format("Warning aptmass>grid(k,i,j),bin,process,dmassdt",5i5,es15.6)
              endif
           enddo
           enddo
           enddo

!        ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!        do ijn=1,gs%n_bin*gs%n_tendpros*L
!          in=(ijn-1)/gs%n_tendpros+1
!          j=ijn-(in-1)*gs%n_tendpros
!          n=(in-1)/gs%n_bin+1
!          i=in-(n-1)*gs%n_bin
           do n = 1, L
           do i = 1, gs%n_bin
           do j = 1, gs%n_tendpros
              if(iupdate_gs(j)==0) cycle
              if(gs%MS(i,n)%dmassdt(imas,j)>1.0e+03.or.&
                 gs%MS(i,n)%dmassdt(imas,j)<-1.0e+03)then
                 write(*,235) KD(n),ID(n),JD(n),i,j,gs%MS(i,n)%dmassdt(imas,j)
235              format("Warning apsmass>grid(k,i,j),bin,process,dmassdt",5i5,es15.6)
              endif
           enddo
           enddo
           enddo
        end if

!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
        do n = 1, L
        do i = 1, gs%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            gs%MS(i,n)%mass(imat)=gs%MS(i,n)%mass(imat)+ndxdt(i,n)*dtd
            fap(i,n)=gs%MS(i,n)%mass(imat)/gs%MS(i,n)%mass(imt)
          endif
        enddo
        enddo
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
        do n = 1, L
        do i = 1, gs%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            if(fap(i,n)>1.0_PS) then
              gs%MS(i,n)%mass(imt)=gs%MS(i,n)%mass(imat)/mx_aprat
              gs%MS(i,n)%con= &
                  gs%MS(i,n)%mass(imt)&
                  /max(brat1*gs%binb(i),min(gs%MS(i,n)%mass(imt)/gs%MS(i,n)%con&
                 ,brat2*gs%binb(i+1)))
            elseif(fap(i,n)<=min_fapt_s) then
              gs%MS(i,n)%mass(imat)=min(max(m_lmt_ap, &
                                   min_fapt_s*gs%MS(i,n)%mass(imt))&
                                  ,mx_aprat*gs%MS(i,n)%mass(imt))
            else
              gs%MS(i,n)%mass(imat)=min(mx_aprat*gs%MS(i,n)%mass(imt), &
                                    max(m_lmt_ap, &
                                    gs%MS(i,n)%mass(imat)))

            endif
            if(gs%MS(i,n)%mass(imt)<1.0e-30.or.gs%MS(i,n)%con<1.0e-30) then
              icond3(i,n)=0
            endif
          endif
        enddo
        enddo
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
        do n = 1, L
        do i = 1, gs%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            gs%MS(i,n)%mass(imas)=gs%MS(i,n)%mass(imas)+ndxdt2(i,n)*dtd
            fap(i,n)=gs%MS(i,n)%mass(imas)/gs%MS(i,n)%mass(imat)
          endif
        enddo
        enddo
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
        do n = 1, L
        do i = 1, gs%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            if(fap(i,n)<min_faps_s) then
              gs%MS(i,n)%mass(imas)=min(max(m_lmt_ap,&
                           min_faps_s*gs%MS(i,n)%mass(imat)) &
                          ,gs%MS(i,n)%mass(imat))
            endif
          endif
        enddo
        enddo

      endif
!dbg      call check_phi('af_massup7')
!
!   constrain particle property variables
!
!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN

!        ierror(in)=0
!        ierror2(in)=0
!        ierror3(in)=0

        if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
          if(gs%MS(i,n)%a_len<=0.0.or.gs%MS(i,n)%c_len<=0.0) then
            call gen_length(gs%MS(i,n)%mass(imc)/gs%MS(i,n)%con,&
                    gs%MS(i,n)%a_len,gs%MS(i,n)%c_len)

            gs%IS(i,n)%d=0.0_PS
            gs%IS(i,n)%ag=0.0_PS
            gs%IS(i,n)%cg=0.0_PS
            gs%IS(i,n)%n_exice=0.0_PS
            if(debug) then
               write(*,*) "length generated:",i,KD(n),ID(n),JD(n),gs%MS(i,n)%con,gs%MS(i,n)%mass,gs%MS(i,n)%a_len &
                     ,gs%MS(i,n)%c_len
               write(*,*) "inigr,growthmode,sh,hab",gs%IS(i,n)%init_growth,gs%IS(i,n)%growth_mode &
                    ,gs%IS(i,n)%sh_type,gs%IS(i,n)%habit
               write(*,*) "T,s_v(2)",ag%TV(n)%T,ag%TV(n)%s_v_n(2)
               write(*,*) "oq_alen",oq_alen(i,n),ndxdt_alen(i,n)
               write(*,*) "oq_clen",oq_clen(i,n),ndxdt_clen(i,n)
               write(*,*) "dcondt:",gs%MS(i,n)%dcondt(:)
               write(*,*) "dmasst:",gs%MS(i,n)%dmassdt(1,:)
               write(*,*) "dvoldt_a:",gs%MS(i,n)%dvoldt(iacr,:)
               write(*,*) "dvoldt_c:",gs%MS(i,n)%dvoldt(iccr,:)
            end if
          else
            !
            ! reality check of axis ratio
            !
             if(gs%MS(i,n)%c_len>gs%MS(i,n)%a_len*phi_max) then
                if(debug) then
                   write(*,*) "axis ratio too big:",i,KD(n),ID(n),JD(n),gs%MS(i,n)%con,gs%MS(i,n)%mass,gs%MS(i,n)%a_len &
                     ,gs%MS(i,n)%c_len,gs%IS(i,n)%d
                   write(*,*) "inigr,growthmode,sh,hab",gs%IS(i,n)%init_growth,gs%IS(i,n)%growth_mode &
                    ,gs%IS(i,n)%sh_type,gs%IS(i,n)%habit
                   write(*,*) "T,s_v(2)",ag%TV(n)%T,ag%TV(n)%s_v_n(2)
                   write(*,*) "oq_alen",oq_alen(i,n),ndxdt_alen(i,n)
                   write(*,*) "oq_clen",oq_clen(i,n),ndxdt_clen(i,n)
                   write(*,*) "dcondt:",gs%MS(i,n)%dcondt(:)
                   write(*,*) "dmasst:",gs%MS(i,n)%dmassdt(1,:)
                   write(*,*) "dvoldt_a:",gs%MS(i,n)%dvoldt(iacr,:)
                   write(*,*) "dvoldt_c:",gs%MS(i,n)%dvoldt(iccr,:)
                end if
                gs%MS(i,n)%c_len=gs%MS(i,n)%a_len*phi_max
             elseif(gs%MS(i,n)%c_len<gs%MS(i,n)%a_len*phi_min) then
                if(debug) then
                   write(*,*) "axis ratio too small:",i,KD(n),ID(n),JD(n),gs%MS(i,n)%con,gs%MS(i,n)%mass,gs%MS(i,n)%a_len &
                          ,gs%MS(i,n)%c_len,gs%IS(i,n)%d
                   write(*,*) "inigr,growthmode,sh,hab",gs%IS(i,n)%init_growth,gs%IS(i,n)%growth_mode &
                          ,gs%IS(i,n)%sh_type,gs%IS(i,n)%habit
                   write(*,*) "T,s_v(2)",ag%TV(n)%T,ag%TV(n)%s_v_n(2)
                   write(*,*) "oq_alen",oq_alen(i,n),ndxdt_alen(i,n)
                   write(*,*) "oq_clen",oq_clen(i,n),ndxdt_clen(i,n)
                   write(*,*) "dcondt:",gs%MS(i,n)%dcondt(:)
                   write(*,*) "dmasst:",gs%MS(i,n)%dmassdt(1,:)
                   write(*,*) "dvoldt_a:",gs%MS(i,n)%dvoldt(iacr,:)
                   write(*,*) "dvoldt_c:",gs%MS(i,n)%dvoldt(iccr,:)
                end if
                gs%MS(i,n)%a_len=gs%MS(i,n)%c_len/phi_min
                gs%IS(i,n)%d=max(0.0_RP,min(0.9_PS*gs%MS(i,n)%a_len,gs%IS(i,n)%d))
            end if

            ! reality check of density
            ! calculate minimum lengths, assuming hexagonal plates
            phi=gs%MS(i,n)%c_len/gs%MS(i,n)%a_len
            psi=gs%IS(i,n)%d/gs%MS(i,n)%a_len
            pag=gs%IS(i,n)%ag/gs%MS(i,n)%a_len
            pcg=gs%IS(i,n)%cg/gs%MS(i,n)%c_len

            a_min=(gs%MS(i,n)%mass(imc)/gs%MS(i,n)%con/&
                        (3.0_PS*sq_three*phi*(1.0_RP-psi)*den_i))**0.33333333

!!c                   a_min = ((gs%MS(i,n)%mass(imc)/gs%MS(i,n)%con/(1.0_PS+gs%IS(i,n)%n_exice))/&
!!c                        (3.0_PS*sq_three*phi*den_i))**0.33333333
            if(phi<1.0e-30_RP.or.psi==1.0_RP) then
               LOG_ERROR("update_group_all",*) "a_min wrong",i,KD(n),ID(n),JD(n),gs%MS(i,n)%con,gs%MS(i,n)%mass,gs%MS(i,n)%a_len &
                    ,gs%MS(i,n)%c_len,gs%IS(i,n)%d
!            write(*,*) "oq_alen",oq_alen(i,n),ndxdt_alen(i,n)
!            write(*,*) "oq_clen",oq_clen(i,n),ndxdt_clen(i,n)
!            write(*,*) "dvoldt",gs%MS(i,n)%dvoldt(iacr,:)
!            write(*,*) "dvoldt",gs%MS(i,n)%dvoldt(iccr,:)
               call PRC_abort
            endif
            c_min=a_min*phi
            d_min=a_min*psi
            ag_min=a_min*pag
            cg_min=c_min*pcg

            if(a_min>gs%MS(i,n)%a_len) then
               if ( debug ) then
                  write(*,*) "alen too small:",i,KD(n),ID(n),JD(n),gs%MS(i,n)%con,gs%MS(i,n)%mass,gs%MS(i,n)%a_len &
                       ,gs%MS(i,n)%c_len,gs%IS(i,n)%d
                  write(*,*) "inigr,growthmode,sh,hab",gs%IS(i,n)%init_growth,gs%IS(i,n)%growth_mode &
                       ,gs%IS(i,n)%sh_type,gs%IS(i,n)%habit
                  write(*,*) "T,s_v(2)",ag%TV(n)%T,ag%TV(n)%s_v_n(2)
                  write(*,*) "oq_alen",oq_alen(i,n),ndxdt_alen(i,n)
                  write(*,*) "oq_clen",oq_clen(i,n),ndxdt_clen(i,n)
                  write(*,*) "dcondt:",gs%MS(i,n)%dcondt(:)
                  write(*,*) "dmasst:",gs%MS(i,n)%dmassdt(1,:)
                  write(*,*) "dvoldt_a:",gs%MS(i,n)%dvoldt(iacr,:)
                  write(*,*) "dvoldt_c:",gs%MS(i,n)%dvoldt(iccr,:)
               end if
              gs%MS(i,n)%c_len=c_min
              gs%MS(i,n)%a_len=a_min
              gs%IS(i,n)%d=d_min
              gs%IS(i,n)%ag=ag_min
              gs%IS(i,n)%cg=cg_min
            elseif(c_min>gs%MS(i,n)%c_len) then
               if ( debug ) then
                  write(*,*) "clen too small:",i,KD(n),ID(n),JD(n),gs%MS(i,n)%con,gs%MS(i,n)%mass,gs%MS(i,n)%a_len &
                       ,gs%MS(i,n)%c_len,gs%IS(i,n)%d
                  write(*,*) "inigr,growthmode,sh,hab",gs%IS(i,n)%init_growth,gs%IS(i,n)%growth_mode &
                       ,gs%IS(i,n)%sh_type,gs%IS(i,n)%habit
                  write(*,*) "T,s_v(2)",ag%TV(n)%T,ag%TV(n)%s_v_n(2)
                  write(*,*) "oq_alen",oq_alen(i,n),ndxdt_alen(i,n)
                  write(*,*) "oq_clen",oq_clen(i,n),ndxdt_clen(i,n)
                  write(*,*) "dcondt:",gs%MS(i,n)%dcondt(:)
                  write(*,*) "dmasst:",gs%MS(i,n)%dmassdt(1,:)
                  write(*,*) "dvoldt_a:",gs%MS(i,n)%dvoldt(iacr,:)
                  write(*,*) "dvoldt_c:",gs%MS(i,n)%dvoldt(iccr,:)
               end if
              gs%MS(i,n)%c_len=c_min
              gs%MS(i,n)%a_len=a_min
              gs%IS(i,n)%d=d_min
              gs%IS(i,n)%ag=ag_min
              gs%IS(i,n)%cg=cg_min
            endif

!            if(phi>1.01) then
!              ierror(i,n)=2
!              write(*,*) "phi>2:",i,KD(n),ID(n),JD(n),gs%MS(i,n)%con,gs%MS(i,n)%mass,gs%MS(i,n)%a_len &
!                   ,gs%MS(i,n)%c_len,gs%IS(i,n)%d
!              write(*,*) "inigr,growthmode,sh,hab",gs%IS(i,n)%init_growth,gs%IS(i,n)%growth_mode &
!                   ,gs%IS(i,n)%sh_type,gs%IS(i,n)%habit
!              write(*,*) "T,s_v(2)",ag%TV(n)%T,ag%TV(n)%s_v_n(2)
!              write(*,*) "old con",oc(i,n)
!              write(*,*) "oq_alen",oq_alen(i,n),ndxdt_alen(i,n)
!              write(*,*) "oq_clen",oq_clen(i,n),ndxdt_clen(i,n)
!              write(*,*) "dcondt:",gs%MS(i,n)%dcondt(:)
!              write(*,*) "dmasst:",gs%MS(i,n)%dmassdt(1,:)
!              write(*,*) "dvoldt_a:",gs%MS(i,n)%dvoldt(iacr,:)
!              write(*,*) "dvoldt_c:",gs%MS(i,n)%dvoldt(iccr,:)
!            endif

            den_ip1=(gs%MS(i,n)%mass(imc)+gs%MS(i,n)%mass(ima)+gs%MS(i,n)%mass(imr))/&
                  gs%MS(i,n)%con/max(1.0e-30_RP,gs%IS(i,n)%V_cs)
!
!!!            if(den_ip1<1.0e-3.and.icond3(i,n)==1) then
!            if(den_ip1<1.0e-2.and.icond3(i,n)==1) then
!               ierror(i,n)=3
!               den_ip1=(gs%MS(i,n)%mass(imc)+gs%MS(i,n)%mass(ima)+gs%MS(i,n)%mass(imr))/&
!                    gs%MS(i,n)%con/gs%IS(i,n)%V_cs
!               write(*,*) "den<1.0e-3:",i,KD(n),ID(n),JD(n),icond3(i,n),gs%MS(i,n)%con,gs%MS(i,n)%mass,gs%MS(i,n)%a_len &
!                    ,gs%MS(i,n)%c_len,gs%IS(i,n)%d,den_ip1,gs%IS(i,n)%v_cs
!               write(*,*) "inigr,growthmode,sh,hab",gs%IS(i,n)%init_growth,gs%IS(i,n)%growth_mode &
!                    ,gs%IS(i,n)%sh_type,gs%IS(i,n)%habit
!               write(*,*) "T,s_v(2)",ag%TV(n)%T,ag%TV(n)%s_v_n(2)
!               write(*,*) "oq_alen",oq_alen(i,n),ndxdt_alen(i,n)
!               write(*,*) "oq_clen",oq_clen(i,n),ndxdt_clen(i,n)
!               write(*,*) "old den",gs%MS(i,n)%den
!               write(*,*) "dcondt:",gs%MS(i,n)%dcondt(:)
!               write(*,*) "dmasst:",gs%MS(i,n)%dmassdt(1,:)
!               write(*,*) "dvoldt_v:",gs%MS(i,n)%dvoldt(ivcs,:)
!               write(*,*) "dvoldt_a:",gs%MS(i,n)%dvoldt(iacr,:)
!               write(*,*) "dvoldt_c:",gs%MS(i,n)%dvoldt(iccr,:)
!!!            endif

          endif
        endif
      enddo
      enddo


!dbg      do in=1,gs%n_bin*L
!dbg        n=(in-1)/gs%N_BIN+1
!dbg        i=in-(n-1)*gs%N_BIN
!dbg        if(gs%MS(i,n)%dmassdt(1,7)>0.0) then
!dbg        write(*,*) "update_group_all ck1:",i,KD(n),ID(n),JD(n),gs%MS(i,n)%con,gs%MS(i,n)%mass,gs%MS(i,n)%a_len &
!dbg                   ,gs%MS(i,n)%c_len,gs%IS(i,n)%d
!dbg        write(*,*) "dcondt_s:",i,gs%MS(i,n)%dcondt(:)
!dbg        write(*,*) "dmasdt_s:",i,gs%MS(i,n)%dmassdt(1,:)
!dbg        endif
!dbg      enddo
!dbg      do in=1,gr%n_bin*L
!dbg        n=(in-1)/gr%N_BIN+1
!dbg        i=in-(n-1)*gr%N_BIN
!dbg        if(gr%MS(i,n)%dmassdt(1,3)>0.0) then
!dbg        write(*,*) "con:",i,gr%MS(i,n)%con
!dbg        write(*,*) "dcondt_r:",i,gr%MS(i,n)%dcondt(:)
!dbg        write(*,*) "dmasst_r:",i,gr%MS(i,n)%dmassdt(1,:)
!dbg        endif
!dbg      enddo
!            write(*,*) "inigr,growthmode,sh,hab",gs%IS(i,n)%init_growth,gs%IS(i,n)%growth_mode &
!               ,gs%IS(i,n)%sh_type,gs%IS(i,n)%habit
!            write(*,*) "T,s_v(2)",ag%TV(n)%T,ag%TV(n)%s_v_n(2)
!!            write(*,*) "old con",oc(i,n)
!            write(*,*) "oq_alen",oq_alen(i,n),ndxdt_alen(i,n)
!!            write(*,*) "oq_clen",oq_clen(i,n),ndxdt_clen(i,n)
!            write(*,*) "dcondt:",gs%MS(i,n)%dcondt(:)
!            write(*,*) "dmasst:",gs%MS(i,n)%dmassdt(1,:)
!       enddo


!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, L
      do i = 1, gs%N_BIN
        if(icond3(i,n)==0.or.icond3(i,n)==2) then
          ! in case of V_cs==0.0
          gs%MS(i,n)%con = 0.0_PS
          gs%MS(i,n)%mass(imt) = 0.0_PS
          gs%MS(i,n)%mass(imc)=0.0_PS
          gs%MS(i,n)%mass(imr)=0.0_PS
          gs%MS(i,n)%mass(ima)=0.0_PS
          gs%MS(i,n)%mass(imw)=0.0_PS
          gs%MS(i,n)%mass(imf)=0.0_PS
          gs%MS(i,n)%mass(imat)=0.0_PS
          gs%MS(i,n)%mass(imas)=0.0_PS
          gs%IS(i,n)%V_cs=0.0_PS
          gs%MS(i,n)%a_len=0.0_PS
          gs%MS(i,n)%c_len=0.0_PS
          gs%IS(i,n)%d=0.0_PS
          gs%IS(i,n)%ag=0.0_PS
          gs%IS(i,n)%cg=0.0_PS
          gs%IS(i,n)%n_exice=0.0_PS
        elseif(icond3(i,n)==3) then
          ! in case of v_cs<v_csmin_hex
          gs%IS(i,n)%V_cs=V_csmin_hex
          gs%MS(i,n)%a_len = 1.0e-4_PS
          gs%MS(i,n)%c_len = 1.0e-4_PS
          gs%IS(i,n)%d = 0.0_PS
          gs%IS(i,n)%ag=0.0_PS
          gs%IS(i,n)%cg=0.0_PS
          gs%IS(i,n)%n_exice=0.0_PS
          gs%MS(i,n)%mass(imr)=0.0_PS
          tmass1=m_icmin*gs%MS(i,n)%con
          gs%MS(i,n)%mass(imc)=tmass1
          gs%MS(i,n)%mass(imt)=tmass1
          tmass2=min(max(m_lmt_ap,min_fapt_s*tmass1),tmass1)
          gs%MS(i,n)%mass(imat)=tmass2
          gs%MS(i,n)%mass(imas)=min(max(m_lmt_ap,min_faps_s*tmass2),tmass2)
          gs%MS(i,n)%mass(imw)=0.0_PS
          gs%MS(i,n)%mass(imf)=0.0_PS

        endif
      enddo
      enddo


      do i=1,gs%n_bin
        do n=1,L
          totmixs(n)=totmixs(n)+max(0.0_PS,gs%MS(i,n)%mass(imt)-gs%MS(i,n)%mass(imat))
        enddo
      end do

    endif

    ! +++ update vapor field +++
    if(update_vapor==1) then

      dtd=max(gr%dt,gs%dt,ga(1)%dt)
      do n=1,L
        ag%tv(n)%rv=ag%tv(n)%rv+ag%tv(n)%dmassdt/ag%tv(n)%den*dtd
        rv2(n)=qtp(n)-(totmixr(n)+totmixs(n))/ag%TV(n)%den
        e_n=ag%TV(n)%P*ag%tv(n)%rv/(Rdvchiarui+ag%tv(n)%rv)
        s_n(n)=e_n/get_sat_vapor_pres_lk(1,ag%TV(n)%T_m,ag%estbar,ag%esitbar)-1.0

!!c        write(*,*) "ck rv",kd(n),id(n),jd(n),ag%TV(n)%rv,ag%TV(n)%dmassdt, &
!!c                    rv2(n),qtp(n),totmixr(n),totmixs(n)

        if(debug .and. s_n(n)>0.10_PS) then
           write(*,*) "update_group:k,i,j,p,qv,qr,qi,t,tn,sv1,svn1,sn,rv2",KD(n),ID(n),JD(n),ag%TV(n)%P,&
               ag%tv(n)%rv,totmixr(n),totmixs(n),ag%TV(n)%T,&
               ag%TV(n)%T_n,ag%TV(n)%s_v(1),ag%TV(n)%s_v_n(1),s_n(n),rv2(n)
        endif
        if(ag%tv(n)%rv<0.0_PS.or.s_n(n)>0.50) then
          if(ag%tv(n)%rv<0.0_PS) then
             ! assume vapor mixing ratio to 10% of ice saturation.
             rv_n=(-0.9_PS+1.0_PS)*ag%tv(n)%e_sat(2)*Rdvchiarui/&
                      (ag%tv(n)%P-(-0.9+1.0_PS)*ag%tv(n)%e_sat(2))

             if ( debug ) then
                write(*,*) "rv is negative"
                write(*,*) "update_group neg", KD(n),ID(n),JD(n),ag%TV(n)%P,&
                     qtp(n),ag%tv(n)%rv,rv_n,totmixr(n),totmixs(n),ag%TV(n)%T,&
                     ag%TV(n)%T_n,ag%TV(n)%s_v(1),ag%TV(n)%s_v_n(1),s_n(n)

                write(*,*) "mes_rc,density of air",mes_rc(n),ag%TV(n)%den
                write(*,*) "rain"
                write(*,*) (gr%MS(i,n)%mass(rmt),i=1,gr%N_BIN)
                do i=1,gr%N_BIN
                   write(*,*) i,(gr%MS(i,n)%dmassdt(rmt,j),j=1,gr%N_tendpros)
                end do

                write(*,*) "ice"
                write(*,*) (gs%MS(i,n)%mass(imt),i=1,gs%N_BIN)
                do i=1,gs%N_BIN
                   write(*,*) i,(gs%MS(i,n)%dmassdt(imt,j),j=1,gs%N_tendpros)
                end do
                write(*,*) ag%TV(n)%s_v(1),ag%TV(n)%s_v(2),&
                     ag%TV(n)%s_v_n(1),ag%TV(n)%s_v_n(2)
             end if

            ag%tv(n)%rv=rv_n
          else
             if ( debug ) then
                write(*,*) "rv(n) is super saturated"
                write(*,*) "update_group su", KD(n),ID(n),JD(n),ag%TV(n)%P,&
                     qtp(n),ag%tv(n)%rv,totmixr(n),totmixs(n),ag%TV(n)%T,&
                     ag%TV(n)%T_n,ag%TV(n)%s_v(1),ag%TV(n)%s_v_n(1),s_n(n)
                write(*,*) "mes_rc,density of air",mes_rc(n),ag%TV(n)%den
                write(*,*) "rain"
                write(*,*) (gr%MS(i,n)%mass(rmt),i=1,gr%N_BIN)
                do i=1,gr%N_BIN
                   write(*,*) i,(gr%MS(i,n)%dmassdt(rmt,j),j=1,gr%N_tendpros)
                end do

                write(*,*) "ice"
                write(*,*) (gs%MS(i,n)%mass(imt),i=1,gs%N_BIN)
                do i=1,gs%N_BIN
                   write(*,*) i,(gs%MS(i,n)%dmassdt(imt,j),j=1,gs%N_tendpros)
                end do
                write(*,*) ag%TV(n)%s_v(1),ag%TV(n)%s_v(2),&
                     ag%TV(n)%s_v_n(1),ag%TV(n)%s_v_n(2)
             end if
          endif
        endif
      enddo

    endif

    !
    ! Aerosol
    !
    if(flagp_a/=0) then

      do ic=1,nacat
        ! +++ for aerosols +++
        if((flagp_a==-1.or.flagp_a==-4).and.ic==2) then
          cycle
        elseif((flagp_a==-2.or.flagp_a==-5).and.ic/=2) then
          cycle
        elseif(flagp_a==-3.or.flagp_a==-6) then
          !write(fid_alog,*) "flaga", MAXVAL(ga(1)%MS(:,:)%con), MAXLOC(ga(1)%MS(:,:)%con), MINVAL(ga(1)%MS(:,:)%con), MINLOC(ga(1)%MS(:,:)%con)
          cycle
        end if

!        do in=1,ga(ic)%n_bin*L
!          n=(in-1)/ga(ic)%N_BIN+1
!          i=in-(n-1)*ga(ic)%N_BIN
        do n = 1, L
        do i = 1, ga(ic)%N_BIN
          ndxdt(i,n) = 0.0d+0
        enddo
        enddo

        do j=1,ga(ic)%n_tendpros
          if(iupdate_ga(j,ic)==0) cycle
!          do in=1,ga(ic)%n_bin*L
!            n=(in-1)/ga(ic)%N_BIN+1
!            i=in-(n-1)*ga(ic)%N_BIN
          do n = 1, L
          do i = 1, ga(ic)%N_BIN
            ndxdt(i,n)=ndxdt(i,n)+ga(ic)%MS(i,n)%dcondt(j)

!!c           if((ga(ic)%MS(i,n)%dcondt(j)<-1.0e-5 .or.ga(ic)%MS(i,n)%dcondt(j)>1.0e-5) ) then
!!c                   if(ic==2.and.ga(ic)%MS(i,n)%dcondt(j)>1.0e-4) then
!!c                      write(*,*) "con_a:ic,process,grid,dcdt",ic,j,n,ga(ic)%MS(i,n)%dcondt(j)
!!c                      write(*,*) "m1,m2,m3",ga(ic)%MS(i,n)%mass
!!c                      write(*,*) "con", ga(ic)%MS(i,n)%con
!!c                      write(*,*) "svw,i",ag%tv(n)%s_v(1),ag%tv(n)%s_v(2)
!!c                      write(*,*) "dsvw,i",ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2)
!!c!!c                      stop
!!c                   end if
          enddo
          enddo
        enddo

!        do in=1,ga(ic)%n_bin*L
!          n=(in-1)/ga(ic)%N_BIN+1
!          i=in-(n-1)*ga(ic)%N_BIN
        do n = 1, L
        do i = 1, ga(ic)%N_BIN
          ga(ic)%MS(i,n)%con=max(0.0_DS,&
                     ga(ic)%MS(i,n)%con+ndxdt(i,n)*ga(ic)%dt)

        enddo
        enddo
        ! check
        icond1(1:ga(ic)%n_bin,1:L)=0
!        do in=1,ga(ic)%n_bin*L
!          n=(in-1)/ga(ic)%N_BIN+1
!          i=in-(n-1)*ga(ic)%N_BIN
        do n = 1, L
        do i = 1, ga(ic)%N_BIN
          if(ga(ic)%MS(i,n)%con>1.0e+9) then
            icond1(i,n)=1
          endif
        enddo
        enddo
        if ( debug ) then
!          do in=1,ga(ic)%n_bin*L
!            n=(in-1)/ga(ic)%N_BIN+1
!            i=in-(n-1)*ga(ic)%N_BIN
           do n = 1, L
           do i = 1, ga(ic)%N_BIN
              if(icond1(i,n)==1) then
                 write(*,*) "update_group:con",KD(n),ID(n),JD(n),ic,ga(ic)%MS(i,n)%con
                 write(*,*) "qr >bin,grid,dcondt", i,ID(n),JD(n),KD(n),(gr%MS(j,n)%dcondt(1),j=1,gr%N_BIN)
                 write(*,*) "qr >bin,grid,dmassdt", i,ID(n),JD(n),KD(n),(gr%MS(j,n)%dmassdt(rmt,1),j=1,gr%N_BIN)
                 write(*,*) "qa >icat,dcondt",ic,(ga(ic)%ms(1,n)%dcondt(j),j=1,ga(ic)%n_tendpros)
                 write(*,*) "qa >icat,dmassdt",ic,(ga(ic)%MS(1,n)%dmassdt(amt,j),j=1,ga(ic)%n_tendpros)
                 write(*,*) "qi >bin,grid,dcondt",i,ID(n),JD(n),KD(n),(gs%MS(j,n)%dcondt(1),j=1,gs%N_BIN)
                 write(*,*) "qi >bin,grid,dmassdt",i,ID(n),JD(n),KD(n),(gs%MS(j,n)%dmassdt(imt,1),j=1,gs%N_BIN)
              end if
!!c             write(*,'("i,n,mten,cten",2i5,10es15.6)') i,n,ndxdt1,ndxdt,&
!!c                  ga(ic)%MS(i,n)%dmassdt(amt,1),ga(ic)%MS(i,n)%dmassdt(amt,9),&
!!c                  ga(ic)%MS(i,n)%dcondt(1),ga(ic)%MS(i,n)%dcondt(9)
           enddo
           enddo
        end if
        if(level>=4) then


          !
          ! mass components
          !
!          do in=1,ga(ic)%n_bin*L
!            n=(in-1)/ga(ic)%N_BIN+1
!            i=in-(n-1)*ga(ic)%N_BIN
           do n = 1, L
           do i = 1, ga(ic)%N_BIN
            ndxdt(i,n) = 0.0d+0
            ndxdt2(i,n) = 0.0d+0
          enddo
          enddo

          do j=1,ga(ic)%n_tendpros
            if(iupdate_ga(j,ic)==0) cycle
!            do in=1,ga(ic)%n_bin*L
!              n=(in-1)/ga(ic)%N_BIN+1
!              i=in-(n-1)*ga(ic)%N_BIN
            do n = 1, L
            do i = 1, ga(ic)%N_BIN
              ! total mass of aerosols (g/g)
              ndxdt(i,n)=ndxdt(i,n)+ga(ic)%MS(i,n)%dmassdt(amt,j)

              ! mass by soluble mass of aerosols (g/g)
              ndxdt2(i,n)=ndxdt2(i,n)+ga(ic)%MS(i,n)%dmassdt(ams,j)
            enddo
            enddo
          enddo

          if ( debug ) then
!          ierror(1:ga(ic)%n_bin*ga(ic)%n_tendpros*L)=0
!          do ijn=1,ga(ic)%n_bin*ga(ic)%n_tendpros*L
!            in=(ijn-1)/ga(ic)%n_tendpros+1
!            j=ijn-(in-1)*ga(ic)%n_tendpros
!            n=(in-1)/ga(ic)%n_bin+1
!            i=in-(n-1)*ga(ic)%n_bin
             do n = 1, L
             do i = 1, ga(ic)%n_bin
             do j = 1, ga(ic)%n_tendpros
                if(iupdate_ga(j,ic)==0) cycle
                if(ga(ic)%MS(i,n)%dmassdt(amt,j)>1.0e-7.or.&
                   ga(ic)%MS(i,n)%dmassdt(amt,j)<-1.0e-7)then
                   write(*,*) "Warning aptm in ap>ga(ic)id,process,grid,id,jd,kd,dmassdt"
                   write(*,*) ic,j,n,ID(n),JD(n),KD(n),ga(ic)%MS(i,n)%dmassdt(amt,j)
                   write(*,*) "m1,m2,m3",ga(ic)%MS(i,n)%mass
                   write(*,*) "con", ga(ic)%MS(i,n)%con
                   write(*,*) "svw,i",ag%tv(n)%s_v(1),ag%tv(n)%s_v(2)
                   write(*,*) "dsvw,i",ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2)
                endif
             enddo
             enddo
             enddo
          end if

          if ( debug ) then
!          ierror(1:ga(ic)%n_bin*ga(ic)%n_tendpros*L)=0
!          do ijn=1,ga(ic)%n_bin*ga(ic)%n_tendpros*L
!            in=(ijn-1)/ga(ic)%n_tendpros+1
!            j=ijn-(in-1)*ga(ic)%n_tendpros
!            n=(in-1)/ga(ic)%n_bin+1
!            i=in-(n-1)*ga(ic)%n_bin
             do n = 1, L
             do i = 1, ga(ic)%n_bin
             do j = 1, ga(ic)%n_tendpros
                if(iupdate_ga(j,ic)==0) cycle
                if(ga(ic)%MS(i,n)%dmassdt(ams,j)>1.0e+08.or.&
                   ga(ic)%MS(i,n)%dmassdt(ams,j)<-1.0e+08)then
                   write(*,*) "Warning apsm in ap>ga(ic)id,process,dmassdt"
                   write(*,*) ic,j,ga(ic)%MS(i,n)%dmassdt(ams,j)
                   write(*,*) "m1,m2,m3",ga(ic)%MS(i,n)%mass
                   write(*,*) "con", ga(ic)%MS(i,n)%con
                   write(*,*) "svw,i",ag%tv(n)%s_v(1),ag%tv(n)%s_v(2)
                   write(*,*) "nmass",ga(ic)%n_masscom
                endif
             enddo
             enddo
             enddo
          end if

!          do in=1,ga(ic)%n_bin*L
!            n=(in-1)/ga(ic)%N_BIN+1
!            i=in-(n-1)*ga(ic)%N_BIN
          do n = 1, L
          do i = 1, ga(ic)%N_BIN
            ga(ic)%MS(i,n)%mass(amt)=max(ga(ic)%MS(i,n)%mass(amt)+&
                           ndxdt(i,n)*ga(ic)%dt,0.0_ds)
          enddo
          enddo

!          do in=1,ga(ic)%n_bin*L
!            n=(in-1)/ga(ic)%N_BIN+1
!            i=in-(n-1)*ga(ic)%N_BIN
          do n = 1, L
          do i = 1, ga(ic)%N_BIN
            ga(ic)%MS(i,n)%mass(ams)=min(max(ga(ic)%MS(i,n)%mass(ams)+&
                           real(ndxdt2(i,n)*ga(ic)%dt,ps),0.0_ps),ga(ic)%MS(i,n)%mass(amt))

          enddo
          enddo

          if(ic==1) then
            !
            ! for CCN accumulation mode (category 1)
            !
!!c                         if(XA(amt_q,i,ic,n)*sep_faps>XA(ams_q,i,ic,n)) then
!!c                            write(*,*) "sol mass is less than limt: ic,n",ic,n
!!c                            write(*,*) XA(amt_q,i,ic,n),XA(acon_q,i,ic,n),XA(ams_q,i,ic,n)
!!c                            write(*,'("dmassdt 1",15ES15.6)') (ga(ic)%MS(i,n)%dmassdt(amt,j),j=1,ga(ic)%N_tendpros)
!!c                            write(*,'("dmassdt 2",15ES15.6)') (ga(ic)%MS(i,n)%dmassdt(ams,j),j=1,ga(ic)%N_tendpros)
!!c                         end if
!            do in=1,ga(ic)%n_bin*L
!              n=(in-1)/ga(ic)%N_BIN+1
!              i=in-(n-1)*ga(ic)%N_BIN
             do n = 1, L
             do i = 1, ga(ic)%N_BIN
              ga(ic)%MS(i,n)%mass(ams)=min(max(ga(ic)%MS(i,n)%mass(amt)*sep_faps,&
                           ga(ic)%MS(i,n)%mass(ams)),ga(ic)%MS(i,n)%mass(amt))
            enddo
            enddo
          elseif(ic==2) then
            !
            ! for IN accumulation mode (category 2)
            !
!            do in=1,ga(ic)%n_bin*L
!              n=(in-1)/ga(ic)%N_BIN+1
!              i=in-(n-1)*ga(ic)%N_BIN
             do n = 1, L
             do i = 1, ga(ic)%N_BIN
              ga(ic)%MS(i,n)%mass(ams)=max(0.0_PS,&
                           min(ga(ic)%MS(i,n)%mass(amt)*0.99_PS*sep_faps,&
                           ga(ic)%MS(i,n)%mass(ams)))

            enddo
            enddo
          elseif(ic==3) then
            !
            ! for CCN nucleation mode (category 3)
            !
!            do in=1,ga(ic)%n_bin*L
!              n=(in-1)/ga(ic)%N_BIN+1
!              i=in-(n-1)*ga(ic)%N_BIN
             do n = 1, L
             do i = 1, ga(ic)%N_BIN
              ga(ic)%MS(i,n)%mass(ams)=min(max(ga(ic)%MS(i,n)%mass(amt)*sep_faps,&
                           ga(ic)%MS(i,n)%mass(ams)),ga(ic)%MS(i,n)%mass(amt))
            enddo
            enddo
          elseif(ic==4) then
            !
            ! for CCN nucleation mode (category 4)
            !
!            do in=1,ga(ic)%n_bin*L
!              n=(in-1)/ga(ic)%N_BIN+1
!              i=in-(n-1)*ga(ic)%N_BIN
             do n = 1, L
             do i = 1, ga(ic)%N_BIN
              ga(ic)%MS(i,n)%mass(ams)=min(max(ga(ic)%MS(i,n)%mass(amt)*sep_faps,&
                           ga(ic)%MS(i,n)%mass(ams)),ga(ic)%MS(i,n)%mass(amt))
            enddo
            enddo
          else
             LOG_ERROR("update_group_all",*) "cal_model_tend_all>this ap cat not defined"
             call PRC_abort
          end if


          m_lmt=coef4pi3*r3_lmt
!test          if(ic==1.or.ic==2.or.ic==4) then
!test            m_lmt=coef4pi3*r3_lmt
!test          elseif(ic==3) then
!test            m_lmt=coef4pi3*r3_lmt_nuc
!test          else
!test            write(*,*) "this category not defined in ini_group_all"
!test            stop
!test          endif

!          do in=1,ga(ic)%n_bin*L
!            n=(in-1)/ga(ic)%N_BIN+1
!            i=in-(n-1)*ga(ic)%N_BIN
          do n = 1, L
          do i = 1, ga(ic)%N_BIN
            if(ga(ic)%MS(i,n)%con<n_lmt_ap) then

!!c                      write(*,'("ap fixed ic",i5,3es15.6)')&
!!c                           ic,XA(acon_q,i,ic,n)*ag%tv(n)%den,n_lmt_ap
!!c                      rmod1=XA(amt_q,i,ic,n)/XA(acon_q,i,ic,n)
!!c                      rmod2=XA(ams_q,i,ic,n)/XA(acon_q,i,ic,n)
              var_n1=n_lmt_ap
              var_n2=var_n1*m_lmt*ga(ic)%MS(i,n)%den
              ga(ic)%MS(i,n)%con=var_n1
              ga(ic)%MS(i,n)%mass(amt)=var_n2
              ga(ic)%MS(i,n)%mass(ams)=var_n2*eps_ap0(ic)
            end if
          enddo
          enddo

!          do in=1,ga(ic)%n_bin*L
!            n=(in-1)/ga(ic)%N_BIN+1
!            i=in-(n-1)*ga(ic)%N_BIN
          do n = 1, L
          do i = 1, ga(ic)%N_BIN
            if(ga(ic)%MS(i,n)%con>n_max_ap) then
!!c                      write(*,'("ap con reached max. so fixed.",i5,3es15.6)')&
!!c                           ic,XA(acon_q,i,ic,n)*ag%tv(n)%den,n_max_ap
              ga(ic)%MS(i,n)%con=n_max_ap
            end if
          enddo
          enddo

!          do in=1,ga(ic)%n_bin*L
!            n=(in-1)/ga(ic)%N_BIN+1
!            i=in-(n-1)*ga(ic)%N_BIN
          do n = 1, L
          do i = 1, ga(ic)%N_BIN
            if(ga(ic)%MS(i,n)%mass(amt)<m_lmt*ga(ic)%MS(i,n)%den*&
               ga(ic)%MS(i,n)%con) then
              var_n1=ga(ic)%MS(i,n)%con*m_lmt*ga(ic)%MS(i,n)%den
              ga(ic)%MS(i,n)%mass(amt)=var_n1
              ga(ic)%MS(i,n)%mass(ams)=var_n1*eps_ap0(ic)
            end if
          enddo
          enddo

!          do in=1,ga(ic)%n_bin*L
!            n=(in-1)/ga(ic)%N_BIN+1
!            i=in-(n-1)*ga(ic)%N_BIN
          do n = 1, L
          do i = 1, ga(ic)%N_BIN
            ga(ic)%MS(i,n)%mass(ams)=&
                min(ga(ic)%MS(i,n)%mass(ams),ga(ic)%MS(i,n)%mass(amt))
          enddo
          enddo
        endif
      enddo
    endif

!!$  contains
!!$    subroutine check_phi(from)
!!$      character*(*) :: from
!!$!      do in=1,gs%n_bin*L
!!$!        n=(in-1)/gs%N_BIN+1
!!$!        i=in-(n-1)*gs%N_BIN
!!$      do n = 1, L
!!$      do i = 1, gs%N_BIN
!!$!        ierror(in)=0
!!$        if(gs%MS(i,n)%a_len>1.0e-30_RP.and.gs%MS(i,n)%c_len/max(1.0e-30_RP,gs%MS(i,n)%a_len)>2.0_RP) then
!!$            write(*,*) "ck_phi detected",trim(from) &
!!$            ,i,KD(n),ID(n),JD(n),gs%MS(i,n)%con,gs%MS(i,n)%mass,gs%MS(i,n)%a_len &
!!$                     ,gs%MS(i,n)%c_len,gs%IS(i,n)%d
!!$            write(*,*) "inigr,growthmode,sh,hab",gs%IS(i,n)%init_growth,gs%IS(i,n)%growth_mode &
!!$               ,gs%IS(i,n)%sh_type,gs%IS(i,n)%habit
!!$        endif
!!$      enddo
!!$      enddo
!!$    end subroutine check_phi
  end subroutine update_group_all

  subroutine ini_tendency(g,iupdate,ilag)
    type (Group), intent(inout)  :: g
    integer,dimension(*),intent(in) :: iupdate
    integer,intent(in) :: ilag
    integer :: i,j,k,m

!    do k=1,g%N_tendpros
!      if(iupdate(k)==1) then
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
    do j = 1, g%L
    do i = 1, g%N_BIN
    do k=1,g%N_tendpros
      if(iupdate(k)==1) then
          g%MS(i,j)%dcondt(k) = 0.0_PS
       end if
    enddo
    enddo
    enddo
!    do m=1,1+g%N_masscom
!      do k=1,g%N_tendpros
!        if(iupdate(k)==1) then
!          do ij=1,g%N_BIN*g%L
!            j=(ij-1)/g%N_BIN+1
!            i=ij-(j-1)*g%N_BIN
    do j = 1, g%L
    do i = 1, g%N_BIN
    do k=1, g%N_tendpros
       if(iupdate(k)==1) then
          do m=1,1+g%N_masscom
             g%MS(i,j)%dmassdt(m,k) = 0.0_PS
          enddo
       endif
    enddo
    enddo
    enddo

!    do m=1,g%N_nonmass
!      do k=1,g%N_tendpros
!        if(iupdate(k)==1) then
!          do ij=1,g%N_BIN*g%L
!            j=(ij-1)/g%N_BIN+1
!            i=ij-(j-1)*g%N_BIN
    do j = 1, g%L
    do i = 1, g%N_BIN
    do k = 1, g%N_tendpros
       if(iupdate(k)==1) then
          do m=1,g%N_nonmass
             g%MS(i,j)%dvoldt(m,k) = 0.0_PS
          enddo
       endif
    enddo
    enddo
    enddo

    ! Lagrangian tendency
    if(ilag==1) then
      if(iupdate(1)==1) then
        ! vapor deposition
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
         do j = 1, g%L
         do i = 1, g%N_BIN
            g%MS(i,j)%Ldmassdt(1) = 0.0_PS
         enddo
         enddo
      endif
      if(iupdate(4)==1) then
        ! riming process
!        do ij=1,g%N_BIN*g%L
!          j=(ij-1)/g%N_BIN+1
!          i=ij-(j-1)*g%N_BIN
         do j = 1, g%L
         do i = 1, g%N_BIN
            g%MS(i,j)%Ldmassdt(2) = 0.0_PS
         enddo
         enddo
      endif
    endif
  end subroutine ini_tendency

  subroutine ini_tendency_ag(ag,ivp)
    ! thermo variable object
    type (AirGroup), intent(inout)  :: ag
    integer,intent(in) :: ivp
    integer :: n

    do n=1,ag%L
      ag%TV(n)%dmassdt_v(ivp)=0.0_PS
    enddo
  end subroutine ini_tendency_ag

  subroutine update_modelvars_all(xc,xr,xs,xa,rv, &
       l,nrtype,nrbin,nrcat,nstype,nsbin,nscat,natype,nabin,nacat,&
       level,ag,gr,gs,ga,&
       flagp_r,flagp_s,flagp_a)
    integer,intent(in)  :: l,nrtype,nrbin,nrcat,nstype,nsbin,nscat,&
         natype,nabin,nacat
    !integer :: ID(*),JD(*),KD(*)
    type (group), intent(in)  :: gr,gs
    type (group), dimension(nacat),intent(in) :: ga
    type (airgroup), intent(in) :: ag
    real(MP_KIND) :: xc(*),XR(nrtype,nrbin,nrcat,*),XS(nstype,nsbin,nscat,*),&
            XA(natype,nabin,nacat,*),rv(*)
    integer,intent(in)  :: flagp_r,flagp_s,flagp_a
    ! level of complexity
    integer, intent(in)  :: level
    !
    ! local vars
    !
    ! new total volume variables
    real(ps),dimension(mxnbin,LMAX)  :: aQ

    ! total mixing ratio of hydrometeors
!tmp    real(PS),dimension(L) :: totmixr,totmixs

!tmp    integer,dimension(mxnbin*LMAX*mxntend) :: ierror
!tmp    integer,dimension(mxnbin*LMAX) :: icond1,icond2

    integer  :: i,k,n,ic!,var_status,j
    integer :: iq1

    ! +++ vapor +++
    do n=1,L
      rv(n)=ag%tv(n)%rv
    enddo
    !
    ! +++ for rain +++
    !
    if( flagp_r > 0 ) then

!      do in=1,gr%n_bin*L
!        n=(in-1)/gr%N_BIN+1
!        i=in-(n-1)*gr%N_BIN
       do n = 1, L
       do i = 1, gr%N_BIN
          XR(rmt_q,i,nrcat,n)=gr%MS(i,n)%mass(rmt)/ag%tv(n)%den
          XR(rcon_q,i,nrcat,n)=gr%MS(i,n)%con/ag%tv(n)%den
       enddo
       enddo
       if(level>=4) then
!        do in=1,gr%n_bin*L
!          n=(in-1)/gr%N_BIN+1
!          i=in-(n-1)*gr%N_BIN
          do n = 1, L
          do i = 1, gr%N_BIN
             XR(rmat_q,i,nrcat,n)=gr%MS(i,n)%mass(rmat)/ag%tv(n)%den
             XR(rmas_q,i,nrcat,n)=gr%MS(i,n)%mass(rmas)/ag%tv(n)%den
          enddo
          enddo
       endif

       ! for cloud droplets
       do n=1,L
          XC(n)=XR(rmt_q,1,nrcat,n)
       enddo
    endif
    !
    ! ice
    !
    if( flagp_s > 0 ) then

!      do in=1,gs%n_bin*L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
       do n = 1, L
       do i = 1, gs%N_BIN
        XS(imt_q,i,nscat,n)=gs%MS(i,n)%mass(imt)/ag%tv(n)%den
        XS(icon_q,i,nscat,n)=gs%MS(i,n)%con/ag%tv(n)%den

        ! mass by ice crystals (g/g)
        XS(imc_q,i,nscat,n)=gs%MS(i,n)%mass(imc)/ag%tv(n)%den
        ! mass by riming process (g/g)
        XS(imr_q,i,nscat,n)=gs%MS(i,n)%mass(imr)/ag%tv(n)%den
        ! mass by agg process (g/g)
        XS(ima_q,i,nscat,n)=gs%MS(i,n)%mass(ima)/ag%tv(n)%den
      enddo
      enddo

      iq1=ivcs_q
      nonmass_loop: do k=1,gs%N_nonmass
        !      1. volume of circumscribing sphere (cm^3/g)
        !      2. con. weighted a-axis length^3 (cm^3/g)
        !      3. con. weighted c-axis length^3 (cm^3/g)
        !      4. con. weighted d-axis length^3 (cm^3/g)
        !      5. con. weighted r-axis length^3 (cm^3/g) (rosetta bulletes)
        !      6. con. weighted e-axis length^3 (cm^3/g) (polycrystals)

        if(k==ivcs) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            aq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%v_cs
          enddo
          enddo
        elseif(k==iacr) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            aq(i,n) = gs%MS(i,n)%con*gs%MS(i,n)%a_len**3.0
          enddo
          enddo
        elseif(k==iccr) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            aq(i,n) = gs%MS(i,n)%con*gs%MS(i,n)%c_len**3.0
          enddo
          enddo
        elseif(k==idcr) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            aq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%d**3.0
          enddo
          enddo
        elseif(k==iag) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            aq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%ag**3.0
          enddo
          enddo
        elseif(k==icg) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            aq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%cg**3.0
          enddo
          enddo
        elseif(k==inex) then
!          do in=1,gs%n_bin*L
!            n=(in-1)/gs%N_BIN+1
!            i=in-(n-1)*gs%N_BIN
           do n = 1, L
           do i = 1, gs%N_BIN
            aq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%n_exice
          enddo
          enddo
        endif

!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
        do n = 1, L
        do i = 1, gs%N_BIN
          XS(iq1-1+k,i,nscat,n)=aq(i,n)/ag%tv(n)%den
        enddo
        enddo
      enddo nonmass_loop

      if(level>=6) then
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
         do n = 1, L
         do i = 1, gs%N_BIN
          ! melt water mass  (g/g)
          XS(imw_q,i,nscat,n)=gs%MS(i,n)%mass(imw)/ag%tv(n)%den
          ! freezing nucleation mass  (g/g)
          XS(imf_q,i,nscat,n)=gs%MS(i,n)%mass(imf)/ag%tv(n)%den
        enddo
        enddo
      endif

      if(level>=4) then
!        do in=1,gs%n_bin*L
!          n=(in-1)/gs%N_BIN+1
!          i=in-(n-1)*gs%N_BIN
         do n = 1, L
         do i = 1, gs%N_BIN
          ! mass by total mass of aerosols (g/g)
          XS(imat_q,i,nscat,n)=gs%MS(i,n)%mass(imat)/ag%tv(n)%den
          ! mass by soluble mass of aerosols (g/g)
          XS(imas_q,i,nscat,n)=gs%MS(i,n)%mass(imas)/ag%tv(n)%den
        enddo
        enddo
      endif

!      totmixs=0.0_PS
!      do i=1,gs%n_bin
!        do n=1,L
!          if(icond1(in)/=1.and.icond2(in)/=1) then
!            totmixs(n)=totmixs(n)+max(0.0_PS,XS(imt_q,i,nscat,n)-XS(imat_q,i,nscat,n))
!          endif
!        enddo
!      end do
!
    endif

    !
    ! Aerosol
    !
    if(flagp_a/=0) then

      do ic=1,nacat
        ! +++ for aerosols +++
        if((flagp_a==-1.or.flagp_a==-4).and.ic==2) then
          cycle
        elseif((flagp_a==-2.or.flagp_a==-5).and.ic/=2) then
          cycle
        elseif(flagp_a==-3.or.flagp_a==-6) then
          cycle
        end if

!        do in=1,ga(ic)%n_bin*L
!          n=(in-1)/ga(ic)%N_BIN+1
!          i=in-(n-1)*ga(ic)%N_BIN
        do n = 1, L
        do i = 1, ga(ic)%N_BIN
          XA(acon_q,i,ic,n)=ga(ic)%MS(i,n)%con/ag%tv(n)%den
          ! total mass of aerosols (g/g)
          XA(amt_q,i,ic,n)=ga(ic)%MS(i,n)%mass(amt)/ag%tv(n)%den
          ! mass by soluble mass of aerosols (g/g)
          XA(ams_q,i,ic,n)=ga(ic)%MS(i,n)%mass(ams)/ag%tv(n)%den
        enddo
        enddo
      enddo
    endif

  end subroutine update_modelvars_all

!!$  subroutine cal_surface_temp( phase, ms, th_var,estbar,esitbar)
!!$    use class_Thermo_Var, only: &
!!$       get_mod_thermal_cond, &
!!$       get_sat_vapor_pres_lk
!!$    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! this routine calculates the surface temperature of hydrometeor,
!!$    ! assuming that the heat tendency is zero (thermal equilibrium).
!!$    !
!!$    ! NOTE: mass tendency by riming process is required, so this should be put
!!$    !       after collection routines.
!!$    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! phase of hydrometeor
!!$    ! 1: liquid, 2: solid
!!$    integer, intent(in)  :: phase
!!$    type (Mass_Bin), intent(inout) :: ms
!!$    type (Thermo_Var), intent(in) :: th_var
!!$    integer :: em
!!$
!!$!!c    real(PS) :: r_e
!!$    ! maximum realisitc super saturation
!!$    real(PS),parameter :: max_Sw=0.001
!!$
!!$    integer             :: i,phase2
!!$    real(PS) :: TS_A1,TS_B1,TS_B11,TS_B12,tmp,tmp_n,es_sfc
!!$
!!$    ! saturation vapor lookup table
!!$    real(DS)  :: estbar(150),esitbar(111)
!!$
!!$
!!$    ! modified conductivity
!!$    real(PS) :: mk_a
!!$
!!$! <<< 2014/12 T. Hashino bug found
!!$    phase2=phase
!!$! >>>> 2014/12 T. Hashino bug found
!!$    if( ms%con > 1.0e-30_PS .and. ms%mass(1) > 1.0e-30_PS ) then
!!$       if( phase == 1 ) then
!!$
!!$!!c       ms%tmp = th_var%T + ( L_e * ms%Ldmassdt(1) )/&
!!$!!c            (4.0_PS*PI*ms%CAP*th_var%k_a*ms%fh)
!!$
!!$          mk_a = get_mod_thermal_cond(ms%a_len,th_var%K_a,th_var%T,th_var%den)
!!$
!!$          ms%tmp=th_var%T +  L_e*(ms%coef(1)*th_var%s_v_n(phase)+ms%coef(2))/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh)
!!$!!c          write(*,'("surfT",10ES15.6)') th_var%s_v_n(phase),th_var%T,ms%tmp,ms%coef(1),ms%coef(2)
!!$
!!$       else if( phase == 2 ) then
!!$
!!$          call cal_coef_Ts(ms,th_var,phase,TS_A1,TS_B11,TS_B12,phase2)
!!$
!!$          if(T_0>th_var%T) then
!!$             if(phase2==1) then
!!$!!c                TS_B1=TS_B11*min(max_Sw,th_var%s_v_n(phase2))+TS_B12
!!$                TS_B1=TS_B11*th_var%s_v_n(phase2)+TS_B12
!!$             else
!!$                TS_B1=TS_B11*th_var%s_v_n(phase2)+TS_B12
!!$             end if
!!$
!!$             em=0
!!$             call cal_tsfc_iter(th_var%T,TS_A1,TS_B1,phase2,estbar,esitbar,tmp,em)
!!$             if(em/=0) then
!!$                write(*,*) "sfc temp error"
!!$                write(*,'(2I5,10Es15.6)') em,phase2,th_var%T,TS_A1,TS_B1,tmp
!!$             endif
!!$             ms%tmp=min(tmp,T_0)
!!$          else
!!$             !
!!$             ! assume that the riming process produce surface water
!!$             ! due to the liquid layer over the hydrometeors:
!!$             ! no fusion heat release by riming process
!!$             !
!!$             ! The liquid layer is turbulent, so the heat received at the
!!$             ! surface is balanced with the one of ice core.
!!$             !
!!$             ! Therefore, the surface temperature is T_0.
!!$! <<< 2015/03 T. Hashino mod for KiD
!!$             ! above assumption seems not working. So set it to the ambient
!!$             ! temperature
!!$             ms%tmp=T_0
!!$!mod             ms%tmp=th_var%T
!!$! >>> 2015/03 T. Hashino
!!$!
!!$          end if
!!$
!!$!!c          ms%e_sat=get_sat_vapor_pres(phase2,ms%tmp)
!!$          ms%e_sat=get_sat_vapor_pres_lk(phase2,ms%tmp,estbar,esitbar)
!!$       end if
!!$    else
!!$       ms%tmp = th_var%T
!!$       ms%e_sat=get_sat_vapor_pres_lk(phase2,ms%tmp,estbar,esitbar)
!!$    end if
!!$  end subroutine cal_surface_temp

!!$  subroutine cal_surface_temp_vec( phase, g,ag,icond1)
!!$    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! this routine calculates the surface temperature of hydrometeor,
!!$    ! assuming that the heat tendency is zero (thermal equilibrium).
!!$    !
!!$    ! NOTE: mass tendency by riming process is required, so this should be put
!!$    !       after collection routines.
!!$    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! phase of hydrometeor
!!$    ! 1: liquid, 2: solid
!!$    integer, intent(in)  :: phase
!!$    type (Group), intent(inout)  :: g
!!$    type (AirGroup), intent(in)   :: ag
!!$    integer,dimension(*),intent(in)    ::  icond1
!!$    integer :: em
!!$
!!$!!c    real(PS) :: r_e
!!$    ! maximum realisitc super saturation
!!$    real(PS),parameter :: max_Sw=0.001
!!$
!!$    integer             :: phase2
!!$    real(PS) :: TS_A1,TS_B1,TS_B11,TS_B12,tmp!,tmp_n,es_sfc
!!$
!!$    ! modified conductivity
!!$    real(PS) :: mk_a
!!$
!!$    integer                     :: i,n,in
!!$
!!$! <<< 2014/12 T. Hashino bug found
!!$    phase2=phase
!!$! >>>> 2014/12 T. Hashino bug found
!!$
!!$    if( phase == 1 ) then
!!$!CDIR NODEP
!!$      do in=1,g%n_bin*g%L
!!$        n=(in-1)/g%N_BIN+1
!!$        i=in-(n-1)*g%N_BIN
!!$        if(icond1(in)==0) then
!!$
!!$!!c       g%MS(i,n)%tmp = ag%TV(n)%T + ( L_e * g%MS(i,n)%Ldmassdt(1) )/&
!!$!!c            (4.0_PS*PI*g%MS(i,n)%CAP*ag%TV(n)%k_a*g%MS(i,n)%fh)
!!$
!!$          mk_a = get_mod_thermal_cond(g%MS(i,n)%a_len,ag%TV(n)%K_a,ag%TV(n)%T,ag%TV(n)%den)
!!$
!!$          g%MS(i,n)%tmp=ag%TV(n)%T +  L_e*(g%MS(i,n)%coef(1)*ag%TV(n)%s_v_n(phase)+g%MS(i,n)%coef(2))/&
!!$               (4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh)
!!$!!c          write(*,'("surfT",10ES15.6)') ag%TV(n)%s_v_n(phase),ag%TV(n)%T,g%MS(i,n)%tmp,g%MS(i,n)%coef(1),g%MS(i,n)%coef(2)
!!$
!!$        endif
!!$      enddo
!!$    else if( phase == 2 ) then
!!$
!!$!CDIR NODEP
!!$      do in=1,g%n_bin*g%L
!!$        n=(in-1)/g%N_BIN+1
!!$        i=in-(n-1)*g%N_BIN
!!$        if(icond1(in)==0) then
!!$
!!$!org          call cal_coef_Ts(g%MS(i,n),ag%TV(n),phase,TS_A1,TS_B11,TS_B12,phase2)
!!$          ! calculate modified thermal conductivity
!!$          mk_a = get_mod_thermal_cond(sqrt(g%MS(i,n)%semi_a**2+g%MS(i,n)%semi_c**2),&
!!$                                 ag%TV(n)%K_a,ag%TV(n)%T,ag%TV(n)%den)
!!$
!!$          if(T_0>ag%TV(n)%T) then
!!$!!c       r_e=ag%TV(n)%e_sat_n(1)/ag%TV(n)%e_sat_n(2)
!!$            if(g%MS(i,n)%inmlt==0) then
!!$              phase2=phase
!!$
!!$              TS_A1=L_s*4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn/&
!!$                   (4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh+C_w*g%MS(i,n)%Ldmassdt(2))
!!$
!!$              TS_B11=L_s*4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn*&
!!$                   ag%TV(n)%e_sat_n(phase2)/ag%TV(n)%T_n/&
!!$                   (4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh+C_w*g%MS(i,n)%Ldmassdt(2))
!!$
!!$              TS_B12=(L_s*4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn*&
!!$               ag%TV(n)%e_sat_n(phase2)/ag%TV(n)%T_n+&
!!$               (L_f+C_w*ag%TV(n)%T_n)*g%MS(i,n)%Ldmassdt(2)-L_f*g%MS(i,n)%dmassdt(1,10)/g%MS(i,n)%con+&
!!$               4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh*ag%TV(n)%T_n)/&
!!$               (4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh+C_w*g%MS(i,n)%Ldmassdt(2))
!!$
!!$
!!$            else
!!$          !
!!$          ! This formulation ignores conduction between the surface of the hydrometeor
!!$          ! and the surface of the ice core.
!!$          !
!!$              phase2=1
!!$
!!$              TS_A1=L_e*4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn/&
!!$                   (4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh+C_w*g%MS(i,n)%Ldmassdt(2))
!!$
!!$              TS_B11=L_e*4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn*&
!!$               ag%TV(n)%e_sat_n(phase2)/ag%TV(n)%T_n/&
!!$               (4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh+C_w*g%MS(i,n)%Ldmassdt(2))
!!$
!!$              TS_B12=(L_e*4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn*&
!!$               ag%TV(n)%e_sat_n(phase2)/ag%TV(n)%T_n+&
!!$               (L_f+C_w*ag%TV(n)%T_n)*g%MS(i,n)%Ldmassdt(2)-L_f*g%MS(i,n)%dmassdt(1,10)/g%MS(i,n)%con+&
!!$               4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh*ag%TV(n)%T_n)/&
!!$               (4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh+C_w*g%MS(i,n)%Ldmassdt(2))
!!$
!!$            end if
!!$
!!$          else
!!$       !
!!$       ! assume that the riming process produce surface water
!!$       ! due to the liquid layer over the hydrometeors:
!!$       ! no fusion heat release by riming process
!!$       !
!!$       ! The liquid layer is turbulent, so the heat received at the
!!$       ! surface is balanced with the one of ice core.
!!$       !
!!$       ! Therefore, the surface temperature is T_0.
!!$            phase2=1
!!$!!c
!!$!!c             A1=L_e*4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn/&
!!$!!c                  (4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh+C_w*g%MS(i,n)%Ldmassdt(2))
!!$!!c
!!$!!c             B1=(L_e*4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn*&
!!$!!c                  ag%TV(n)%e_sat_n(phase2)/ag%TV(n)%T_n*&
!!$!!c                  (ag%TV(n)%s_v_n(phase2)+1.0_PS)+&
!!$!!c                  (C_w*ag%TV(n)%T_n)*g%MS(i,n)%Ldmassdt(2)-L_f*g%MS(i,n)%dmassdt(1,10)/g%MS(i,n)%con+&
!!$!!c                  4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh*ag%TV(n)%T_n)/&
!!$!!c                  (4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh+C_w*g%MS(i,n)%Ldmassdt(2))
!!$          end if
!!$
!!$          if(T_0>ag%TV(n)%T) then
!!$            if(phase2==1) then
!!$!!c                TS_B1=TS_B11*min(max_Sw,ag%TV(n)%s_v_n(phase2))+TS_B12
!!$              TS_B1=TS_B11*ag%TV(n)%s_v_n(phase2)+TS_B12
!!$            else
!!$              TS_B1=TS_B11*ag%TV(n)%s_v_n(phase2)+TS_B12
!!$            end if
!!$
!!$            em=0
!!$            call cal_tsfc_iter(ag%TV(n)%T,TS_A1,TS_B1,phase2,ag%estbar,ag%esitbar,tmp,em)
!!$            if(em/=0) then
!!$              write(*,*) "sfc temp error"
!!$              write(*,'(2I5,10Es15.6)') em,phase2,ag%TV(n)%T,TS_A1,TS_B1,tmp
!!$            endif
!!$            g%MS(i,n)%tmp=min(tmp,T_0)
!!$          else
!!$            !
!!$            ! assume that the riming process produce surface water
!!$            ! due to the liquid layer over the hydrometeors:
!!$            ! no fusion heat release by riming process
!!$            !
!!$            ! The liquid layer is turbulent, so the heat received at the
!!$            ! surface is balanced with the one of ice core.
!!$            !
!!$            ! Therefore, the surface temperature is T_0.
!!$! <<< 2015/03 T. Hashino mod for KiD
!!$            ! above assumption seems not working. So set it to the ambient
!!$            ! temperature
!!$            g%MS(i,n)%tmp=T_0
!!$!mod             g%MS(i,n)%tmp=ag%TV(n)%T
!!$! >>> 2015/03 T. Hashino
!!$!
!!$          end if
!!$
!!$!!c          g%MS(i,n)%e_sat=get_sat_vapor_pres(phase2,g%MS(i,n)%tmp)
!!$          g%MS(i,n)%e_sat=get_sat_vapor_pres_lk(phase2,g%MS(i,n)%tmp,ag%estbar,ag%esitbar)
!!$        end if
!!$      enddo
!!$    end if
!!$  end subroutine cal_surface_temp_vec

  subroutine cal_surface_temp2_vec( phase, g,ag,icond1)
    use class_Mass_Bin, only: &
       cal_coef_Ts3
    use class_Thermo_Var, only: &
       get_mod_thermal_cond, &
       get_sat_vapor_pres_lk
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! this routine calculates the surface temperature of hydrometeor,
    ! assuming that the heat tendency is zero (thermal equilibrium).
    !
    ! NOTE: mass tendency by riming process is required, so this should be put
    !       after collection routines.
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! phase of hydrometeor
    ! 1: liquid, 2: solid
    integer, intent(in)  :: phase
    type (Group), intent(inout)  :: g
    type (AirGroup), intent(in)   :: ag
    integer,dimension(g%N_BIN,g%L),intent(in)    ::  icond1
    !integer :: em

!!c    real(PS) :: r_e
    ! maximum realisitc super saturation
    real(PS),parameter :: max_Sw=0.001

    integer             :: phase2

    ! modified conductivity
    real(PS) :: mk_a

    real(PS) :: x0,x1,x2,gx0,gx1,gx2
    real(PS),dimension(mxnbin,LMAX) :: tsfc,tsfc_b
    real(8) :: TS_A1,TS_B1,TS_D1,TS_D2,TS_B11,TS_B12,TS_B13,w0,w1,w2 &
              ,aL,bL,dL
    real(PS) :: Tmax
    real(PS),dimension(mxnbin,LMAX) :: diff
    integer,dimension(mxnbin,LMAX) :: icond2
    real(PS),dimension(5) :: dT_w

    integer                     :: i,n,iter

! <<< 2014/12 T. Hashino bug found
    phase2=phase
! >>>> 2014/12 T. Hashino bug found

    if( phase == 1 ) then
!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then

!!c       g%MS(i,n)%tmp = ag%TV(n)%T + ( L_e * g%MS(i,n)%Ldmassdt(1) )/&
!!c            (4.0_PS*PI*g%MS(i,n)%CAP*ag%TV(n)%k_a*g%MS(i,n)%fh)
          mk_a = get_mod_thermal_cond(g%MS(i,n)%a_len,ag%TV(n)%K_a,ag%TV(n)%T,ag%TV(n)%den)

          g%MS(i,n)%tmp=ag%TV(n)%T +  L_e*(g%MS(i,n)%coef(1)*ag%TV(n)%s_v_n(phase)+g%MS(i,n)%coef(2))/&
               (4.0_PS*PI*g%MS(i,n)%CAP*mk_a*g%MS(i,n)%fh)
!!c          write(*,'("surfT",10ES15.6)') ag%TV(n)%s_v_n(phase),ag%TV(n)%T,g%MS(i,n)%tmp,g%MS(i,n)%coef(1),g%MS(i,n)%coef(2)

        endif
      enddo
      enddo
    else if( phase == 2 ) then

!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN

        diff(i,n)=0.0
        tsfc(i,n)=ag%TV(n)%T_n
        if(T_0<=ag%TV(n)%T) then
          icond2(i,n)=1
       !
       ! assume that the riming process produce surface water
       ! due to the liquid layer over the hydrometeors:
       ! no fusion heat release by riming process
       !
       ! The liquid layer is turbulent, so the heat received at the
       ! surface is balanced with the one of ice core.
       !
       ! Therefore, the surface temperature is T_0.
          phase2=1
          g%MS(i,n)%tmp=T_0
          g%MS(i,n)%e_sat=get_sat_vapor_pres_lk(phase2,g%MS(i,n)%tmp,ag%estbar,ag%esitbar)

        else
          icond2(i,n)=icond1(i,n)
        endif
      enddo
      enddo

      dT_w(1:5)=(/20.0,10.0,5.0,1.0,0.5/)
      do iter=1,5

!        do in=1,g%n_bin*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
         do n = 1, g%L
         do i = 1, g%N_BIN

          if(icond2(i,n)==0) then

            mk_a = get_mod_thermal_cond(sqrt(g%MS(i,n)%semi_a**2+g%MS(i,n)%semi_c**2),&
                                        ag%TV(n)%K_a,ag%TV(n)%T,ag%TV(n)%den)
            call cal_coef_Ts3(g%MS(i,n),ag%TV(n),g%token &
                             ,TS_A1,TS_B11,TS_B12,TS_B13,TS_D1 &
                             ,phase2,Tmax,mk_a)


            TS_B1=TS_B11*ag%TV(n)%e_sat_n(phase2)/ag%TV(n)%T_n * &
                   (ag%TV(n)%s_v_n(phase2) + 1.0_PS)  + &
            !
                  TS_B12+TS_B13*ag%TV(n)%T_n

            ! locally fit a parabora using Lagrange polynomial
            ! the fit equation is
            !   g(x)=esi(Tsfc)/Tsfc
            !
            ! fit the polynomical
            x2=min(Tmax,tsfc(i,n)+dT_w(iter))
            x1=x2-dT_w(iter)
            x0=x1-dT_w(iter)
            gx0=get_sat_vapor_pres_lk(phase2,x0,ag%estbar,ag%esitbar)/x0
            gx1=get_sat_vapor_pres_lk(phase2,x1,ag%estbar,ag%esitbar)/x1
            gx2=get_sat_vapor_pres_lk(phase2,x2,ag%estbar,ag%esitbar)/x2

            w0=gx0/((x0-x1)*(x0-x2))
            w1=gx1/((x1-x0)*(x1-x2))
            w2=gx2/((x2-x0)*(x2-x1))
            aL=w0+w1+w2
            bL=-w0*(x1+x2)-w1*(x0+x2)-w2*(x0+x1)
            dL=w0*x1*x2+w1*x0*x2+w2*x0*x1

            tsfc_b(i,n)=tsfc(i,n)

            TS_D2=(TS_D1*bL+TS_A1)*(TS_D1*bL+TS_A1)-&
                          4.0*TS_D1*aL*(TS_D1*dL-TS_B1)


            if(TS_D2>=0.0d+0) then
              tsfc(i,n)=(-(TS_D1*bL+TS_A1)+sqrt(TS_D2)  )/ &
                         (2.0*TS_D1*aL)
            else
              tsfc(i,n)=tsfc_b(i,n)
            endif

            diff(i,n)=(tsfc(i,n)-tsfc_b(i,n))/tsfc_b(i,n)

!tmp            write(*,'("tsfc1",3I5,6ES15.6)') i,n,iter,x0,x1,x2,ag%TV(n)%T_n,tsfc(i,n),diff(i,n)
!tmp            write(*,'("tsfc1-1",9ES15.6)') gx0,gx1,gx2,w0,w1,w2,aL,bL,dL
!tmp            write(*,'("tsfc1-2",4ES15.6)') TS_A1,TS_B1,TS_D1,(TS_D1*bL+TS_A1)*(TS_D1*bL+TS_A1)-4.0*TS_D1*aL*(TS_D1*dL-TS_B1)
!tmp            write(*,'("tsfc1-3",6ES15.6)') g%MS(i,n)%Ldmassdt(2),g%MS(i,n)%dmassdt(1,10),&
!tmp                          ag%TV(n)%s_v_n(phase2),g%MS(i,n)%CAP,g%MS(i,n)%fv,g%MS(i,n)%fh
!tmp            write(*,'("tsfc1-4",3ES15.6)') mk_a,ag%TV(n)%k_a,&
!tmp                     sqrt(g%MS(i,n)%semi_a**2+g%MS(i,n)%semi_c**2)
!tmp            write(*,*) "ck e/T", aL*tsfc(i,n)**2+bL*tsfc(i,n)+dL, &
!tmp                   get_sat_vapor_pres_lk(phase2,tsfc(i,n),ag%estbar,ag%esitbar)/tsfc(i,n)

            if(tsfc(i,n).ge.Tmax) then
              tsfc(i,n)=min(Tmax,tsfc(i,n))
              icond2(i,n)=1
              diff(i,n)=0.0
            elseif(abs(diff(i,n))<1.0e-4) then
              ! convergence check
              icond2(i,n)=1
              g%MS(i,n)%tmp=min(tsfc(i,n),Tmax)
              g%MS(i,n)%e_sat=get_sat_vapor_pres_lk(phase2,g%MS(i,n)%tmp,ag%estbar,ag%esitbar)
            elseif(iter<5) then

            else
              icond2(i,n)=1
              g%MS(i,n)%tmp=min(tsfc(i,n),Tmax)
              g%MS(i,n)%e_sat=get_sat_vapor_pres_lk(phase2,g%MS(i,n)%tmp,ag%estbar,ag%esitbar)
            endif
          end if
        enddo
        enddo
      enddo

      if ( debug ) then
!      if(any(abs(diff(1:g%n_bin*g%L))>=1.0e-4)) then
!        do in=1,g%n_bin*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
         do n = 1, g%L
         do i = 1, g%N_BIN
            if(abs(diff(i,n))>1.0e-4) then
               write(*,'("cal_surface_temp2_vec> large difference",2I5,5ES15.6)') &
                    i,n,diff(i,n),ag%TV(n)%T_n,g%MS(i,n)%tmp,tsfc(i,n),tsfc_b(i,n)
            endif
         enddo
         enddo
!      endif
      end if

    end if

  end subroutine cal_surface_temp2_vec

!!$  subroutine cal_model_tendency_all(level,out_type,ag,mes_rc,&
!!$       gr,gs,ga,&
!!$       l,nrtype,nrbin,nrcat,nstype,nsbin,nscat,natype,nabin,nacat,&
!!$       qtp,rv,xc,xr,xs,xa,flagp_r,flagp_s,flagp_a,eps_ap0,ID,JD,KD)
!!$!!c  subroutine cal_model_tendency( level, out_type, g, th_var, mes_rc,x, dxdt)
!!$    integer,intent(in)  :: l,nrtype,nrbin,nrcat,nstype,nsbin,nscat,&
!!$         natype,nabin,nacat
!!$    integer :: ID(*),JD(*),KD(*)
!!$    type (group), intent(in)  :: gr,gs
!!$    type (group), dimension(nacat) :: ga
!!$    type (airgroup), intent(in) :: ag
!!$    ! message from reality-check
!!$!tmp    integer,pointer,dimension(:)  :: mes_rc
!!$    integer,dimension(*)   :: mes_rc
!!$    real(PS) :: qtp(*),rv(*),xc(*),XR(nrtype,nrbin,nrcat,*),XS(nstype,nsbin,nscat,*),&
!!$            XA(natype,nabin,nacat,*)
!!$    real(PS),dimension(*)  :: eps_ap0
!!$    integer,intent(in)  :: flagp_r,flagp_s,flagp_a
!!$    real (ps),dimension(mxnbin,LMAX)     :: ndxdt,ndxdt2
!!$    !real(PS) :: ndxdt3
!!$    real (ps)     :: var_n1,var_n2
!!$    ! level of complexity
!!$    integer, intent(in)  :: level
!!$    integer,intent(in)  :: out_type
!!$    ! original total volume variables
!!$!!!    real(ps),dimension(:,:,:),allocatable  :: oQ
!!$    real(ps),dimension(mxnbin,LMAX)  :: oQ
!!$!!!    real(ps),dimension(mxnbin)  :: oQ
!!$
!!$    ! mass fraction of aerosol particles to total mass of hydrometeors and
!!$    ! mass fraction of soluble aerosol particles to total aerosol mass
!!$    ! those are used to maintain positive aerosol masses, and use fraction
!!$    ! of smaller bin.
!!$    real(ps),dimension(max(gr%n_bin,gs%n_bin),L) :: fap
!!$    !real(ps) :: fapt_r,faps_r,fapt_s,faps_s, fap3
!!$    ! minimum possible fractions
!!$    real(ps),parameter :: min_fapt_r=1.0e-18_ps,min_faps_r=1.0e-5_ps&
!!$!!c    real(ps),parameter :: min_fapt_r=1.0e-18_ps,min_faps_r=0.0_ps&
!!$         ,min_fapt_s=1.0e-18_ps,min_faps_s=0.0_ps,sep_faps=1.0e-5_PS
!!$    ! possible ratio of mean mass to bin boundaries
!!$    real(ps),parameter :: brat1=1.001,brat2=0.999
!!$    real(ps),parameter :: mlmt=1.0e-30,nlmt=1.0e-30,m_lmt_ap=1.0e-25
!!$    !     the minimum possible background concentration for large particle
!!$    !     is assumed to be 1.0e-5 cm^-3
!!$    real(ps),parameter :: n_lmt_ap=1.0e-5,n_max_ap=1.0e+4
!!$!parcel model    real(ps),parameter :: n_lmt_ap=1.0e-15,n_max_ap=1.0e+4
!!$    !     The minimum radius possible for accumulation particles
!!$    real(ps),parameter :: r3_lmt=1.0e-18
!!$    !     The minimum radius possible for nucleation particles
!!$    real(PS),parameter :: r3_lmt_nuc=1.0e-21
!!$
!!$
!!$    real(ps),parameter :: mx_aprat=0.99999
!!$    real(ps) :: m_lmt
!!$    ! realistic bounds for length predictions
!!$    real(ps),parameter :: min_hexlen3=1.0e-12,max_hexlen3=27.0,&
!!$                          max_roslen3=27.0,max_irrlen3=27.0,max_exice=1.0
!!$    ! realistic bounds for crystal mass
!!$    real(ps),parameter :: m_icmin=4.763209003e-12
!!$
!!$    ! accuracy limit to prevent dubious crystal features by numerical errors
!!$    real(PS),parameter :: accuracy_lmt=1.0e-6
!!$
!!$    ! vapor specific humidity cacluated with activation scheme
!!$    real(PS) :: rv_n
!!$    ! total mixing ratio of hydrometeors
!!$    real(PS),dimension(L) :: totmixr,totmixs
!!$    real(PS),dimension(L) :: s_n,rv2
!!$    real(PS) :: e_n
!!$
!!$    integer,dimension(mxnbin*LMAX*mxntend) :: ierror
!!$    integer,dimension(mxnbin*LMAX) :: icond1,icond2
!!$
!!$    integer  :: i,j,k,n,ic,in,ijn!,ij,var_status
!!$    !integer :: ie
!!$    !integer :: ic1,ic2
!!$    integer :: iq1
!!$
!!$    !integer :: imark
!!$
!!$!!c    real(PS) :: alen_n,clen_n,den_n
!!$
!!$!    allocate(&
!!$!      ndxdt(max(maxval(ga(1:nacat)%N_bin),gr%N_bin,gs%N_bin),L), &
!!$!      ndxdt2(max(maxval(ga(1:nacat)%N_bin),gr%N_bin,gs%N_bin),L),stat=var_status)
!!$!    if(var_status/=0) then
!!$!      write(*,*) "cal_model_tendecy_all: allocation failed 1"
!!$!      stop
!!$!    endif
!!$!    allocate(&
!!$!      ierror(max(maxval(ga(1:nacat)%N_bin),gr%N_bin,gs%N_bin)* &
!!$!             max(maxval(ga(1:nacat)%n_tendpros),gr%n_tendpros,gs%n_tendpros)*L), &
!!$!      icond1(max(maxval(ga(1:nacat)%N_bin),gr%N_bin,gs%N_bin)*L), &
!!$!      icond2(max(maxval(ga(1:nacat)%N_bin),gr%N_bin,gs%N_bin)*L), &
!!$!!             stat=var_status)
!!$!    if(var_status/=0) then
!!$!      write(*,*) "cal_model_tendecy_all: allocation failed 2"
!!$!      stop
!!$!    endif
!!$!
!!$!    allocate(&
!!$!      oQ(gs%N_bin,L),stat=var_status)
!!$!!!!      oQ(gs%N_bin,L,gs%N_nonmass),stat=var_status)
!!$!    if(var_status/=0) then
!!$!      write(*,*) "cal_model_tendecy_all: allocation failed 3"
!!$!      stop
!!$!    endif
!!$
!!$
!!$!!c       do i=1,gs%n_bin
!!$!!c          write(*,'("XS out 0>",4I4,20ES14.6)') KD(n),ID(n),JD(n),i,XS(1:mxnpi,i,NSCAT,n)
!!$!!c       enddo
!!$
!!$    !
!!$    ! +++ for rain +++
!!$    !
!!$    if( flagp_r > 0 ) then
!!$      do in=1,gr%n_bin*L
!!$        n=(in-1)/gr%N_BIN+1
!!$        i=in-(n-1)*gr%N_BIN
!!$
!!$        ndxdt(i,n) = 0.0_ps
!!$        ndxdt2(i,n) = 0.0_ps
!!$      enddo
!!$
!!$      do j=1,gr%n_tendpros
!!$!CDIR NODEP
!!$        do in=1,gr%n_bin*L
!!$          n=(in-1)/gr%N_BIN+1
!!$          i=in-(n-1)*gr%N_BIN
!!$
!!$          ndxdt(i,n) = ndxdt(i,n)+gr%MS(i,n)%dmassdt(rmt,j)
!!$          ndxdt2(i,n) = ndxdt2(i,n)+gr%MS(i,n)%dcondt(j)
!!$        enddo
!!$      enddo
!!$
!!$      ierror(1:gr%n_bin*gr%n_tendpros*L)=0
!!$      do ijn=1,gr%n_bin*gr%n_tendpros*L
!!$        in=(ijn-1)/gr%n_tendpros+1
!!$        j=ijn-(in-1)*gr%n_tendpros
!!$        n=(in-1)/gr%n_bin+1
!!$        i=in-(n-1)*gr%n_bin
!!$
!!$        if(gr%MS(i,n)%dmassdt(rmt,j)>1.0e+03.or.&
!!$           gr%MS(i,n)%dmassdt(rmt,j)<-1.0e+03) then
!!$          ierror(ijn)=1
!!$        endif
!!$      enddo
!!$      if(any(ierror(1:gr%n_bin*gr%n_tendpros*L)>0)) then
!!$        do ijn=1,gr%n_bin*gr%n_tendpros*L
!!$          in=(ijn-1)/gr%n_tendpros+1
!!$          j=ijn-(in-1)*gr%n_tendpros
!!$          n=(in-1)/gr%n_bin+1
!!$          i=in-(n-1)*gr%n_bin
!!$          if(ierror(ijn)==1) then
!!$            write(*,201) n,i,j,gr%MS(i,n)%dmassdt(rmt,j)
!!$201       format("Warning mt_rain>grid,bin,process,dmassdt",3i5,es15.6)
!!$!!c           stop
!!$          endif
!!$        enddo
!!$      endif
!!$
!!$      ierror(1:gr%n_bin*gr%n_tendpros*L)=0
!!$      do ijn=1,gr%n_bin*gr%n_tendpros*L
!!$        in=(ijn-1)/gr%n_tendpros+1
!!$        j=ijn-(in-1)*gr%n_tendpros
!!$        n=(in-1)/gr%n_bin+1
!!$        i=in-(n-1)*gr%n_bin
!!$
!!$        if(gr%MS(i,n)%dcondt(j)<-1.0e+8 .or. &
!!$           gr%MS(i,n)%dcondt(j)>1.0e+8 ) then
!!$          ierror(ijn)=1
!!$        endif
!!$      enddo
!!$      if(any(ierror(1:gr%n_bin*gr%n_tendpros*L)>0)) then
!!$        do ijn=1,gr%n_bin*gr%n_tendpros*L
!!$          in=(ijn-1)/gr%n_tendpros+1
!!$          j=ijn-(in-1)*gr%n_tendpros
!!$          n=(in-1)/gr%n_bin+1
!!$          i=in-(n-1)*gr%n_bin
!!$          if(ierror(ijn)==1) then
!!$            write(*,*) "Warning con_r:bin,grid,process",i,n,j,gr%MS(i,n)%dcondt(j)
!!$          endif
!!$        enddo
!!$      endif
!!$
!!$      if( out_type == 2 ) then
!!$!CDIR NODEP
!!$        do in=1,gr%n_bin*L
!!$          n=(in-1)/gr%N_BIN+1
!!$          i=in-(n-1)*gr%N_BIN
!!$
!!$          XR(rmt_q,i,nrcat,n)=(gr%MS(i,n)%mass(rmt)+ndxdt(i,n)*gr%dt)/ag%tv(n)%den
!!$          XR(rcon_q,i,nrcat,n)=(gr%MS(i,n)%con+ndxdt2(i,n)*gr%dt)/ag%tv(n)%den
!!$        enddo
!!$        !
!!$        ! check
!!$        !
!!$        icond1(1:gr%n_bin*L)=0
!!$        icond2(1:gr%n_bin*L)=0
!!$        do in=1,gr%n_bin*L
!!$          n=(in-1)/gr%N_BIN+1
!!$          i=in-(n-1)*gr%N_BIN
!!$
!!$          if(XR(rmt_q,i,nrcat,n)*ag%TV(n)%den<gr%MS(i,n)%mass(rmt)*accuracy_lmt.or.&
!!$                   XR(rmt_q,i,nrcat,n)<1.0e-30) then
!!$            icond1(in)=1
!!$          elseif(XR(rmt_q,i,nrcat,n)>1.0e-1) then
!!$            icond1(in)=2
!!$          end if
!!$
!!$          if(XR(rcon_q,i,nrcat,n)*ag%TV(n)%den<gr%MS(i,n)%con*accuracy_lmt.or.&
!!$                     XR(rcon_q,i,nrcat,n)<1.0e-30) then
!!$            icond2(in)=1
!!$          else
!!$            icond2(in)=3
!!$
!!$          endif
!!$        enddo
!!$
!!$        do ijn=1,gr%n_bin*nrtype*L
!!$          in=(ijn-1)/nrtype+1
!!$          j=ijn-(in-1)*nrtype
!!$          n=(in-1)/gr%n_bin+1
!!$          i=in-(n-1)*gr%n_bin
!!$          if(icond1(in)==1.or.icond2(in)==1) then
!!$            XR(j,i,nrcat,n)=0.0_ps
!!$          endif
!!$        end do
!!$!CDIR NODEP
!!$        do in=1,gr%n_bin*L
!!$          n=(in-1)/gr%N_BIN+1
!!$          i=in-(n-1)*gr%N_BIN
!!$          if(icond1(in)/=1.and.icond2(in)==3) then
!!$!            write(*,*) "xr check",i,n,xr(rcon_q,i,nrcat,n),xr(rmt_q,i,nrcat,n)&
!!$!                   ,gr%binb(i),gr%binb(i+1)
!!$            XR(rcon_q,i,nrcat,n)=&
!!$                  XR(rmt_q,i,nrcat,n)&
!!$                  /max(brat1*gr%binb(i),min(XR(rmt_q,i,nrcat,n)/XR(rcon_q,i,nrcat,n)&
!!$                 ,brat2*gr%binb(i+1)))
!!$          endif
!!$        enddo
!!$        if(any(icond1(1:gr%n_bin*L)==2)) then
!!$          do in=1,gr%n_bin*L
!!$            n=(in-1)/gr%N_BIN+1
!!$            i=in-(n-1)*gr%N_BIN
!!$            if(icond1(in)==2) then
!!$              write(*,251) i,n,XR(rmt_q,i,nrcat,n)
!!$251        format("warning:qr is large at bin,grid:",2i5,2es15.6)
!!$            endif
!!$          enddo
!!$        endif
!!$      endif
!!$      if(level>=4) then
!!$        do in=1,gr%n_bin*L
!!$          n=(in-1)/gr%N_BIN+1
!!$          i=in-(n-1)*gr%N_BIN
!!$          ndxdt(i,n) = 0.0_ps
!!$          ndxdt2(i,n) = 0.0_ps
!!$        enddo
!!$        do j=1,gr%n_tendpros
!!$!CDIR NODEP
!!$          do in=1,gr%n_bin*L
!!$            n=(in-1)/gr%N_BIN+1
!!$            i=in-(n-1)*gr%N_BIN
!!$
!!$            ! mass by total mass of aerosols (g/g)
!!$            ndxdt(i,n)=ndxdt(i,n)+gr%MS(i,n)%dmassdt(rmat,j)
!!$
!!$            ! mass by soluble mass of aerosols (g/g)
!!$            ndxdt2(i,n)=ndxdt2(i,n)+gr%MS(i,n)%dmassdt(rmas,j)
!!$          enddo
!!$        enddo
!!$
!!$        ierror(1:gr%n_bin*gr%n_tendpros*L)=0
!!$        do ijn=1,gr%n_bin*gr%n_tendpros*L
!!$          in=(ijn-1)/gr%n_tendpros+1
!!$          j=ijn-(in-1)*gr%n_tendpros
!!$          n=(in-1)/gr%n_bin+1
!!$          i=in-(n-1)*gr%n_bin
!!$          if(gr%MS(i,n)%dmassdt(rmat,j)>1.0e+08.or.&
!!$             gr%MS(i,n)%dmassdt(rmat,j)<-1.0e+08)then
!!$            ierror(ijn)=1
!!$          endif
!!$        enddo
!!$        if(any(ierror(1:gr%n_bin*gr%n_tendpros*L)>0)) then
!!$          do ijn=1,gr%n_bin*gr%n_tendpros*L
!!$            in=(ijn-1)/gr%n_tendpros+1
!!$            j=ijn-(in-1)*gr%n_tendpros
!!$            n=(in-1)/gr%n_bin+1
!!$            i=in-(n-1)*gr%n_bin
!!$            if(ierror(ijn)==1) then
!!$              write(*,236) i,n,j,gr%MS(i,n)%dmassdt(rmat,j)
!!$236     format("Warning aptm in rain>bin,grid,process,dmassdt",3i5,es15.6)
!!$            endif
!!$          enddo
!!$        endif
!!$
!!$        ierror(1:gr%n_bin*gr%n_tendpros*L)=0
!!$        do ijn=1,gr%n_bin*gr%n_tendpros*L
!!$          in=(ijn-1)/gr%n_tendpros+1
!!$          j=ijn-(in-1)*gr%n_tendpros
!!$          n=(in-1)/gr%n_bin+1
!!$          i=in-(n-1)*gr%n_bin
!!$          if(gr%MS(i,n)%dmassdt(rmas,j)>1.0e+08.or.&
!!$             gr%MS(i,n)%dmassdt(rmas,j)<-1.0e+08)then
!!$            ierror(ijn)=1
!!$          endif
!!$        enddo
!!$        if(any(ierror(1:gr%n_bin*gr%n_tendpros*L)>0)) then
!!$          do ijn=1,gr%n_bin*gr%n_tendpros*L
!!$            in=(ijn-1)/gr%n_tendpros+1
!!$            j=ijn-(in-1)*gr%n_tendpros
!!$            n=(in-1)/gr%n_bin+1
!!$            i=in-(n-1)*gr%n_bin
!!$            if(ierror(ijn)==1) then
!!$              write(*,237) i,n,j,gr%MS(i,n)%dmassdt(rmas,j)
!!$237       format("Warning apsm in solid>bin,grid,process,dmassdt",3i5,es15.6)
!!$            endif
!!$          enddo
!!$        endif
!!$
!!$        if( out_type == 2 ) then
!!$!CDIR NODEP
!!$          do in=1,gr%n_bin*L
!!$            n=(in-1)/gr%N_BIN+1
!!$            i=in-(n-1)*gr%N_BIN
!!$            if(icond1(in)/=1.and.icond2(in)/=1) then
!!$              var_n1=(gr%MS(i,n)%mass(rmat)+ndxdt(i,n)*gr%dt)/ag%tv(n)%den
!!$              XR(rmat_q,i,nrcat,n)=var_n1
!!$              fap(i,n)=var_n1/XR(rmt_q,i,nrcat,n)
!!$!            else
!!$!              fap(i,n)=1.0
!!$            endif
!!$          enddo
!!$!CDIR NODEP
!!$          do in=1,gr%n_bin*L
!!$            n=(in-1)/gr%N_BIN+1
!!$            i=in-(n-1)*gr%N_BIN
!!$            if(icond1(in)/=1.and.icond2(in)/=1) then
!!$              if(fap(i,n)>1.0_PS) then
!!$                XR(rmt_q,i,nrcat,n)=XR(rmat_q,i,nrcat,n)
!!$              elseif(fap(i,n)<=min_fapt_r) then
!!$                XR(rmat_q,i,nrcat,n)=min(max(m_lmt_ap/ag%tv(n)%den,min_fapt_r*XR(rmt_q,i,nrcat,n))&
!!$                                  ,mx_aprat*XR(rmt_q,i,nrcat,n))
!!$              endif
!!$            endif
!!$          enddo
!!$!CDIR NODEP
!!$          do in=1,gr%n_bin*L
!!$            n=(in-1)/gr%N_BIN+1
!!$            i=in-(n-1)*gr%N_BIN
!!$            if(icond1(in)/=1.and.icond2(in)/=1) then
!!$              var_n1=(gr%MS(i,n)%mass(rmas)+ndxdt2(i,n)*gr%dt)/ag%tv(n)%den
!!$              XR(rmas_q,i,nrcat,n)=var_n1
!!$              fap(i,n)=var_n1/XR(rmat_q,i,nrcat,n)
!!$!            else
!!$!              fap(i,n)=1.0
!!$            endif
!!$          enddo
!!$!CDIR NODEP
!!$          do in=1,gr%n_bin*L
!!$            n=(in-1)/gr%N_BIN+1
!!$            i=in-(n-1)*gr%N_BIN
!!$            if(icond1(in)/=1.and.icond2(in)/=1) then
!!$              if(fap(i,n)<min_faps_r) then
!!$                XR(rmas_q,i,nrcat,n)=min(max(m_lmt_ap/ag%tv(n)%den,&
!!$                                   min_faps_r*XR(rmat_q,i,nrcat,n)),XR(rmat_q,i,nrcat,n))
!!$              endif
!!$            endif
!!$          enddo
!!$        endif
!!$      endif
!!$
!!$      totmixr=0.0_PS
!!$      do i=1,gr%n_bin
!!$        do n=1,L
!!$          if(icond1(in)/=1.and.icond2(in)/=1) then
!!$            totmixr(n)=totmixr(n)+max(0.0_PS,XR(rmt_q,i,nrcat,n)-XR(rmat_q,i,nrcat,n))
!!$          endif
!!$        enddo
!!$      end do
!!$
!!$
!!$      if(debug) then
!!$        do n=1,L
!!$!         if(sum(XR(rcon_q,1:gr%n_bin,nrcat,n))*ag%tv(n)%den>100.0) then
!!$         if(sum(XR(rcon_q,1:gr%n_bin,nrcat,n))*ag%tv(n)%den>0.0) then
!!$         write(*,'("cal_model_tendency:bin,i,j,k,total rcon, others",4i5,100es15.6)') n,ID(n),JD(n),KD(n) &
!!$                      ,sum(XR(rcon_q,1:gr%n_bin,nrcat,n))*ag%tv(n)%den &
!!$                      ,XR(rcon_q,1:gr%n_bin,nrcat,n)*ag%tv(n)%den
!!$
!!$         write(*,'("cal_model_tendency:bin,i,j,k,total act, total vap mass tendency",4i5,100es15.6)') n,ID(n),JD(n),KD(n) &
!!$                      ,sum(gr%MS(1:gr%n_bin,n)%dmassdt(rmt,3)) &  ! activation
!!$                      ,sum(gr%MS(1:gr%n_bin,n)%dmassdt(rmt,1))    ! vapor
!!$         write(*,*) "cal_model_tendency:con-bin,i,j,k,total rim, total melt, total mossop, total coal"
!!$         write(*,'(4I5,10ES15.6)') n,ID(n),JD(n),KD(n) &
!!$                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(4)) &  ! riming
!!$                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(10)) &  ! melting_shedding
!!$                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(9)) &  ! mossop and hallet
!!$                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(2))    ! collision-coal
!!$         endif
!!$       enddo
!!$      endif
!!$
!!$      ! for cloud droplets
!!$      do n=1,L
!!$        XC(n)=XR(rmt_q,1,nrcat,n)
!!$      enddo
!!$    endif
!!$
!!$    !
!!$    ! ice
!!$    !
!!$    if( flagp_s > 0 ) then
!!$      do in=1,gs%n_bin*L
!!$        n=(in-1)/gs%N_BIN+1
!!$        i=in-(n-1)*gs%N_BIN
!!$
!!$        ndxdt(i,n) = 0.0_ps
!!$        ndxdt2(i,n) = 0.0_ps
!!$      enddo
!!$      do j=1,gs%n_tendpros
!!$!CDIR NODEP
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$          ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(imt,j)
!!$          ndxdt2(i,n)=ndxdt2(i,n)+gs%MS(i,n)%dcondt(j)
!!$        enddo
!!$      enddo
!!$
!!$      ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!!$      do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$        in=(ijn-1)/gs%n_tendpros+1
!!$        j=ijn-(in-1)*gs%n_tendpros
!!$        n=(in-1)/gs%n_bin+1
!!$        i=in-(n-1)*gs%n_bin
!!$        if(gs%MS(i,n)%dmassdt(imt,j)>1.0e+08.or.&
!!$           gs%MS(i,n)%dmassdt(imt,j)<-1.0e+08)then
!!$          ierror(ijn)=1
!!$        endif
!!$      enddo
!!$      if(any(ierror(1:gs%n_bin*gs%n_tendpros*L)>0)) then
!!$        do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$          in=(ijn-1)/gs%n_tendpros+1
!!$          j=ijn-(in-1)*gs%n_tendpros
!!$          n=(in-1)/gs%n_bin+1
!!$          i=in-(n-1)*gs%n_bin
!!$          if(ierror(ijn)==1) then
!!$            write(*,202) n,i,j,gs%MS(i,n)%dmassdt(imt,j)
!!$202     format("Warning mt_s>grid,bin,process,dmassdt",3i5,es15.6)
!!$          endif
!!$        enddo
!!$      endif
!!$
!!$      ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!!$      do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$        in=(ijn-1)/gs%n_tendpros+1
!!$        j=ijn-(in-1)*gs%n_tendpros
!!$        n=(in-1)/gs%n_bin+1
!!$        i=in-(n-1)*gs%n_bin
!!$        if(gs%MS(i,n)%dcondt(j)<-1.0e+8 .or. &
!!$           gs%MS(i,n)%dcondt(j)>1.0e+8 ) then
!!$          ierror(ijn)=1
!!$        endif
!!$      enddo
!!$      if(any(ierror(1:gs%n_bin*gs%n_tendpros*L)>0)) then
!!$        do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$          in=(ijn-1)/gs%n_tendpros+1
!!$          j=ijn-(in-1)*gs%n_tendpros
!!$          n=(in-1)/gs%n_bin+1
!!$          i=in-(n-1)*gs%n_bin
!!$          if(ierror(ijn)==1) then
!!$            write(*,*) "Warning con_s:bin,grid,process",i,n,j,gs%MS(i,n)%dcondt(j)
!!$          endif
!!$        enddo
!!$      endif
!!$
!!$      if( out_type == 2 ) then
!!$!CDIR NODEP
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$          XS(imt_q,i,nscat,n)=(gs%MS(i,n)%mass(imt)+ndxdt(i,n)*gs%dt)/ag%tv(n)%den
!!$          XS(icon_q,i,nscat,n)=(gs%MS(i,n)%con+ndxdt2(i,n)*gs%dt)/ag%tv(n)%den
!!$        enddo
!!$        !
!!$        ! check
!!$        !
!!$        icond1(1:gs%n_bin*L)=0
!!$        icond2(1:gs%n_bin*L)=0
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$          if(XS(imt_q,i,nscat,n)*ag%TV(n)%den<gs%MS(i,n)%mass(imt)*accuracy_lmt.or.&
!!$                     XS(imt_q,i,nscat,n)<1.0e-30) then
!!$            icond1(in)=1
!!$          elseif(XS(imt_q,i,nscat,n)>1.0e-1) then
!!$            icond1(in)=2
!!$          endif
!!$          if(XS(icon_q,i,nscat,n)*ag%TV(n)%den<gs%MS(i,n)%con*accuracy_lmt.or.&
!!$                     XS(icon_q,i,nscat,n)<1.0e-30) then
!!$            icond2(in)=1
!!$          else
!!$            icond2(in)=3
!!$
!!$          endif
!!$        enddo
!!$        do ijn=1,gs%n_bin*nstype*L
!!$          in=(ijn-1)/nstype+1
!!$          j=ijn-(in-1)*nstype
!!$          n=(in-1)/gs%n_bin+1
!!$          i=in-(n-1)*gs%n_bin
!!$          if(icond1(in)==1.or.icond2(in)==1) then
!!$            XS(j,i,nscat,n)=0.0_ps
!!$          endif
!!$        end do
!!$!CDIR NODEP
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$          if(icond1(in)/=1.and.icond2(in)==3) then
!!$            XS(icon_q,i,nscat,n)=&
!!$                 XS(imt_q,i,nscat,n)&
!!$                  /max(brat1*gs%binb(i),min(XS(imt_q,i,nscat,n)/XS(icon_q,i,nscat,n)&
!!$                 ,brat2*gs%binb(i+1)))
!!$          endif
!!$        enddo
!!$        if(any(icond1(1:gs%n_bin*L)==2)) then
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            if(icond1(in)==2) then
!!$              write(*,252) i,n,XS(imt_q,i,nscat,n)
!!$252      format("warning:qs is large at bin,grid:",2i5,2es15.6)
!!$            endif
!!$          enddo
!!$        endif
!!$      endif
!!$
!!$      ! volume of circumscribing sphere (cm^3/g) and
!!$      ! volume of axis length production (cm^3/g)
!!$      iq1=ivcs_q
!!$      nonmass_loop: do k=1,gs%N_nonmass
!!$        !      1. volume of circumscribing sphere (cm^3/g)
!!$        !      2. con. weighted a-axis length^3 (cm^3/g)
!!$        !      3. con. weighted c-axis length^3 (cm^3/g)
!!$        !      4. con. weighted d-axis length^3 (cm^3/g)
!!$        !      5. con. weighted r-axis length^3 (cm^3/g) (rosetta bulletes)
!!$        !      6. con. weighted e-axis length^3 (cm^3/g) (polycrystals)
!!$
!!$        if(k==ivcs) then
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            oq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%v_cs
!!$          enddo
!!$        elseif(k==iacr) then
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            oq(i,n) = gs%MS(i,n)%con*gs%MS(i,n)%a_len**3.0
!!$          enddo
!!$        elseif(k==iccr) then
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            oq(i,n) = gs%MS(i,n)%con*gs%MS(i,n)%c_len**3.0
!!$          enddo
!!$        elseif(k==idcr) then
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            oq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%d**3.0
!!$          enddo
!!$        elseif(k==iag) then
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            oq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%ag**3.0
!!$          enddo
!!$        elseif(k==icg) then
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            oq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%cg**3.0
!!$          enddo
!!$        elseif(k==inex) then
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            oq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%n_exice
!!$          enddo
!!$        endif
!!$
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$
!!$          ndxdt(i,n) = 0.0_ps
!!$        enddo
!!$
!!$        do j=1,gs%n_tendpros
!!$!CDIR NODEP
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dvoldt(k,j)
!!$          enddo
!!$        enddo
!!$
!!$        ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!!$        do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$          in=(ijn-1)/gs%n_tendpros+1
!!$          j=ijn-(in-1)*gs%n_tendpros
!!$          n=(in-1)/gs%n_bin+1
!!$          i=in-(n-1)*gs%n_bin
!!$          if(gs%MS(i,n)%dvoldt(k,j)>1.0e+20.or.&
!!$             gs%MS(i,n)%dvoldt(k,j)<-1.0e+20)then
!!$            ierror(ijn)=1
!!$          endif
!!$        enddo
!!$        if(any(ierror(1:gs%n_bin*gs%n_tendpros*L)>0)) then
!!$          do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$            in=(ijn-1)/gs%n_tendpros+1
!!$            j=ijn-(in-1)*gs%n_tendpros
!!$            n=(in-1)/gs%n_bin+1
!!$            i=in-(n-1)*gs%n_bin
!!$
!!$            if(ierror(ijn)==1) then
!!$              write(*,231) n,k,j,gs%MS(i,n)%dvoldt(k,j)
!!$231         format("axis>grid,axis#,process,dvoldt",3i5,es15.6)
!!$              write(*,'("phi_cs,den",2es15.6)') gs%IS(i,n)%phi_cs,gs%MS(i,n)%den
!!$              write(*,'("psi_ic,phi_ic,alen",3es15.6)') gs%IS(i,n)%psi_ic,gs%IS(i,n)%phi_ic,gs%MS(i,n)%a_len
!!$              write(*,'("i,ht,sh",3i4)') i,gs%IS(i,n)%habit,gs%IS(i,n)%sh_type
!!$
!!$            endif
!!$          enddo
!!$        endif
!!$
!!$        if( out_type == 2 ) then
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            if(icond1(in)/=1.and.icond2(in)/=1) then
!!$              XS(iq1-1+k,i,nscat,n)=max(0.0_PS,(oq(i,n)+ndxdt(i,n)*gs%dt)/ag%tv(n)%den)
!!$            endif
!!$          enddo
!!$        endif
!!$
!!$      enddo nonmass_loop
!!$
!!$      ! bound with maximum and minium possible length.
!!$!CDIR NODEP
!!$      do in=1,gs%n_bin*L
!!$        n=(in-1)/gs%N_BIN+1
!!$        i=in-(n-1)*gs%N_BIN
!!$
!!$        if(icond1(in)/=1.and.icond2(in)/=1) then
!!$          XS(iacr_q,i,nscat,n)=max(min_hexlen3*XS(icon_q,i,nscat,n),&
!!$                  min(max_hexlen3*XS(icon_q,i,nscat,n),XS(iacr_q,i,nscat,n)))
!!$          XS(iccr_q,i,nscat,n)=max(min_hexlen3*XS(icon_q,i,nscat,n),&
!!$                  min(max_hexlen3*XS(icon_q,i,nscat,n),XS(iccr_q,i,nscat,n)))
!!$
!!$          XS(idcr_q,i,nscat,n)=max(0.0_RP,min(XS(iacr_q,i,nscat,n),XS(idcr_q,i,nscat,n)))
!!$
!!$          XS(iag_q,i,nscat,n)=max(0.0_RP,min(max_roslen3*XS(icon_q,i,nscat,n),XS(iag_q,i,nscat,n)))
!!$          XS(icg_q,i,nscat,n)=max(0.0_RP,min(max_irrlen3*XS(icon_q,i,nscat,n),XS(icg_q,i,nscat,n)))
!!$          XS(inex_q,i,nscat,n)=max(0.0_RP,min(max_exice*XS(icon_q,i,nscat,n),XS(inex_q,i,nscat,n)))
!!$        endif
!!$      enddo
!!$
!!$
!!$      ! another constraints on d-axis length and a and c coordinates according to
!!$      ! the a and c axis lengths
!!$!CDIR NODEP
!!$      do in=1,gs%n_bin*L
!!$        n=(in-1)/gs%N_BIN+1
!!$        i=in-(n-1)*gs%N_BIN
!!$        if(XS(iacr_q,i,nscat,n)<=1.0e-9_PS*XS(icon_q,i,nscat,n).or.&
!!$           XS(iccr_q,i,nscat,n)<=1.0e-12_PS*XS(icon_q,i,nscat,n)) then
!!$          XS(idcr_q,i,nscat,n)=0.0_PS
!!$        end if
!!$        if(XS(iacr_q,i,nscat,n)<=min_hexlen3*XS(icon_q,i,nscat,n).or.&
!!$           XS(iccr_q,i,nscat,n)<=min_hexlen3*XS(icon_q,i,nscat,n)) then
!!$          XS(iag_q,i,nscat,n)=0.0_PS
!!$          XS(icg_q,i,nscat,n)=0.0_PS
!!$          XS(inex_q,i,nscat,n)=0.0_PS
!!$        end if
!!$      enddo
!!$
!!$      !
!!$      ! mass by ice crystals (g/g)
!!$      !
!!$      do in=1,gs%n_bin*L
!!$        n=(in-1)/gs%N_BIN+1
!!$        i=in-(n-1)*gs%N_BIN
!!$
!!$        ndxdt(i,n) = 0.0_ps
!!$      enddo
!!$      do j=1,gs%n_tendpros
!!$!CDIR NODEP
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$          ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(imc,j)
!!$        enddo
!!$      enddo
!!$
!!$      ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!!$      do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$        in=(ijn-1)/gs%n_tendpros+1
!!$        j=ijn-(in-1)*gs%n_tendpros
!!$        n=(in-1)/gs%n_bin+1
!!$        i=in-(n-1)*gs%n_bin
!!$        if(gs%MS(i,n)%dmassdt(imc,j)>1.0e+08.or.&
!!$           gs%MS(i,n)%dmassdt(imc,j)<-1.0e+08)then
!!$          ierror(ijn)=1
!!$        endif
!!$      enddo
!!$      if(any(ierror(1:gs%n_bin*gs%n_tendpros*L)>0)) then
!!$        do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$          in=(ijn-1)/gs%n_tendpros+1
!!$          j=ijn-(in-1)*gs%n_tendpros
!!$          n=(in-1)/gs%n_bin+1
!!$          i=in-(n-1)*gs%n_bin
!!$          if(ierror(ijn)==1) then
!!$            write(*,233) n,j,gs%MS(i,n)%dmassdt(imc,j)
!!$233     format("Warning icemass>grid,process,dmassdt",2i5,es15.6)
!!$          endif
!!$        enddo
!!$      endif
!!$
!!$      if( out_type == 2 ) then
!!$!CDIR NODEP
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$          if(icond1(in)/=1.and.icond2(in)/=1) then
!!$            XS(imc_q,i,nscat,n)=min(XS(imt_q,i,NSCAT,n),&
!!$                     max((gs%MS(i,n)%mass(imc)+ndxdt(i,n)*gs%dt)/ag%tv(n)%den,m_icmin*XS(icon_q,i,NSCAT,n)))
!!$          endif
!!$        enddo
!!$      end if
!!$
!!$      !
!!$      ! mass by riming process (g/g)
!!$      !
!!$      do in=1,gs%n_bin*L
!!$        n=(in-1)/gs%N_BIN+1
!!$        i=in-(n-1)*gs%N_BIN
!!$
!!$        ndxdt(i,n) = 0.0_ps
!!$      enddo
!!$      do j=1,gs%n_tendpros
!!$!CDIR NODEP
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$          ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(imr,j)
!!$        enddo
!!$      enddo
!!$
!!$      ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!!$      do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$        in=(ijn-1)/gs%n_tendpros+1
!!$        j=ijn-(in-1)*gs%n_tendpros
!!$        n=(in-1)/gs%n_bin+1
!!$        i=in-(n-1)*gs%n_bin
!!$        if(gs%MS(i,n)%dmassdt(imr,j)>1.0e+03.or.&
!!$           gs%MS(i,n)%dmassdt(imr,j)<-1.0e+03)then
!!$          ierror(ijn)=1
!!$        endif
!!$      enddo
!!$      if(any(ierror(1:gs%n_bin*gs%n_tendpros*L)>0)) then
!!$        do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$          in=(ijn-1)/gs%n_tendpros+1
!!$          j=ijn-(in-1)*gs%n_tendpros
!!$          n=(in-1)/gs%n_bin+1
!!$          i=in-(n-1)*gs%n_bin
!!$          if(ierror(ijn)==1) then
!!$            write(*,232) n,j,gs%MS(i,n)%dmassdt(imr,j)
!!$232     format("Warning rimemass>grid,process,dmassdt",2i5,es15.6)
!!$          endif
!!$        enddo
!!$      endif
!!$
!!$      if( out_type == 2 ) then
!!$!CDIR NODEP
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$          if(icond1(in)/=1.and.icond2(in)/=1) then
!!$            XS(imr_q,i,nscat,n)=min(XS(imt_q,i,NSCAT,n)-XS(imc_q,i,NSCAT,n),&
!!$                     max((gs%MS(i,n)%mass(imr)+ndxdt(i,n)*gs%dt)/ag%tv(n)%den,0.0_ps))
!!$          endif
!!$        enddo
!!$      end if
!!$      !
!!$      ! mass by agg process (g/g)
!!$      !
!!$      do in=1,gs%n_bin*L
!!$        n=(in-1)/gs%N_BIN+1
!!$        i=in-(n-1)*gs%N_BIN
!!$
!!$        ndxdt(i,n) = 0.0_ps
!!$      enddo
!!$      do j=1,gs%n_tendpros
!!$!CDIR NODEP
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$          ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(ima,j)
!!$        enddo
!!$      enddo
!!$      ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!!$      do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$        in=(ijn-1)/gs%n_tendpros+1
!!$        j=ijn-(in-1)*gs%n_tendpros
!!$        n=(in-1)/gs%n_bin+1
!!$        i=in-(n-1)*gs%n_bin
!!$        if(gs%MS(i,n)%dmassdt(ima,j)>1.0e+03.or.&
!!$           gs%MS(i,n)%dmassdt(ima,j)<-1.0e+03)then
!!$          ierror(ijn)=1
!!$        endif
!!$      enddo
!!$      if(any(ierror(1:gs%n_bin*gs%n_tendpros*L)>0)) then
!!$        do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$          in=(ijn-1)/gs%n_tendpros+1
!!$          j=ijn-(in-1)*gs%n_tendpros
!!$          n=(in-1)/gs%n_bin+1
!!$          i=in-(n-1)*gs%n_bin
!!$          if(ierror(ijn)==1) then
!!$            write(*,247) n,j,gs%MS(i,n)%dmassdt(ima,j)
!!$247     format("Warning aggmass>grid,process,dmassdt",2i5,es15.6)
!!$          endif
!!$        enddo
!!$      endif
!!$
!!$      if( out_type == 2 ) then
!!$!CDIR NODEP
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$          if(icond1(in)/=1.and.icond2(in)/=1) then
!!$            XS(ima_q,i,nscat,n)=min(XS(imt_q,i,NSCAT,n)-XS(imc_q,i,NSCAT,n),&
!!$                     max((gs%MS(i,n)%mass(ima)+ndxdt(i,n)*gs%dt)/ag%tv(n)%den,0.0_ps))
!!$          endif
!!$        enddo
!!$      end if
!!$
!!$      if(level>=6) then
!!$        !
!!$        ! melt water mass  (g/g)
!!$        !
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$
!!$          ndxdt(i,n) = 0.0_ps
!!$        enddo
!!$
!!$        do j=1,gs%n_tendpros
!!$!CDIR NODEP
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(imw,j)
!!$          enddo
!!$        enddo
!!$
!!$        ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!!$        do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$          in=(ijn-1)/gs%n_tendpros+1
!!$          j=ijn-(in-1)*gs%n_tendpros
!!$          n=(in-1)/gs%n_bin+1
!!$          i=in-(n-1)*gs%n_bin
!!$          if(gs%MS(i,n)%dmassdt(imw,j)>1.0e+08.or.&
!!$             gs%MS(i,n)%dmassdt(imw,j)<-1.0e+08)then
!!$            ierror(ijn)=1
!!$          endif
!!$        enddo
!!$        if(any(ierror(1:gs%n_bin*gs%n_tendpros*L)>0)) then
!!$          do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$            in=(ijn-1)/gs%n_tendpros+1
!!$            j=ijn-(in-1)*gs%n_tendpros
!!$            n=(in-1)/gs%n_bin+1
!!$            i=in-(n-1)*gs%n_bin
!!$            if(ierror(ijn)==1) then
!!$              write(*,239) n,j,gs%MS(i,n)%dmassdt(imw,j)
!!$239     format("Warning meltmass>grid,process,dmassdt",2i5,es15.6)
!!$            endif
!!$          enddo
!!$        endif
!!$
!!$        if( out_type == 2 ) then
!!$!CDIR NODEP
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            if(icond1(in)/=1.and.icond2(in)/=1) then
!!$              XS(imw_q,i,nscat,n)=min(max(XS(imt_q,i,NSCAT,n)-XS(imc_q,i,NSCAT,n)-&
!!$                        XS(imr_q,i,NSCAT,N),0.0_PS),&
!!$                        max((gs%MS(i,n)%mass(imw)+ndxdt(i,n)*gs%dt)/ag%tv(n)%den,0.0_ps))
!!$            endif
!!$          enddo
!!$        end if
!!$        !
!!$        ! freezing nucleation mass  (g/g)
!!$        !
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$
!!$          ndxdt(i,n) = 0.0_ps
!!$        enddo
!!$        do j=1,gs%n_tendpros
!!$!CDIR NODEP
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(imf,j)
!!$          enddo
!!$        enddo
!!$
!!$        ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!!$        do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$          in=(ijn-1)/gs%n_tendpros+1
!!$          j=ijn-(in-1)*gs%n_tendpros
!!$          n=(in-1)/gs%n_bin+1
!!$          i=in-(n-1)*gs%n_bin
!!$          if(gs%MS(i,n)%dmassdt(imf,j)>1.0e+08.or.&
!!$             gs%MS(i,n)%dmassdt(imf,j)<-1.0e+08)then
!!$            ierror(ijn)=1
!!$          endif
!!$        enddo
!!$        if(any(ierror(1:gs%n_bin*gs%n_tendpros*L)>0)) then
!!$          do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$            in=(ijn-1)/gs%n_tendpros+1
!!$            j=ijn-(in-1)*gs%n_tendpros
!!$            n=(in-1)/gs%n_bin+1
!!$            i=in-(n-1)*gs%n_bin
!!$            if(ierror(ijn)==1) then
!!$              write(*,249) n,j,gs%MS(i,n)%dmassdt(imf,j)
!!$249          format("Warning freezmass>grid,process,dmassdt",2i5,es15.6)
!!$            endif
!!$          enddo
!!$        endif
!!$
!!$        if( out_type == 2 ) then
!!$!CDIR NODEP
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            if(icond1(in)/=1.and.icond2(in)/=1) then
!!$              XS(imf_q,i,nscat,n)=min(XS(imc_q,i,NSCAT,n),&
!!$                        max((gs%MS(i,n)%mass(imf)+ndxdt(i,n)*gs%dt)/ag%tv(n)%den,0.0_ps))
!!$            endif
!!$          enddo
!!$        end if
!!$      endif
!!$
!!$      if(level>=4) then
!!$        do in=1,gs%n_bin*L
!!$          n=(in-1)/gs%N_BIN+1
!!$          i=in-(n-1)*gs%N_BIN
!!$          ndxdt(i,n) = 0.0_ps
!!$          ndxdt2(i,n) = 0.0_ps
!!$        enddo
!!$        do j=1,gs%n_tendpros
!!$!CDIR NODEP
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$
!!$            ! mass by total mass of aerosols (g/g)
!!$            ndxdt(i,n)=ndxdt(i,n)+gs%MS(i,n)%dmassdt(imat,j)
!!$
!!$            ! mass by soluble mass of aerosols (g/g)
!!$            ndxdt2(i,n)=ndxdt2(i,n)+gs%MS(i,n)%dmassdt(imas,j)
!!$          enddo
!!$        enddo
!!$
!!$        ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!!$        do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$          in=(ijn-1)/gs%n_tendpros+1
!!$          j=ijn-(in-1)*gs%n_tendpros
!!$          n=(in-1)/gs%n_bin+1
!!$          i=in-(n-1)*gs%n_bin
!!$          if(gs%MS(i,n)%dmassdt(imat,j)>1.0e+03.or.&
!!$             gs%MS(i,n)%dmassdt(imat,j)<-1.0e+03)then
!!$            ierror(ijn)=1
!!$          endif
!!$        enddo
!!$        if(any(ierror(1:gs%n_bin*gs%n_tendpros*L)>0)) then
!!$          do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$            in=(ijn-1)/gs%n_tendpros+1
!!$            j=ijn-(in-1)*gs%n_tendpros
!!$            n=(in-1)/gs%n_bin+1
!!$            i=in-(n-1)*gs%n_bin
!!$            if(ierror(ijn)==1) then
!!$              write(*,234) KD(n),ID(n),JD(n),i,j,gs%MS(i,n)%dmassdt(imat,j)
!!$234       format("Warning aptmass>grid(k,i,j),bin,process,dmassdt",5i5,es15.6)
!!$            endif
!!$          enddo
!!$        endif
!!$
!!$        ierror(1:gs%n_bin*gs%n_tendpros*L)=0
!!$        do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$          in=(ijn-1)/gs%n_tendpros+1
!!$          j=ijn-(in-1)*gs%n_tendpros
!!$          n=(in-1)/gs%n_bin+1
!!$          i=in-(n-1)*gs%n_bin
!!$          if(gs%MS(i,n)%dmassdt(imas,j)>1.0e+03.or.&
!!$             gs%MS(i,n)%dmassdt(imas,j)<-1.0e+03)then
!!$            ierror(ijn)=1
!!$          endif
!!$        enddo
!!$        if(any(ierror(1:gs%n_bin*gs%n_tendpros*L)>0)) then
!!$          do ijn=1,gs%n_bin*gs%n_tendpros*L
!!$            in=(ijn-1)/gs%n_tendpros+1
!!$            j=ijn-(in-1)*gs%n_tendpros
!!$            n=(in-1)/gs%n_bin+1
!!$            i=in-(n-1)*gs%n_bin
!!$            if(ierror(ijn)==1) then
!!$              write(*,235) KD(n),ID(n),JD(n),i,j,gs%MS(i,n)%dmassdt(imas,j)
!!$235       format("Warning apsmass>grid(k,i,j),bin,process,dmassdt",5i5,es15.6)
!!$            endif
!!$          enddo
!!$        endif
!!$
!!$        if( out_type == 2 ) then
!!$!CDIR NODEP
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            if(icond1(in)/=1.and.icond2(in)/=1) then
!!$              var_n1=(gs%MS(i,n)%mass(imat)+ndxdt(i,n)*gs%dt)/ag%tv(n)%den
!!$              XS(imat_q,i,nscat,n)=var_n1
!!$              fap(i,n)=var_n1/XS(imt_q,i,nscat,n)
!!$!            else
!!$!              fap(i,n)=1.0
!!$            endif
!!$          enddo
!!$!CDIR NODEP
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            if(icond1(in)/=1.and.icond2(in)/=1) then
!!$              if(fap(i,n)>1.0_PS) then
!!$                XS(imt_q,i,nscat,n)=XS(imat_q,i,nscat,n)
!!$              elseif(fap(i,n)<=min_fapt_s) then
!!$                XS(imat_q,i,nscat,n)=min(max(m_lmt_ap/ag%tv(n)%den,min_fapt_s*XS(imt_q,i,nscat,n))&
!!$                                  ,mx_aprat*XS(imt_q,i,nscat,n))
!!$              endif
!!$            endif
!!$          enddo
!!$!CDIR NODEP
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            if(icond1(in)/=1.and.icond2(in)/=1) then
!!$              var_n1=(gs%MS(i,n)%mass(imas)+ndxdt2(i,n)*gs%dt)/ag%tv(n)%den
!!$              XS(imas_q,i,nscat,n)=var_n1
!!$              fap(i,n)=var_n1/XS(imat_q,i,nscat,n)
!!$!            else
!!$!              fap(i,n)=1.0
!!$            endif
!!$          enddo
!!$!CDIR NODEP
!!$          do in=1,gs%n_bin*L
!!$            n=(in-1)/gs%N_BIN+1
!!$            i=in-(n-1)*gs%N_BIN
!!$            if(icond1(in)/=1.and.icond2(in)/=1) then
!!$              if(fap(i,n)<min_faps_s) then
!!$                XS(imas_q,i,nscat,n)=min(max(m_lmt_ap/ag%tv(n)%den,&
!!$                           min_faps_s*XS(imat_q,i,nscat,n)),XS(imat_q,i,nscat,n))
!!$              endif
!!$            endif
!!$          enddo
!!$        endif
!!$      endif
!!$
!!$!      totmixs=0.0_PS
!!$!      do n=1,L
!!$!          do i=1,gs%n_bin
!!$!             in=i+(n-1)*gs%n_bin
!!$!             if(icond1(in)==1.or.icond2(in)==1) cycle
!!$!
!!$!             totmixs(n)=totmixs(n)+max(0.0_PS,XS(imt_q,i,nscat,n)-XS(imat_q,i,nscat,n))
!!$!          end do
!!$!
!!$!      enddo
!!$      totmixs=0.0_PS
!!$      do i=1,gs%n_bin
!!$        do n=1,L
!!$          if(icond1(in)/=1.and.icond2(in)/=1) then
!!$            totmixs(n)=totmixs(n)+max(0.0_PS,XS(imt_q,i,nscat,n)-XS(imat_q,i,nscat,n))
!!$          endif
!!$        enddo
!!$      end do
!!$
!!$    endif
!!$
!!$    ! +++ update vapor field +++
!!$    if( out_type == 2 ) then
!!$      icond1(1:L)=0
!!$      icond2(1:L)=0
!!$      do n=1,L
!!$        rv(n)=ag%tv(n)%rv+ag%tv(n)%dmassdt/ag%tv(n)%den*gr%dt
!!$        rv2(n)=qtp(n)-totmixr(n)-totmixs(n)
!!$        e_n=ag%TV(n)%P*rv(n)/(Rdvchiarui+rv(n))
!!$        s_n(n)=e_n/get_sat_vapor_pres_lk(1,ag%TV(n)%T_m,ag%estbar,ag%esitbar)-1.0_PS
!!$
!!$        if(s_n(n)>0.10_PS) then
!!$          icond1(n)=1
!!$        endif
!!$        if(rv(n)<0.0_PS.or.s_n(n)>0.50) then
!!$          if(rv(n)<0.0_PS) then
!!$            icond2(n)=2
!!$          else
!!$            icond2(n)=5
!!$          endif
!!$        endif
!!$      enddo
!!$      if(any(icond1(1:L)>0).or.any(icond2(1:L)>0)) then
!!$        do n=1,L
!!$          if(icond1(n)==1) then
!!$            write(*,*) "cal_model:k,i,j,p,qv,qr,qi,t,tn,sv1,svn1,sn,rv2",KD(n),ID(n),JD(n),ag%TV(n)%P,&
!!$               rv(n),totmixr(n),totmixs(n),ag%TV(n)%T,&
!!$               ag%TV(n)%T_n,ag%TV(n)%s_v(1),ag%TV(n)%s_v_n(1),s_n(n),rv2(n)
!!$          endif
!!$          if(icond2(n)==2) then
!!$            write(*,*) "rv is negative"
!!$
!!$            ! assume vapor mixing ratio to 10% of ice saturation.
!!$            rv_n=(-0.9_PS+1.0_PS)*ag%tv(n)%e_sat(2)*Rdvchiarui/&
!!$                      (ag%tv(n)%P-(-0.9+1.0_PS)*ag%tv(n)%e_sat(2))
!!$
!!$            write(*,*) "cal_model neg", KD(n),ID(n),JD(n),ag%TV(n)%P,&
!!$                  qtp(n),rv(n),rv_n,totmixr(n),totmixs(n),ag%TV(n)%T,&
!!$                  ag%TV(n)%T_n,ag%TV(n)%s_v(1),ag%TV(n)%s_v_n(1),s_n(n)
!!$
!!$            rv(n)=rv_n
!!$          elseif(icond2(n)==5) then
!!$            write(*,*) "rv(n) is super saturated"
!!$            write(*,*) "cal_model su", KD(n),ID(n),JD(n),ag%TV(n)%P,&
!!$                 qtp(n),rv(n),totmixr(n),totmixs(n),ag%TV(n)%T,&
!!$                 ag%TV(n)%T_n,ag%TV(n)%s_v(1),ag%TV(n)%s_v_n(1),s_n(n)
!!$          endif
!!$          if(icond2(n)>0) then
!!$            write(*,*) "mes_rc,density of air",mes_rc(n),ag%TV(n)%den
!!$            write(*,*) "rain"
!!$            write(*,*) (gr%MS(i,n)%mass(rmt),i=1,gr%N_BIN)
!!$            do i=1,gr%N_BIN
!!$              write(*,*) i,(gr%MS(i,n)%dmassdt(rmt,j),j=1,gr%N_tendpros)
!!$            end do
!!$
!!$            write(*,*) "ice"
!!$            write(*,*) (gs%MS(i,n)%mass(imt),i=1,gs%N_BIN)
!!$            do i=1,gs%N_BIN
!!$              write(*,*) i,(gs%MS(i,n)%dmassdt(imt,j),j=1,gs%N_tendpros)
!!$            end do
!!$            write(*,*) ag%TV(n)%s_v(1),ag%TV(n)%s_v(2),&
!!$                  ag%TV(n)%s_v_n(1),ag%TV(n)%s_v_n(2)
!!$
!!$          endif
!!$        enddo
!!$      endif
!!$
!!$    endif
!!$
!!$    !
!!$    ! Aerosol
!!$    !
!!$    if(flagp_a/=0) then
!!$
!!$      do ic=1,nacat
!!$        ! +++ for aerosols +++
!!$        if((flagp_a==-1.or.flagp_a==-4).and.ic==2) then
!!$          cycle
!!$        elseif((flagp_a==-2.or.flagp_a==-5).and.ic/=2) then
!!$          cycle
!!$        elseif(flagp_a==-3.or.flagp_a==-6) then
!!$          cycle
!!$        end if
!!$
!!$        do in=1,ga(ic)%n_bin*L
!!$          n=(in-1)/ga(ic)%N_BIN+1
!!$          i=in-(n-1)*ga(ic)%N_BIN
!!$          ndxdt(i,n) = 0.0_ps
!!$        enddo
!!$
!!$        do j=1,ga(ic)%n_tendpros
!!$!CDIR NODEP
!!$          do in=1,ga(ic)%n_bin*L
!!$            n=(in-1)/ga(ic)%N_BIN+1
!!$            i=in-(n-1)*ga(ic)%N_BIN
!!$
!!$            ndxdt(i,n)=ndxdt(i,n)+ga(ic)%MS(i,n)%dcondt(j)
!!$
!!$!!c           if((ga(ic)%MS(i,n)%dcondt(j)<-1.0e-5 .or.ga(ic)%MS(i,n)%dcondt(j)>1.0e-5) ) then
!!$!!c                   if(ic==2.and.ga(ic)%MS(i,n)%dcondt(j)>1.0e-4) then
!!$!!c                      write(*,*) "con_a:ic,process,grid,dcdt",ic,j,n,ga(ic)%MS(i,n)%dcondt(j)
!!$!!c                      write(*,*) "m1,m2,m3",ga(ic)%MS(i,n)%mass
!!$!!c                      write(*,*) "con", ga(ic)%MS(i,n)%con
!!$!!c                      write(*,*) "svw,i",ag%tv(n)%s_v(1),ag%tv(n)%s_v(2)
!!$!!c                      write(*,*) "dsvw,i",ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2)
!!$!!c!!c                      stop
!!$!!c                   end if
!!$          enddo
!!$        enddo
!!$        if( out_type == 2 ) then
!!$          do in=1,ga(ic)%n_bin*L
!!$            n=(in-1)/ga(ic)%N_BIN+1
!!$            i=in-(n-1)*ga(ic)%N_BIN
!!$            XA(acon_q,i,ic,n)=max(0.0_PS,&
!!$                     (ga(ic)%MS(i,n)%con+ndxdt(i,n)*ga(ic)%dt))/ag%tv(n)%den
!!$
!!$          enddo
!!$          ! check
!!$          icond1(1:ga(ic)%n_bin*L)=0
!!$          do in=1,ga(ic)%n_bin*L
!!$            n=(in-1)/ga(ic)%N_BIN+1
!!$            i=in-(n-1)*ga(ic)%N_BIN
!!$            if(XA(acon_q,i,ic,n)>1.0e+9) then
!!$              icond1(in)=1
!!$            endif
!!$          enddo
!!$          if(any(icond1(1:ga(ic)%n_bin*L)>0)) then
!!$            do in=1,ga(ic)%n_bin*L
!!$              n=(in-1)/ga(ic)%N_BIN+1
!!$              i=in-(n-1)*ga(ic)%N_BIN
!!$              if(icond1(in)==1) then
!!$                write(*,*) "cal_model_out:xa,con,connew",KD(n),ID(n),JD(n),ic,&
!!$                        XA(acon_q,i,ic,n),ga(ic)%MS(i,n)%con,XA(acon_q,i,ic,n)*ag%tv(n)%den
!!$                write(*,*) "qr >bin,grid,dcondt", i,ID(n),JD(n),KD(n),(gr%MS(j,n)%dcondt(1),j=1,gr%N_BIN)
!!$                write(*,*) "qr >bin,grid,dmassdt", i,ID(n),JD(n),KD(n),(gr%MS(j,n)%dmassdt(rmt,1),j=1,gr%N_BIN)
!!$                write(*,*) "qa >icat,dcondt",ic,(ga(ic)%ms(1,n)%dcondt(j),j=1,ga(ic)%n_tendpros)
!!$                write(*,*) "qa >icat,dmassdt",ic,(ga(ic)%MS(1,n)%dmassdt(amt,j),j=1,ga(ic)%n_tendpros)
!!$                write(*,*) "qi >bin,grid,dcondt",i,ID(n),JD(n),KD(n),(gs%MS(j,n)%dcondt(1),j=1,gs%N_BIN)
!!$                write(*,*) "qi >bin,grid,dmassdt",i,ID(n),JD(n),KD(n),(gs%MS(j,n)%dmassdt(imt,1),j=1,gs%N_BIN)
!!$              end if
!!$!!c             write(*,'("i,n,mten,cten",2i5,10es15.6)') i,n,ndxdt1,ndxdt,&
!!$!!c                  ga(ic)%MS(i,n)%dmassdt(amt,1),ga(ic)%MS(i,n)%dmassdt(amt,9),&
!!$!!c                  ga(ic)%MS(i,n)%dcondt(1),ga(ic)%MS(i,n)%dcondt(9)
!!$            enddo
!!$          endif
!!$        endif
!!$        if(level>=4) then
!!$
!!$          do in=1,ga(ic)%n_bin*L
!!$            n=(in-1)/ga(ic)%N_BIN+1
!!$            i=in-(n-1)*ga(ic)%N_BIN
!!$            ndxdt(i,n) = 0.0_ps
!!$            ndxdt2(i,n) = 0.0_ps
!!$          enddo
!!$
!!$          do j=1,ga(ic)%n_tendpros
!!$!CDIR NODEP
!!$            do in=1,ga(ic)%n_bin*L
!!$              n=(in-1)/ga(ic)%N_BIN+1
!!$              i=in-(n-1)*ga(ic)%N_BIN
!!$
!!$              ! total mass of aerosols (g/g)
!!$              ndxdt(i,n)=ndxdt(i,n)+ga(ic)%MS(i,n)%dmassdt(amt,j)
!!$
!!$              ! mass by soluble mass of aerosols (g/g)
!!$              ndxdt2(i,n)=ndxdt2(i,n)+ga(ic)%MS(i,n)%dmassdt(ams,j)
!!$            enddo
!!$          enddo
!!$
!!$          ierror(1:ga(ic)%n_bin*ga(ic)%n_tendpros*L)=0
!!$          do ijn=1,ga(ic)%n_bin*ga(ic)%n_tendpros*L
!!$            in=(ijn-1)/ga(ic)%n_tendpros+1
!!$            j=ijn-(in-1)*ga(ic)%n_tendpros
!!$            n=(in-1)/ga(ic)%n_bin+1
!!$            i=in-(n-1)*ga(ic)%n_bin
!!$            if(ga(ic)%MS(i,n)%dmassdt(amt,j)>1.0e-7.or.&
!!$               ga(ic)%MS(i,n)%dmassdt(amt,j)<-1.0e-7)then
!!$              ierror(ijn)=1
!!$            endif
!!$          enddo
!!$          if(any(ierror(1:ga(ic)%n_bin*ga(ic)%n_tendpros*L)>0)) then
!!$            do ijn=1,ga(ic)%n_bin*ga(ic)%n_tendpros*L
!!$              in=(ijn-1)/ga(ic)%n_tendpros+1
!!$              j=ijn-(in-1)*ga(ic)%n_tendpros
!!$              n=(in-1)/ga(ic)%n_bin+1
!!$              i=in-(n-1)*ga(ic)%n_bin
!!$              if(ierror(ijn)==1) then
!!$                write(*,*) "Warning aptm in ap>ga(ic)id,process,grid,id,jd,kd,dmassdt"
!!$                write(*,*) ic,j,n,ID(n),JD(n),KD(n),ga(ic)%MS(i,n)%dmassdt(amt,j)
!!$                write(*,*) "m1,m2,m3",ga(ic)%MS(i,n)%mass
!!$                write(*,*) "con", ga(ic)%MS(i,n)%con
!!$                write(*,*) "svw,i",ag%tv(n)%s_v(1),ag%tv(n)%s_v(2)
!!$                write(*,*) "dsvw,i",ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2)
!!$              endif
!!$            enddo
!!$          endif
!!$
!!$          ierror(1:ga(ic)%n_bin*ga(ic)%n_tendpros*L)=0
!!$          do ijn=1,ga(ic)%n_bin*ga(ic)%n_tendpros*L
!!$            in=(ijn-1)/ga(ic)%n_tendpros+1
!!$            j=ijn-(in-1)*ga(ic)%n_tendpros
!!$            n=(in-1)/ga(ic)%n_bin+1
!!$            i=in-(n-1)*ga(ic)%n_bin
!!$            if(ga(ic)%MS(i,n)%dmassdt(ams,j)>1.0e+08.or.&
!!$               ga(ic)%MS(i,n)%dmassdt(ams,j)<-1.0e+08)then
!!$              ierror(ijn)=1
!!$            endif
!!$          enddo
!!$          if(any(ierror(1:ga(ic)%n_bin*ga(ic)%n_tendpros*L)>0)) then
!!$            do ijn=1,ga(ic)%n_bin*ga(ic)%n_tendpros*L
!!$              in=(ijn-1)/ga(ic)%n_tendpros+1
!!$              j=ijn-(in-1)*ga(ic)%n_tendpros
!!$              n=(in-1)/ga(ic)%n_bin+1
!!$              i=in-(n-1)*ga(ic)%n_bin
!!$              if(ierror(ijn)==1) then
!!$                  write(*,*) "Warning apsm in ap>ga(ic)id,process,dmassdt"
!!$                  write(*,*) ic,j,ga(ic)%MS(i,n)%dmassdt(ams,j)
!!$                  write(*,*) "m1,m2,m3",ga(ic)%MS(i,n)%mass
!!$                  write(*,*) "con", ga(ic)%MS(i,n)%con
!!$                  write(*,*) "svw,i",ag%tv(n)%s_v(1),ag%tv(n)%s_v(2)
!!$                  write(*,*) "nmass",ga(ic)%n_masscom
!!$              endif
!!$            enddo
!!$          endif
!!$
!!$          if( out_type == 2 ) then
!!$            do in=1,ga(ic)%n_bin*L
!!$              n=(in-1)/ga(ic)%N_BIN+1
!!$              i=in-(n-1)*ga(ic)%N_BIN
!!$              XA(amt_q,i,ic,n)=max((ga(ic)%MS(i,n)%mass(amt)+&
!!$                           ndxdt(i,n)*ga(ic)%dt)/ag%tv(n)%den,0.0_ps)
!!$            enddo
!!$
!!$!CDIR NODEP
!!$            do in=1,ga(ic)%n_bin*L
!!$              n=(in-1)/ga(ic)%N_BIN+1
!!$              i=in-(n-1)*ga(ic)%N_BIN
!!$              XA(ams_q,i,ic,n)=min(max((ga(ic)%MS(i,n)%mass(ams)+&
!!$                           ndxdt2(i,n)*ga(ic)%dt)/ag%tv(n)%den,0.0_ps),XA(amt_q,i,ic,n))
!!$
!!$            enddo
!!$
!!$            if(ic==1) then
!!$              !
!!$              ! for CCN accumulation mode (category 1)
!!$              !
!!$!!c                         if(XA(amt_q,i,ic,n)*sep_faps>XA(ams_q,i,ic,n)) then
!!$!!c                            write(*,*) "sol mass is less than limt: ic,n",ic,n
!!$!!c                            write(*,*) XA(amt_q,i,ic,n),XA(acon_q,i,ic,n),XA(ams_q,i,ic,n)
!!$!!c                            write(*,'("dmassdt 1",15ES15.6)') (ga(ic)%MS(i,n)%dmassdt(amt,j),j=1,ga(ic)%N_tendpros)
!!$!!c                            write(*,'("dmassdt 2",15ES15.6)') (ga(ic)%MS(i,n)%dmassdt(ams,j),j=1,ga(ic)%N_tendpros)
!!$!!c                         end if
!!$!CDIR NODEP
!!$              do in=1,ga(ic)%n_bin*L
!!$                n=(in-1)/ga(ic)%N_BIN+1
!!$                i=in-(n-1)*ga(ic)%N_BIN
!!$                XA(ams_q,i,ic,n)=min(max(XA(amt_q,i,ic,n)*sep_faps,&
!!$                              XA(ams_q,i,ic,n)),XA(amt_q,i,ic,n))
!!$              enddo
!!$            elseif(ic==2) then
!!$              !
!!$              ! for IN accumulation mode (category 2)
!!$              !
!!$!CDIR NODEP
!!$              do in=1,ga(ic)%n_bin*L
!!$                n=(in-1)/ga(ic)%N_BIN+1
!!$                i=in-(n-1)*ga(ic)%N_BIN
!!$                XA(ams_q,i,ic,n)=max(0.0_PS,min(XA(amt_q,i,ic,n)*0.99_PS*sep_faps,&
!!$                              XA(ams_q,i,ic,n)))
!!$              enddo
!!$            elseif(ic==3) then
!!$              !
!!$              ! for CCN nucleation mode (category 3)
!!$              !
!!$!CDIR NODEP
!!$              do in=1,ga(ic)%n_bin*L
!!$                n=(in-1)/ga(ic)%N_BIN+1
!!$                i=in-(n-1)*ga(ic)%N_BIN
!!$                XA(ams_q,i,ic,n)=min(max(XA(amt_q,i,ic,n)*sep_faps,&
!!$                              XA(ams_q,i,ic,n)),XA(amt_q,i,ic,n))
!!$              enddo
!!$            elseif(ic==4) then
!!$              !
!!$              ! for CCN nucleation mode (category 4)
!!$              !
!!$!CDIR NODEP
!!$              do in=1,ga(ic)%n_bin*L
!!$                n=(in-1)/ga(ic)%N_BIN+1
!!$                i=in-(n-1)*ga(ic)%N_BIN
!!$                XA(ams_q,i,ic,n)=min(max(XA(amt_q,i,ic,n)*sep_faps,&
!!$                              XA(ams_q,i,ic,n)),XA(amt_q,i,ic,n))
!!$              enddo
!!$            else
!!$              write(*,*) "cal_model_tend_all>this ap cat not defined"
!!$              stop
!!$            end if
!!$          endif
!!$
!!$          if(ic==1.or.ic==2.or.ic==4) then
!!$            m_lmt=coef4pi3*r3_lmt
!!$          elseif(ic==3) then
!!$            m_lmt=coef4pi3*r3_lmt_nuc
!!$          else
!!$            write(*,*) "this category not defined in ini_group_all"
!!$            stop
!!$          endif
!!$
!!$!CDIR NODEP
!!$          do in=1,ga(ic)%n_bin*L
!!$            n=(in-1)/ga(ic)%N_BIN+1
!!$            i=in-(n-1)*ga(ic)%N_BIN
!!$
!!$            if(XA(acon_q,i,ic,n)<n_lmt_ap/ag%tv(n)%den) then
!!$
!!$!!c                      write(*,'("ap fixed ic",i5,3es15.6)')&
!!$!!c                           ic,XA(acon_q,i,ic,n)*ag%tv(n)%den,n_lmt_ap
!!$!!c                      rmod1=XA(amt_q,i,ic,n)/XA(acon_q,i,ic,n)
!!$!!c                      rmod2=XA(ams_q,i,ic,n)/XA(acon_q,i,ic,n)
!!$              var_n1=n_lmt_ap/ag%tv(n)%den
!!$              var_n2=var_n1*m_lmt*ga(ic)%MS(i,n)%den
!!$              XA(acon_q,i,ic,n)=var_n1
!!$              XA(amt_q,i,ic,n)=var_n2
!!$              XA(ams_q,i,ic,n)=var_n2*eps_ap0(ic)
!!$!!c                      XA(amt_q,i,ic,n)=XA(acon_q,i,ic,n)*rmod1
!!$!!c                      XA(ams_q,i,ic,n)=XA(acon_q,i,ic,n)*rmod2
!!$            end if
!!$          enddo
!!$
!!$!CDIR NODEP
!!$          do in=1,ga(ic)%n_bin*L
!!$            n=(in-1)/ga(ic)%N_BIN+1
!!$            i=in-(n-1)*ga(ic)%N_BIN
!!$            if(XA(acon_q,i,ic,n)>n_max_ap/ag%tv(n)%den) then
!!$!!c                      write(*,'("ap con reached max. so fixed.",i5,3es15.6)')&
!!$!!c                           ic,XA(acon_q,i,ic,n)*ag%tv(n)%den,n_max_ap
!!$              XA(acon_q,i,ic,n)=n_max_ap/ag%tv(n)%den
!!$            end if
!!$          enddo
!!$
!!$!CDIR NODEP
!!$          do in=1,ga(ic)%n_bin*L
!!$            n=(in-1)/ga(ic)%N_BIN+1
!!$            i=in-(n-1)*ga(ic)%N_BIN
!!$            if(XA(amt_q,i,ic,n)<m_lmt*ga(ic)%MS(i,n)%den*XA(acon_q,i,ic,n)) then
!!$!!c                      write(*,'("ap reached min:cat,xa1,min,eps0",i5,10es15.6)')&
!!$!!c                           ic,XA(amt_q,i,ic,n),m_lmt*ga(ic)%MS(i,n)%den*XA(acon_q,i,ic,n),eps_ap0(ic)
!!$              var_n1=XA(acon_q,i,ic,n)*m_lmt*ga(ic)%MS(i,n)%den
!!$              XA(amt_q,i,ic,n)=var_n1
!!$              XA(ams_q,i,ic,n)=var_n1*eps_ap0(ic)
!!$            end if
!!$          enddo
!!$
!!$!CDIR NODEP
!!$          do in=1,ga(ic)%n_bin*L
!!$            n=(in-1)/ga(ic)%N_BIN+1
!!$            i=in-(n-1)*ga(ic)%N_BIN
!!$            XA(ams_q,i,ic,n)=min(XA(ams_q,i,ic,n),XA(amt_q,i,ic,n))
!!$          enddo
!!$!!c                      if(XA(ams_q,i,ic,n)>0.0_ps.and.XA(amt_q,i,ic,n)>0.0_ps) then
!!$!!c                         write(*,'("ap fixed2, ic,mmass,min_mass",i5,3es15.6)')&
!!$!!c                              ic,XA(amt_q,i,ic,n)/XA(acon_q,i,ic,n),coef4pi3*r3_lmt*ga(ic)%MS(i,n)%den
!!$!!c                         rmod1=XA(ams_q,i,ic,n)/XA(amt_q,i,ic,n)
!!$!!c                         XA(amt_q,i,ic,n)=XA(acon_q,i,ic,n)*coef4pi3*r3_lmt*ga(ic)%MS(i,n)%den
!!$!!c                         XA(ams_q,i,ic,n)=XA(amt_q,i,ic,n)*rmod1
!!$!!c                      else
!!$!!c                         write(*,'("ap fixed3, ic,mmass,min_mass",i5,3es15.6)')&
!!$!!c                              ic,XA(amt_q,i,ic,n)/XA(acon_q,i,ic,n),coef4pi3*r3_lmt*ga(ic)%MS(i,n)%den
!!$!!c                         XA(amt_q,i,ic,n)=XA(acon_q,i,ic,n)*coef4pi3*r3_lmt*ga(ic)%MS(i,n)%den
!!$!!c                         XA(ams_q,i,ic,n)=XA(amt_q,i,ic,n)*ga(ic)%MS(i,n)%eps_map
!!$!!c                      end if
!!$!!c                   end if
!!$!!c                    if(XA(amt_q,i,ic,n)>1.0e-7) then
!!$!!c                       write(*,'("ap mass is large: ic,n,id,jd,kd",6I5)') ic,n,ID(n),JD(n),KD(n)
!!$!!c                       write(*,'("xa1,xa2,xa3",5ES15.6)') XA(amt_q,i,ic,n)*ag%TV(n)%den,XA(acon_q,i,ic,n)*ag%TV(n)%den&
!!$!!c                                                     ,XA(ams_q,i,ic,n)*ag%TV(n)%den
!!$!!c                       write(*,'("density,con,mass_t,mass_sol",6ES15.6)') ga(ic)%MS(i,n)%den&
!!$!!c                                      ,ga(ic)%MS(i,n)%con,ga(ic)%MS(i,n)%mass(amt),ga(ic)%MS(i,n)%mass(ams)
!!$!!c                       write(*,'("dcondt 1",15ES15.6)') (ga(ic)%MS(i,n)%dcondt(j),j=1,ga(ic)%N_tendpros)
!!$!!c                       write(*,'("dmassdt 1",15ES15.6)') (ga(ic)%MS(i,n)%dmassdt(amt,j),j=1,ga(ic)%N_tendpros)
!!$!!c                       write(*,'("dmassdt 2",15EES15.6)') (ga(ic)%MS(i,n)%dmassdt(ams,j),j=1,ga(ic)%N_tendpros)
!!$!!c                    end if
!!$        endif
!!$      enddo
!!$    endif
!!$
!!$!    deallocate(oq, &
!!$!             icond1,icond2,ierror,&
!!$!             ndxdt2,ndxdt)
!!$
!!$
!!$!!c    call qiprint2(xs,nstype,nsbin,nscat,l,"af mten")
!!$!!c    call qiprint3(gs,l,"af mten")
!!$  end subroutine cal_model_tendency_all

!!$  subroutine cal_model_tendency_all_scl(level,out_type,ag,mes_rc,&
!!$       gr,gs,ga,&
!!$       l,nrtype,nrbin,nrcat,nstype,nsbin,nscat,natype,nabin,nacat,&
!!$       qtp,rv,xc,xr,xs,xa,flagp_r,flagp_s,flagp_a,eps_ap0,ID,JD,KD)
!!$!!c  subroutine cal_model_tendency( level, out_type, g, th_var, mes_rc,x, dxdt)
!!$    integer,intent(in)  :: l,nrtype,nrbin,nrcat,nstype,nsbin,nscat,&
!!$         natype,nabin,nacat
!!$    integer :: ID(*),JD(*),KD(*)
!!$    type (group), intent(in)  :: gr,gs
!!$    type (group), dimension(nacat) :: ga
!!$    type (airgroup), intent(in) :: ag
!!$    ! message from reality-check
!!$!tmp    integer,pointer,dimension(:)  :: mes_rc
!!$    integer,dimension(*)   :: mes_rc
!!$    real(PS) :: qtp(*),rv(*),xc(*),XR(nrtype,nrbin,nrcat,*),XS(nstype,nsbin,nscat,*),&
!!$            XA(natype,nabin,nacat,*)
!!$    real(PS),dimension(*)  :: eps_ap0
!!$    integer,intent(in)  :: flagp_r,flagp_s,flagp_a
!!$    real (ps)     :: ndxdt,ndxdt1
!!$    ! level of complexity
!!$    integer, intent(in)  :: level
!!$    integer  :: out_type
!!$    integer  :: i,j,k,n,ic!,var_status
!!$    ! original total volume variables
!!$    real(ps),dimension(gs%N_nonmass)  :: oQ
!!$    !         1. volume of circumscribing sphere (cm^3/g)
!!$    !         2. con. weighted a-axis length^3 (cm^3/g)
!!$    !         3. con. weighted c-axis length^3 (cm^3/g)
!!$    !         4. con. weighted d-axis length^3 (cm^3/g)
!!$    !         5. con. weighted r-axis length^3 (cm^3/g) (rosetta bulletes)
!!$    !         6. con. weighted e-axis length^3 (cm^3/g) (polycrystals)
!!$
!!$    integer,dimension(gs%N_nonmass)  :: ick
!!$
!!$!!c    ! mass fraction of ice crystal to total mass
!!$!!c    real(ps) :: fic,fic_s
!!$!!c    ! minimum possible fractions
!!$!!c    real(ps),parameter :: min_fic_s=
!!$
!!$    ! mass fraction of aerosol particles to total mass of hydrometeors and
!!$    ! mass fraction of soluble aerosol particles to total aerosol mass
!!$    ! those are used to maintain positive aerosol masses, and use fraction
!!$    ! of smaller bin.
!!$    real(ps) :: fap
!!$    real(ps) :: fapt_r,faps_r,fapt_s,faps_s
!!$    ! minimum possible fractions
!!$    real(ps),parameter :: min_fapt_r=1.0e-18_ps,min_faps_r=1.0e-5_ps&
!!$!!c    real(ps),parameter :: min_fapt_r=1.0e-18_ps,min_faps_r=0.0_ps&
!!$         ,min_fapt_s=1.0e-18_ps,min_faps_s=0.0_ps,sep_faps=1.0e-5_PS
!!$    ! possible ratio of mean mass to bin boundaries
!!$    real(ps),parameter :: brat1=1.001,brat2=0.999
!!$    real(ps),parameter :: mlmt=1.0e-30,nlmt=1.0e-30,m_lmt_ap=1.0e-25
!!$    !     the minimum possible background concentration for large particle
!!$    !     is assumed to be 1.0e-5 cm^-3
!!$    real(ps),parameter :: n_lmt_ap=1.0e-5,n_max_ap=1.0e+4
!!$!parcel model    real(ps),parameter :: n_lmt_ap=1.0e-15,n_max_ap=1.0e+4
!!$    !     The minimum radius possible for accumulation particles
!!$    real(ps),parameter :: r3_lmt=1.0e-18
!!$    !     The minimum radius possible for nucleation particles
!!$    real(PS),parameter :: r3_lmt_nuc=1.0e-21
!!$
!!$
!!$    real(ps),parameter :: mx_aprat=0.99999
!!$    real(ps) :: m_lmt
!!$    ! realistic bounds for length predictions
!!$    real(ps),parameter :: min_hexlen3=1.0e-12,max_hexlen3=27.0,&
!!$                          max_roslen3=27.0,max_irrlen3=27.0,max_exice=1.0
!!$    ! realistic bounds for crystal mass
!!$    real(ps),parameter :: m_icmin=4.763209003e-12
!!$
!!$    ! accuracy limit to prevent dubious crystal features by numerical errors
!!$    real(PS),parameter :: accuracy_lmt=1.0e-6
!!$
!!$    ! vapor specific humidity cacluated with activation scheme
!!$    real(PS) :: rv_n
!!$    ! total mixing ratio of hydrometeors
!!$    real(PS) :: totmixr,totmixs
!!$    real(PS) :: s_n,e_n,rv2
!!$
!!$    integer :: ie
!!$    !integer :: ic1,ic2
!!$    integer :: iq1
!!$
!!$    !integer :: imark
!!$
!!$!!c    real(PS) :: alen_n,clen_n,den_n
!!$
!!$
!!$    do n=1,L
!!$
!!$       ! +++ loop over grids +++
!!$       if( mes_rc(n) == 0 ) cycle
!!$
!!$!!c       do i=1,gs%n_bin
!!$!!c          write(*,'("XS out 0>",4I4,20ES14.6)') KD(n),ID(n),JD(n),i,XS(1:mxnpi,i,NSCAT,n)
!!$!!c       enddo
!!$
!!$
!!$       totmixr=0.0_PS
!!$       totmixs=0.0_PS
!!$
!!$       ! initialize the mass fractions
!!$       fapt_r=1.0_ps
!!$       faps_r=1.0_ps
!!$       fapt_s=1.0_ps
!!$       faps_s=0.0_ps
!!$
!!$       ! +++ for rain +++
!!$       if( flagp_r > 0 ) then
!!$          do i=1,gr%n_bin
!!$             ndxdt = 0.0_ps
!!$             do j = 1, gr%n_tendpros
!!$                if(gr%MS(i,n)%dmassdt(rmt,j)>1.0e+03.or.&
!!$                     gr%MS(i,n)%dmassdt(rmt,j)<-1.0e+03)then
!!$                   write(*,201) n,i,j,gr%MS(i,n)%dmassdt(rmt,j)
!!$201                format("Warning mt_rain>grid,bin,process,dmassdt",3i5,es15.6)
!!$!!c                   stop
!!$                end if
!!$
!!$                ndxdt = ndxdt+gr%MS(i,n)%dmassdt(rmt,j)
!!$             end do
!!$             if( out_type == 2 ) then
!!$                XR(rmt_q,i,nrcat,n)=(gr%MS(i,n)%mass(rmt)+ndxdt*gr%dt)/ag%tv(n)%den
!!$                ! check
!!$!!c                if(XR(rmt_q,i,nrcat,n)<=mlmt/ag%tv(n)%den) then
!!$                if(XR(rmt_q,i,nrcat,n)*ag%TV(n)%den<gr%MS(i,n)%mass(rmt)*accuracy_lmt.or.&
!!$                     XR(rmt_q,i,nrcat,n)<1.0e-30) then
!!$                   do j=1,nrtype
!!$                      XR(j,i,nrcat,n)=0.0_ps
!!$                   end do
!!$                   cycle
!!$                else if(XR(rmt_q,i,nrcat,n)>1.0e-1) then
!!$                   write(*,251) i,n,XR(rmt_q,i,nrcat,n)
!!$251                format("warning:qr is large at bin,grid:",2i5,2es15.6)
!!$                end if
!!$             end if
!!$!!c          dxdt(1,i,n)=dxdt(1,i,n)+ndxdt/ag%tv(n)%den
!!$             ndxdt1=ndxdt
!!$
!!$             ndxdt = 0.0_ps
!!$             do j = 1, gr%n_tendpros
!!$                if(gr%MS(i,n)%dcondt(j)<-1.0e+8 .or.gr%MS(i,n)%dcondt(j)>1.0e+8 ) then
!!$                   write(*,*) "Warning con_r:bin,grid,process",i,n,j,gr%MS(i,n)%dcondt(j)
!!$!!c                   stop
!!$                end if
!!$                ndxdt = ndxdt+gr%MS(i,n)%dcondt(j)
!!$             end do
!!$             if( out_type == 2 ) then
!!$                XR(rcon_q,i,nrcat,n)=(gr%MS(i,n)%con+ndxdt*gr%dt)/ag%tv(n)%den
!!$!!c                if(XR(rcon_q,i,nrcat,n)>nlmt/ag%tv(n)%den) then
!!$                if(XR(rcon_q,i,nrcat,n)*ag%TV(n)%den<gr%MS(i,n)%con*accuracy_lmt.or.&
!!$                     XR(rcon_q,i,nrcat,n)<1.0e-30) then
!!$                   do j=1,nrtype
!!$                      XR(j,i,nrcat,n)=0.0_ps
!!$                   end do
!!$                   cycle
!!$                else
!!$                   XR(rcon_q,i,nrcat,n)=&
!!$                        XR(rmt_q,i,nrcat,n)&
!!$                        /max(brat1*gr%binb(i),min(XR(rmt_q,i,nrcat,n)/XR(rcon_q,i,nrcat,n)&
!!$                        ,brat2*gr%binb(i+1)))
!!$                   ie=0
!!$!!c                else
!!$!!c                   XR(rcon_q,i,nrcat,n)=&
!!$!!c                        XR(rmt_q,i,nrcat,n)/(0.5_ps*(gr%binb(i)+gr%binb(i+1)))
!!$!!c                 ie=1
!!$                end if
!!$!!c                if(XR(rcon_q,i,nrcat,n)>1.0e+6) then
!!$!!c                   write(*,'("cal_model_out,ie,xr2,con",I3,3ES15.6)') ie,XR(rcon_q,i,nrcat,n),gr%MS(i,n)%con
!!$!!c                   write(*,771) i,ID(n),JD(n),KD(n),(gr%MS(i,n)%dcondt(j),j=1,gr%n_tendpros)
!!$!!c771                   format("qr >bin,grid,dcondt",4i5,30es15.6)
!!$!!c                   write(*,772) i,ID(n),JD(n),KD(n),(gr%MS(i,n)%dmassdt(rmt,j),j=1,gr%n_tendpros)
!!$!!c772                   format("qr >bin,grid,dmassdt",4i5,30es15.6)
!!$!!c                   write(*,*) "con of CCN,IN",ga(1)%MS(1,n)%con,ga(2)%MS(1,n)%con
!!$!!c                   write(*,773) 1,(ga(1)%ms(1,n)%dcondt(j),j=1,ga(1)%n_tendpros)
!!$!!c773                   format("qa >icat,dcondt",i5,30es15.6)
!!$!!c                   write(*,774) 1,(ga(1)%ms(1,n)%dmassdt(amt,j),j=1,ga(1)%n_tendpros)
!!$!!c774                   format("qa >icat,dmassdt",i5,30es15.6)
!!$
!!$!!c                      stop
!!$!!c                end if
!!$
!!$             end if
!!$
!!$!!c             write(*,'("cal_model_tendency:bin,i,j,k,mten,cten",4i5,10es15.6)') i,ID(n),JD(n),KD(n) &
!!$!!c                 ,ndxdt1,ndxdt &
!!$!!c                 ,gr%MS(i,n)%dmassdt(rmt,1),gr%MS(i,n)%dmassdt(rmt,3) &
!!$!!c                 ,gr%MS(i,n)%dcondt(3),gr%MS(i,n)%dcondt(3)
!!$
!!$!!c          dxdt(2,i,n)=dxdt(2,i,n)+ndxdt/ag%tv(n)%den
!!$
!!$
!!$             if(level>=4) then
!!$                ! mass by total mass of aerosols (g/g)
!!$                ndxdt=0.0_ps
!!$                do j=1,gr%n_tendpros
!!$                   if(gr%MS(i,n)%dmassdt(rmat,j)>1.0e+08.or.&
!!$                        gr%MS(i,n)%dmassdt(rmat,j)<-1.0e+08)then
!!$                      write(*,236) i,n,j,gr%MS(i,n)%dmassdt(rmat,j)
!!$236                   format("Warning aptm in rain>bin,grid,process,dmassdt",3i5,es15.6)
!!$!!c                      stop
!!$                   end if
!!$
!!$                   ndxdt=ndxdt+gr%MS(i,n)%dmassdt(rmat,j)
!!$                end do
!!$                if( out_type == 2 ) then
!!$                   ndxdt=(gr%MS(i,n)%mass(rmat)+ndxdt*gr%dt)/ag%tv(n)%den
!!$                   fap=ndxdt/XR(rmt_q,i,nrcat,n)
!!$                   if(fap>1.0_PS) then
!!$!!c                      write(*,*) "total mass is less than ap total",XR(rmt_q,i,nrcat,n),XR(rmat_q,i,nrcat,n)
!!$                      XR(rmat_q,i,nrcat,n)=ndxdt
!!$                      XR(rmt_q,i,nrcat,n)=XR(rmat_q,i,nrcat,n)
!!$                   elseif(fap>min_fapt_r) then
!!$                      XR(rmat_q,i,nrcat,n)=ndxdt
!!$                      fapt_r=fap
!!$                   else
!!$!!c                      XR(rmat_q,i,nrcat,n)=fapt_r*XR(rmt_q,i,nrcat,n)
!!$                      XR(rmat_q,i,nrcat,n)=min(max(m_lmt_ap/ag%tv(n)%den,min_fapt_r*XR(rmt_q,i,nrcat,n))&
!!$                                         ,mx_aprat*XR(rmt_q,i,nrcat,n))
!!$                   end if
!!$                end if
!!$!!c                if(gr%MS(i,n)%mass(rmat)+ndxdt*gr%dt<=0.0_ps) then
!!$!!c                   write(*,*) "here"
!!$!!c                   do j=1,gr%n_tendpros
!!$!!c                      write(*,236) i,n,j,gr%MS(i,n)%dmassdt(rmat,j)
!!$!!c                   end do
!!$!!c                end if
!!$                ! mass by soluble mass of aerosols (g/g)
!!$                ndxdt=0.0_ps
!!$                do j=1,gr%n_tendpros
!!$                   if(gr%MS(i,n)%dmassdt(rmas,j)>1.0e+08.or.&
!!$                        gr%MS(i,n)%dmassdt(rmas,j)<-1.0e+08)then
!!$                      write(*,237) i,n,j,gr%MS(i,n)%dmassdt(rmas,j)
!!$237                   format("Warning apsm in solid>bin,grid,process,dmassdt",3i5,es15.6)
!!$!!c                      stop
!!$                   end if
!!$
!!$                   ndxdt=ndxdt+gr%MS(i,n)%dmassdt(rmas,j)
!!$                end do
!!$                if( out_type == 2 ) then
!!$                   ndxdt=(gr%MS(i,n)%mass(rmas)+ndxdt*gr%dt)/ag%tv(n)%den
!!$                   fap=ndxdt/XR(rmat_q,i,nrcat,n)
!!$                   if(fap>=min_faps_r) then
!!$                      XR(rmas_q,i,nrcat,n)=ndxdt
!!$                      faps_r=fap
!!$                   else
!!$!!c                      XR(rmas_q,i,nrcat,n)=faps_r*XR(rmat_q,i,nrcat,n)
!!$                      XR(rmas_q,i,nrcat,n)=min(max(m_lmt_ap/ag%tv(n)%den,&
!!$                           min_faps_r*XR(rmat_q,i,nrcat,n)),XR(rmat_q,i,nrcat,n))
!!$!!c                      XR(rmas_q,i,nrcat,n)=min(eps_ap0(1)*XR(rmat_q,i,nrcat,n),XR(rmat_q,i,nrcat,n))
!!$                   end if
!!$                end if
!!$             end if
!!$
!!$             totmixr=totmixr+max(0.0_PS,XR(rmt_q,i,nrcat,n)-XR(rmat_q,i,nrcat,n))
!!$          end do
!!$
!!$       end if
!!$
!!$       if(debug) then
!!$!         if(sum(XR(rcon_q,1:gr%n_bin,nrcat,n))*ag%tv(n)%den>100.0) then
!!$         if(sum(XR(rcon_q,1:gr%n_bin,nrcat,n))*ag%tv(n)%den>0.0) then
!!$         write(*,'("cal_model_tendency:bin,i,j,k,total rcon, others",4i5,100es15.6)') n,ID(n),JD(n),KD(n) &
!!$                      ,sum(XR(rcon_q,1:gr%n_bin,nrcat,n))*ag%tv(n)%den &
!!$                      ,XR(rcon_q,1:gr%n_bin,nrcat,n)*ag%tv(n)%den
!!$
!!$         write(*,'("cal_model_tendency:bin,i,j,k,total act, total vap mass tendency",4i5,100es15.6)') n,ID(n),JD(n),KD(n) &
!!$                      ,sum(gr%MS(1:gr%n_bin,n)%dmassdt(rmt,3)) &  ! activation
!!$                      ,sum(gr%MS(1:gr%n_bin,n)%dmassdt(rmt,1))    ! vapor
!!$         write(*,*) "cal_model_tendency:con-bin,i,j,k,total rim, total melt, total mossop, total coal"
!!$         write(*,'(4I5,10ES15.6)') n,ID(n),JD(n),KD(n) &
!!$                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(4)) &  ! riming
!!$                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(10)) &  ! melting_shedding
!!$                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(9)) &  ! mossop and hallet
!!$                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(2))    ! collision-coal
!!$         endif
!!$       endif
!!$
!!$       ! for cloud droplets
!!$       XC(n)=XR(rmt_q,1,nrcat,n)
!!$
!!$       ! ice
!!$       if( flagp_s > 0 ) then
!!$          do i=1,gs%n_bin
!!$             ndxdt = 0.0_ps
!!$             do j = 1, gs%n_tendpros
!!$                if(gs%MS(i,n)%dmassdt(imt,j)>1.0e+08.or.&
!!$                     gs%MS(i,n)%dmassdt(imt,j)<-1.0e+08)then
!!$                   write(*,202) n,i,j,gs%MS(i,n)%dmassdt(imt,j)
!!$202                format("Warning mt_s>grid,bin,process,dmassdt",3i5,es15.6)
!!$!!c                   stop
!!$                end if
!!$                ndxdt = ndxdt+gs%MS(i,n)%dmassdt(imt,j)
!!$             end do
!!$             if( out_type == 2 ) then
!!$                XS(imt_q,i,nscat,n)=(gs%MS(i,n)%mass(imt)+ndxdt*gs%dt)/ag%tv(n)%den
!!$!!c                if(XS(imt_q,i,nscat,n)<=mlmt/ag%tv(n)%den) then
!!$                if(XS(imt_q,i,nscat,n)*ag%TV(n)%den<gs%MS(i,n)%mass(imt)*accuracy_lmt.or.&
!!$                     XS(imt_q,i,nscat,n)<1.0e-30) then
!!$                   do j=1,nstype
!!$                      XS(j,i,nscat,n)=0.0_ps
!!$                   end do
!!$                   cycle
!!$                else if(XS(imt_q,i,nscat,n)>1.0e-1) then
!!$                   write(*,252) i,n,XS(imt_q,i,nscat,n)
!!$252                format("warning:qs is large at bin,grid:",2i5,2es15.6)
!!$
!!$                end if
!!$             end if
!!$!!c          dxdt(1,i,n)=dxdt(1,i,n)+ndxdt/ag%tv(n)%den
!!$
!!$             ndxdt = 0.0_ps
!!$             do j = 1, gs%n_tendpros
!!$                if(gs%MS(i,n)%dcondt(j)<-1.0e+8 .or.gs%MS(i,n)%dcondt(j)>1.0e+8 ) then
!!$                   write(*,*) "Warning con_s:bin,grid,process",i,n,j,gs%MS(i,n)%dcondt(j)
!!$!!c                   stop
!!$                end if
!!$                ndxdt = ndxdt+gs%MS(i,n)%dcondt(j)
!!$             end do
!!$             if( out_type == 2 ) then
!!$                XS(icon_q,i,nscat,n)=(gs%MS(i,n)%con+ndxdt*gs%dt)/ag%tv(n)%den
!!$
!!$                if(XS(icon_q,i,nscat,n)*ag%TV(n)%den<gs%MS(i,n)%con*accuracy_lmt.or.&
!!$                     XS(icon_q,i,nscat,n)<1.0e-30) then
!!$                   do j=1,nstype
!!$                      XS(j,i,nscat,n)=0.0_ps
!!$                   end do
!!$                   cycle
!!$                else
!!$!!c                if(XS(icon_q,i,nscat,n)>nlmt/ag%tv(n)%den) then
!!$                   XS(icon_q,i,nscat,n)=&
!!$                        XS(imt_q,i,nscat,n)&
!!$                        /max(brat1*gs%binb(i),min(XS(imt_q,i,nscat,n)/XS(icon_q,i,nscat,n)&
!!$                        ,brat2*gs%binb(i+1)))
!!$!!c                else
!!$!!c                   XS(icon_q,i,nscat,n)=&
!!$!!c                        XS(imt_q,i,nscat,n)/(0.5_ps*(gs%binb(i)+gs%binb(i+1)))
!!$                end if
!!$             end if
!!$
!!$             iq1=ivcs_q
!!$             oq(ivcs) = gs%MS(i,n)%con*gs%IS(i,n)%v_cs
!!$             oq(iacr) = gs%MS(i,n)%con*gs%MS(i,n)%a_len**3.0
!!$             oq(iccr) = gs%MS(i,n)%con*gs%MS(i,n)%c_len**3.0
!!$             oq(idcr) = gs%MS(i,n)%con*gs%IS(i,n)%d**3.0
!!$             oq(iag) = gs%MS(i,n)%con*gs%IS(i,n)%ag**3.0
!!$             oq(icg) = gs%MS(i,n)%con*gs%IS(i,n)%cg**3.0
!!$             oq(inex) = gs%MS(i,n)%con*gs%IS(i,n)%n_exice
!!$
!!$
!!$             ick=0
!!$             ! volume of circumscribing sphere (cm^3/g) and
!!$             ! volume of axis length production (cm^3/g)
!!$             do k=1,gs%N_nonmass
!!$                ndxdt=0.0_ps
!!$                do j=1,gs%n_tendpros
!!$                   if(gs%MS(i,n)%dvoldt(k,j)>1.0e+20.or.&
!!$                        gs%MS(i,n)%dvoldt(k,j)<-1.0e+20)then
!!$                      write(*,231) n,k,j,gs%MS(i,n)%dvoldt(k,j)
!!$231                   format("axis>grid,axis#,process,dvoldt",3i5,es15.6)
!!$                      write(*,'("phi_cs,den",2es15.6)') gs%IS(i,n)%phi_cs,gs%MS(i,n)%den
!!$                      write(*,'("psi_ic,phi_ic,alen",3es15.6)') gs%IS(i,n)%psi_ic,gs%IS(i,n)%phi_ic,gs%MS(i,n)%a_len
!!$                      write(*,'("i,ht,sh",3i4)') i,gs%IS(i,n)%habit,gs%IS(i,n)%sh_type
!!$!!c                      stop
!!$                   end if
!!$                   ndxdt=ndxdt+gs%MS(i,n)%dvoldt(k,j)
!!$
!!$!!c                   if(ag%tv(n)%t<260.15.and.ag%tv(n)%t>250.15) then
!!$!!c                      if(k==5.and.gs%MS(i,n)%dvoldt(k,j)>1.0e-20) then
!!$!!c                         write(*,235) n,k,j,gs%MS(i,n)%dvoldt(k,j)
!!$!!c235                      format("bad ros len: axis>grid,axis#,process,dvoldt",3i5,es15.6)
!!$!!c                         write(*,'("dvoldt for all axis:",15es15.6)') (gs%MS(i,n)%dvoldt(k,j),k=1,1+gs%n_axis)
!!$!!c                         write(*,'("i,ht,sh",3i4)') i,gs%IS(i,n)%habit,gs%IS(i,n)%sh_type
!!$!!c                         write(*,'("a,c,d,e,r",5es15.6)') gs%MS(i,n)%a_len,gs%MS(i,n)%c_len,&
!!$!!c                              gs%IS(i,n)%d,gs%IS(i,n)%e,gs%IS(i,n)%r
!!$!!c                         if(i>=2) then
!!$!!c                            write(*,'("s: a,c,d,e,r",5es15.6)') gs%ms(i-1,n)%a_len,gs%ms(i-1,n)%c_len,&
!!$!!c                                 gs%is(i-1,n)%d,gs%is(i-1,n)%e,gs%is(i-1,n)%r
!!$!!c                         end if
!!$!!c                         stop
!!$!!c                      end if
!!$!!c                   end if
!!$                end do
!!$!!c                if(k==1) then
!!$!!c                   write(*,'("grid,axis,bin,dvoldt:",3i5,12es15.5)') n,k,i,(gs%MS(i,n)%dvoldt(k,j),j=1,gs%n_tendpros)
!!$!!c                end if
!!$
!!$                if( out_type == 2 ) then
!!$                   ndxdt=(oq(k)+ndxdt*gs%dt)/ag%tv(n)%den
!!$                   if(ndxdt>0.0_ps) then
!!$                      ick(k)=0
!!$                      XS(iq1-1+k,i,nscat,n)=ndxdt
!!$                   elseif(ndxdt.eq.0.0_ps) then
!!$                      ick(k)=1
!!$                      XS(iq1-1+k,i,nscat,n)=0.0_ps
!!$                   elseif(ndxdt.lt.0.0_ps) then
!!$!!c                      write(*,551) k,(oq(k)+ndxdt*gs%dt)/ag%tv(n)%den
!!$551                   format("new q is neg at axis",i5,es15.6)
!!$!!c                      write(*,552) i,n,(gs%MS(i,n)%dvoldt(k,j),j=1,gs%n_tendpros)
!!$552                   format("i,n,dvoldt",2i5,15es15.6)
!!$                      ick(k)=2
!!$                      XS(iq1-1+k,i,nscat,n)=0.0_ps
!!$                   end if
!!$
!!$!!c                   if((ick(2)>=1.or.ick(3)>=1).and.ick(k)==0) then
!!$!!c                      write(*,591) k,i,n,XS(iacr_q,i,nscat,n),XS(iccr_q,i,nscat,n),XS(2+k,i,nscat,n)
!!$!!c591                   format("a or c is zero (axis,bin,grid), alen,clen,xlen",3i5,3es15.6)
!!$!!c                      XS(2+k,i,nscat,n)=0.0_ps
!!$!!c                   end if
!!$                end if
!!$             end do
!!$
!!$             ! bound with maximum and minium possible length.
!!$             XS(iacr_q,i,nscat,n)=max(min_hexlen3*XS(icon_q,i,nscat,n),&
!!$                  min(max_hexlen3*XS(icon_q,i,nscat,n),XS(iacr_q,i,nscat,n)))
!!$             XS(iccr_q,i,nscat,n)=max(min_hexlen3*XS(icon_q,i,nscat,n),&
!!$                  min(max_hexlen3*XS(icon_q,i,nscat,n),XS(iccr_q,i,nscat,n)))
!!$             XS(idcr_q,i,nscat,n)=max(0.0_RP,min(XS(iacr_q,i,nscat,n),XS(idcr_q,i,nscat,n)))
!!$!!c             if((XS(idcr_q,i,nscat,n)/XS(icon_q,i,nscat,n))**(1.0/3.0)>1.0e-10.and.ag%TV(n)%T<253.16) then
!!$!!c                write(*,*) "something wrong here",i,n
!!$!!c                write(*,*) (gs%MS(i,n)%dvoldt(idcr,j),j=1,gs%N_tendpros)
!!$!!c             end if
!!$             XS(iag_q,i,nscat,n)=max(0.0_RP,min(max_roslen3*XS(icon_q,i,nscat,n),XS(iag_q,i,nscat,n)))
!!$             XS(icg_q,i,nscat,n)=max(0.0_RP,min(max_irrlen3*XS(icon_q,i,nscat,n),XS(icg_q,i,nscat,n)))
!!$
!!$!!c             k=gs%N_nonmass
!!$!!c             ndxdt=0.0_ps
!!$!!c             do j=1,gs%n_tendpros
!!$!!c                if(gs%MS(i,n)%dvoldt(k,j)>1.0e+20.or.&
!!$!!c                     gs%MS(i,n)%dvoldt(k,j)<-1.0e+20)then
!!$!!c                   write(*,231) n,k,j,gs%MS(i,n)%dvoldt(k,j)
!!$!!c                   write(*,'("phi_cs,den",2es15.6)') gs%IS(i,n)%phi_cs,gs%MS(i,n)%den
!!$!!c                   write(*,'("psi_ic,phi_ic,alen",3es15.6)') gs%IS(i,n)%psi_ic,gs%IS(i,n)%phi_ic,gs%MS(i,n)%a_len
!!$!!c                   write(*,'("i,ht,sh",3i4)') i,gs%IS(i,n)%habit,gs%IS(i,n)%sh_type
!!$!!c                   stop
!!$!!c                end if
!!$!!c
!!$!!c                ndxdt=ndxdt+gs%MS(i,n)%dvoldt(k,j)
!!$!!c
!!$!!c             end do
!!$!!c             XS(2+k,i,nscat,n)=(oq(k)+ndxdt*gs%dt)/ag%tv(n)%den
!!$!!c!!c
!!$!!c
!!$!!c             if(XS(imt_q,i,nscat,n)>0.0_PS.and.max_exice*XS(icon_q,i,nscat,n)*1.01<XS(inex_q,i,nscat,n)) then
!!$!!c               write(*,*) "XS2,XS9,nice",XS(icon_q,i,nscat,n),XS(inex_q,i,nscat,n),XS(inex_q,i,nscat,n)/XS(icon_q,i,nscat,n)
!!$!!c               write(*,'("dvoldt",20ES15.6)') (gs%MS(i,n)%dvoldt(7,j),j=1,gs%n_tendpros)
!!$!!c               write(*,'("dvdt/dcdt",20ES15.6)') (gs%MS(i,n)%dvoldt(7,j)/gs%MS(i,n)%dcondt(j),j=1,gs%n_tendpros)
!!$!!c               write(*,'("phi_cs,den",2es15.6)') gs%IS(i,n)%phi_cs,gs%MS(i,n)%den
!!$!!c               write(*,'("psi_ic,phi_ic,alen",3es15.6)') gs%IS(i,n)%psi_ic,gs%IS(i,n)%phi_ic,gs%MS(i,n)%a_len
!!$!!c               write(*,'("i,ht,sh,N_nonmass",5i4)') i,gs%IS(i,n)%habit,gs%IS(i,n)%sh_type,gs%N_nonmass
!!$!!c               write(*,'("ag,cg,nexice",4ES15.6)') gs%IS(i,n)%ag,gs%IS(i,n)%cg,gs%IS(i,n)%n_exice
!!$!!c               write(*,*) "con,new_con,rat",gs%MS(i,n)%con,XS(icon_q,i,nscat,n)*ag%TV(n)%den,XS(icon_q,i,nscat,n)*ag%TV(n)%den/gs%MS(i,n)%con
!!$!!c               write(*,*) "mass,new_mass,rat",gs%MS(i,n)%mass(imt),XS(imt_q,i,nscat,n)*ag%TV(n)%den,XS(imt_q,i,nscat,n)*ag%TV(n)%den/gs%MS(i,n)%mass(imt)
!!$!!c
!!$!!c             end if
!!$
!!$
!!$!!c                   if(ag%tv(n)%t<260.15.and.ag%tv(n)%t>250.15) then
!!$!!c                      if(k==5.and.gs%MS(i,n)%dvoldt(k,j)>1.0e-20) then
!!$!!c                         write(*,235) n,k,j,gs%MS(i,n)%dvoldt(k,j)
!!$!!c235                      format("bad ros len: axis>grid,axis#,process,dvoldt",3i5,es15.6)
!!$!!c                         write(*,'("dvoldt for all axis:",15es15.6)') (gs%MS(i,n)%dvoldt(k,j),k=1,1+gs%n_axis)
!!$!!c                         write(*,'("i,ht,sh",3i4)') i,gs%IS(i,n)%habit,gs%IS(i,n)%sh_type
!!$!!c                         write(*,'("a,c,d,e,r",5es15.6)') gs%MS(i,n)%a_len,gs%MS(i,n)%c_len,&
!!$!!c                              gs%IS(i,n)%d,gs%IS(i,n)%e,gs%IS(i,n)%r
!!$
!!$             XS(inex_q,i,nscat,n)=max(0.0_RP,min(max_exice*XS(icon_q,i,nscat,n),XS(inex_q,i,nscat,n)))
!!$
!!$
!!$             ! another constraints on d-axis length and a and c coordinates according to
!!$             ! the a and c axis lengths
!!$             if(XS(iacr_q,i,nscat,n)<=1.0e-9_PS*XS(icon_q,i,nscat,n).or.&
!!$                  XS(iccr_q,i,nscat,n)<=1.0e-12_PS*XS(icon_q,i,nscat,n)) then
!!$                XS(idcr_q,i,nscat,n)=0.0_PS
!!$             end if
!!$             if(XS(iacr_q,i,nscat,n)<=min_hexlen3*XS(icon_q,i,nscat,n).or.&
!!$                  XS(iccr_q,i,nscat,n)<=min_hexlen3*XS(icon_q,i,nscat,n)) then
!!$                XS(iag_q,i,nscat,n)=0.0_PS
!!$                XS(icg_q,i,nscat,n)=0.0_PS
!!$                XS(inex_q,i,nscat,n)=0.0_PS
!!$             end if
!!$
!!$!!c             call check_rosrel("calmodel",i,n,XS(imt_q,i,nscat,n)/XS(icon_q,i,nscat,n),&
!!$!!c                               (XS(iacr_q,i,nscat,n)/XS(icon_q,i,nscat,n))**(1.0/3.0),&
!!$!!c                               (XS(iccr_q,i,nscat,n)/XS(icon_q,i,nscat,n))**(1.0/3.0),&
!!$!!c                               XS(inex_q,i,nscat,n)/XS(icon_q,i,nscat,n),&
!!$!!c                               (XS(iag_q,i,nscat,n)/XS(icon_q,i,nscat,n))**(1.0/3.0),&
!!$!!c                               (XS(icg_q,i,nscat,n)/XS(icon_q,i,nscat,n))**(1.0/3.0),imark)
!!$
!!$!!c             if(imark==1) then
!!$!!c               write(*,'("i,n,alen0,clen0,nex0,ag0,cg0,mmassIN",2I5,10ES15.6)') i,n,gs%MS(i,n)%a_len,gs%MS(i,n)%c_len,&
!!$!!c                        gs%IS(i,n)%n_exice,gs%IS(i,n)%ag,gs%IS(i,n)%cg,ga(2)%MS(1,n)%mean_mass
!!$!!c               write(*,'("dmdt1",20ES15.6)') (gs%MS(i,n)%dmassdt(imt,j),j=1,gs%n_tendpros)
!!$!!c               write(*,'("dcdt1",20ES15.6)') (gs%MS(i,n)%dcondt(j),j=1,gs%n_tendpros)
!!$!!c               write(*,'("dvoldt1",20ES15.6)') (gs%MS(i,n)%dvoldt(ivcs,j),j=1,gs%n_tendpros)
!!$!!c               write(*,'("dvoldt2",20ES15.6)') (gs%MS(i,n)%dvoldt(iacr,j),j=1,gs%n_tendpros)
!!$!!c               write(*,'("dvoldt3",20ES15.6)') (gs%MS(i,n)%dvoldt(iccr,j),j=1,gs%n_tendpros)
!!$!!c             end if
!!$
!!$!!c             if(XS(imt_q,i,nscat,n)>0.0_PS.and.XS(inex_q,i,nscat,n)/XS(icon_q,i,nscat,n)<0.5.and.&
!!$!!c                ag%TV(n)%T<273.15-22) then
!!$!!c             if(XS(imt_q,i,nscat,n)>0.0_PS.and.(XS(iccr_q,i,nscat,n)/XS(iacr_q,i,nscat,n))**(1.0/3.0)>10.0) then
!!$!!c               write(*,'("Large axis ratio: i,n,ht,sh,N_nonmass",5i4)') i,n,gs%IS(i,n)%habit,gs%IS(i,n)%sh_type,gs%N_nonmass
!!$!!c               write(*,'("Single ice: i,n,ht,sh,N_nonmass",5i4)') i,n,gs%IS(i,n)%habit,gs%IS(i,n)%sh_type,gs%N_nonmass
!!$!!c               write(*,'("dvoldt2",20ES15.6)') (gs%MS(i,n)%dvoldt(iacr,j),j=1,gs%n_tendpros)
!!$!!c               write(*,'("dvoldt3",20ES15.6)') (gs%MS(i,n)%dvoldt(iccr,j),j=1,gs%n_tendpros)
!!$!!c               write(*,'("dmdt1",20ES15.6)') (gs%MS(i,n)%dmassdt(imt,j),j=1,gs%n_tendpros)
!!$!!c               write(*,'("n_mmass,n_a,n_c,n_phi,n_nexice",20ES15.6)') XS(imt_q,i,nscat,n)/XS(icon_q,i,nscat,n),&
!!$!!c                          (XS(iacr_q,i,nscat,n)/XS(icon_q,i,nscat,n))**(1.0/3.0),(XS(iccr_q,i,nscat,n)/XS(icon_q,i,nscat,n))**(1.0/3.0),&
!!$!!c                          (XS(iccr_q,i,nscat,n)/XS(iacr_q,i,nscat,n))**(1.0/3.0), XS(inex_q,i,nscat,n)/XS(icon_q,i,nscat,n)
!!$!!c               write(*,'("phi_cs,den",2es15.6)') gs%IS(i,n)%phi_cs,gs%MS(i,n)%den
!!$!!c               write(*,'("psi_ic,phi_ic,alen,clen",4es15.6)') gs%IS(i,n)%psi_ic,gs%IS(i,n)%phi_ic,gs%MS(i,n)%a_len,gs%MS(i,n)%c_len
!!$!!c               write(*,'("ag,cg,nexice",4ES15.6)') gs%IS(i,n)%ag,gs%IS(i,n)%cg,gs%IS(i,n)%n_exice
!!$!!c               write(*,*) "con,new_con,rat",gs%MS(i,n)%con,XS(icon_q,i,nscat,n)*ag%TV(n)%den,XS(icon_q,i,nscat,n)*ag%TV(n)%den/gs%MS(i,n)%con
!!$!!c               write(*,*) "mass,new_mass,rat",gs%MS(i,n)%mass(imt),XS(imt_q,i,nscat,n)*ag%TV(n)%den,XS(imt_q,i,nscat,n)*ag%TV(n)%den/gs%MS(i,n)%mass(imt)
!!$!!c               write(*,*) "T,svw,svi",ag%TV(n)%T,ag%TV(n)%s_v(1),ag%TV(n)%s_v(2)
!!$!!c               write(*,'("alen3,clen3",3es15.6)') gs%MS(i,n)%a_len**3.0,gs%MS(i,n)%c_len**3.0
!!$!!c               write(*,*) " "
!!$!!c
!!$!!c             end if
!!$
!!$!!c             if( (XS(idcr_q,i,nscat,n)/XS(iacr_q,i,nscat,n))**(1.0/3.0)>2.0/3.0.and.&
!!$!!c                 (XS(iccr_q,i,nscat,n)/XS(iacr_q,i,nscat,n))**(1.0/3.0)<1.0) then
!!$!!c
!!$!!c                write(*,'("3 dendrite forming at bin,grid,ID,JD,KD,T",5I5,ES15.6)') &
!!$!!c                      i,n,ID(n),JD(n),KD(n),ag%TV(n)%T-273.16
!!$!!c                write(*,*) "alen,clen,d",(XS(iacr_q,i,nscat,n)/XS(icon_q,i,nscat,n))**(1.0/3.0),(XS(iccr_q,i,nscat,n)/XS(icon_q,i,nscat,n))**(1.0/3.0),(XS(idcr_q,i,nscat,n)/XS(icon_q,i,nscat,n))**(1.0/3.0)
!!$!!c                write(*,*) "con,mass",XS(icon_q,i,nscat,n),XS(imt_q,i,nscat,n)
!!$!!c             end if
!!$
!!$             ! mass by ice crystals (g/g)
!!$             ndxdt=0.0_ps
!!$             do j=1,gs%n_tendpros
!!$                if(gs%MS(i,n)%dmassdt(imc,j)>1.0e+08.or.&
!!$                     gs%MS(i,n)%dmassdt(imc,j)<-1.0e+08)then
!!$                   write(*,233) n,j,gs%MS(i,n)%dmassdt(imc,j)
!!$233                format("Warning icemass>grid,process,dmassdt",2i5,es15.6)
!!$!!c                   stop
!!$                end if
!!$
!!$                ndxdt=ndxdt+gs%MS(i,n)%dmassdt(imc,j)
!!$             end do
!!$             if( out_type == 2 ) then
!!$                XS(imc_q,i,nscat,n)=min(XS(imt_q,i,NSCAT,n),&
!!$!tmp                     max((gs%MS(i,n)%mass(imc)+ndxdt*gs%dt)/ag%tv(n)%den,0.0_ps))
!!$                     max((gs%MS(i,n)%mass(imc)+ndxdt*gs%dt)/ag%tv(n)%den,m_icmin*XS(icon_q,i,NSCAT,n)))
!!$             end if
!!$
!!$
!!$             ! mass by riming process (g/g)
!!$             ndxdt=0.0_ps
!!$             do j=1,gs%n_tendpros
!!$                if(gs%MS(i,n)%dmassdt(imr,j)>1.0e+03.or.&
!!$                     gs%MS(i,n)%dmassdt(imr,j)<-1.0e+03)then
!!$                   write(*,232) n,j,gs%MS(i,n)%dmassdt(imr,j)
!!$232                format("Warning rimemass>grid,process,dmassdt",2i5,es15.6)
!!$!!c                   stop
!!$                end if
!!$
!!$                ndxdt=ndxdt+gs%MS(i,n)%dmassdt(imr,j)
!!$             end do
!!$             if( out_type == 2 ) then
!!$                XS(imr_q,i,nscat,n)=min(XS(imt_q,i,NSCAT,n)-XS(imc_q,i,NSCAT,n),&
!!$                     max((gs%MS(i,n)%mass(imr)+ndxdt*gs%dt)/ag%tv(n)%den,0.0_ps))
!!$             end if
!!$
!!$             ! mass by agg process (g/g)
!!$             ndxdt=0.0_ps
!!$             do j=1,gs%n_tendpros
!!$                if(gs%MS(i,n)%dmassdt(ima,j)>1.0e+03.or.&
!!$                     gs%MS(i,n)%dmassdt(ima,j)<-1.0e+03)then
!!$                   write(*,247) n,j,gs%MS(i,n)%dmassdt(ima,j)
!!$247                format("Warning aggmass>grid,process,dmassdt",2i5,es15.6)
!!$!!c                   stop
!!$                end if
!!$
!!$                ndxdt=ndxdt+gs%MS(i,n)%dmassdt(ima,j)
!!$             end do
!!$             if( out_type == 2 ) then
!!$                XS(ima_q,i,nscat,n)=min(XS(imt_q,i,NSCAT,n)-XS(imc_q,i,NSCAT,n),&
!!$                     max((gs%MS(i,n)%mass(ima)+ndxdt*gs%dt)/ag%tv(n)%den,0.0_ps))
!!$             end if
!!$!!c             if(ndxdt>0.0) then
!!$!!c                write(*,'("cal_model_out> agg exist",2I5,5ES15.6)') i,n,XS(imt_q,i,nscat,n),XS(imc_q,i,nscat,n),XS(ima_q,i,nscat,n),ndxdt
!!$!!c             end if
!!$
!!$             if(level>=6) then
!!$                ! melt water mass  (g/g)
!!$                ndxdt=0.0_ps
!!$                do j=1,gs%n_tendpros
!!$                   if(gs%MS(i,n)%dmassdt(imw,j)>1.0e+08.or.&
!!$                        gs%MS(i,n)%dmassdt(imw,j)<-1.0e+08)then
!!$                      write(*,239) n,j,gs%MS(i,n)%dmassdt(imw,j)
!!$239                   format("Warning meltmass>grid,process,dmassdt",2i5,es15.6)
!!$!!c                      stop
!!$                   end if
!!$
!!$                   ndxdt=ndxdt+gs%MS(i,n)%dmassdt(imw,j)
!!$                end do
!!$                if( out_type == 2 ) then
!!$                   XS(imw_q,i,nscat,n)=min(max(XS(imt_q,i,NSCAT,n)-XS(imc_q,i,NSCAT,n)-&
!!$                        XS(imr_q,i,NSCAT,N),0.0_PS),&
!!$                        max((gs%MS(i,n)%mass(imw)+ndxdt*gs%dt)/ag%tv(n)%den,0.0_ps))
!!$                end if
!!$             end if
!!$
!!$             if(level>=6) then
!!$                ! freezing nucleation mass  (g/g)
!!$                ndxdt=0.0_ps
!!$                do j=1,gs%n_tendpros
!!$                   if(gs%MS(i,n)%dmassdt(imf,j)>1.0e+08.or.&
!!$                        gs%MS(i,n)%dmassdt(imf,j)<-1.0e+08)then
!!$                      write(*,249) n,j,gs%MS(i,n)%dmassdt(imf,j)
!!$249                   format("Warning meltmass>grid,process,dmassdt",2i5,es15.6)
!!$!!c                      stop
!!$                   end if
!!$
!!$                   ndxdt=ndxdt+gs%MS(i,n)%dmassdt(imf,j)
!!$                end do
!!$                if( out_type == 2 ) then
!!$                   XS(imf_q,i,nscat,n)=min(XS(imc_q,i,NSCAT,n),&
!!$                        max((gs%MS(i,n)%mass(imf)+ndxdt*gs%dt)/ag%tv(n)%den,0.0_ps))
!!$                end if
!!$             end if
!!$
!!$             if(level>=4) then
!!$                ! mass by total mass of aerosols (g/g)
!!$                ndxdt=0.0_ps
!!$                do j=1,gs%n_tendpros
!!$                   if(gs%MS(i,n)%dmassdt(imat,j)>1.0e+03.or.&
!!$                        gs%MS(i,n)%dmassdt(imat,j)<-1.0e+03)then
!!$                      write(*,234) KD(n),ID(n),JD(n),i,j,gs%MS(i,n)%dmassdt(imat,j)
!!$234                   format("Warning aptmass>grid(k,i,j),bin,process,dmassdt",5i5,es15.6)
!!$!!c                      stop
!!$                   end if
!!$
!!$                   ndxdt=ndxdt+gs%MS(i,n)%dmassdt(imat,j)
!!$                end do
!!$                if( out_type == 2 ) then
!!$                   ndxdt=(gs%MS(i,n)%mass(imat)+ndxdt*gs%dt)/ag%tv(n)%den
!!$                   fap=ndxdt/XS(imt_q,i,nscat,n)
!!$                   if(fap>1.0_PS) then
!!$!!c                      write(*,*) "total mass is less than ap total 2",XS(imt_q,i,nscat,n),XS(imat_q,i,nscat,n)
!!$                      XS(imat_q,i,nscat,n)=ndxdt
!!$                      XS(imt_q,i,nscat,n)=XS(imat_q,i,nscat,n)
!!$                   elseif(fap>min_fapt_s) then
!!$                      XS(imat_q,i,nscat,n)=ndxdt
!!$                      fapt_s=fap
!!$                   else
!!$!!c                      XS(imat_q,i,nscat,n)=fapt_s*XS(imt_q,i,nscat,n)
!!$                      XS(imat_q,i,nscat,n)=min(max(m_lmt_ap/ag%tv(n)%den,min_fapt_s*XS(imt_q,i,nscat,n))&
!!$                                         ,mx_aprat*XS(imt_q,i,nscat,n))
!!$                   end if
!!$                end if
!!$                ! mass by soluble mass of aerosols (g/g)
!!$                ndxdt=0.0_ps
!!$                do j=1,gs%n_tendpros
!!$                   if(gs%MS(i,n)%dmassdt(imas,j)>1.0e+03.or.&
!!$                        gs%MS(i,n)%dmassdt(imas,j)<-1.0e+03)then
!!$                      write(*,235) KD(n),ID(n),JD(n),i,j,gs%MS(i,n)%dmassdt(imas,j)
!!$235                   format("Warning apsmass>grid(k,i,j),bin,process,dmassdt",5i5,es15.6)
!!$!!c                      stop
!!$                   end if
!!$
!!$                   ndxdt=ndxdt+gs%MS(i,n)%dmassdt(imas,j)
!!$                end do
!!$                if( out_type == 2 ) then
!!$                   ndxdt=(gs%MS(i,n)%mass(imas)+ndxdt*gs%dt)/ag%tv(n)%den
!!$                   fap=XS(imas_q,i,nscat,n)/XS(imat_q,i,nscat,n)
!!$                   if(fap>=min_faps_s) then
!!$                      XS(imas_q,i,nscat,n)=ndxdt
!!$                      faps_s=fap
!!$                   else
!!$!!c                      XS(imas_q,i,nscat,n)=faps_s*XS(imat_q,i,nscat,n)
!!$                      XS(imas_q,i,nscat,n)=min(max(m_lmt_ap/ag%tv(n)%den,&
!!$                           min_faps_s*XS(imat_q,i,nscat,n)),XS(imat_q,i,nscat,n))
!!$!!c                      XS(imas_q,i,nscat,n)=min(eps_ap0(2)*XS(imat_q,i,nscat,n),XS(imat_q,i,nscat,n))
!!$                   end if
!!$                end if
!!$             end if
!!$
!!$!!c             write(*,'("XS out>",4I4,20ES14.6)') KD(n),ID(n),JD(n),i,XS(1:mxnpi,i,NSCAT,n)
!!$!!c             write(*,'("tend out>",4I4,20ES14.6)') KD(n),ID(n),JD(n),i,(gs%MS(i,n)%dmassdt(imt,j),j=1,gs%n_tendpros)
!!$!!c             write(*,'("tend out vol>",4I4,20ES14.6)') KD(n),ID(n),JD(n),i,gs%MS(i,n)%dvoldt(1:7,1)
!!$
!!$!!c             if(JD(n)>50.and.JD(n)<150.and.KD(n)>=37) then
!!$!!c                 write(*,'("cal_model:ibin,ngrid,j,k,con,mass4",4I5,2ES15.6)') &
!!$!!c                               i,n,JD(n),KD(n),XS(icon_q,i,nscat,n)*ag%TV(n)%den,XS(imc_q,i,nscat,n)*ag%TV(n)%den
!!$!!c             end if
!!$
!!$
!!$!!c             alen_n=(XS(iacr_q,i,nscat,n)/XS(icon_q,i,nscat,n))**(1.0/3.0)
!!$!!c             clen_n=(XS(iccr_q,i,nscat,n)/XS(icon_q,i,nscat,n))**(1.0/3.0)
!!$!!c             den_n=XS(imc_q,i,nscat,n)/XS(icon_q,i,nscat,n)/(4.0_PS*PI/3.0_PS*(alen_n**2+clen_n**2)**1.5)
!!$
!!$!!c             if(den_n<1.0e-3) then
!!$!!c                write(*,'("cal_model_tend>small den:bin,n,i,j,k,mass_comp",5i5,10ES15.6)') i,n,ID(n),JD(n),KD(n),&
!!$!!c                     XS(imt_q,i,nscat,n),XS(imr_q,i,nscat,n),XS(ima_q,i,nscat,n),&
!!$!!c                     XS(imc_q,i,nscat,n),XS(imw_q,i,nscat,n)/XS(imt_q,i,nscat,n)
!!$!!c                write(*,*) "p,temp,sv1",ag%TV(n)%P,ag%TV(n)%T,ag%TV(n)%s_v_n(1)
!!$!!c                write(*,'("dmassdt mass1",20ES15.6)') (gs%MS(i,n)%dmassdt(imt,j),j=1,gs%N_tendpros)
!!$!!c                write(*,'("dmassdt mass4",20ES15.6)') (gs%MS(i,n)%dmassdt(imc,j),j=1,gs%N_tendpros)
!!$!!c                write(*,'("dmassdt mass5",20ES15.6)') (gs%MS(i,n)%dmassdt(imat,j),j=1,gs%N_tendpros)
!!$!!c                write(*,'("dvoldt alen",20ES15.6)') (gs%MS(i,n)%dvoldt(2,j),j=1,gs%N_tendpros)
!!$!!c                write(*,'("dvoldt clen",20ES15.6)') (gs%MS(i,n)%dvoldt(3,j),j=1,gs%N_tendpros)
!!$!!c                write(*,'("con,den,denip,denic,mm",10ES15.6)') gs%MS(i,n)%con,gs%MS(i,n)%den,gs%IS(i,n)%den_ip,gs%IS(i,n)%den_ic,&
!!$!!c                                gs%MS(i,n)%mean_mass
!!$!!c                write(*,*) "mass comp",gs%MS(i,n)%mass
!!$!!c                write(*,*) "sh,ht",gs%IS(i,n)%sh_type,gs%IS(i,n)%habit
!!$!!c                write(*,*) "alen,clen,dlen,rlen,elen",gs%MS(i,n)%a_len,gs%MS(i,n)%c_len,gs%IS(i,n)%d,gs%IS(i,n)%r,gs%IS(i,n)%e
!!$!!c                write(*,*) "new alen,clen,dlen,rlen,elen",(XS(iacr_q,i,NSCAT,N)/XS(icon_q,i,nscat,n))**(1.0/3.0)&
!!$!!c                           ,(XS(iccr_q,i,NSCAT,N)/XS(icon_q,i,nscat,n))**(1.0/3.0),(XS(idcr_q,i,NSCAT,N)/XS(icon_q,i,nscat,n))**(1.0/3.0)&
!!$!!c                           ,(XS(iag_q,i,NSCAT,N)/XS(icon_q,i,nscat,n))**(1.0/3.0),(XS(icg_q,i,NSCAT,N)/XS(icon_q,i,nscat,n))**(1.0/3.0)
!!$!!c                write(*,'("new con,den_n,mm_n",5ES15.6)') XS(icon_q,i,NSCAT,N),den_n,XS(imt_q,i,nscat,n)/XS(icon_q,i,NSCAT,n)
!!$!!c             end if
!!$!!c             if(XS(imr_q,i,nscat,n)<0.0) then
!!$!!c                write(*,'("mass xs",2i5,3es15.6)') i,n,XS(imt_q,i,nscat,n),XS(inex_q,i,nscat,n),XS(imr_q,i,nscat,n)
!!$!!c             end if
!!$
!!$             ! volume by riming process (cm^3/g)
!!$!!c          ndxdt=0.0_ps
!!$!!c          do j=1,gs%n_tendpros
!!$!!c             ndxdt=ndxdt+gs%MS(i,n)%dvoldt(icg,j)
!!$!!c          end do
!!$!!c          if( out_type == 2 ) then
!!$!!c             XS(i,10)=(gs%MS(i,n)%vol(1)*gs%MS(i,n)%con+ndxdt*gs%dt)/ag%tv(n)%den
!!$!!c          end if
!!$!!c          dxdt(i,10)=dxdt(i,10)+ndxdt/ag%tv(n)%den
!!$
!!$!!c             do j=1,gs%n_tendpros
!!$!!c                if(abs(gs%MS(i,n)%dmassdt(imt,j))>0.0_ps) then
!!$!!c                   if((gs%MS(i,n)%dmassdt(imc,j)/gs%MS(i,n)%dmassdt(imt,j)<0.9999).or.&
!!$!!c                      (gs%MS(i,n)%dmassdt(imc,j)/gs%MS(i,n)%dmassdt(imt,j)>1.0001)) then
!!$!!c                      write(*,25) n,j,gs%MS(i,n)%dmassdt(imt,j),gs%MS(i,n)%dmassdt(imc,j)
!!$!!c25                format("ratcheck>grid,process,dmassdt1,4",2i5,2es15.6)
!!$!!c                   stop
!!$!!c                   end if
!!$!!c                end if
!!$!!c             end do
!!$
!!$             totmixs=totmixs+max(0.0_PS,XS(imt_q,i,nscat,n)-XS(imat_q,i,nscat,n))
!!$          end do
!!$       end if
!!$
!!$
!!$       ! +++ update vapor field +++
!!$       if( out_type == 2 ) then
!!$          rv(n)=ag%tv(n)%rv+ag%tv(n)%dmassdt/ag%tv(n)%den*gr%dt
!!$          rv2=qtp(n)-totmixr-totmixs
!!$          e_n=ag%TV(n)%P*rv(n)/(Rdvchiarui+rv(n))
!!$          s_n=e_n/get_sat_vapor_pres_lk(1,ag%TV(n)%T_m,ag%estbar,ag%esitbar)-1.0
!!$
!!$!kid          if(s_n>0.10_PS) then
!!$          if(s_n>0.50_PS) then
!!$             write(*,'("cal_model:k,i,j,p,qv,qr,qi,t,tn,sv1,svn1,sn,rv2",3I5,20ES15.6)') KD(n),ID(n),JD(n),ag%TV(n)%P,&
!!$                  rv(n),totmixr,totmixs,ag%TV(n)%T,&
!!$                  ag%TV(n)%T_n,ag%TV(n)%s_v(1),ag%TV(n)%s_v_n(1),s_n,rv2
!!$          endif
!!$
!!$!kid          if(rv(n)<0.0_PS.or.s_n>0.10) then
!!$          if(rv(n)<0.0_PS.or.s_n>0.50) then
!!$             if(rv(n)<0.0_PS) then
!!$                write(*,*) "rv is negative"
!!$
!!$                ! assume vapor mixing ratio to 10% of ice saturation.
!!$                rv_n=(-0.9+1.0_PS)*ag%tv(n)%e_sat(2)*Rdvchiarui/&
!!$                      (ag%tv(n)%P-(-0.9+1.0_PS)*ag%tv(n)%e_sat(2))
!!$
!!$                write(*,'("cal_model neg",3I5,20ES15.6)') KD(n),ID(n),JD(n),ag%TV(n)%P,&
!!$                  qtp(n),rv(n),rv_n,totmixr,totmixs,ag%TV(n)%T,&
!!$                  ag%TV(n)%T_n,ag%TV(n)%s_v(1),ag%TV(n)%s_v_n(1),s_n
!!$
!!$                rv(n)=rv_n
!!$
!!$             else
!!$                write(*,*) "rv(n) is super saturated"
!!$                write(*,'("cal_model su",3I5,20ES15.6)') KD(n),ID(n),JD(n),ag%TV(n)%P,&
!!$                  qtp(n),rv(n),totmixr,totmixs,ag%TV(n)%T,&
!!$                  ag%TV(n)%T_n,ag%TV(n)%s_v(1),ag%TV(n)%s_v_n(1),s_n
!!$             endif
!!$
!!$             write(*,*) "mes_rc,density of air",mes_rc(n),ag%TV(n)%den
!!$             write(*,*) "rain"
!!$             write(*,'(30ES15.6)') (gr%MS(i,n)%mass(rmt),i=1,gr%N_BIN)
!!$             do i=1,gr%N_BIN
!!$                write(*,'(I5,15ES15.6)') i,(gr%MS(i,n)%dmassdt(rmt,j),j=1,gr%N_tendpros)
!!$             end do
!!$
!!$             write(*,*) "ice"
!!$             write(*,'(30ES15.6)') (gs%MS(i,n)%mass(imt),i=1,gs%N_BIN)
!!$             do i=1,gs%N_BIN
!!$                write(*,'(I5,15ES15.6)') i,(gs%MS(i,n)%dmassdt(imt,j),j=1,gr%N_tendpros)
!!$             end do
!!$             write(*,'("sv1,sv2",5ES15.6)') ag%TV(n)%s_v(1),ag%TV(n)%s_v(2),&
!!$                   ag%TV(n)%s_v_n(1),ag%TV(n)%s_v_n(2)
!!$!!c             stop
!!$          endif
!!$!!c          rv_n=(ag%tv(n)%s_v_n(1)+1.0_PS)*ag%tv(n)%e_sat_n(1)*0.622/&
!!$!!c               (ag%tv(n)%P-(ag%tv(n)%s_v_n(1)+1.0_PS)*ag%tv(n)%e_sat_n(1))
!!$!!c
!!$!!c          write(*,*) "rv(n),rv_n ",KD(n),ID(n),JD(n),qtp(n),rv(n),totmix,rv_n
!!$!!c          if(abs(rv(n)-rv_n)/rv_n>1.0e-4_PS) then
!!$!!c             write(*,*) "rv(n),rv_n dont match!",KD(n),ID(n),JD(n),rv(n),rv_n
!!$!!c             write(*,*) ag%tv(n)%rv,ag%tv(n)%dmassdt,ag%tv(n)%s_v_n(1),ag%tv(n)%T,ag%tv(n)%P
!!$!!c!!c             stop
!!$!!c          end if
!!$       end if
!!$
!!$
!!$    end do
!!$
!!$    if( flagp_a/=0) then
!!$
!!$       do ic=1,nacat
!!$          ! +++ for aerosols +++
!!$          if((flagp_a==-1.or.flagp_a==-4).and.ic==2) then
!!$             cycle
!!$          elseif((flagp_a==-2.or.flagp_a==-5).and.ic/=2) then
!!$             cycle
!!$          elseif(flagp_a==-3.or.flagp_a==-6) then
!!$             cycle
!!$          end if
!!$
!!$          do n=1,L
!!$
!!$             do i=1,ga(ic)%n_bin
!!$
!!$!!                ndxdt = 0.0_ps
!!$!!                do j = 1, ga(ic)%n_tendpros
!!$!!                   if(ga(ic)%MS(i,n)%dmassdt(amt,j)>1.0e+03.or.&
!!$!!                        ga(ic)%MS(i,n)%dmassdt(amt,j)<-1.0e+03)then
!!$!!                      write(*,201) KD(n),ID(n),JD(n),i,j,ga(ic)%MS(i,n)%dmassdt(amt,j)
!!$!!201                   format("mt_rain>grid,bin,process,dmassdt",5i5,es15.6)
!!$!!                      stop
!!$!!                   end if
!!$!!
!!$!!                   ndxdt = ndxdt+ga(ic)%MS(i,n)%dmassdt(amt,j)
!!$!!                end do
!!$!!                if( out_type == 2 ) then
!!$!!                   XA(amt_q,i,ic,n)=max(0.0_ps,&
!!$!!                        (ga(ic)%MS(i,n)%mass(amt)+ndxdt*ga(ic)%dt)/ag%tv(n)%den)
!!$!!                   ! check
!!$!!                   if(XA(amt_q,i,ic,n)==0.0_ps) then
!!$!!                      XA(acon_q,i,ic,n)=0.0_ps
!!$!!                      cycle
!!$!!                   else if(XA(amt_q,i,ic,n)>1.0e-1) then
!!$!!                      write(*,251) KD(n),ID(n),JD(n),i,XA(amt_q,i,ic,n)
!!$!!251                   format("warning:qa is large at bin,gaid:",4i5,2es15.6)
!!$!!                   end if
!!$!!                end if
!!$
!!$                ndxdt = 0.0_ps
!!$                do j = 1, ga(ic)%n_tendpros
!!$                  ! added for kwajex
!!$                  ! skip evap feed back.
!!$                  if(j==1) cycle
!!$                  ! end added for kwajex
!!$
!!$!!c                   if((ga(ic)%MS(i,n)%dcondt(j)<-1.0e-5 .or.ga(ic)%MS(i,n)%dcondt(j)>1.0e-5) ) then
!!$!!c                   if(ic==2.and.ga(ic)%MS(i,n)%dcondt(j)>1.0e-4) then
!!$!!c                      write(*,*) "con_a:ic,process,grid,dcdt",ic,j,n,ga(ic)%MS(i,n)%dcondt(j)
!!$!!c                      write(*,*) "m1,m2,m3",ga(ic)%MS(i,n)%mass
!!$!!c                      write(*,*) "con", ga(ic)%MS(i,n)%con
!!$!!c                      write(*,*) "svw,i",ag%tv(n)%s_v(1),ag%tv(n)%s_v(2)
!!$!!c                      write(*,*) "dsvw,i",ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2)
!!$!!c!!c                      stop
!!$!!c                   end if
!!$                   ndxdt = ndxdt+ga(ic)%MS(i,n)%dcondt(j)
!!$                end do
!!$                if( out_type == 2 ) XA(acon_q,i,ic,n)=max(0.0_PS,&
!!$                     (ga(ic)%MS(i,n)%con+ndxdt*ga(ic)%dt))/ag%tv(n)%den
!!$                if(XA(acon_q,i,ic,n)>1.0e+9) then
!!$                   write(*,'("cal_model_out:xa,con,connew",4I5,4ES15.6)') KD(n),ID(n),JD(n),ic,&
!!$                            XA(acon_q,i,ic,n),ga(ic)%MS(i,n)%con,XA(acon_q,i,ic,n)*ag%tv(n)%den
!!$
!!$                   write(*,971) i,ID(n),JD(n),KD(n),(gr%MS(j,n)%dcondt(1),j=1,gr%N_BIN)
!!$971                   format("qr >bin,grid,dcondt",4i5,50es15.6)
!!$                   write(*,972) i,ID(n),JD(n),KD(n),(gr%MS(j,n)%dmassdt(rmt,1),j=1,gr%N_BIN)
!!$972                   format("qr >bin,grid,dmassdt",4i5,50es15.6)
!!$                   write(*,973) ic,(ga(ic)%ms(1,n)%dcondt(j),j=1,ga(ic)%n_tendpros)
!!$973                   format("qa >icat,dcondt",i5,50es15.6)
!!$                   write(*,974) ic,(ga(ic)%MS(1,n)%dmassdt(amt,j),j=1,ga(ic)%n_tendpros)
!!$974                   format("qa >icat,dmassdt",i5,50es15.6)
!!$                   write(*,975) i,ID(n),JD(n),KD(n),(gs%MS(j,n)%dcondt(1),j=1,gs%N_BIN)
!!$975                   format("qi >bin,grid,dcondt",4i5,50es15.6)
!!$                   write(*,976) i,ID(n),JD(n),KD(n),(gs%MS(j,n)%dmassdt(imt,1),j=1,gs%N_BIN)
!!$976                   format("qi >bin,grid,dmassdt",4i5,50es15.6)
!!$                end if
!!$!!c             write(*,'("i,n,mten,cten",2i5,10es15.6)') i,n,ndxdt1,ndxdt,&
!!$!!c                  ga(ic)%MS(i,n)%dmassdt(amt,1),ga(ic)%MS(i,n)%dmassdt(amt,9),&
!!$!!c                  ga(ic)%MS(i,n)%dcondt(1),ga(ic)%MS(i,n)%dcondt(9)
!!$!!c          dxdt(2,i,n)=dxdt(2,i,n)+ndxdt/ag%tv(n)%den
!!$
!!$                if(level>=4) then
!!$                   ! total mass of aerosols (g/g)
!!$                   ndxdt=0.0_ps
!!$                   do j=1,ga(ic)%n_tendpros
!!$                      ! added for kwajex
!!$                      ! skip evap feed back.
!!$                      if(j==1) cycle
!!$                      ! end added for kwajex
!!$
!!$                      if(ga(ic)%MS(i,n)%dmassdt(amt,j)>1.0e-7.or.&
!!$                           ga(ic)%MS(i,n)%dmassdt(amt,j)<-1.0e-7)then
!!$                         write(*,288) ic,j,n,ID(n),JD(n),KD(n),ga(ic)%MS(i,n)%dmassdt(amt,j)
!!$288                      format("Warning apsm in ap>ga(ic)id,process,grid,id,jd,kd,dmassdt",6i5,es15.6)
!!$                         write(*,*) "m1,m2,m3",ga(ic)%MS(i,n)%mass
!!$                         write(*,*) "con", ga(ic)%MS(i,n)%con
!!$                         write(*,*) "svw,i",ag%tv(n)%s_v(1),ag%tv(n)%s_v(2)
!!$                         write(*,*) "dsvw,i",ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2)
!!$!!c                         stop
!!$                      end if
!!$!!c                      if(ic==2.and.abs(ga(ic)%MS(i,n)%dmassdt(amt,j))>1.0e-30) then
!!$!!c                         write(*,288) ic,j,ga(ic)%MS(i,n)%dmassdt(amt,j)
!!$!!c                         write(*,*) "m1,m2,m3",ga(ic)%MS(i,n)%mass
!!$!!c                         write(*,*) "con", ga(ic)%MS(i,n)%con
!!$!!c                         write(*,*) "svw,i",ag%tv(n)%s_v(1),ag%tv(n)%s_v(2)
!!$!!c                         write(*,*) "nmass",ga(ic)%n_masscom
!!$!!c                         stop
!!$!!c                      end if
!!$
!!$                      ndxdt=ndxdt+ga(ic)%MS(i,n)%dmassdt(amt,j)
!!$                   end do
!!$                   if( out_type == 2 ) then
!!$                      XA(amt_q,i,ic,n)=max((ga(ic)%MS(i,n)%mass(amt)+&
!!$                           ndxdt*ga(ic)%dt)/ag%tv(n)%den,0.0_ps)
!!$                   end if
!!$
!!$                   ! mass by soluble mass of aerosols (g/g)
!!$                   ndxdt=0.0_ps
!!$                   do j=1,ga(ic)%n_tendpros
!!$
!!$                      ! added for kwajex
!!$                      ! skip evap feed back.
!!$                      if(j==1) cycle
!!$                      ! end added for kwajex
!!$
!!$                      if(ga(ic)%MS(i,n)%dmassdt(ams,j)>1.0e+08.or.&
!!$                           ga(ic)%MS(i,n)%dmassdt(ams,j)<-1.0e+08)then
!!$                         write(*,287) ic,j,ga(ic)%MS(i,n)%dmassdt(ams,j)
!!$287                      format("Warning apsm in ap>ga(ic)id,process,dmassdt",2i5,es15.6)
!!$                         write(*,*) "m1,m2,m3",ga(ic)%MS(i,n)%mass
!!$                         write(*,*) "con", ga(ic)%MS(i,n)%con
!!$                         write(*,*) "svw,i",ag%tv(n)%s_v(1),ag%tv(n)%s_v(2)
!!$                         write(*,*) "nmass",ga(ic)%n_masscom
!!$!!c                         stop
!!$                      end if
!!$
!!$                      ndxdt=ndxdt+ga(ic)%MS(i,n)%dmassdt(ams,j)
!!$                   end do
!!$                   if( out_type == 2 ) then
!!$                      XA(ams_q,i,ic,n)=min(max((ga(ic)%MS(i,n)%mass(ams)+&
!!$                           ndxdt*ga(ic)%dt)/ag%tv(n)%den,0.0_ps),XA(amt_q,i,ic,n))
!!$                      if(ic==1) then
!!$                         ! for CCN accumulation mode (category 1)
!!$!!c                         if(XA(amt_q,i,ic,n)*sep_faps>XA(ams_q,i,ic,n)) then
!!$!!c                            write(*,*) "sol mass is less than limt: ic,n",ic,n
!!$!!c                            write(*,*) XA(amt_q,i,ic,n),XA(acon_q,i,ic,n),XA(ams_q,i,ic,n)
!!$!!c                            write(*,'("dmassdt 1",15ES15.6)') (ga(ic)%MS(i,n)%dmassdt(amt,j),j=1,ga(ic)%N_tendpros)
!!$!!c                            write(*,'("dmassdt 2",15ES15.6)') (ga(ic)%MS(i,n)%dmassdt(ams,j),j=1,ga(ic)%N_tendpros)
!!$!!c                         end if
!!$                         XA(ams_q,i,ic,n)=min(max(XA(amt_q,i,ic,n)*sep_faps,&
!!$                              XA(ams_q,i,ic,n)),XA(amt_q,i,ic,n))
!!$                      elseif(ic==2) then
!!$                         ! for IN accumulation mode (category 2)
!!$                         XA(ams_q,i,ic,n)=max(0.0_PS,min(XA(amt_q,i,ic,n)*0.99_PS*sep_faps,&
!!$                              XA(ams_q,i,ic,n)))
!!$                      elseif(ic==3) then
!!$                         ! for CCN nucleation mode (category 3)
!!$                         XA(ams_q,i,ic,n)=min(max(XA(amt_q,i,ic,n)*sep_faps,&
!!$                              XA(ams_q,i,ic,n)),XA(amt_q,i,ic,n))
!!$                      elseif(ic==4) then
!!$                         ! for CCN nucleation mode (category 4)
!!$                         XA(ams_q,i,ic,n)=min(max(XA(amt_q,i,ic,n)*sep_faps,&
!!$                              XA(ams_q,i,ic,n)),XA(amt_q,i,ic,n))
!!$                      else
!!$                         write(*,*) "cal_model_tend_all>this ap cat not defined"
!!$                         stop
!!$                      end if
!!$                   end if
!!$
!!$                   if(ic==1.or.ic==2.or.ic==4) then
!!$                      m_lmt=coef4pi3*r3_lmt*ga(ic)%MS(i,n)%den
!!$                   elseif(ic==3) then
!!$                      m_lmt=coef4pi3*r3_lmt_nuc*ga(ic)%MS(i,n)%den
!!$                   else
!!$                      write(*,*) "this category not defined in ini_group_all"
!!$                      stop
!!$                   endif
!!$
!!$                   if(XA(acon_q,i,ic,n)<n_lmt_ap/ag%tv(n)%den) then
!!$!!c                      write(*,'("ap fixed ic",i5,3es15.6)')&
!!$!!c                           ic,XA(acon_q,i,ic,n)*ag%tv(n)%den,n_lmt_ap
!!$!!c                      rmod1=XA(amt_q,i,ic,n)/XA(acon_q,i,ic,n)
!!$!!c                      rmod2=XA(ams_q,i,ic,n)/XA(acon_q,i,ic,n)
!!$                      XA(acon_q,i,ic,n)=n_lmt_ap/ag%tv(n)%den
!!$                      XA(amt_q,i,ic,n)=XA(acon_q,i,ic,n)*m_lmt
!!$                      XA(ams_q,i,ic,n)=XA(amt_q,i,ic,n)*eps_ap0(ic)
!!$!!c                      XA(amt_q,i,ic,n)=XA(acon_q,i,ic,n)*rmod1
!!$!!c                      XA(ams_q,i,ic,n)=XA(acon_q,i,ic,n)*rmod2
!!$                   end if
!!$                   if(XA(acon_q,i,ic,n)>n_max_ap/ag%tv(n)%den) then
!!$!!c                      write(*,'("ap con reached max. so fixed.",i5,3es15.6)')&
!!$!!c                           ic,XA(acon_q,i,ic,n)*ag%tv(n)%den,n_max_ap
!!$                      XA(acon_q,i,ic,n)=n_max_ap/ag%tv(n)%den
!!$                   end if
!!$                   if(XA(amt_q,i,ic,n)<m_lmt*XA(acon_q,i,ic,n)) then
!!$!!c                      write(*,'("ap reached min:cat,xa1,min,eps0",i5,10es15.6)')&
!!$!!c                           ic,XA(amt_q,i,ic,n),m_lmt*XA(acon_q,i,ic,n),eps_ap0(ic)
!!$                      XA(amt_q,i,ic,n)=XA(acon_q,i,ic,n)*m_lmt
!!$                      XA(ams_q,i,ic,n)=XA(amt_q,i,ic,n)*eps_ap0(ic)
!!$                   end if
!!$                   XA(ams_q,i,ic,n)=min(XA(ams_q,i,ic,n),XA(amt_q,i,ic,n))
!!$!!c                      if(XA(ams_q,i,ic,n)>0.0_ps.and.XA(amt_q,i,ic,n)>0.0_ps) then
!!$!!c                         write(*,'("ap fixed2, ic,mmass,min_mass",i5,3es15.6)')&
!!$!!c                              ic,XA(amt_q,i,ic,n)/XA(acon_q,i,ic,n),coef4pi3*r3_lmt*ga(ic)%MS(i,n)%den
!!$!!c                         rmod1=XA(ams_q,i,ic,n)/XA(amt_q,i,ic,n)
!!$!!c                         XA(amt_q,i,ic,n)=XA(acon_q,i,ic,n)*coef4pi3*r3_lmt*ga(ic)%MS(i,n)%den
!!$!!c                         XA(ams_q,i,ic,n)=XA(amt_q,i,ic,n)*rmod1
!!$!!c                      else
!!$!!c                         write(*,'("ap fixed3, ic,mmass,min_mass",i5,3es15.6)')&
!!$!!c                              ic,XA(amt_q,i,ic,n)/XA(acon_q,i,ic,n),coef4pi3*r3_lmt*ga(ic)%MS(i,n)%den
!!$!!c                         XA(amt_q,i,ic,n)=XA(acon_q,i,ic,n)*coef4pi3*r3_lmt*ga(ic)%MS(i,n)%den
!!$!!c                         XA(ams_q,i,ic,n)=XA(amt_q,i,ic,n)*ga(ic)%MS(i,n)%eps_map
!!$!!c                      end if
!!$!!c                   end if
!!$!!c                    if(XA(amt_q,i,ic,n)>1.0e-7) then
!!$!!c                       write(*,'("ap mass is large: ic,n,id,jd,kd",6I5)') ic,n,ID(n),JD(n),KD(n)
!!$!!c                       write(*,'("xa1,xa2,xa3",5ES15.6)') XA(amt_q,i,ic,n)*ag%TV(n)%den,XA(acon_q,i,ic,n)*ag%TV(n)%den&
!!$!!c                                                     ,XA(ams_q,i,ic,n)*ag%TV(n)%den
!!$!!c                       write(*,'("density,con,mass_t,mass_sol",6ES15.6)') ga(ic)%MS(i,n)%den&
!!$!!c                                      ,ga(ic)%MS(i,n)%con,ga(ic)%MS(i,n)%mass(amt),ga(ic)%MS(i,n)%mass(ams)
!!$!!c                       write(*,'("dcondt 1",15ES15.6)') (ga(ic)%MS(i,n)%dcondt(j),j=1,ga(ic)%N_tendpros)
!!$!!c                       write(*,'("dmassdt 1",15ES15.6)') (ga(ic)%MS(i,n)%dmassdt(amt,j),j=1,ga(ic)%N_tendpros)
!!$!!c                       write(*,'("dmassdt 2",15EES15.6)') (ga(ic)%MS(i,n)%dmassdt(ams,j),j=1,ga(ic)%N_tendpros)
!!$!!c                    end if
!!$
!!$                end if
!!$             end do
!!$          end do
!!$       end do
!!$    end if
!!$
!!$
!!$
!!$
!!$!!c    call qiprint2(xs,nstype,nsbin,nscat,l,"af mten")
!!$!!c    call qiprint3(gs,l,"af mten")
!!$  end subroutine cal_model_tendency_all_scl

!!c  subroutine cal_dmdt_ratio(g,dmdt_ratio1,dmdt_ratio2)
!!c    type (group), intent(in)         :: g
!!c    real(ps),pointer,dimension(:) :: dmdt_ratio1
!!c    real(ps),pointer,dimension(:,:) :: dmdt_ratio2
!!c    ! dmdt_ratio
!!c    real(ps) :: dmdt1
!!c    integer :: i,j,k
!!c    dmdt1=0.0_ps
!!c    dmdt_ratio1=0.0_ps
!!c    dmdt_ratio2=0.0_ps
!!c    do i = 1, g%n_bin
!!c       do j = 1, g%n_tendpros
!!c          dmdt1=dmdt1+g%ms(i)%dmassdt(1,j)
!!c       end do
!!c       do k = 1, g%n_masscom
!!c          do j = 1, g%n_tendpros
!!c             dmdt_ratio2(k,j)=g%ms(i)%dmassdt(1+k,j)/dmdt1
!!c             dmdt_ratio1(k)=dmdt_ratio1(k)+g%ms(i)%dmassdt(1+k,j)
!!c          end do
!!c          dmdt_ratio1(k)=dmdt_ratio1(k)/dmdt1
!!c       end do
!!c    end do
!!c  end subroutine cal_dmdt_ratio

!!$  subroutine cal_model_variable(g,x,dxdt,nxtype)
!!$    type (group), intent(in) :: g
!!$    real(PS), dimension(:,:) :: x
!!$    real(PS), dimension(:,:) :: dxdt
!!$    integer, intent(in)    :: nxtype
!!$    integer  :: i,j
!!$    do i=1,g%n_bin
!!$       do j=1,nxtype
!!$          x(i,j) = x(i,j) + dxdt(i,j)*g%dt
!!$       end do
!!$    end do
!!$   end subroutine cal_model_variable

!!$  subroutine update_surface_temp(level, g, ag,nu_aps,phi_aps,m_aps)
!!$    integer, intent(in)         :: level
!!$    type (group), intent(inout)         :: g
!!$    type (airgroup),intent(in)    :: ag
!!$    ! message from reality-check
!!$!tmp    integer,pointer,dimension(:)  :: mes_rc
!!$    !integer,dimension(*)   :: mes_rc
!!$    real(PS),dimension(*),intent(in) :: nu_aps,phi_aps,m_aps
!!$
!!$    integer    :: i,n
!!$    do n=1,g%L
!!$       ! +++ loop over grids +++
!!$       if( g%mark_cm(n) == 3 ) cycle
!!$       do i=1,g%n_bin
!!$
!!$          call cal_coef_vapdep2(2,level,g%MS(i,n),ag%TV(n),&
!!$               nu_aps,phi_aps,m_aps,g%IS(i,n))
!!$
!!$          call cal_surface_temp( 2, g%MS(i,n), ag%tv(n),ag%estbar,ag%esitbar)
!!$!!c          call cal_coef_vapdep(2,level,g%MS(i,n),ag%TV(n),g%IS(i,n))
!!$
!!$       end do
!!$    end do
!!$  end subroutine update_surface_temp

!!$  subroutine cal_capacitance(level,phase,th_var, ms, is)
!!$    ! +++ Calculate capacitance based on the spheroid assumption,
!!$    !     for example, see Chen and Lamb (1994a).                     +++
!!$    ! phase of hydrometeor
!!$    ! 1: liquid, 2: solid
!!$    integer, intent(in)  :: phase
!!$    integer, intent(in)  :: level
!!$    type (Thermo_Var), intent(in) :: th_var
!!$    type (Mass_Bin), intent(inout) :: ms
!!$    type (Ice_Shape), optional, intent(in)    :: is
!!$    real(PS)                    :: eps, d, out, phi
!!$    integer                     :: i
!!$    ! coefficient for capacitance by Chiruta and Wang (2003)
!!$!!c    real(PS),parameter   :: CAP_ros=0.434*N_lob**0.257
!!$    real(PS),parameter   :: CAP_ros=0.619753727
!!$    ! temperature of freezing
!!$    real(PS), parameter           :: T_0 = 273.16
!!$
!!$    ! initialize
!!$    ms%CAP=0.0_PS
!!$    ms%CAP_hex=0.0_PS
!!$    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return
!!$
!!$    if( phase == 1 ) then
!!$       ! --- in case of liquid phase ---
!!$       ms%CAP = ms%len/2.0_PS
!!$    else if( phase == 2 ) then
!!$       ! --- in case of solid phase ---
!!$       ! +++ first caclulate for a hexagonal monocrystal +++
!!$       if( is%phi_ic < 1.0_PS ) then
!!$          ! --- in case of oblate spheroids ---
!!$          d = sqrt( ms%a_len**2.0 - ms%c_len**2.0)
!!$          eps = sqrt( 1.0_PS - is%phi_ic**2.0)
!!$          ms%CAP_hex = d/asin(eps)
!!$       else if( is%phi_ic> 1.0_PS ) then
!!$          ! --- in case of prolate spheroids ---
!!$          d = sqrt( ms%c_len**2.0 - ms%a_len**2.0)
!!$          eps = sqrt( 1.0_PS - is%phi_ic**(-2.0))
!!$          ms%CAP_hex = d/log((1.0_PS + eps)*is%phi_ic )
!!$!!c             ms%CAP_hex = ( 0.708_PS + 0.615_PS * (is%phi_ic**(-0.76)))*ms%c_len
!!$       else
!!$          ! --- in case of sphere ---
!!$          ms%CAP_hex = ms%a_len
!!$       end if
!!$
!!$
!!$       ! ++++ second calculate the capacitance for polycrystals ++++
!!$       if( is%sh_type <= 2 ) then
!!$          ! --- case of ice crystals ---
!!$          if(is%is_mod(2)==1) then
!!$             if(is%habit<=3) then
!!$                ms%CAP=ms%CAP_hex
!!$             elseif(is%habit==4) then
!!$                ms%CAP = CAP_ros*is%r
!!$             elseif(is%habit==5) then
!!$                ! assume phi=0.25 with a spheroid
!!$                eps=0.968245837_PS
!!$                d=is%e*eps
!!$                ms%CAP=d/asin(eps)
!!$             else
!!$                ! assume phi=0.25 with a spheroid
!!$                eps=0.968245837_PS
!!$                d=max(is%e,is%r)*eps
!!$                ms%CAP=d/asin(eps)
!!$   !!c             ms%CAP = max(is%r,is%e)
!!$             end if
!!$          else
!!$             ! --- in case of sphere ---
!!$             ms%CAP = ms%semi_a
!!$          endif
!!$       elseif( is%sh_type <= 3 ) then
!!$          ! --- case of aggregates and rimed aggregates ---
!!$          ! use Capacitance/D_max = 0.25 derived by
!!$          !      Westbrook et al. (2008), JAS
!!$          ms%CAP=(is%V_cs/coefpi6)**(1.0/3.0)*0.25_PS
!!$       else
!!$          ! --- case of graupels ---
!!$          phi = is%phi_cs
!!$          if( phi < 1.0_PS ) then
!!$             ! --- in case of oblate spheroids ---
!!$             d = sqrt( ms%semi_a**2.0 - ms%semi_c**2.0)
!!$             eps = sqrt( 1.0_PS - phi**2.0)
!!$             ms%CAP = d/asin(eps)
!!$          else if( phi > 1.0_PS ) then
!!$             ! --- in case of prolate spheroids ---
!!$             d = sqrt( ms%semi_c**2.0 - ms%semi_a**2.0)
!!$             eps = sqrt( 1.0_PS - phi**(-2.0))
!!$             ms%CAP = d/log((1.0_PS + eps)*phi )
!!$!!c             ms%CAP = ( 0.708_PS + 0.615_PS * (phi**(-0.76)))*ms%c_len
!!$          else
!!$             ! --- in case of sphere ---
!!$             ms%CAP = ms%semi_a
!!$          end if
!!$       end if
!!$       if(level<=5.and.th_var%T>=T_0) then
!!$          ms%CAP=max(ms%CAP,(ms%mean_mass/coef4pi3)**(1.0/3.0))
!!$       elseif(th_var%T>=T_0) then
!!$          ms%CAP=max(ms%CAP,(ms%mass(imw)/ms%con/coef4pi3)**(1.0/3.0))
!!$       end if
!!$
!!$       ! +++ additional factor to represent the shape +++
!!$       ! under construction
!!$    end if
!!$  end subroutine cal_capacitance

  subroutine cal_capacitance_vec(level,phase,g,ag,icond1)
    ! +++ Calculate capacitance based on the spheroid assumption,
    !     for example, see Chen and Lamb (1994a).                     +++
    ! phase of hydrometeor
    ! 1: liquid, 2: solid
    integer, intent(in)  :: phase
    integer, intent(in)  :: level
    type (Group), intent(inout)  :: g
    type (AirGroup), intent(in)   :: ag
    integer,dimension(g%N_BIN,g%L),intent(in)    ::  icond1
    !
    real(PS)                    :: eps, d, phi!, out
    ! coefficient for capacitance by Chiruta and Wang (2003)
!!c    real(PS),parameter   :: CAP_ros=0.434*N_lob**0.257
    real(PS),parameter   :: CAP_ros=0.619753727

    integer                     :: i,n

    if( phase == 1 ) then
       ! --- in case of liquid phase ---
!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then
          g%MS(i,n)%CAP = 0.5_PS*g%MS(i,n)%len
        endif
      enddo
      enddo
    else if( phase == 2 ) then
      ! --- in case of solid phase ---
      ! +++ first caclulate for a hexagonal monocrystal +++
!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then

          if( g%IS(i,n)%phi_ic < 1.0_PS ) then
            ! --- in case of oblate spheroids ---
            d = sqrt( max(g%MS(i,n)%a_len**2 - g%MS(i,n)%c_len**2, 0.0_PS) )
            eps = sqrt( 1.0_PS - g%IS(i,n)%phi_ic**2)
            g%MS(i,n)%CAP_hex = d/asin(eps)
          else if( g%IS(i,n)%phi_ic> 1.0_PS ) then
            ! --- in case of prolate spheroids ---
            d = sqrt( max(g%MS(i,n)%c_len**2 - g%MS(i,n)%a_len**2, 0.0_PS) )
            eps = sqrt( 1.0_PS - 1.0_PS/g%IS(i,n)%phi_ic**2)
            g%MS(i,n)%CAP_hex = d/log((1.0_PS + eps)*g%IS(i,n)%phi_ic )
!!c             g%MS(i,n)%CAP_hex = ( 0.708_PS + 0.615_PS * (g%IS(i,n)%phi_ic**(-0.76)))*g%MS(i,n)%c_len
          else
            ! --- in case of sphere ---
            g%MS(i,n)%CAP_hex = g%MS(i,n)%a_len
          end if


          ! ++++ second calculate the capacitance for polycrystals ++++
          if( g%IS(i,n)%sh_type <= 2 ) then
            ! --- case of ice crystals ---
            if(g%IS(i,n)%is_mod(2)==1) then
              if(g%IS(i,n)%habit<=3) then
                g%MS(i,n)%CAP=g%MS(i,n)%CAP_hex
              elseif(g%IS(i,n)%habit==4) then
                g%MS(i,n)%CAP = CAP_ros*g%IS(i,n)%r
              elseif(g%IS(i,n)%habit==5) then
                ! assume phi=0.25 with a spheroid
                eps=0.968245837_PS
                d=g%IS(i,n)%e*eps
                g%MS(i,n)%CAP=d/asin(eps)
              else
                ! assume phi=0.25 with a spheroid
                eps=0.968245837_PS
                d=max(g%IS(i,n)%e,g%IS(i,n)%r)*eps
                g%MS(i,n)%CAP=d/asin(eps)
   !!c             g%MS(i,n)%CAP = max(g%IS(i,n)%r,g%IS(i,n)%e)
              end if
            else
              ! --- in case of sphere ---
              g%MS(i,n)%CAP = g%MS(i,n)%semi_a
            endif
          elseif( g%IS(i,n)%sh_type <= 3 ) then
            ! --- case of aggregates and rimed aggregates ---
            ! use Capacitance/D_max = 0.25 derived by
            !      Westbrook et al. (2008), JAS
            g%MS(i,n)%CAP=(g%IS(i,n)%V_cs/coefpi6)**(1.0/3.0)*0.25_PS
          else
            ! --- case of graupels ---
            phi = g%IS(i,n)%phi_cs
            if( phi < 1.0_PS ) then
              ! --- in case of oblate spheroids ---
              d = sqrt( g%MS(i,n)%semi_a**2.0 - g%MS(i,n)%semi_c**2.0)
              eps = sqrt( 1.0_PS - phi**2.0)
              g%MS(i,n)%CAP = d/asin(eps)
            else if( phi > 1.0_PS ) then
              ! --- in case of prolate spheroids ---
              d = sqrt( g%MS(i,n)%semi_c**2.0 - g%MS(i,n)%semi_a**2.0)
              eps = sqrt( 1.0_PS - phi**(-2.0))
              g%MS(i,n)%CAP = d/log((1.0_PS + eps)*phi )
!!c             g%MS(i,n)%CAP = ( 0.708_PS + 0.615_PS * (phi**(-0.76)))*g%MS(i,n)%c_len
            else
              ! --- in case of sphere ---
              g%MS(i,n)%CAP = g%MS(i,n)%semi_a
            end if
          end if
          if(level<=5.and.ag%TV(n)%T>=T_0) then
            g%MS(i,n)%CAP=max(g%MS(i,n)%CAP,(g%MS(i,n)%mean_mass/coef4pi3)**(1.0/3.0))
          else
            if(ag%TV(n)%T>=T_0) then
              g%MS(i,n)%CAP=max(g%MS(i,n)%CAP,(g%MS(i,n)%mass(imw)/g%MS(i,n)%con/coef4pi3)**(1.0/3.0))
            endif
          end if

          ! +++ additional factor to represent the shape +++
          ! under construction
        endif
      enddo
      enddo
    end if
  end subroutine cal_capacitance_vec

!!$  subroutine cal_ventilation_coef( phase, ms, th_var,ishape)
!!$    use class_Thermo_Var, only: &
!!$       get_fkn
!!$    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! Calculate ventilation coefficient
!!$    ! based on empirical expressions by Hall and Pruppacher (1976).
!!$    !
!!$    ! Note:
!!$    ! * Reynolds number is recalculated using a perimeter for
!!$    !   characteristic length.
!!$    ! * Characteristic length is needed to be defined beforehand.
!!$    ! * Heat capacity is assumed to be constant.
!!$    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! phase of hydrometeor
!!$    ! 1: liquid, 2: solid
!!$    integer, intent(in)  :: phase
!!$    type (Mass_Bin), intent(inout) :: ms
!!$    type (Thermo_Var), intent(in) :: th_var
!!$    type (Ice_Shape), optional       :: ishape
!!$
!!$    ! Schmidt number
!!$    real(PS)                   :: N_sc
!!$    ! Nusselt number
!!$    real(PS)                   :: N_ns
!!$    real(PS)                   :: Xv, Xh
!!$    ! major semiaxis
!!$    real(PS)                   :: r_m
!!$
!!$    ! initialize
!!$    ms%Nre=0.0_PS
!!$    ms%fv=1.0_PS
!!$    ms%fh=1.0_PS
!!$    ms%fkn=1.0_PS
!!$    ms%fac=1.0_PS
!!$    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return
!!$
!!$    ! +++ calculate the Schmidt number +++
!!$    N_sc = th_var%d_vis/th_var%den/th_var%D_v
!!$
!!$    ! +++ calculate the Nusselt or Prandtl number +++
!!$    ! +++
!!$!!c  need to check
!!$!!c    N_ns = th_var%d_vis/C_pa/th_var%k_a
!!$    N_ns = th_var%d_vis/th_var%den/th_var%k_a
!!$!!c    N_ns = 0.72_PS
!!$
!!$    if( phase == 1 ) then
!!$       ! --- in case of liquid phase ---
!!$       ! assume the spherical water drops
!!$       ! +++ calculate Reynolds number +++
!!$       ms%Nre = ms%len*ms%vtm*th_var%den&
!!$            /th_var%d_vis
!!$
!!$       ! +++ calculate dimensionless number X +++
!!$       Xv = (N_sc**(1.0/3.0))*sqrt(ms%Nre)
!!$       Xh = (N_ns**(1.0/3.0))*sqrt(ms%Nre)
!!$
!!$       if( Xv < 1.4_PS ) then
!!$          ms%fv = 1.0_PS + 0.108_PS*(Xv**2.0)
!!$       else if( 1.4_PS <= Xv .and. Xv <= 51.4_PS ) then
!!$          ms%fv = 0.78_PS + 0.308_PS*Xv
!!$       else
!!$          ms%fv = 0.78_PS + 0.308_PS*51.4_PS
!!$       end if
!!$
!!$!!c       ms%fv=1.0_PS
!!$
!!$       if( Xh < 1.4_PS ) then
!!$          ms%fh = 1.0_PS + 0.108_PS*(Xh**2.0)
!!$       else if( 1.4_PS <= Xh .and. Xh <= 51.4_PS ) then
!!$          ms%fh = 0.78_PS + 0.308_PS*Xh
!!$       else
!!$          ms%fh = 0.78_PS + 0.308_PS*51.4_PS
!!$       end if
!!$
!!$
!!$       ! +++ calculate kinetic effect +++
!!$       r_m=ms%len*0.5_PS
!!$       ms%fkn = get_fkn(th_var,phase,r_m)
!!$
!!$    else if( phase == 2 ) then
!!$       ! --- in case of ice phase ---
!!$       ! assume the ice crystal model
!!$       if(.not.present(ishape)) then
!!$          write(*,*) " cal_ventilation > Need ice shape object!"
!!$          stop
!!$       end if
!!$       ! +++ calculate Reynolds number +++
!!$       ms%Nre = ms%len*ms%vtm*th_var%den&
!!$            /th_var%d_vis
!!$
!!$       ! +++ calculate dimensionless number X +++
!!$       Xv = (N_sc**(1.0/3.0))*sqrt(ms%Nre)
!!$       Xh = (N_ns**(1.0/3.0))*sqrt(ms%Nre)
!!$
!!$       if(ishape%sh_type<=2) then
!!$          if( Xv < 1.0_PS ) then
!!$             ms%fac = (1.0_PS + 0.14_PS*(Xv**2.0)*&
!!$                  sqrt(ms%c_len/ms%len/2.0_PS))/&
!!$                  (1.0_PS + 0.14_PS*(Xv**2.0)*&
!!$                  sqrt(ms%a_len/ms%len/2.0_PS))
!!$
!!$          else if( 1.0_PS <= Xv ) then
!!$             ms%fac = (0.86_PS + 0.28_PS*Xv*&
!!$                  sqrt(ms%c_len/ms%len/2.0_PS))/&
!!$                  (0.86_PS + 0.28_PS*Xv*&
!!$                  sqrt(ms%a_len/ms%len/2.0_PS))
!!$          end if
!!$       else
!!$          ! assume crystals within aggregates, rimed aggregates, graupel
!!$          ! won't have ventilation effect on the shape.
!!$          ms%fac = 1.0_PS
!!$!tmp          if( Xv < 1.0_PS ) then
!!$!tmp             ms%fac = (1.0_PS + 0.14_PS*(Xv**2.0)*&
!!$!tmp                  sqrt(ms%semi_c/ms%len/2.0_PS))/&
!!$!tmp                  (1.0_PS + 0.14_PS*(Xv**2.0)*&
!!$!tmp                  sqrt(ms%semi_a/ms%len/2.0_PS))
!!$!tmp
!!$!tmp          else if( 1.0_PS <= Xv ) then
!!$!tmp             ms%fac = (0.86_PS + 0.28_PS*Xv*&
!!$!tmp                  sqrt(ms%semi_c/ms%len/2.0_PS))/&
!!$!tmp                  (0.86_PS + 0.28_PS*Xv*&
!!$!tmp                  sqrt(ms%semi_a/ms%len/2.0_PS))
!!$!tmp          end if
!!$       end if
!!$
!!$       if(ishape%sh_type<=2.and.ishape%habit==4) then
!!$          ! if the ice crystal is bullet rosettes
!!$          ! use Lie et al (2003)'s parameterization.
!!$          ms%fv = 1.0_PS + 0.3005_PS*Xv - 0.0022_PS*Xv*Xv
!!$          ms%fh = 1.0_PS + 0.3005_PS*Xh - 0.0022_PS*Xh*Xh
!!$       else
!!$          ! if the ice crystal is hexagonal plates, dendrites, or columns,
!!$          ! if it is polycrystal,
!!$          ! or if it is aggregate,
!!$          if( Xv < 1.0_PS ) then
!!$             ms%fv = 1.0_PS + 0.14_PS*(Xv**2.0)
!!$          else if( 1.0_PS <= Xv ) then
!!$             ms%fv = 0.86_PS + 0.28_PS*Xv
!!$          end if
!!$          if( Xh < 1.0_PS ) then
!!$             ms%fh = 1.0_PS + 0.14_PS*(Xh**2.0)
!!$          else if( 1.0_PS <= Xh ) then
!!$             ms%fh = 0.86_PS + 0.28_PS*Xh
!!$          end if
!!$       end if
!!$
!!$       ! +++ calculate kinetic effect +++
!!$       ms%fkn=1.0_PS
!!$
!!$       r_m = sqrt(ms%semi_a**2+ms%semi_c**2)
!!$       ms%fkn=get_fkn(th_var,phase,r_m)
!!$
!!$       if(ms%inmlt==1.and.th_var%T>=T_0) then
!!$          if( Xh < 1.4_PS ) then
!!$             ms%fh = max(ms%fh,1.0_PS + 0.108_PS*(Xh**2.0))
!!$          else if( 1.4_PS <= Xh .and. Xh <= 51.4_PS ) then
!!$             ms%fh = max(ms%fh,0.78_PS + 0.308_PS*Xh)
!!$          else
!!$             ms%fh = max(ms%fh,0.78_PS + 0.308_PS*51.4_PS)
!!$          end if
!!$       end if
!!$    end if
!!$
!!$  END subroutine cal_ventilation_coef

  subroutine cal_ventilation_coef_vec(phase,g,ag,icond1)
    use class_Thermo_Var, only: &
       get_fkn
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Calculate ventilation coefficient
    ! based on empirical expressions by Hall and Pruppacher (1976).
    !
    ! Note:
    ! * Reynolds number is recalculated using a perimeter for
    !   characteristic length.
    ! * Characteristic length is needed to be defined beforehand.
    ! * Heat capacity is assumed to be constant.
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! phase of hydrometeor
    ! 1: liquid, 2: solid
    integer, intent(in)  :: phase
    type (Group), intent(inout)  :: g
    type (AirGroup), intent(in)   :: ag
    integer,dimension(g%N_BIN,g%L),intent(in)    ::  icond1

    ! Schmidt number
    real(PS)                   :: N_sc
    ! Nusselt number
    real(PS)                   :: N_ns
    real(PS)                   :: Xv, Xh
    ! major semiaxis
    real(PS)                   :: r_m

    integer :: i,n

    if( phase == 1 ) then
      ! --- in case of liquid phase ---

!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then

          ! +++ calculate the Schmidt number +++
          N_sc = ag%TV(n)%d_vis/ag%TV(n)%den/ag%TV(n)%D_v

          ! +++ calculate the Nusselt or Prandtl number +++
          ! +++
          N_ns = ag%TV(n)%d_vis/ag%TV(n)%den/ag%TV(n)%k_a

          ! assume the spherical water drops
          ! +++ calculate Reynolds number +++
          g%MS(i,n)%Nre = g%MS(i,n)%len*g%MS(i,n)%vtm*ag%TV(n)%den&
                          /ag%TV(n)%d_vis

          ! +++ calculate dimensionless number X +++
          Xv = (N_sc**(1.0/3.0))*sqrt(g%MS(i,n)%Nre)
          Xh = (N_ns**(1.0/3.0))*sqrt(g%MS(i,n)%Nre)

          if( Xv < 1.4_PS ) then
            g%MS(i,n)%fv = 1.0_PS + 0.108_PS*(Xv**2.0)
          else if( 1.4_PS <= Xv .and. Xv <= 51.4_PS ) then
            g%MS(i,n)%fv = 0.78_PS + 0.308_PS*Xv
          else
            g%MS(i,n)%fv = 0.78_PS + 0.308_PS*51.4_PS
          end if

          if( Xh < 1.4_PS ) then
            g%MS(i,n)%fh = 1.0_PS + 0.108_PS*(Xh**2.0)
          else if( 1.4_PS <= Xh .and. Xh <= 51.4_PS ) then
            g%MS(i,n)%fh = 0.78_PS + 0.308_PS*Xh
          else
            g%MS(i,n)%fh = 0.78_PS + 0.308_PS*51.4_PS
          end if


          ! +++ calculate kinetic effect +++
          r_m=g%MS(i,n)%len*0.5_PS
          g%MS(i,n)%fkn = get_fkn(ag%TV(n),phase,r_m)
        endif
      enddo
      enddo

    else if( phase == 2 ) then

      ! --- in case of ice phase ---
!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then

          ! +++ calculate the Schmidt number +++
          N_sc = ag%TV(n)%d_vis/ag%TV(n)%den/ag%TV(n)%D_v

          ! +++ calculate the Nusselt or Prandtl number +++
          ! +++
          N_ns = ag%TV(n)%d_vis/ag%TV(n)%den/ag%TV(n)%k_a

          ! +++ calculate Reynolds number +++
          g%MS(i,n)%Nre = g%MS(i,n)%len*g%MS(i,n)%vtm*ag%TV(n)%den&
                         /ag%TV(n)%d_vis

          ! +++ calculate dimensionless number X +++
          Xv = (N_sc**(1.0/3.0))*sqrt(g%MS(i,n)%Nre)
          Xh = (N_ns**(1.0/3.0))*sqrt(g%MS(i,n)%Nre)

          if(g%IS(i,n)%sh_type<=2) then
            if( Xv < 1.0_PS ) then
              g%MS(i,n)%fac = (1.0_PS + 0.14_PS*(Xv**2.0)*&
                  sqrt(g%MS(i,n)%c_len/g%MS(i,n)%len/2.0_PS))/&
                  (1.0_PS + 0.14_PS*(Xv**2.0)*&
                  sqrt(g%MS(i,n)%a_len/g%MS(i,n)%len/2.0_PS))

            else
              g%MS(i,n)%fac = (0.86_PS + 0.28_PS*Xv*&
                  sqrt(g%MS(i,n)%c_len/g%MS(i,n)%len/2.0_PS))/&
                  (0.86_PS + 0.28_PS*Xv*&
                  sqrt(g%MS(i,n)%a_len/g%MS(i,n)%len/2.0_PS))
            end if
          else
            ! assume crystals within aggregates, rimed aggregates, graupel
            ! won't have ventilation effect on the shape.
            g%MS(i,n)%fac = 1.0_PS
!tmp          if( Xv < 1.0_PS ) then
!tmp             g%MS(i,n)%fac = (1.0_PS + 0.14_PS*(Xv**2.0)*&
!tmp                  sqrt(g%MS(i,n)%semi_c/g%MS(i,n)%len/2.0_PS))/&
!tmp                  (1.0_PS + 0.14_PS*(Xv**2.0)*&
!tmp                  sqrt(g%MS(i,n)%semi_a/g%MS(i,n)%len/2.0_PS))
!tmp
!tmp          else if( 1.0_PS <= Xv ) then
!tmp             g%MS(i,n)%fac = (0.86_PS + 0.28_PS*Xv*&
!tmp                  sqrt(g%MS(i,n)%semi_c/g%MS(i,n)%len/2.0_PS))/&
!tmp                  (0.86_PS + 0.28_PS*Xv*&
!tmp                  sqrt(g%MS(i,n)%semi_a/g%MS(i,n)%len/2.0_PS))
!tmp          end if
          end if

          if(g%IS(i,n)%sh_type<=2.and.g%IS(i,n)%habit==4) then
            ! if the ice crystal is bullet rosettes
            ! use Lie et al (2003)'s parameterization.
            g%MS(i,n)%fv = 1.0_PS + 0.3005_PS*Xv - 0.0022_PS*Xv*Xv
            g%MS(i,n)%fh = 1.0_PS + 0.3005_PS*Xh - 0.0022_PS*Xh*Xh
          else
            ! if the ice crystal is hexagonal plates, dendrites, or columns,
            ! if it is polycrystal,
            ! or if it is aggregate,
            if( Xv < 1.0_PS ) then
              g%MS(i,n)%fv = 1.0_PS + 0.14_PS*(Xv**2.0)
            else
              g%MS(i,n)%fv = 0.86_PS + 0.28_PS*Xv
            end if
            if( Xh < 1.0_PS ) then
              g%MS(i,n)%fh = 1.0_PS + 0.14_PS*(Xh**2.0)
            else
              g%MS(i,n)%fh = 0.86_PS + 0.28_PS*Xh
            end if
          end if

          ! +++ calculate kinetic effect +++
          g%MS(i,n)%fkn=1.0_PS

          r_m = sqrt(g%MS(i,n)%semi_a**2+g%MS(i,n)%semi_c**2)
          g%MS(i,n)%fkn=get_fkn(ag%TV(n),phase,r_m)

          if(g%MS(i,n)%inmlt==1.and.ag%TV(n)%T>=T_0) then
            if( Xh < 1.4_PS ) then
              g%MS(i,n)%fh = max(g%MS(i,n)%fh,1.0_PS + 0.108_PS*(Xh**2.0))
            else if( 1.4_PS <= Xh .and. Xh <= 51.4_PS ) then
              g%MS(i,n)%fh = max(g%MS(i,n)%fh,0.78_PS + 0.308_PS*Xh)
            else
              g%MS(i,n)%fh = max(g%MS(i,n)%fh,0.78_PS + 0.308_PS*51.4_PS)
            end if
          end if
        endif
      enddo
      enddo
    end if

  END subroutine cal_ventilation_coef_vec

!!$  subroutine cal_terminal_vel( phase, ms, th_var, mode, level,em, ishape)
!!$    implicit none
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! Calculate
!!$    ! 1. Characteristic length based on total surface area and
!!$    !    perimeter normal to the flow direction,
!!$    ! 2. Reynolds number and then,
!!$    ! 3. terminal velocity based on the spheroid assumption,
!!$    ! see Bohm (1992).
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    !
!!$    ! phase of hydrometeor
!!$    ! 1: liquid, 2: solid
!!$    integer, intent(in)  :: phase
!!$    type (Mass_Bin), intent(inout)   :: ms
!!$    type (Ice_Shape), optional       :: ishape
!!$    type (Thermo_Var), intent(in)    :: th_var
!!$    integer                          :: mode
!!$    integer, intent(in)              :: level
!!$
!!$    ! effective cross-sectional area
!!$    real(PS)                    :: A_e
!!$    ! circumscribed cross-sectional area
!!$    real(PS)                    :: A_c
!!$    ! the ratio of the effective over the circumscribed cross-sectional area, Ae/A
!!$    real(PS)                    :: q
!!$    ! charasteristic length defined for Bohm (1992)
!!$    real(PS)                    :: char_len
!!$
!!$    ! axial ratio, c/(2*(a+d))
!!$    real(PS)                    :: alpha
!!$    ! shape factor
!!$    real(PS)                    :: k
!!$    ! the Davies or Best number
!!$    real(PS)                    :: X, X_0, X_old
!!$    ! pressure drag coefficient
!!$    real(PS)                    :: C_DP
!!$    !
!!$    real(PS)                    :: beta, gamma, xx
!!$    ! Oseen-type solution for the drag coefficient
!!$    real(PS)                    :: C_DO
!!$    ! Reynolds number
!!$    real(PS)                    :: N_re_0
!!$
!!$    ! total surface area of the body
!!$    real(PS)                    :: OMEGA
!!$    ! Perimeter
!!$    real(PS)                    :: P
!!$
!!$    ! Stoke's terminal velocity
!!$    real(PS) :: U_S
!!$    ! radius of a drop
!!$    real(PS) :: rad
!!$    ! Bond number and physical property number
!!$    real(PS) :: NBO,NP
!!$    real(PS) :: Y,lambda_a
!!$
!!$    real(PS)  :: beta0
!!$    integer                    :: i
!!$    integer,intent(inout) :: em
!!$
!!$    ! initialize
!!$    ms%vtm=0.0_PS
!!$    ms%Nre=0.0_PS
!!$    if(phase==2) then
!!$       ms%len=0.0_PS
!!$    endif
!!$    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return
!!$
!!$    if( phase == 1 ) then
!!$       ! +++ calculate the mass component for one particle +++
!!$       rad=ms%len*0.5_PS
!!$
!!$
!!$       if(rad<0.5e-4_PS) then  ! bug found 2016/05/30
!!$          ms%vtm=0.0_PS
!!$       elseif(0.5e-4_PS<=rad.and.rad<10.0e-4_PS) then
!!$          ! calculate Stokes terminal velocity
!!$          U_S=rad**2.0*gg*(den_w-th_var%den_a)/4.5_PS/th_var%d_vis
!!$          lambda_a=6.6e-6_PS*(th_var%d_vis/1.818e-4_PS)*(1013250.0_PS/th_var%P)*(th_var%T/293.15_PS)
!!$          ms%vtm=(1.0_PS+1.26_PS*lambda_a/rad)*U_S
!!$
!!$       elseif(rad<535.0e-4_PS) then
!!$          X=log(32.0_PS*rad**3.0*(den_w-th_var%den_a)*th_var%den_a*gg/3.0_PS/th_var%d_vis**2.0)
!!$
!!$
!!$!          Y=-0.318657e+1_PS+0.992696_PS*X-0.153193e-2_PS*X**2.0-0.987059e-3_PS*X**3.0&
!!$!               -0.578878e-3_PS*X**4.0+0.855176e-4*X**5.0-0.327815e-5*X**6.0
!!$
!!$          Y=-0.318657e+1_PS+0.992696_PS*X-0.153193e-2_PS*X*X-0.987059e-3_PS*X*X*X &
!!$               -0.578878e-3_PS*X*X*X*X+0.855176e-4*X*X*X*X*X &
!!$               -0.327815e-5*X*X*X*X*X*X
!!$          ms%Nre=exp(Y)
!!$          ! +++ calculate the terminal velocity +++
!!$          ms%vtm=ms%Nre*th_var%d_vis/(2.0_PS*rad*th_var%den_a)
!!$
!!$       else
!!$          rad=min(rad,3500.0e-4_RP)
!!$
!!$          NBO=gg*(den_w-th_var%den_a)*rad**2/th_var%sig_wa
!!$          NP=th_var%sig_wa**3*th_var%den_a**2/th_var%d_vis**4/gg
!!$
!!$          X=log(NBO*NP**(1.0_RP/6.0_RP)*16.0_PS/3.0_PS)
!!$!          Y=-0.500015e+1_PS+0.523778e+1_PS*X-0.204914e+1_PS*X**2.0+0.475294_PS*X**3.0&
!!$!               -0.542819e-1_PS*X**4.0+0.238449e-2_PS*X**5.0
!!$          Y=-0.500015e+1_PS+0.523778e+1_PS*X-0.204914e+1_PS*X*X+0.475294_PS*X*X*X &
!!$               -0.542819e-1_PS*X*X*X*X+0.238449e-2_PS*X*X*X*X*X
!!$
!!$
!!$          ms%Nre=NP**(1.0_RP/6.0_RP)*exp(Y)
!!$          ! +++ calculate the terminal velocity +++
!!$          ms%vtm=ms%Nre*th_var%d_vis/(2.0_PS*rad*th_var%den_a)
!!$
!!$
!!$       end if
!!$
!!$    else if( phase == 2 ) then
!!$       ! --- in case of ice phase ---
!!$       if(.not.present(ishape)) then
!!$          write(*,*) " cal_terminal_vel > Need ice shape object!"
!!$          stop
!!$       end if
!!$
!!$       ! +++ calculate some important variables, depending on level +++
!!$       call cal_some_vel(ms,ishape,mode,alpha,OMEGA,A_e,A_c,P,q,char_len)
!!$       ms%len = OMEGA/P
!!$       ishape%q_e=q
!!$
!!$!!c       if(ishape%sh_type>=3) then
!!$!!c          open( unit=21, file="check_svtm.dat",POSITION='APPEND')
!!$!!c          write(21,'(10ES15.6)') ms%mean_mass,ms%semi_a,ms%semi_c,alpha,q,char_len,ms%Nre
!!$!!c          close(21)
!!$!!c       end if
!!$
!!$
!!$
!!$       ! +++ calculate the shape factor +++
!!$       if( alpha <= 0.16_PS ) then
!!$          k = 0.85_PS
!!$       else if ( 0.16_PS < alpha .and. alpha <= 1.0_PS ) then
!!$          k = 0.82_PS + 0.18_PS*alpha
!!$       else if ( 1.0_PS < alpha .and. alpha <= 3.0_PS ) then
!!$          k = 0.37_PS + 0.63_PS/alpha
!!$       else if ( 3.0_PS < alpha ) then
!!$          k = 1.33_PS/( log(alpha) + 1.19_PS)
!!$       end if
!!$
!!$       ! +++ calculate the Davies or Best number +++
!!$       if( q <= 1.0_PS ) then
!!$          X = 8.0_PS*ms%mean_mass*gg*th_var%den_a/&
!!$               (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*(q**0.25))
!!$       else if( q > 1.0_PS ) then
!!$          X = 8.0_PS*ms%mean_mass*gg*th_var%den_a/&
!!$               (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*q)
!!$       end if
!!$
!!$       ! +++ calculate the pressure drag C_DP +++
!!$       if( 0.3_PS < alpha .and. alpha < 1.0_PS ) then
!!$          gamma = 3.76_PS - 8.41*alpha + 9.18*alpha**2.0 - 3.53*alpha**3.0
!!$          C_DP = 0.292*k*gamma
!!$       else if( alpha <= 0.3_PS ) then
!!$          gamma = 1.98
!!$          C_DP = 0.292*k*gamma
!!$!!c          else if( 1.0_PS <= alpha .and. 1.0_PS < q ) then
!!$       else if( 1.0_PS <= alpha ) then
!!$          if( 1.0_PS <= q ) then
!!$             C_DP = (1.46*q-0.46)*(0.492-0.200/sqrt(alpha))
!!$          else
!!$!!c             gamma = 3.76_PS - 8.41*alpha + 9.18*alpha**2.0 - 3.53*alpha**3.0
!!$             gamma = 1.0
!!$             C_DP = 0.292*k*gamma
!!$          end if
!!$       end if
!!$!!c       write(*,'(I3,5ES14.5)') i,k,gamma,alpha,q,C_DP
!!$
!!$       ! +++ calculate the Oseen-type drag coefficient, C_DO +++
!!$       C_DO = 4.5*(k**2.0)*max(alpha, 1.0_PS)
!!$
!!$       ! +++ calculate the Reynolds number +++
!!$       beta = sqrt(1.0_PS + (C_DP/(6.0_PS*k))*sqrt(X/C_DP)) - 1.0_PS
!!$       beta0=beta
!!$       N_re_0 = 6.0*k*beta**2.0/C_DP
!!$
!!$       ! +++ correction for transition to turbulent flow +++
!!$       if( 1000.0 <= N_re_0 ) then
!!$          X_0 = 2.8e+06
!!$          C_DP=C_DP*(1.0_PS + 1.6_PS*(X/X_0)**2.0)/(1.0_PS + (X/X_0)**2.0)
!!$
!!$          beta = sqrt(1.0_PS + (C_DP/(6.0_PS*k))*sqrt(X/C_DP)) - 1.0_PS
!!$          N_re_0 = 6.0*k*beta**2.0/C_DP
!!$
!!$       end if
!!$
!!$       ! +++ matching theory +++
!!$       gamma = (C_DO-C_DP)/(4.0_PS*C_DP)
!!$       ms%Nre = N_re_0*(1.0_PS+2.0_PS*beta*exp(-beta*gamma)/&
!!$            ((2.0_PS+beta)*(1.0_PS+beta)))
!!$       if( ms%Nre < 0.0_PS ) then
!!$          write(*,*) "Reynolds number is negative!!"
!!$          write(*,'(10ES14.6)') ms%mass,ms%con,ms%a_len,ms%c_len,q,k,gamma
!!$          write(*,'(10ES14.6)') N_re_0,exp(-beta*gamma),beta,X,C_DP
!!$          write(*,'(2I5,10ES14.6)') ishape%sh_type,ishape%habit,ms%semi_a,ms%semi_c,ishape%phi_cs
!!$          write(*,'(10ES14.6)') ishape%semi_aip,ishape%semi_cip,ishape%V_cs
!!$          em=2
!!$          return
!!$       end if
!!$
!!$       ! +++ calculate the terminal velocity +++
!!$       ms%vtm = ms%Nre*th_var%d_vis/(char_len*th_var%den_a)
!!$       if(ms%vtm>5.0e+3.and.ms%den>0.8_PS) then
!!$          write(*,*) "terminal_vel > The ice vel is set to 50.0m/s"
!!$          write(*,*) "old vel",ms%vtm
!!$          write(*,'("mass,con,alen,clen,q,k,gamma",10ES14.6)') ms%mass(1),ms%con,ms%a_len,ms%c_len,q,k,gamma
!!$          write(*,'("N_re0,ex,beta,X,C_DP",10ES14.6)') N_re_0,exp(-beta*gamma),beta,X,C_DP
!!$          write(*,'("semi_a,semi_c,mean_mass,beta0",10ES14.6)') ms%semi_a,ms%semi_c,ms%mean_mass,beta0
!!$          write(*,'("char_len,d_vis,den,Nre,A_c",10ES14.6)') char_len,th_var%d_vis,th_var%den,ms%Nre,A_c
!!$          write(*,'("pressure,T",10ES14.6)') th_var%P,th_var%T
!!$          ms%vtm=5.0e+3_PS
!!$       elseif(ms%vtm>2.0e+3) then
!!$          ms%vtm=2.0e+3_PS
!!$       elseif(ms%vtm<0.0) then
!!$          write(*,*) "terminal_vel > The ice vel is set to 0.0m/s"
!!$          write(*,*) "old vel",ms%vtm
!!$          write(*,'("mass,con,alen,clen,q,k,gamma",10ES14.6)') ms%mass(1),ms%con,ms%a_len,ms%c_len,q,k,gamma
!!$          write(*,'("N_re0,ex,beta,X,C_DP",10ES14.6)') N_re_0,exp(-beta*gamma),beta,X,C_DP
!!$          write(*,'("semi_a,semi_c,mean_mass,beta0",10ES14.6)') ms%semi_a,ms%semi_c,ms%mean_mass,beta0
!!$          write(*,'("char_len,d_vis,den,Nre,A_c",10ES14.6)') char_len,th_var%d_vis,th_var%den,ms%Nre,A_c
!!$          write(*,'("pressure,T",10ES14.6)') th_var%P,th_var%T
!!$          ms%vtm=0.0_PS
!!$       endif
!!$    end if
!!$
!!$  end subroutine cal_terminal_vel

  subroutine cal_terminal_vel_vec(phase,g,ag,icond1,mode)
    use scale_prc, only: &
       PRC_abort
    use class_Ice_Shape, only: &
       get_total_sfc_area, &
       get_effect_area, &
       get_circum_area, &
       get_perimeter, &
       get_vip
    use class_Mass_Bin, only: &
       get_aspect_ratio_sol, &
       get_area_ratio_sol
    implicit none
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Calculate
    ! 1. Characteristic length based on total surface area and
    !    perimeter normal to the flow direction,
    ! 2. Reynolds number and then,
    ! 3. terminal velocity based on the spheroid assumption,
    ! see Bohm (1992).
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    integer, intent(in)  :: phase
    type (Group), intent(inout)  :: g
    type (AirGroup), intent(in)   :: ag
    integer,dimension(g%N_BIN,g%L),intent(in)    ::  icond1
    integer                          :: mode

    ! effective cross-sectional area
    real(PS)                    :: A_e
    ! circumscribed cross-sectional area
    real(PS)                    :: A_c
    ! effective cross-sectional area for ice crystals inside
    real(PS)                    :: A_e_ic
    ! circumscribed cross-sectional area
    real(PS)                    :: A_c_ic
    ! the ratio of the effective over the circumscribed cross-sectional area, Ae/A
    real(PS)                    :: q
    real(PS)                    :: q_cs
    ! charasteristic length defined for Bohm (1992)
    real(PS)                    :: char_len

    ! axial ratio, c/(2*(a+d))
    real(PS)                    :: alpha
    ! shape factor
    real(PS)                    :: k
    ! the Davies or Best number
    real(PS)                    :: X, X_0!, X_old
    ! pressure drag coefficient
    real(PS)                    :: C_DP
    !
    real(PS)                    :: beta, gamma!, xx
    ! Oseen-type solution for the drag coefficient
    real(PS)                    :: C_DO
    ! Reynolds number
    real(PS)                    :: N_re_0

    ! total surface area of the body
    real(PS)                    :: OMEGA
    ! Perimeter
    real(PS)                    :: P

    ! Stoke's terminal velocity
    real(PS) :: U_S
    ! radius of a drop
    real(PS) :: rad
    ! Bond number and physical property number
    real(PS) :: NBO,NP
    real(PS) :: Y,lambda_a

    ! artificial reduction of terminal velocity
    real(PS) :: artificial_tv

    integer                    :: i,n
!    integer,dimension(g%N_bin*g%L) :: ierror


!    ierror=0

    artificial_tv = 0.95_PS

    if( phase == 1 ) then

!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then

          ! +++ calculate the mass component for one particle +++
          rad=g%MS(i,n)%len*0.5_PS


          if(rad<0.5e-4_PS) then  ! bug found 2016/05/30
            g%MS(i,n)%vtm=0.0_PS
          elseif(0.5e-4_PS<=rad.and.rad<10.0e-4_PS) then
            ! calculate Stokes terminal velocity
            U_S=rad**2.0*gg*(den_w-ag%TV(n)%den_a)/4.5_PS/ag%TV(n)%d_vis
            lambda_a=6.6e-6_PS*(ag%TV(n)%d_vis/1.818e-4_PS)*(1013250.0_PS/ag%TV(n)%P)*(ag%TV(n)%T/293.15_PS)
            g%MS(i,n)%vtm=(1.0_PS+1.26_PS*lambda_a/rad)*U_S

          elseif(rad<535.0e-4_PS) then
            X=log(32.0_PS*rad**3.0*(den_w-ag%TV(n)%den_a)*ag%TV(n)%den_a*gg/3.0_PS/ag%TV(n)%d_vis**2.0)

            Y=-0.318657e+1_PS+0.992696_PS*X-0.153193e-2_PS*X*X-0.987059e-3_PS*X*X*X &
               -0.578878e-3_PS*X*X*X*X+0.855176e-4*X*X*X*X*X &
               -0.327815e-5*X*X*X*X*X*X
            g%MS(i,n)%Nre=exp(Y)
            ! +++ calculate the terminal velocity +++
            g%MS(i,n)%vtm=g%MS(i,n)%Nre*ag%TV(n)%d_vis/(2.0_PS*rad*ag%TV(n)%den_a)

          else
            rad=min(rad,3500.0e-4_RP)

            NBO=gg*(den_w-ag%TV(n)%den_a)*rad**2.0/ag%TV(n)%sig_wa
            NP=ag%TV(n)%sig_wa**3.0*ag%TV(n)%den_a**2.0/ag%TV(n)%d_vis**4.0/gg

            X=log(NBO*NP**(1.0/6.0)*16.0_PS/3.0_PS)
            Y=-0.500015e+1_PS+0.523778e+1_PS*X-0.204914e+1_PS*X*X+0.475294_PS*X*X*X &
               -0.542819e-1_PS*X*X*X*X+0.238449e-2_PS*X*X*X*X*X


            g%MS(i,n)%Nre=NP**(1.0/6.0)*exp(Y)
            ! +++ calculate the terminal velocity +++
            g%MS(i,n)%vtm=g%MS(i,n)%Nre*ag%TV(n)%d_vis/(2.0_PS*rad*ag%TV(n)%den_a)


          end if
        endif
      enddo
      enddo
    else if( phase == 2 ) then
       ! --- in case of ice phase ---
!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then

          ! +++ get aspect ratio +++
          alpha=get_aspect_ratio_sol(g%MS(i,n),g%IS(i,n))

          ! +++ calculate the total surface area +++
          OMEGA=get_total_sfc_area( g%IS(i,n), g%MS(i,n)%semi_a, g%MS(i,n)%semi_c, mode)

          ! +++ calculate the area ratio +++
          ! ++++++ calculate the effective cross section ++++++
          A_e = get_effect_area( g%IS(i,n),g%IS(i,n)%sh_type, &
                                 g%MS(i,n)%semi_a, g%MS(i,n)%semi_c, mode)

          A_e_ic = get_effect_area( g%IS(i,n),1, &
                                    g%MS(i,n)%a_len, g%MS(i,n)%c_len, mode)

          ! ++++++ calculate the circumscribed cross-sectional area ++++++
          A_c = get_circum_area( g%IS(i,n),g%IS(i,n)%sh_type, &
                              g%MS(i,n)%semi_a, g%MS(i,n)%semi_c, mode)

          A_c_ic = get_circum_area( g%IS(i,n),1, &
                              g%MS(i,n)%a_len, g%MS(i,n)%c_len, mode)

          ! +++ calculate the characteristic length based on
          !     perimeter and total surface area +++
          P = get_perimeter( g%IS(i,n), g%MS(i,n)%semi_a, g%MS(i,n)%semi_c, mode)

          ! +++ calculate characteristic length +++
          char_len = 2.0*sqrt(A_c/PI)

          ! +++ calculate area ratio +++
          q_cs=(g%MS(i,n)%mass(imc)+g%MS(i,n)%mass(ima)+g%MS(i,n)%mass(imr))&
                     /g%MS(i,n)%con/den_i/&
                 get_vip(g%IS(i,n)%is_mod(2),g%IS(i,n)%phi_cs,g%IS(i,n)%semi_aip)

!!!!CDIR NOIEXPAND
          q = get_area_ratio_sol(g%MS(i,n),g%IS(i,n),mode,A_e,A_c,A_e_ic,A_c_ic,q_cs)

          g%MS(i,n)%len = OMEGA/P
          g%IS(i,n)%q_e=q

!!c       if(g%IS(i,n)%sh_type>=3) then
!!c          open( unit=21, file="check_svtm.dat",POSITION='APPEND')
!!c          write(21,'(10ES15.6)') g%MS(i,n)%mean_mass,g%MS(i,n)%semi_a,g%MS(i,n)%semi_c,alpha,q,char_len,g%MS(i,n)%Nre
!!c          close(21)
!!c       end if

          ! +++ calculate the shape factor +++
          if( alpha <= 0.16_PS ) then
            k = 0.85_PS
          else if ( 0.16_PS < alpha .and. alpha <= 1.0_PS ) then
            k = 0.82_PS + 0.18_PS*alpha
          else if ( 1.0_PS < alpha .and. alpha <= 3.0_PS ) then
            k = 0.37_PS + 0.63_PS/alpha
          else
            k = 1.33_PS/( log(alpha) + 1.19_PS)
          end if

          ! +++ calculate the Davies or Best number +++
          if( q <= 1.0_PS ) then
            X = 8.0_PS*g%MS(i,n)%mean_mass*gg*ag%TV(n)%den_a/&
               (PI*(ag%TV(n)%d_vis**2)*max(alpha, 1.0_PS)*(q**0.25_PS))
          else
            X = 8.0_PS*g%MS(i,n)%mean_mass*gg*ag%TV(n)%den_a/&
               (PI*(ag%TV(n)%d_vis**2)*max(alpha, 1.0_PS)*q)
          end if

          ! +++ calculate the pressure drag C_DP +++
          if( alpha <= 0.3_PS ) then
            gamma = 1.98_PS
            C_DP = 0.292*k*gamma
          elseif( 0.3_PS < alpha .and. alpha < 1.0_PS ) then
            gamma = 3.76_PS - 8.41_PS*alpha + 9.18_PS*alpha**2 - 3.53_PS*alpha**3
            C_DP = 0.292*k*gamma
          else
            if( 1.0_PS <= q ) then
              gamma = 1.0_PS
              C_DP = (1.46_PS*q-0.46_PS)*(0.492_PS-0.200_PS/sqrt(alpha))
            else
              gamma = 1.0_PS
              C_DP = 0.292_PS*k*gamma
            end if
          end if


          ! +++ calculate the Oseen-type drag coefficient, C_DO +++
          C_DO = 4.5_PS*k*k*max(alpha, 1.0_PS)

          ! +++ calculate the Reynolds number +++
          beta = sqrt(1.0_PS + (C_DP/(6.0_PS*k))*sqrt(X/C_DP)) - 1.0_PS
          N_re_0 = 6.0_PS*k*beta*beta/C_DP

          ! +++ correction for transition to turbulent flow +++
          if( 1000.0_PS <= N_re_0 ) then

            X_0 = 2.8e+06_PS
            C_DP=C_DP*(1.0_PS + 1.6_PS*(X/X_0)**2)/(1.0_PS + (X/X_0)**2)

            beta = sqrt(1.0_PS + (C_DP/(6.0_PS*k))*sqrt(X/C_DP)) - 1.0_PS
            N_re_0 = 6.0_PS*k*beta*beta/C_DP
          endif


          ! +++ matching theory +++
          gamma = (C_DO-C_DP)/(4.0_PS*C_DP)
          g%MS(i,n)%Nre = N_re_0*(1.0_PS+2.0_PS*beta*exp(-beta*gamma)/&
               ((2.0_PS+beta)*(1.0_PS+beta)))

          if( g%MS(i,n)%Nre < 0.0_PS ) then
            g%MS(i,n)%Nre=0.0_PS
            LOG_ERROR("cal_terminal_vel_vec",*) "Reynolds number is negative!!"
            LOG_ERROR_CONT('(10ES14.6)') g%MS(i,n)%mass,g%MS(i,n)%con,g%MS(i,n)%a_len,g%MS(i,n)%c_len
            LOG_ERROR_CONT('(2I5,10ES14.6)') g%IS(i,n)%sh_type,g%IS(i,n)%habit,g%MS(i,n)%semi_a,g%MS(i,n)%semi_c,g%IS(i,n)%phi_cs
            LOG_ERROR_CONT('(10ES14.6)') g%IS(i,n)%semi_aip,g%IS(i,n)%semi_cip,g%IS(i,n)%V_cs
            call PRC_abort
          end if

          ! +++ calculate the terminal velocity +++
          g%MS(i,n)%vtm = g%MS(i,n)%Nre*ag%TV(n)%d_vis/(char_len*ag%TV(n)%den_a)

          if(g%MS(i,n)%vtm>5.0e+3.and.g%MS(i,n)%den>0.8_PS) then
            g%MS(i,n)%vtm=5.0e+3_PS
            LOG_ERROR("cal_terminal_vel_vec",*) "terminal_vel > The ice vel is set to 50.0m/s"
            LOG_ERROR_CONT(*) "old vel",g%MS(i,n)%vtm
            LOG_ERROR_CONT('("mass,con,alen,clen,q,k,gamma",10ES14.6)') g%MS(i,n)%mass(1),g%MS(i,n)%con,g%MS(i,n)%a_len,g%MS(i,n)%c_len,q,k,gamma
            LOG_ERROR_CONT('("N_re0,ex,beta,X,C_DP",10ES14.6)') N_re_0,exp(-beta*gamma),beta,X,C_DP
            LOG_ERROR_CONT('("semi_a,semi_c,mean_mass",10ES14.6)') g%MS(i,n)%semi_a,g%MS(i,n)%semi_c,g%MS(i,n)%mean_mass
            LOG_ERROR_CONT('("char_len,d_vis,den,Nre,A_c",10ES14.6)') char_len,ag%TV(n)%d_vis,ag%TV(n)%den,g%MS(i,n)%Nre,A_c
            LOG_ERROR_CONT('("pressure,T",10ES14.6)') ag%TV(n)%P,ag%TV(n)%T
            call PRC_abort
          elseif(g%MS(i,n)%vtm>2.0e+3) then
            g%MS(i,n)%vtm=2.0e+3_PS

          elseif(g%MS(i,n)%vtm<0.0) then
            g%MS(i,n)%vtm=0.0_PS
            LOG_ERROR("cal_terminal_vel_vec",*) "terminal_vel > The ice vel is set to 0.0m/s"
!tmp          LOG_ERROR_CONT(*) "old vel",g%MS(i,n)%vtm
!tmp          LOG_ERROR_CONT('("mass,con,alen,clen,q,k,gamma",10ES14.6)') g%MS(i,n)%mass(1),g%MS(i,n)%con,g%MS(i,n)%a_len,g%MS(i,n)%c_len,q,k,gamma
!tmp          LOG_ERROR_CONT('("N_re0,ex,beta,X,C_DP",10ES14.6)') N_re_0,exp(-beta*gamma),beta,X,C_DP
!tmp          LOG_ERROR_CONT('("semi_a,semi_c,mean_mass,beta0",10ES14.6)') g%MS(i,n)%semi_a,g%MS(i,n)%semi_c,g%MS(i,n)%mean_mass,beta0
!tmp          LOG_ERROR_CONT('("char_len,d_vis,den,Nre,A_c",10ES14.6)') char_len,ag%TV(n)%d_vis,ag%TV(n)%den,g%MS(i,n)%Nre,A_c
!tmp          LOG_ERROR_CONT('("pressure,T",10ES14.6)') ag%TV(n)%P,ag%TV(n)%T
            call PRC_abort
          endif

        endif
      enddo
      enddo

    else
!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then
          g%MS(i,n)%vtm = 0.0_PS
        endif
      enddo
      enddo
    end if

  end subroutine cal_terminal_vel_vec

!!$  subroutine cal_wterm_vel_v3(phase, ms, th_var, mode, level, &
!!$       a, m1, m2, WVT,ishape)
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! calculate the terminal velocity weighted by concentration and mass.
!!$    ! use the approximation of NRe= a*X**b fitted to results by Bohm equation.
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! phase of hydrometeor
!!$    integer, intent(in)  :: phase
!!$    type (Mass_Bin), intent(in)   :: ms
!!$    type (Thermo_Var), intent(in)    :: th_var
!!$
!!$    ! parameter of linear distribution
!!$    real(DS),dimension(*)     :: a
!!$    real(DS), dimension(*)  :: WVT
!!$    type (Ice_Shape), optional       :: ishape
!!$    !
!!$    real(DS),intent(in) :: m1,m2
!!$    real(PS) :: alpha,char_len,OMEGA,P,&
!!$         A_e,A_c,q,X1,X2
!!$    real(DS) :: A1,B1,D1,Q1,M,N,MT,NT
!!$    real(PS), dimension(4) :: a_rex,b_rex
!!$
!!$    ! parameters to define NRe=a*X**b
!!$    ! for liquid
!!$!!c    real(PS),parameter :: &
!!$!!c         la1=4.072942E-02,lb1=9.953352E-01,&
!!$!!c         la2=7.810803E-02,lb2=7.935192E-01,&
!!$!!c         la3=3.443135E-01,lb3=6.126759E-01,&
!!$!!c         la4=2.372899E+00,lb4=4.645688E-01,&
!!$!!c         LXlim1=2.492336E+01,LXlim2=3.439475E+03,LXlim3=4.490854E+05
!!$
!!$    ! does not include diameter > 5.0mm
!!$    real(PS),parameter :: &
!!$         la1=4.049863E-02,lb1=9.939921E-01,&
!!$         la2=7.395895E-02,lb2=8.073323E-01,&
!!$         la3=2.557850E-01,lb3=6.412501E-01,&
!!$         la4=1.033221E+00,lb4=5.235578E-01,&
!!$         LXlim1=2.517541E+01,LXlim2=1.700411E+03,LXlim3=1.356160E+05
!!$
!!$    real(DS),parameter :: D_max=0.45,m_max=0.047712938426395
!!$    real(DS) :: v_max
!!$
!!$    ! for solid
!!$    ! for ylim   1.000000E+00   2.000000E+01   2.000000E+02
!!$    real(PS),parameter :: &
!!$         sa1=4.742587E-02,sb1=9.910232E-01,&
!!$         sa2=9.447403E-02,sb2=7.722376E-01,&
!!$         sa3=2.549053E-01,sb3=6.299768E-01,&
!!$         sa4=7.880099E-01,sb4=5.236878E-01,&
!!$         SXlim1=2.167592E+01,SXlim2=1.027213E+03,SXlim3=3.934172E+04
!!$    ! for ylim   1.000000E+00   1.000000E+01   2.000000E+02
!!$!!c    real(PS),parameter :: &
!!$!!c         sa1=4.688108E-02,sb1=9.868769E-01,&
!!$!!c         sa2=8.315829E-02,sb2=8.049257E-01,&
!!$!!c         sa3=2.125120E-01,sb3=6.480806E-01,&
!!$!!c         sa4=8.059516E-01,sb4=5.225413E-01,&
!!$!!c         SXlim1=2.221647E+01,SXlim2=3.838900E+02,SXlim3=3.876058E+04
!!$
!!$    ! define terminal velocity regimes by masses
!!$    real(PS),dimension(5) :: X_b
!!$    ! masses and diameter devided by terminal velocity regimes.
!!$    real(DS) :: m_sub1,m_sub2,D_1,D_2
!!$
!!$    ! +++ mode of shape to calculate the ventilation coefficient
!!$    !     and terminal velocity +++
!!$    ! 1 : spheroid with a_len and c_len
!!$    ! 2 : ice crystal model
!!$    integer,intent(in)   :: mode
!!$
!!$    ! level of complexity
!!$    integer,intent(in)   :: level
!!$    integer  :: i, var_Status
!!$    integer :: em
!!$
!!$    real(PS), parameter :: chiarui_1=0.0_PS, chiarui_2=1.0e+30
!!$
!!$
!!$    WVT(1:2)=0.0_DS
!!$    NT=0.0_DS
!!$    MT=0.0_DS
!!$
!!$
!!$    if( phase == 1 ) then
!!$       ! --- in case of liquid phase ---
!!$       a_rex=(/la1,la2,la3,la4/)
!!$       b_rex=(/lb1,lb2,lb3,lb4/)
!!$       X_b=(/chiarui_1,LXlim1,LXlim2,LXlim3,chiarui_2/)
!!$
!!$       ! calculate Best number for the bin limits
!!$       ! +++ calculate the axial ratio +++
!!$       alpha=ms%c_len/ms%a_len
!!$
!!$       ! +++ calculate the effective cross section +++
!!$       A_e = PI * ms%a_len*ms%a_len
!!$
!!$       ! +++ calculate the circumscribed cross-sectional area +++
!!$       A_c = A_e
!!$
!!$       q = 1.0_PS
!!$
!!$       char_len=2.0_PS*ms%a_len
!!$
!!$       ! +++ calculate the Davies or Best number +++
!!$       D1=8.0_DS*gg*th_var%den/&
!!$            (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*(q**0.25))
!!$
!!$       Q1=q**0.25
!!$
!!$       X1=D1*m1
!!$       X2=D1*m2
!!$
!!$
!!$       do i=1,4
!!$
!!$          if(X_b(i+1)<X1) cycle
!!$          if(X_b(i)>X2) exit
!!$
!!$          m_sub1=max(m1,real(X_b(i)/D1,DS))
!!$          m_sub2=min(m2,real(X_b(i+1)/D1,DS))
!!$
!!$          D_1=(m_sub1/coedpi6)**(1.0_RP/3.0_RP)
!!$          D_2=(m_sub2/coedpi6)**(1.0_RP/3.0_RP)
!!$
!!$          N=(m_sub2-m_sub1)*(a(1)-a(3)*(a(2)-(m_sub2+m_sub1)/2.0))
!!$          M=a(1)*(m_sub2**2.0_RP-m_sub1**2.0_RP)/2.0_RP &
!!$               +a(3)*((m_sub2**3-m_sub1**3)/3.0_RP-a(2)*(m_sub2**2-m_sub1**2)/2.0_RP)
!!$
!!$
!!$          A1=th_var%d_vis**(1.0_RP-2.0_RP*b_rex(i))*th_var%den_a**(b_rex(i)-1.0_RP)*&
!!$               (8.0_PS*gg/PI)**b_rex(i)*(Q1*max(alpha,1.0_RP))**(-b_rex(i))*&
!!$               ms%mean_mass**(1.0_RP/3.0_RP)/char_len*a_rex(i)
!!$
!!$          B1=b_rex(i)-1.0_DS/3.0_DS
!!$
!!$          if(D_2<D_max) then
!!$             ! 1. mass-weighted terminal velocity
!!$             WVT(1)=WVT(1)+max(0.0_DS,&
!!$                  A1*( (a(1)-a(3)*a(2))/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS))+&
!!$                  a(3)/(B1+3.0_DS)*(m_sub2**(B1+3.0_DS)-m_sub1**(B1+3.0_DS))))
!!$
!!$             ! 2. con-weighted terminal velocity
!!$             WVT(2)=WVT(2)+max(0.0_DS,&
!!$                  A1*( (a(1)-a(3)*a(2))/(B1+1.0_DS)*(m_sub2**(B1+1.0_DS)-m_sub1**(B1+1.0_DS))+&
!!$                  a(3)/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS))))
!!$          elseif(D_1<D_max) then
!!$
!!$             v_max=A1*m_max**B1
!!$
!!$             ! 1. mass-weighted terminal velocity
!!$             WVT(1)=WVT(1)+max(0.0_DS,&
!!$                  A1*( (a(1)-a(3)*a(2))/(B1+2.0_DS)*(m_max**(B1+2.0_DS)-m_sub1**(B1+2.0_DS))+&
!!$                  a(3)/(B1+3.0_DS)*(m_max**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)))+&
!!$                  v_max*( (a(1)-a(3)*a(2))/2.0_DS*(m_sub2**2.0-m_max**2.0)+&
!!$                  a(3)/3.0_DS*(m_sub2**3.0-m_max**3.0)))
!!$
!!$             ! 2. con-weighted terminal velocity
!!$             WVT(2)=WVT(2)+max(0.0_DS,&
!!$                  A1*( (a(1)-a(3)*a(2))/(B1+1.0_DS)*(m_max**(B1+1.0_DS)-m_sub1**(B1+1.0_DS))+&
!!$                  a(3)/(B1+2.0_DS)*(m_max**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)))+&
!!$                  v_max*( (a(1)-a(3)*a(2))*(m_sub2-m_max)+&
!!$                  a(3)/2.0_DS*(m_sub2**2.0-m_max**2.0)))
!!$          else
!!$             v_max=A1*m_max**B1
!!$
!!$             ! 1. mass-weighted terminal velocity
!!$             WVT(1)=WVT(1)+max(0.0_DS,&
!!$                  v_max*( (a(1)-a(3)*a(2))/2.0_DS*(m_sub2**2.0-m_sub1**2.0)+&
!!$                  a(3)/3.0_DS*(m_sub2**3.0-m_sub1**3.0)))
!!$
!!$             ! 2. con-weighted terminal velocity
!!$             WVT(2)=WVT(2)+max(0.0_DS,&
!!$                  v_max*( (a(1)-a(3)*a(2))*(m_sub2-m_sub1)+&
!!$                  a(3)/2.0_DS*(m_sub2**2.0-m_sub1**2.0)))
!!$
!!$          end if
!!$          NT=NT+N
!!$          MT=MT+M
!!$       end do
!!$       WVT(1)=WVT(1)/MT
!!$       WVT(2)=WVT(2)/NT
!!$!!c       WVT(1)=max(min(WVT(1)/MT,ms%vtm*2.0),ms%vtm*0.5)
!!$!!c       WVT(2)=max(min(WVT(2)/NT,ms%vtm*2.0),ms%vtm*0.5)
!!$
!!$    else if( phase == 2 ) then
!!$       ! --- in case of liquid phase ---
!!$       a_rex=(/sa1,sa2,sa3,sa4/)
!!$       b_rex=(/sb1,sb2,sb3,sb4/)
!!$       X_b=(/chiarui_1,SXlim1,SXlim2,SXlim3,chiarui_2/)
!!$
!!$       ! calculate Best number for the bin limits
!!$       call cal_some_vel(ms,ishape,mode,alpha,OMEGA,A_e,A_c,P,q,char_len)
!!$
!!$       ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$       ! coef. from X
!!$       if( q <= 1.0_PS ) then
!!$          D1 = 8.0_PS*gg*th_var%den/&
!!$               (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*(q**0.25))
!!$          Q1=q**0.25
!!$       else if( q > 1.0_PS ) then
!!$          D1 = 8.0_PS*gg*th_var%den/&
!!$               (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*q)
!!$          Q1=q
!!$       end if
!!$
!!$       X1=D1*m1
!!$       X2=D1*m2
!!$
!!$       do i=1,4
!!$
!!$          if(X_b(i+1)<X1) cycle
!!$          if(X_b(i)>X2) exit
!!$
!!$          m_sub1=max(m1,real(X_b(i)/D1,DS))
!!$          m_sub2=min(m2,real(X_b(i+1)/D1,DS))
!!$
!!$          N=(m_sub2-m_sub1)*(a(1)-a(3)*(a(2)-(m_sub2+m_sub1)/2.0_RP))
!!$          M=a(1)*(m_sub2**2-m_sub1**2.0)/2.0_RP &
!!$               +a(3)*((m_sub2**3-m_sub1**3)/3.0_RP-a(2)*(m_sub2**2-m_sub1**2)/2.0_RP)
!!$
!!$
!!$          A1=th_var%d_vis**(1.0_RP-2.0_RP*b_rex(i))*th_var%den_a**(b_rex(i)-1.0_RP)*&
!!$               (8.0_DS*gg/PI)**b_rex(i)*(Q1*max(alpha,1.0_RP))**(-b_rex(i))*&
!!$               ms%mean_mass**(1.0_RP/3.0_RP)/char_len*a_rex(i)
!!$
!!$          B1=b_rex(i)-1.0_DS/3.0_DS
!!$
!!$          ! 1. mass-weighted terminal velocity
!!$          WVT(1)=WVT(1)+max(0.0_DS,&
!!$               A1*( (a(1)-a(3)*a(2))/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS))+&
!!$               a(3)/(B1+3.0_DS)*(m_sub2**(B1+3.0_DS)-m_sub1**(B1+3.0_DS))))
!!$
!!$          ! 2. con-weighted terminal velocity
!!$          WVT(2)=WVT(2)+max(0.0_DS,&
!!$               A1*( (a(1)-a(3)*a(2))/(B1+1.0_DS)*(m_sub2**(B1+1.0_DS)-m_sub1**(B1+1.0_DS))+&
!!$               a(3)/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS))))
!!$
!!$          NT=NT+N
!!$          MT=MT+M
!!$       end do
!!$       WVT(1)=WVT(1)/MT
!!$       WVT(2)=WVT(2)/NT
!!$    end  if
!!$
!!$
!!$  end subroutine cal_wterm_vel_v3

  subroutine cal_wterm_vel_v3_vec(g,ag, mode, &
       a2d, m1, m2, error_number, WVT)
    use class_Mass_Bin, only: &
       get_aspect_ratio_sol, &
       get_area_ratio_sol, &
       cal_best_number
    use class_Ice_Shape, only: &
       get_total_sfc_area, &
       get_effect_area, &
       get_circum_area, &
       get_perimeter, &
       get_vip
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate the terminal velocity weighted by concentration and mass.
    ! use the approximation of NRe= a*X**b fitted to results by Bohm equation.
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    type (Group), intent(in)   :: g
    ! thermo variable object
    type (AirGroup), intent(in)  :: ag

    ! paarameter of linear distribution
    real(8),dimension(g%n_bin,g%L,4)     :: a2d
    real(DS), dimension(2,g%n_bin,*)  :: WVT

    real(8),dimension(g%n_bin,*),intent(in) :: m1,m2
    integer,dimension(g%n_bin,*),intent(in) :: error_number


    real(PS), dimension(4) :: a_rex,b_rex

    real(PS) :: alpha,char_len,OMEGA,P,&
         A_e,A_c,A_e_ic,A_c_ic,q_cs,q,X1,X2
    real(DS) :: A1,B1,D1,Q1,aM,aN
    real(DS),dimension(g%n_bin,g%L) :: aMT,aNT
    real(8),dimension(4)     :: a
    real(8) :: xA,xB,xC,xD,aN_abmx,aM_abmx

    ! parameters to define NRe=a*X**b
    ! for liquid
!!c    real(PS),parameter :: &
!!c         la1=4.072942E-02,lb1=9.953352E-01,&
!!c         la2=7.810803E-02,lb2=7.935192E-01,&
!!c         la3=3.443135E-01,lb3=6.126759E-01,&
!!c         la4=2.372899E+00,lb4=4.645688E-01,&
!!c         LXlim1=2.492336E+01,LXlim2=3.439475E+03,LXlim3=4.490854E+05

    ! does not include diameter > 5.0mm
    real(PS),parameter :: &
         la1=4.049863E-02,lb1=9.939921E-01,&
         la2=7.395895E-02,lb2=8.073323E-01,&
         la3=2.557850E-01,lb3=6.412501E-01,&
         la4=1.033221E+00,lb4=5.235578E-01,&
         LXlim1=2.517541E+01,LXlim2=1.700411E+03,LXlim3=1.356160E+05

    real(DS),parameter :: D_max=0.45,m_max=0.047712938426395
    real(DS) :: v_max

    ! for solid
    ! for ylim   1.000000E+00   2.000000E+01   2.000000E+02
    real(PS),parameter :: &
         sa1=4.742587E-02,sb1=9.910232E-01,&
         sa2=9.447403E-02,sb2=7.722376E-01,&
         sa3=2.549053E-01,sb3=6.299768E-01,&
         sa4=7.880099E-01,sb4=5.236878E-01,&
         SXlim1=2.167592E+01,SXlim2=1.027213E+03,SXlim3=3.934172E+04
    ! for ylim   1.000000E+00   1.000000E+01   2.000000E+02
!!c    real(PS),parameter :: &
!!c         sa1=4.688108E-02,sb1=9.868769E-01,&
!!c         sa2=8.315829E-02,sb2=8.049257E-01,&
!!c         sa3=2.125120E-01,sb3=6.480806E-01,&
!!c         sa4=8.059516E-01,sb4=5.225413E-01,&
!!c         SXlim1=2.221647E+01,SXlim2=3.838900E+02,SXlim3=3.876058E+04

    ! define terminal velocity regimes by masses
    real(PS),dimension(5) :: X_b
    ! masses and diameter devided by terminal velocity regimes.
    real(DS) :: m_sub1,m_sub2,D_1,D_2

    ! +++ mode of shape to calculate the ventilation coefficient
    !     and terminal velocity +++
    ! 1 : spheroid with a_len and c_len
    ! 2 : ice crystal model
    integer,intent(in)   :: mode

    integer :: i_den_gt0p5

    ! level of complexity
    !integer,intent(in)   :: level
    integer  :: i,j,n!,var_Status
    !integer :: em

    real(PS),parameter :: chiarui_1=0.0_PS, chiarui_2=1.0e+30

    if( g%token == 1 ) then
      ! --- in case of liquid phase ---
      a_rex=(/la1,la2,la3,la4/)
      b_rex=(/lb1,lb2,lb3,lb4/)
      X_b=(/chiarui_1,LXlim1,LXlim2,LXlim3,chiarui_2/)


!      do in=1,g%N_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
      do n = 1, g%L
      do i = 1, g%N_BIN
        WVT(1,i,n)=0.0_DS
        WVT(2,i,n)=0.0_DS
        aNT(i,n)=0.0_DS
        aMT(i,n)=0.0_DS
      enddo
      enddo

      do j=1,4
!        do in=1,g%N_bin*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
         do n = 1, g%L
         do i = 1, g%N_BIN

      no_error_if1: if(error_number(i,n)==0) then

          a(1)=a2d(i,n,1)
          a(2)=a2d(i,n,2)
          a(3)=a2d(i,n,3)
          a(4)=a2d(i,n,4)

          ! calculate Best number for the bin limits
          alpha=g%MS(i,n)%c_len/g%MS(i,n)%a_len

          q = 1.0_PS
!tmp          q = 0.99999_PS

          char_len=2.0_PS*g%MS(i,n)%a_len

          ! ++ calculate Best number for the bin limits +++
          ! coef. from X
          call cal_best_number(D1,Q1,ag%TV(n),q,alpha)

          m_sub1=max(m1(i,n),real(X_b(j)/D1,DS))
          m_sub2=min(m2(i,n),real(X_b(j+1)/D1,DS))

          D_1=(m_sub1/coedpi6)**(1.0/3.0)
          D_2=(m_sub2/coedpi6)**(1.0/3.0)

          X1=D1*m1(i,n)
          X2=D1*m2(i,n)

          xA=m_sub2+m_sub1
          xB=m_sub2-m_sub1
          xC=m_sub2*m_sub1
          xD=m_sub2*m_sub2+m_sub1*m_sub1

          aN=a(1)*xB &
                +0.5d+0*a(2)*xA*xB &
                +a(3)*xB*(xD+xC)/3.0_RP &
                +0.25_RP*a(4)*xD*xA*xB

          aM=0.5d+0*a(1)*xA*xB &
                   +a(2)*xB*(xD+xC)/3.0d+0 &
                   +0.25d+0*a(3)*xD*xA*xB &
                   +0.2d+0*a(4)*xB*(xA*xA*(xD-xC)+xC*xC)



          !
          ! Vt=d_vis/den_a * Nre/d
          !   =d_vis/den_a * (1/d_bar) * (mass_bar/m)**(1/3) * a*X^b
          !   =A1*mass^B1
          !
          A1=ag%TV(n)%d_vis**(1.0_RP-2.0_RP*b_rex(j))*ag%TV(n)%den_a**(b_rex(j)-1.0_RP)*&
             (8.0_PS*gg/PI)**b_rex(j)*(Q1*max(alpha,1.0_RP))**(-b_rex(j))*&
             g%MS(i,n)%mean_mass**(1.0_RP/3.0_RP)/char_len*a_rex(j)

          B1=b_rex(j)-1.0_DS/3.0_DS

          if(X_b(j+1)<X1.or.X_b(j)>X2) then

          elseif(D_2<D_max) then
            !
            ! in case of all the drops are within D_max
            !
            ! 1. mass-weighted terminal velocity
            WVT(1,i,n)=WVT(1,i,n)+max(0.0_DS,&
                  A1*( a(1)/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)) &
                     + a(2)/(B1+3.0_DS)*(m_sub2**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)) &
                     + a(3)/(B1+4.0_DS)*(m_sub2**(B1+4.0_DS)-m_sub1**(B1+4.0_DS)) &
                     + a(4)/(B1+5.0_DS)*(m_sub2**(B1+5.0_DS)-m_sub1**(B1+5.0_DS)) &
                       ))

            ! 2. con-weighted terminal velocity
            WVT(2,i,n)=WVT(2,i,n)+max(0.0_DS,&
                  A1*( a(1)/(B1+1.0_DS)*(m_sub2**(B1+1.0_DS)-m_sub1**(B1+1.0_DS)) &
                     + a(2)/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)) &
                     + a(3)/(B1+3.0_DS)*(m_sub2**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)) &
                     + a(4)/(B1+4.0_DS)*(m_sub2**(B1+4.0_DS)-m_sub1**(B1+4.0_DS)) &
                       ))

            aNT(i,n)=aNT(i,n)+aN
            aMT(i,n)=aMT(i,n)+aM
          elseif(D_1<D_max) then
            !
            ! in case of D_2 boundary is above D_max
            !
            !   above D_max, the velocity becomes constant.
            !

            v_max=A1*m_max**B1

            xA=m_sub2+m_max
            xB=m_sub2-m_max
            xC=m_sub2*m_max
            xD=m_sub2*m_sub2+m_max*m_max

            aN_abmx=a(1)*xB &
                  +0.5d+0*a(2)*xA*xB &
                  +a(3)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a(4)*xD*xA*xB

            aM_abmx=0.5d+0*a(1)*xA*xB &
                     +a(2)*xB*(xD+xC)/3.0d+0 &
                     +0.25d+0*a(3)*xD*xA*xB &
                     +0.2d+0*a(4)*xB*(xA*xA*(xD-xC)+xC*xC)

            ! 1. mass-weighted terminal velocity
            WVT(1,i,n)=WVT(1,i,n)+max(0.0_DS,&
                  A1*( a(1)/(B1+2.0_DS)*(m_max**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)) &
                     + a(2)/(B1+3.0_DS)*(m_max**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)) &
                     + a(3)/(B1+4.0_DS)*(m_max**(B1+4.0_DS)-m_sub1**(B1+4.0_DS)) &
                     + a(4)/(B1+5.0_DS)*(m_max**(B1+5.0_DS)-m_sub1**(B1+5.0_DS)) &
                      ) &
                  + &
                  v_max*aM_abmx)

            ! 2. con-weighted terminal velocity
            WVT(2,i,n)=WVT(2,i,n)+max(0.0_DS,&
                  A1*( a(1)/(B1+1.0_DS)*(m_max**(B1+1.0_DS)-m_sub1**(B1+1.0_DS)) &
                     + a(2)/(B1+2.0_DS)*(m_max**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)) &
                     + a(3)/(B1+3.0_DS)*(m_max**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)) &
                     + a(4)/(B1+4.0_DS)*(m_max**(B1+4.0_DS)-m_sub1**(B1+4.0_DS)) &
                       ) &
                  + &
                  v_max*aN_abmx)

            aNT(i,n)=aNT(i,n)+aN
            aMT(i,n)=aMT(i,n)+aM

          else
            !
            ! in case of all the drops are larger than D_max.
            !
            !   above D_max, the velocity becomes constant.
            !
            v_max=A1*m_max**B1

            xA=m_sub2+m_sub1
            xB=m_sub2-m_sub1
            xC=m_sub2*m_sub1
            xD=m_sub2*m_sub2+m_sub1*m_sub1

            aN_abmx=a(1)*xB &
                  +0.5d+0*a(2)*xA*xB &
                  +a(3)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a(4)*xD*xA*xB

            aM_abmx=0.5d+0*a(1)*xA*xB &
                     +a(2)*xB*(xD+xC)/3.0d+0 &
                     +0.25d+0*a(3)*xD*xA*xB &
                     +0.2d+0*a(4)*xB*(xA*xA*(xD-xC)+xC*xC)

            ! 1. mass-weighted terminal velocity
            WVT(1,i,n)=WVT(1,i,n)+max(0.0_DS,&
                       v_max*aM_abmx)

            ! 2. con-weighted terminal velocity
            WVT(2,i,n)=WVT(2,i,n)+max(0.0_DS,&
                       v_max*aN_abmx)


            aNT(i,n)=aNT(i,n)+aN
            aMT(i,n)=aMT(i,n)+aM
          end if
         endif no_error_if1

        end do
        end do
      end do

!      do in=1,g%N_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
      do n = 1, g%L
      do i = 1, g%N_BIN
        if(error_number(i,n)==0) then
          WVT(1,i,n)=min(max(0.0_DS,WVT(1,i,n)/aMT(i,n)),1.5e+3_DS)
          WVT(2,i,n)=min(max(0.0_DS,WVT(2,i,n)/aNT(i,n)),1.5e+3_DS)
        elseif(error_number(i,n)<10) then
          WVT(1,i,n)=real(min(max(0.0_ps,g%MS(i,n)%vtm),1.5e+3_ps),8)
          WVT(2,i,n)=real(min(max(0.0_ps,g%MS(i,n)%vtm),1.5e+3_ps),8)
        endif
      enddo
      enddo

    elseif( g%token == 2 ) then
      ! --- in case of liquid phase ---
      a_rex=(/sa1,sa2,sa3,sa4/)
      b_rex=(/sb1,sb2,sb3,sb4/)
      X_b=(/chiarui_1,SXlim1,SXlim2,SXlim3,chiarui_2/)

!      do in=1,g%N_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
      do n = 1, g%L
      do i = 1, g%N_BIN
        WVT(1,i,n)=0.0_DS
        WVT(2,i,n)=0.0_DS
        aNT(i,n)=0.0_DS
        aMT(i,n)=0.0_DS
      enddo
      enddo

      do j=1,4
!        do in=1,g%N_bin*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
         do n = 1, g%L
         do i = 1, g%N_BIN

      no_error_if2: if(error_number(i,n)==0) then

          a(1)=a2d(i,n,1)
          a(2)=a2d(i,n,2)
          a(3)=a2d(i,n,3)
          a(4)=a2d(i,n,4)

          ! +++ get aspect ratio +++
          alpha=get_aspect_ratio_sol(g%MS(i,n),g%IS(i,n))

          ! +++ calculate the total surface area +++
          OMEGA=get_total_sfc_area( g%IS(i,n), g%MS(i,n)%semi_a, g%MS(i,n)%semi_c, mode)

          ! +++ calculate the area ratio +++
          ! ++++++ calculate the effective cross section ++++++
          A_e = get_effect_area( g%IS(i,n),g%IS(i,n)%sh_type, &
                           g%MS(i,n)%semi_a, g%MS(i,n)%semi_c, mode)

          A_e_ic = get_effect_area( g%IS(i,n),1, &
                           g%MS(i,n)%a_len, g%MS(i,n)%c_len, mode)

          ! ++++++ calculate the circumscribed cross-sectional area ++++++
          A_c = get_circum_area( g%IS(i,n),g%IS(i,n)%sh_type, &
                           g%MS(i,n)%semi_a, g%MS(i,n)%semi_c, mode)

          A_c_ic = get_circum_area( g%IS(i,n),1, &
                           g%MS(i,n)%a_len, g%MS(i,n)%c_len, mode)

          ! +++ calculate the characteristic length based on
          !     perimeter and total surface area +++
          P = get_perimeter( g%IS(i,n), g%MS(i,n)%semi_a, g%MS(i,n)%semi_c, mode)

          ! +++ calculate characteristic length +++
          char_len = 2.0*sqrt(A_c/PI)

          ! +++ calculate area ratio +++
          q_cs=(g%MS(i,n)%mass(imc)+g%MS(i,n)%mass(ima)+g%MS(i,n)%mass(imr))&
                     /g%MS(i,n)%con/den_i/&
              get_vip(g%IS(i,n)%is_mod(2),g%IS(i,n)%phi_cs,g%IS(i,n)%semi_aip)

          q = get_area_ratio_sol(g%MS(i,n),g%IS(i,n),mode,A_e,A_c,A_e_ic,A_c_ic,q_cs)


          ! ++ calculate Best number for the bin limits +++
          ! coef. from X
          call cal_best_number(D1,Q1,ag%TV(n),q,alpha)

          X1=D1*m1(i,n)
          X2=D1*m2(i,n)


          if(X_b(j+1)<X1.or.X_b(j)>X2) then

          else

            m_sub1=max(m1(i,n),real(X_b(j)/D1,DS))
            m_sub2=min(m2(i,n),real(X_b(j+1)/D1,DS))

            xA=m_sub2+m_sub1
            xB=m_sub2-m_sub1
            xC=m_sub2*m_sub1
            xD=m_sub2*m_sub2+m_sub1*m_sub1

            aN=a(1)*xB &
                  +0.5_RP*a(2)*xA*xB &
                  +a(3)*xB*(xD+xC)/3.0_RP &
                  +0.25_RP*a(4)*xD*xA*xB

            aM=0.5_RP*a(1)*xA*xB &
                     +a(2)*xB*(xD+xC)/3.0_RP &
                     +0.25_RP*a(3)*xD*xA*xB &
                     +0.2_RP*a(4)*xB*(xA*xA*(xD-xC)+xC*xC)


            !
            ! Vt=d_vis/den_a * Nre/d
            !   =d_vis/den_a * (1/d_bar) * (mass_bar/m)**(1/3) * a*X^b
            !   =A1*mass^B1
            !
            A1=ag%TV(n)%d_vis**(1.0_RP-2.0_RP*b_rex(j))*ag%TV(n)%den_a**(b_rex(j)-1.0)*&
                 (8.0_DS*gg/PI)**b_rex(j)*(Q1*max(alpha,1.0_RP))**(-b_rex(j))*&
                 g%MS(i,n)%mean_mass**(1.0_RP/3.0_RP)/char_len*a_rex(j)

            B1=b_rex(j)-1.0_DS/3.0_DS

            ! 1. mass-weighted terminal velocity
            WVT(1,i,n)=WVT(1,i,n)+max(0.0_DS,&
                  A1*( a(1)/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)) &
                     + a(2)/(B1+3.0_DS)*(m_sub2**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)) &
                     + a(3)/(B1+4.0_DS)*(m_sub2**(B1+4.0_DS)-m_sub1**(B1+4.0_DS)) &
                     + a(4)/(B1+5.0_DS)*(m_sub2**(B1+5.0_DS)-m_sub1**(B1+5.0_DS)) &
                       ))

            ! 2. con-weighted terminal velocity
            WVT(2,i,n)=WVT(2,i,n)+max(0.0_DS,&
                  A1*( a(1)/(B1+1.0_DS)*(m_sub2**(B1+1.0_DS)-m_sub1**(B1+1.0_DS)) &
                     + a(2)/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)) &
                     + a(3)/(B1+3.0_DS)*(m_sub2**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)) &
                     + a(4)/(B1+4.0_DS)*(m_sub2**(B1+4.0_DS)-m_sub1**(B1+4.0_DS)) &
                       ))

            aNT(i,n)=aNT(i,n)+aN
            aMT(i,n)=aMT(i,n)+aM
          endif

         endif no_error_if2
        end do
        end do
      enddo


!      do in=1,g%N_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
      do n = 1, g%L
      do i = 1, g%N_BIN
        if(error_number(i,n)==0) then
          i_den_gt0p5=0.5*(1.0-sign(1.0_PS,0.5_PS-g%MS(i,n)%den))
          v_max=real(i_den_gt0p5,PS_KIND)*5.0e+3_DS+&
                  (1.0-real(i_den_gt0p5,PS_KIND))*2.0e+3_DS

          WVT(1,i,n)=min(max(0.0_DS,WVT(1,i,n)/aMT(i,n)),v_max)
          WVT(2,i,n)=min(max(0.0_DS,WVT(2,i,n)/aNT(i,n)),v_max)
        elseif(error_number(i,n)<10) then
          WVT(1,i,n)=real(min(max(0.0_ps,g%MS(i,n)%vtm),1.5e+3_ps),8)
          WVT(2,i,n)=real(min(max(0.0_ps,g%MS(i,n)%vtm),1.5e+3_ps),8)
        endif
      enddo
      enddo

!dbg!CDIR NOVECTOR
!dbg      do in=1,g%N_bin*g%L
!dbg        n=(in-1)/g%N_BIN+1
!dbg        i=in-(n-1)*g%N_BIN
!dbg        if(error_number(i,n)==0) then
!dbg          if(WVT(1,i,n).gt.1.0e-2.and.WVT(2,i,n)<1.0e-5) then
!dbg            write(*,*) "wvt2 is 0",i,n,g%MS(i,n)%vtm,WVT(1:2,i,n) &
!dbg                ,aMT(i,n),aNT(i,n),a2d(i,n,1:4)
!dbg          endif
!dbg        endif
!dbg      enddo

    end  if


  end subroutine cal_wterm_vel_v3_vec

!!$  subroutine cal_growth_mode(level,mes_rc,ms,is, th_var)
!!$!tmp    use mod_amps_utility, only: get_growth_mode
!!$    use mod_amps_utility, only: get_growth_mode_max
!!$    integer,intent(in) :: level,mes_rc
!!$    type (Mass_Bin), intent(in) :: ms
!!$    type (Ice_Shape), intent(inout)  :: is
!!$    type (Thermo_Var), intent(in)  :: th_var
!!$    real(PS) :: svi
!!$
!!$    ! initialize
!!$    IS%growth_mode=0
!!$    IS%init_growth=0
!!$    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return
!!$
!!$    if(level==3.or.level==5.or.level==7) then
!!$!!c       if(max(IS%r,IS%e,MS%a_len,MS%c_len)<10.0e-4_PS.and.th_var%T<253.16) then
!!$       if(max(IS%r,IS%e,sqrt(MS%a_len**2.0+MS%c_len**2.0))<10.0e-4_PS.and.th_var%T<253.16_PS) then
!!$
!!$          svi=th_var%s_v(2)
!!$          if(mes_rc==2.or.mes_rc==4) then
!!$             ! if liquid co exist
!!$             svi=th_var%e_sat(1)/th_var%e_sat(2)-1.0_PS
!!$          end if
!!$!tmp          IS%growth_mode=get_growth_mode(th_var%T,svi)
!!$          IS%growth_mode=get_growth_mode_max(th_var%T,svi)
!!$          IS%init_growth=1
!!$
!!$!!c          IS%growth_mode=1
!!$
!!$       else
!!$          IS%init_growth=0
!!$          if(IS%habit==1.or.IS%habit==2) then
!!$             IS%growth_mode=2
!!$          else if(IS%habit==3) then
!!$             IS%growth_mode=3
!!$          else if(IS%habit==4) then
!!$             IS%growth_mode=4
!!$          else if(IS%habit==5) then
!!$             IS%growth_mode=5
!!$          else if(IS%habit>=6) then
!!$             IS%growth_mode=1
!!$          end if
!!$       end if
!!$    else
!!$       if(max(IS%r,IS%e,MS%a_len,MS%c_len)<10.0e-4_PS.and.th_var%T<253.16_PS) then
!!$          IS%init_growth=1
!!$       else
!!$          IS%init_growth=0
!!$       end if
!!$       if(IS%habit==1.or.IS%habit==2) then
!!$          IS%growth_mode=2
!!$       else if(IS%habit==3) then
!!$          IS%growth_mode=3
!!$       else
!!$          write(*,*) "habit is not right:level,habit",level,IS%habit
!!$          stop
!!$       end if
!!$    end if
!!$  end subroutine cal_growth_mode

  subroutine cal_growth_mode_vec(level,mes_rc,g,ag,icond1,rdsd,ihabit_gm_random)
    use scale_prc, only: &
       PRC_abort
    use mod_amps_utility, only: &
       cal_growth_mode_inl_vec, &
       random_genvar
    integer,intent(in) :: level
    integer,dimension(*),intent(in) :: mes_rc
    type (Group), intent(inout)  :: g
    type (AirGroup), intent(in)   :: ag
    integer,dimension(g%N_BIN,g%L),intent(in)    ::  icond1

    type(random_genvar),intent(inout) :: rdsd
    ! random generaion: 1, max frequency: 0
    integer, intent(in)           :: ihabit_gm_random

    real(PS),dimension(g%L) :: svi
    integer,dimension(g%N_bin,g%L) :: igm
!    integer,dimension(g%N_bin*g%L) :: ierror
    integer                     :: i,n

!    ierror=0
    if(level==3.or.level==5.or.level==7) then

!!c      write(*,*) "rdsd3:",rdsd
      do n=1,g%L
        svi(n)=ag%TV(n)%s_v(2)
        if(mes_rc(n)==2.or.mes_rc(n)==4) then
          ! if liquid co exist
          svi(n)=ag%TV(n)%e_sat(1)/ag%TV(n)%e_sat(2)-1.0_PS
        end if
      enddo
      call cal_growth_mode_inl_vec(igm,g%N_bin,g%L,ihabit_gm_random &
                          ,ag%TV(1:g%L)%T,svi,rdsd)
!!c      write(*,*) "rdsd4:",rdsd

!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
      do n = 1, g%L
      do i = 1, g%N_BIN
        if(icond1(i,n)==0) then
!!c         if(max(g%IS(i,n)%r,g%IS(i,n)%e,g%MS(i,n)%a_len,g%MS(i,n)%c_len)<10.0e-4_PS.and.ag%TV(n)%T<253.16) then
          if(max(g%IS(i,n)%r,g%IS(i,n)%e,sqrt(g%MS(i,n)%a_len**2.0+g%MS(i,n)%c_len**2.0))<10.0e-4_PS.and.ag%TV(n)%T<253.16_PS) then

            g%IS(i,n)%growth_mode=igm(i,n)
            g%IS(i,n)%init_growth=1


!!c            write(*,*) "ck growth_mode1:",i,n,mes_rc(n),g%IS(i,n)%init_growth,g%IS(i,n)%growth_mode,&
!!c                  ag%TV(n)%T,svi(n)

          else
            g%IS(i,n)%init_growth=0
            if(g%IS(i,n)%habit==1.or.g%IS(i,n)%habit==2) then
              g%IS(i,n)%growth_mode=2
            else if(g%IS(i,n)%habit==3) then
              g%IS(i,n)%growth_mode=3
            else if(g%IS(i,n)%habit==4) then
              g%IS(i,n)%growth_mode=4
            else if(g%IS(i,n)%habit==5) then
              g%IS(i,n)%growth_mode=5
            else if(g%IS(i,n)%habit>=6) then
              g%IS(i,n)%growth_mode=1
            end if
!!c            write(*,*) "ck growth_mode2:",i,n,mes_rc(n),g%IS(i,n)%init_growth,g%IS(i,n)%growth_mode,&
!!c                  ag%TV(n)%T,ag%TV(n)%s_v(2)
          end if
        endif
      enddo
      enddo
    else
!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then

          if(max(g%IS(i,n)%r,g%IS(i,n)%e,g%MS(i,n)%a_len,g%MS(i,n)%c_len)<10.0e-4_PS.and.ag%TV(n)%T<253.16_PS) then
            g%IS(i,n)%init_growth=1
          else
            g%IS(i,n)%init_growth=0
          end if
          if(g%IS(i,n)%habit==1.or.g%IS(i,n)%habit==2) then
            g%IS(i,n)%growth_mode=2
          else if(g%IS(i,n)%habit==3) then
            g%IS(i,n)%growth_mode=3
          else
             LOG_ERROR("cal_growth_mode_vec",*) "habit is not right:level,habit",i,n,level,g%IS(i,n)%habit
             call PRC_abort
          end if
        end if
      enddo
      enddo
    endif
  end subroutine cal_growth_mode_vec

!!$  subroutine cal_coef_vapdep(phase,level,ms,th_var,is)
!!$    integer, intent(in)  :: phase,level
!!$    type (Mass_Bin), intent(inout) :: ms
!!$    type (Thermo_Var), intent(in) :: th_var
!!$    type (Ice_Shape), optional       :: is
!!$
!!$    ! initialize
!!$    ms%coef=0.0_PS
!!$    if( ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return
!!$
!!$
!!$    if(phase==1) then
!!$       ms%coef(1)=4.0_PS*PI*th_var%GTP(phase)*ms%CAP*ms%fv*ms%fkn
!!$       ms%coef(2)=0.0_PS
!!$    elseif(phase==2) then
!!$       if(th_var%T<T_0.and.ms%inmlt==0) then
!!$          ms%coef(1)=  4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$               th_var%e_sat(phase)/th_var%T
!!$          ms%coef(2)= -4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$               ms%e_sat/ms%tmp
!!$       else
!!$          if(level>=6) then
!!$             if(ms%mass(imw)/ms%con>1.0e-15_PS) then
!!$                ms%coef(1)=  4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$                     th_var%e_sat(1)/th_var%T
!!$                ms%coef(2)= -4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$                     ms%e_sat/ms%tmp
!!$             else
!!$                ms%coef=0.0_PS
!!$             end if
!!$          else
!!$             ms%coef=0.0_PS
!!$          end if
!!$       end if
!!$    end if
!!$  end subroutine cal_coef_vapdep

  subroutine cal_coef_vapdep2(phase,level,ms,th_var,nu_aps,phi_aps,m_aps,is)
    use class_Mass_Bin, only: &
       get_osm, &
       get_critrad_itr, &
       get_hazerad_itr
    use class_Thermo_Var, only: &
       get_mod_diffusivity, &
       get_mod_thermal_cond
    integer, intent(in)  :: phase,level
    type (Mass_Bin), intent(inout) :: ms
    type (Thermo_Var), intent(in) :: th_var
    type (Ice_Shape), optional       :: is
    real(PS),dimension(*),intent(in) :: phi_aps,nu_aps,m_aps

    ! density, radius, and cubic radius of mixed aerosol
    real(PS) :: den_ap, r_n, r_n3
    ! coefficients for Kohler curve
    real(PS) :: AA,BB
    ! molality
    real(PS) :: molality
    ! modified diffusivity and conductivity
    real(PS) :: mD_v,mk_a
    ! ambient saturation vapor density and its first derivative wrt T
    real(PS) :: rho_s, rho_sp

    ! supersaturation ratio for droplet r at the ambient temperature
    real(PS) :: S_r

    ! moles of soluble fraction
    real(PS) :: n_aps

    ! saturation ratio when haze droplets go back to dry CCN
    real(PS) :: SRW=0.35_PS
    real(PS) :: rd_c

    ! initialize
    ms%coef=0.0_PS
    ms%r_crt=0.0_PS
    ms%r_act=0.0_PS
    if( ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return


    if(phase==1) then

!!$       ms%coef(1)=4.0_PS*PI*th_var%GTP(phase)*ms%CAP*ms%fv*ms%fkn
!!$       ms%coef(2)=0.0_PS
!!$       return

       ! This formulation is based on Davies (1985).
       !
       ! assuming that category 1 is the soluble material.

       ! calculate density of mixed aerosols
       den_ap=ms%den_ai/(1.0_PS-ms%eps_map*(1.0_PS-ms%den_ai/ms%den_as))
       r_n3=(ms%mass(rmat)/ms%con)/coef4pi3/den_ap
       r_n=r_n3**(1.0/3.0)

       ! unit is now in cm
       AA=2.0_PS*th_var%sig_wa/(R_v*th_var%T*den_w)

       ! unitless
       BB=nu_aps(1)*ms%eps_map*M_W*den_ap/(M_aps(1)*den_w)

       ! moles of soluble mass ( mols )
       ! used to get molality
       molality=1.0e+3_PS*ms%mass(rmas)/(ms%mass(rmt)-ms%mass(rmat))/M_aps(1)

       ! get supersaturation ratio
       if(r_n<1.0e-7) then
         ! in case of dry radius less than 1.0e-3 micron, ignore the solute effect
         S_r=exp(AA/ms%a_len)

       else
         S_r=exp(AA/ms%a_len-BB*get_osm(1,molality)/&
                max(ms%a_len**3/max(1.0e-25_PS,r_n3)-1.0_PS,1.0e-25_PS))

       endif
       mD_v = get_mod_diffusivity(phase,ms%a_len,th_var%D_v,th_var%T)
       mk_a = get_mod_thermal_cond(ms%a_len,th_var%K_a,th_var%T,th_var%den)

       rho_s=th_var%e_sat(1)/(th_var%T*R_v)
       rho_sp=rho_s/th_var%T*(L_e/(R_v*th_var%T)-1.0_PS)

       ms%coef(1)=mD_v*rho_sp*4.0_PS*PI*ms%CAP/(mk_a+L_e*mD_v*S_r*rho_sp)*&
            (mk_a*rho_s/rho_sp)*ms%fv

       ! you can add radiative heating/cooling term as coef(2) right here
       !          ms%coef(2)=ms%coef(1)*(1.0_PS-S_r)
       !         -(S_r*mD_v*rho_sp/(mk_a+L_e*mD_v*S_r*rho_sp)*Qr)

       ms%coef(2)=ms%coef(1)*(1.0_PS-S_r)


!!c       if(ms%a_len>1.0e-6_PS.and.ms%a_len<50.0e-4_PS) then
       if(ms%a_len>1.0e-7_PS.and.ms%a_len<50.0e-4_PS.and.r_n>1.0e-7_PS) then
          ! moles of soluble material (in mol um^3/cm^3)
          n_aps=ms%mass(rmas)/ms%con/M_aps(1)*1.0e+12

          ! calculate activation radius (input and output in micron meter)
!!c          ms%r_act=1.0e-4*get_critrad_itr(AA*1.0e+4,BB,n_aps,(r_n3**(1.0/3.0))*1.0e+4)
          rd_c=get_critrad_itr(AA*1.0e+4_PS,BB,n_aps,(r_n3**(1.0_PS/3.0_PS))*1.0e+4_PS)
          ms%r_crt=1.0e-4_PS*rd_c

          if(debug .and. rd_c<(r_n3**(1.0_PS/3.0_PS))*1.0e+4_PS) then
            write(*,*) "rd_c is smaller than dry radius",rd_c,(r_n3**(1.0_PS/3.0_PS))*1.0e+4_PS
          endif

          ! calcluate haze size with SRW.   This is used for evaporation only.
          ms%r_act=1.0e-4_PS*get_hazerad_itr(AA*1.0e+4_PS,BB,SRW,n_aps &
                  ,(r_n3**(1.0_PS/3.0_PS))*1.0e+4_PS,rd_c,0)

       else
          ! according to Chen 1992, JAS
          ms%r_crt=1.0e-4_PS*r_n
          ms%r_act=1.0e-4_PS*r_n

       endif

    elseif(phase==2) then
       if(th_var%T<T_0.and.ms%inmlt==0) then
          ! dry growth
          ms%coef(1)=  4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn
          ms%coef(2)= -ms%coef(1)
       else
          ! wet growth
          ms%coef(1)=  4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn
          ms%coef(2)= -ms%coef(1)
!!c          if(level>=6) then
!!c             if(ms%mass(imw)/ms%con>1.0e-15_PS) then
!!c                ms%coef(1)=  4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn
!!c                ms%coef(2)= -ms%coef(1)
!!c             else
!!c                ms%coef=0.0_PS
!!c             end if
!!c          else
!!c             ms%coef=0.0_PS
!!c          end if
       end if
    end if
  end subroutine cal_coef_vapdep2

  subroutine cal_coef_vapdep2_vec(phase,g,ag,icond1,nu_aps,phi_aps,m_aps)
    use class_Mass_Bin, only: &
       get_critrad_anal, &
       get_hazerad_anal
    integer, intent(in)  :: phase!,level
    type (Group), intent(inout)  :: g
    type (AirGroup), intent(in)   :: ag
    integer,dimension(g%N_BIN,g%L),intent(in)    ::  icond1
    real(PS),dimension(*),intent(in) :: phi_aps
    real(PS),dimension(*),intent(in) :: nu_aps,m_aps

    ! density and radius of mixed aerosol
    real(PS) :: den_ap, r_n
    ! cubic radius of mixed aerosol
    real(PS) :: r_n3

    ! coefficients for Kohler curve
    real(PS) :: AA,sb,beta!,BB
    ! molality
    !real(PS) :: molality
    ! modified diffusivity and conductivity
    !real(PS) :: mD_v,mk_a
    ! ambient saturation vapor density and its first derivative wrt T
    real(PS) :: rho_s!, rho_sp

    ! supersaturation ratio for droplet r at the ambient temperature
    !real(PS) :: S_r

    ! moles of soluble fraction
    !real(PS) :: n_aps

    ! saturation ratio when haze droplets go back to dry CCN
    real(PS) :: SRW=0.99_PS
    real(PS) :: rd_c

    ! supersaturation reduction due to salt effect
    real(PS) :: s_salt
    ! mean thermal velocity of vapor
    real(PS) :: vw
    real(PS) :: buzai_con
    ! psychrometric correction due to latent heat release
    real(PS) :: gamma_w

    !real(PS) :: dum1
    integer                     :: i,n

    if(phase==1) then
!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then

          ! This formulation is based on Khvorostyanov and Curry (2014).
          !
          ! assuming that category 1 is the soluble material.
          !
          ! Kohler's equation is expressed with constants A and B defined
          ! by (6.2.10) of Khvorostyanov and Curry (2014) for consistency

          ! calculate density of mixed aerosols
          den_ap=g%MS(i,n)%den_ai/(1.0_PS-g%MS(i,n)%eps_map*(1.0_PS-g%MS(i,n)%den_ai/g%MS(i,n)%den_as))
          r_n3=(g%MS(i,n)%mass(rmat)/g%MS(i,n)%con)/coef4pi3/den_ap
          r_n=r_n3**(1.0_PS/3.0_PS)

          ! unit is now in cm
          AA=2.0_PS*ag%TV(n)%sig_wa/(R_v*ag%TV(n)%T*den_w)

          ! unitless for soluble mass proportional to volume,      (beta=0.5)
          !       cm for soluble mass proportional to surface area (beta=0.0)
          sb=nu_aps(1)*g%MS(i,n)%eps_map*M_W*den_ap/(M_aps(1)*den_w)*phi_aps(1)
          beta=0.5_PS

          s_salt=  AA/g%MS(i,n)%a_len-sb*r_n**(2.0_PS*(1.0_PS+beta))/&
                   (g%MS(i,n)%a_len**3-r_n3)

          vw=sqrt(8.0_PS/PI*R_v*ag%TV(n)%T)
          buzai_con=4.0_PS*ag%TV(n)%D_v/(vw*a_cliq)

          rho_s=ag%TV(n)%e_sat(1)/(ag%TV(n)%T*R_v)

          gamma_w=1.0_PS+L_e*rho_s/(c_pa*ag%TV(n)%T*ag%TV(n)%den)* &
               (L_e/R_v/ag%TV(n)%T-1.0)

          g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%D_v*rho_s/gamma_w * &
                g%MS(i,n)%CAP*g%MS(i,n)%CAP/(g%MS(i,n)%CAP+buzai_con) *&
                g%MS(i,n)%fv

          g%MS(i,n)%coef(2)=-g%MS(i,n)%coef(1)*s_salt


!!c       if(g%MS(i,n)%a_len>1.0e-6_PS.and.g%MS(i,n)%a_len<50.0e-4_PS) then
          if(g%MS(i,n)%a_len>1.0e-7_PS.and.r_n>1.0e-7_PS) then
!org          if(g%MS(i,n)%a_len>1.0e-7_PS.and.g%MS(i,n)%a_len<50.0e-4_PS.and.r_n>1.0e-7_PS) then
            !
            !   calculation of ciritical radius and haze at SRW is done with
            !    formula by Khvorostyanov and Curry (2014)
            !
            ! moles of soluble material (in mol um^3/cm^3)
!tmp            n_aps=g%MS(i,n)%mass(rmas)/g%MS(i,n)%con/M_aps(1)*1.0e+12

            ! calculate activation radius (input and output in micron meter for itr)
!!c          g%MS(i,n)%r_act=1.0e-4*get_critrad_itr(AA*1.0e+4,BB,n_aps,(r_n3**(1.0/3.0))*1.0e+4)
!            rd_c=get_critrad_itr(AA*1.0e+4,BB,n_aps,(r_n3**(1.0/3.0))*1.0e+4)

            rd_c=get_critrad_anal(AA,sb,beta,r_n)  ! [cm]

            g%MS(i,n)%r_crt=rd_c

!tmp           if(rd_c<(r_n3**(1.0/3.0))*1.0e+4) then
!tmp              write(*,*) "rd_c is smaller than dry radius",rd_c,(r_n3**(1.0/3.0))*1.0e+4
!tmp            endif

            ! calcluate haze size with SRW.   This is used for evaporation only.
!!c            g%MS(i,n)%r_act=1.0e-4_PS*get_hazerad_itr(AA*1.0e+4_PS,BB,SRW,n_aps &
!!c                  ,(r_n3**(1.0/3.0))*1.0e+4,rd_c,0)
            g%MS(i,n)%r_act=get_hazerad_anal(AA,sb,beta,SRW,r_n)

          else
            ! according to Chen 1992, JAS
            g%MS(i,n)%r_crt=r_n
            g%MS(i,n)%r_act=r_n

          endif
        endif
      enddo
      enddo
    elseif(phase==2) then
!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then
          if(ag%TV(n)%T<T_0.and.g%MS(i,n)%inmlt==0) then
            ! dry growth
            g%MS(i,n)%coef(1)=  4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn
            g%MS(i,n)%coef(2)= -g%MS(i,n)%coef(1)
          else
            ! wet growth
            g%MS(i,n)%coef(1)=  4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn
            g%MS(i,n)%coef(2)= -g%MS(i,n)%coef(1)
!!c          if(level>=6) then
!!c             if(g%MS(i,n)%mass(imw)/g%MS(i,n)%con>1.0e-15_PS) then
!!c                g%MS(i,n)%coef(1)=  4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn
!!c                g%MS(i,n)%coef(2)= -g%MS(i,n)%coef(1)
!!c             else
!!c                g%MS(i,n)%coef=0.0_PS
!!c             end if
!!c          else
!!c             g%MS(i,n)%coef=0.0_PS
!!c          end if
          end if
        end if
      enddo
      enddo
    end if
  end subroutine cal_coef_vapdep2_vec

!!$  subroutine cal_coef_vapdep_ap(level,ica,ms,th_var)
!!$    use class_Thermo_Var, only: &
!!$       get_fkn
!!$    integer, intent(in)  :: level,ica
!!$    type (Mass_Bin), intent(inout) :: ms
!!$    type (Thermo_Var), intent(in) :: th_var
!!$    integer :: i,n
!!$    ! coefficient for capacitance by Chiruta and Wang (2003)
!!$!!c    real(PS),parameter   :: CAP_ros=0.434*N_lob**0.257
!!$    real(PS),parameter   :: CAP_ros=0.619753727
!!$    real(PS) :: fkn
!!$
!!$    ! initialize
!!$    ms%coef=0.0_PS
!!$    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return
!!$
!!$    if(ica==2) then
!!$       ! +++ calculate kinetic effect +++
!!$       fkn=get_fkn(th_var,2,ms%a_len)
!!$       if(level<=2.or.level==4.or.level==6) then
!!$          ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*ms%a_len*fkn
!!$       elseif(level==3.or.level==5.or.level==7) then
!!$          if(th_var%T>-20.0+273.16) then
!!$             ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*ms%a_len*fkn
!!$          else
!!$             if(th_var%nuc_gmode==2.or.th_var%nuc_gmode==3) then
!!$                ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*ms%a_len*fkn
!!$             elseif(th_var%nuc_gmode==4) then
!!$                ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*ms%a_len*fkn
!!$!!c                ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*CAP_ros*ms%a_len*fkn
!!$!!c             if(th_var%nuc_gmode==4) then
!!$!!c                ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*CAP_ros*ms%a_len*fkn
!!$             else
!!$                ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*ms%a_len*fkn
!!$!!c                ms%coef(1)=get_dmdt0(th_var%nuc_gmode,th_var%T)
!!$             end if
!!$          end if
!!$       end if
!!$    end if
!!$  end subroutine cal_coef_vapdep_ap

  subroutine cal_coef_vapdep_ap_vec(level,ica,mes_rc,g,ag,icond1)
    use class_Thermo_Var, only: &
       get_fkn
    use mod_amps_utility, only: &
       get_growth_mode
!tmp    use mod_amps_utility, only: get_growth_mode_max
    integer,intent(in) :: level,ica
    integer,dimension(*),intent(in) :: mes_rc
    type (Group), intent(inout)  :: g
    type (AirGroup), intent(inout)   :: ag
    integer,dimension(g%N_BIN,g%L),intent(in)    ::  icond1

    ! coefficient for capacitance by Chiruta and Wang (2003)
!!c    real(PS),parameter   :: CAP_ros=0.434*N_lob**0.257
    real(PS),parameter   :: CAP_ros=0.619753727
    real(PS) :: fkn

    integer                     :: i,n

    if(level==3.or.level==5.or.level==7) then
!      do in=1,g%n_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then

          if((mes_rc(n)==2.or.mes_rc(n)==4).and.ag%TV(n)%T<253.16_PS) then
          ! if liquid co exist
            ag%TV(n)%nuc_gmode=get_growth_mode(ag%TV(n)%T,ag%TV(n)%e_sat(1)/ag%TV(n)%e_sat(2)-1.0_PS)
!tmp            ag%TV(n)%nuc_gmode=get_growth_mode_max(ag%TV(n)%T,ag%TV(n)%e_sat(1)/ag%TV(n)%e_sat(2)-1.0_PS)

!!c          ag%TV(n)%nuc_gmode=1
          endif

          if(ica==2) then
          ! +++ calculate kinetic effect +++
            fkn=get_fkn(ag%TV(n),2,g%MS(i,n)%a_len)
            if(ag%TV(n)%T>-20.0+273.16) then
              g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%gtp(2)*g%MS(i,n)%a_len*fkn
            else
              if(ag%TV(n)%nuc_gmode==2.or.ag%TV(n)%nuc_gmode==3) then
                g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%gtp(2)*g%MS(i,n)%a_len*fkn
              elseif(ag%TV(n)%nuc_gmode==4) then
                g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%gtp(2)*g%MS(i,n)%a_len*fkn
!!c                g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%gtp(2)*CAP_ros*g%MS(i,n)%a_len*fkn
!!c             if(ag%TV(n)%nuc_gmode==4) then
!!c                g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%gtp(2)*CAP_ros*g%MS(i,n)%a_len*fkn
              else
                g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%gtp(2)*g%MS(i,n)%a_len*fkn
!!c                g%MS(i,n)%coef(1)=get_dmdt0(ag%TV(n)%nuc_gmode,ag%TV(n)%T)
              end if
            end if
          endif
        end if
      enddo
      enddo
    elseif(level<=2.or.level==4.or.level==6) then
      if(ica==2) then
!        do in=1,g%n_bin*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
         do n = 1, g%L
         do i = 1, g%N_BIN
          if(icond1(i,n)==0) then
            ! +++ calculate kinetic effect +++
            fkn=get_fkn(ag%TV(n),2,g%MS(i,n)%a_len)

            g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%gtp(2)*g%MS(i,n)%a_len*fkn
          end if
        enddo
        enddo
      endif
    end if
  end subroutine cal_coef_vapdep_ap_vec

  subroutine cal_needgive(g_1, g_2, ngrid, m_t1, m_t2, need_m, give_m)
    ! group 1 is the higher group, and group 2 is lower one.
    ! (lower) - cloud drop - rain - solid_hydro (higher)
    type (group), intent(inout)         :: g_1, g_2
    ! grid number
    integer, intent(in)       :: ngrid
    ! bugdet (total) of mass of the group
    real(ps), intent(inout)      :: m_t1, m_t2
    ! new bugdet (total) of mass of the group
    real(ps)      :: nm_t2
    ! necessary mass transfer from vapor
    real(ps), intent(inout)   :: need_m, give_m
    ! message from cal_needgiv subroutine
    ! 0. variables were not changed.
    ! 1. mass variables were changed.
    ! 2. mass and concentration were changed.
    ! 3. no concentration and mass defined for whole the group
    ! 4. no mass tendency obtained for the bin
    !    (this is defined in cal_mass_budget)
    ! 5. no con tendency obtained for the bin
    !    (this is defined in cal_mass_budget)
    ! 6. no mass and con tendency obtained for the bin
    !    (this is defined in cal_mass_budget)

    !integer                   :: mes_g1, mes_g2
    real(ps), parameter           :: min_mt = 1.0e-20

    ! initial mass of cloud
    real(ps), parameter       :: ini_mcloud = (4.0_ps*pi/3.0_ps)*1.0e-9

    real(ps)   :: mod_dum, total_pos, total_neg
    real(ps),dimension(g_1%n_bin) :: old_mass,modrg
    integer    :: i,j,k
    integer,dimension(g_1%n_bin) :: iswitch


    ! calculate total positive and negative mass
    total_pos = 0.0_ps
    total_neg = 0.0_ps
    do i = 1, g_1%n_bin
      if( g_1%ms(i,ngrid)%mass(1) < 0.0_ps) then
        total_neg = total_neg + (-g_1%ms(i,ngrid)%mass(1))
      else
        total_pos = total_pos + g_1%ms(i,ngrid)%mass(1)
      end if
    end do

    if( total_neg == 0.0_ps .and. total_pos == 0.0_ps ) then
      ! --- case of no hydrometeors in the group ---
      g_1%mark_cm(ngrid) = 3
    else if( total_neg == 0.0_ps .and. total_pos > 0.0_ps ) then
      ! check wheather the mean mass is confined by the
      ! bin limits
      modrg=1.0
      call check_con(g_1,ngrid,modrg)
    else if( total_neg > 0.0_ps .and. total_pos > 0.0_ps ) then

      ! check each bin from mass
      if( m_t1 >= min_mt ) then
        ! --- cases where the group can be fixed within the group ---

        iswitch=0
        do i=1,g_1%n_bin
          old_mass(i)=g_1%ms(i,ngrid)%mass(1)
        enddo
        do i=1,g_1%n_bin
          if( g_1%ms(i,ngrid)%mass(1) < 0.0 ) then
            g_1%ms(i,ngrid)%mark = 1

            mod_dum = (total_pos-(-g_1%ms(i,ngrid)%mass(1)))/total_pos
            total_pos = total_pos-(-g_1%ms(i,ngrid)%mass(1))
            do j = 1, g_1%n_bin
              if( i == j ) then
                iswitch(i)=1
                g_1%ms(j,ngrid)%mass(1)=0.0_PS
              else
                if( g_1%ms(j,ngrid)%mass(1) > 0.0_ps ) then
                  iswitch(i)=2
                  g_1%ms(j,ngrid)%mass(1)=g_1%ms(j,ngrid)%mass(1)*mod_dum
                  g_1%ms(j,ngrid)%mark = 1
!!c                  write(*,*) "cal_needgive: mass fixed 1",mod_dum
                end if
              end if
            end do
          end if
        end do

        do i=1,g_1%n_bin
          modrg(i)=1.0_PS
          if(old_mass(i)>1.0e-30) then
            modrg(i)=g_1%ms(i,ngrid)%mass(1)/old_mass(i)
          endif
        enddo


!        do ik=1,g_1%N_masscom*g_1%n_bin
!          k=(ik-1)/g_1%n_bin+1
!          i=ik-(k-1)*g_1%N_BIN
        do i = 1, g_1%N_BIN
          if(iswitch(i)>0) then
             do k = 1, g_1%N_masscom
                g_1%ms(i,ngrid)%mass(1+k)=g_1%ms(i,ngrid)%mass(1+k)&
                     *modrg(i)
             enddo
          endif
        enddo

        ! check
        ! calculate total positive and negative mass
        total_pos = 0.0_ps
        total_neg = 0.0_ps
        do i = 1, g_1%n_bin
          if( g_1%ms(i,ngrid)%mass(1) < 0.0_ps) then
            total_neg = total_neg + (-g_1%ms(i,ngrid)%mass(1))
          else
            total_pos = total_pos + g_1%ms(i,ngrid)%mass(1)
          end if
        end do
        ! change the concentration depending on the new total mass
        call check_con(g_1,ngrid,modrg)

      else if( m_t1 < 0.0_ps ) then
        ! --- cases where the group need to borrow mass from other groups ---
        iswitch=0
        do i = 1, g_1%n_bin
          modrg(i)=1.0
          if( g_1%ms(i,ngrid)%mass(1) < 0.0_ps) then
            iswitch(i)=1
            g_1%ms(i,ngrid)%mass(1) = 0.0_ps
            g_1%ms(i,ngrid)%con = 0.0_ps
            g_1%ms(i,ngrid)%mark = 2
            modrg(i)=0.0
          end if
        end do
!        do ik=1,g_1%N_masscom*g_1%n_bin
!          k=(ik-1)/g_1%n_bin+1
!          i=ik-(k-1)*g_1%N_BIN
        do i = 1, g_1%N_BIN
           if(iswitch(i)>0) then
              do k = 1, g_1%N_masscom
                 g_1%ms(i,ngrid)%mass(1+k)=0.0_PS
              enddo
          endif
        enddo
        call check_con(g_1,ngrid,modrg)

        if( m_t2 > min_mt ) then
          ! borrow from the second group
          if( total_neg <= m_t2 ) then
            nm_t2 = m_t2 - total_neg
          else if( total_neg > m_t2 ) then
            need_m = need_m + (total_neg - (m_t2-min_mt) )
            nm_t2 = min_mt
          end if
          mod_dum = nm_t2/m_t2
          m_t2 = nm_t2
          modrg=mod_dum

!          do ik=1,(1+g_2%N_masscom)*g_2%n_bin
!            k=(ik-1)/g_2%n_bin+1
!            i=ik-(k-1)*g_2%N_BIN
          do i = 1, g_2%N_BIN
          do k = 1, 1+g_2%N_masscom
            g_2%ms(i,ngrid)%mass(k)=g_2%ms(i,ngrid)%mass(k)*mod_dum
          end do
          end do
          ! change the concentration depending on the new total mass
          call check_con(g_2,ngrid,modrg)
        else
          ! get mass from vapor
          need_m = need_m + total_neg
        end if
      else
        do i=1,g_1%n_bin
          g_1%ms(i,ngrid)%con=0.0_ps
          modrg(i)=0.0_PS
        end do
!        do ik=1,(1+g_1%N_masscom)*g_1%n_bin
!          k=(ik-1)/g_1%n_bin+1
!          i=ik-(k-1)*g_1%N_BIN
        do i = 1, g_1%N_BIN
        do k = 1, (1+g_1%N_masscom)
          g_1%ms(i,ngrid)%mass(k)=0.0_PS
        enddo
        enddo
        call check_con(g_1,ngrid,modrg)

        g_1%mark_cm(ngrid) = 3
        ! give the mass to vapor
        give_m = give_m + m_t1
      end if
    end if
  end subroutine cal_needgive


  subroutine check_con( g, j, modrg)
    use scale_prc, only: &
       PRC_abort
    implicit none
    type (group), intent(inout)         :: g
    integer,intent(in)            :: j
    real(PS),dimension(*) :: modrg
    ! mean mass and middle mass
    real(ps)       :: mean_mass, mid_mass, old_con,new_con
    real(ps)       :: a_new,c_new
    integer        :: i,k

    ! minimum mass of ice crystal with 1 um radius
    real(ps),parameter :: m_icmin=4.763209003e-12

    ! realistic bounds for length predictions
    real(ps),parameter :: max_exice=1.0

    ! possible ratio of mean mass to bin boundaries
    real(ps),parameter :: brat1=1.001,brat2=0.999

    real(PS),dimension(g%n_bin) :: mod_ratio
    integer,dimension(g%n_bin) :: ierror,iswitch

    do i=1,g%n_bin
      old_con = g%ms(i,j)%con
      g%ms(i,j)%con=g%MS(i,j)%con*modrg(i)

      if(g%token==2) then
        g%is(i,j)%v_cs=g%is(i,j)%v_cs*modrg(i)
        a_new=max(1.0e-4_PS,min(5.0_PS,((g%ms(i,j)%a_len**3)*modrg(i))**0.33333333))
        g%ms(i,j)%a_len=a_new
        c_new=max(1.0e-4_PS,min(5.0_PS,((g%ms(i,j)%c_len**3)*modrg(i))**0.33333333))
        g%ms(i,j)%c_len=c_new

        g%is(i,j)%d=max(0.0_PS,min(0.9_PS*a_new,((g%is(i,j)%d**3)*modrg(i))**0.33333333))
        g%is(i,j)%ag=max(0.0_PS,min(a_new,((g%is(i,j)%ag**3)*modrg(i))**0.33333333))
        g%is(i,j)%cg=max(0.0_PS,min(c_new,((g%is(i,j)%cg**3)*modrg(i))**0.33333333))
        g%is(i,j)%n_exice=max(0.0_PS,min(max_exice,g%IS(i,j)%n_exice*modrg(i)))

      endif
    enddo

    ierror=0
    iswitch=0
    do i=1,g%n_bin
      old_con = g%ms(i,j)%con
      if( g%ms(i,j)%con > 0.0_ps .and. g%ms(i,j)%mass(1) > 0.0_ps) then
        mean_mass = g%ms(i,j)%mass(1)/g%ms(i,j)%con
        if( mean_mass < brat1*g%binb(i) .or. mean_mass > brat2*g%binb(i+1) ) then
!!c             mid_mass =(g%binb(i)+g%binb(i+1))/2.0_ps
          iswitch(i)=1

          mid_mass =max(brat1*g%binb(i),min(mean_mass,brat2*g%binb(i+1)))
!tmp          if(mid_mass<1.0e-30) then
!tmp            write(*,*) "ckcon",mean_mass,g%binb(i),g%binb(i+1),g%ms(i,j)%mass,g%ms(i,j)%con
!tmp          endif
          new_con=g%ms(i,j)%mass(1)/mid_mass
          g%ms(i,j)%con=new_con
          mod_ratio(i)=new_con/max(old_con,1.0e-30_RP)

          if(g%token==2) then
            if(g%ms(i,j)%mass(imc)<new_con*m_icmin) then
              g%ms(i,j)%mass(imc)=new_con*m_icmin
!!c                   write(*,*) "check_con> ice_mass fixed"
            end if
            g%is(i,j)%v_cs=g%is(i,j)%v_cs*mod_ratio(i)
            a_new=max(1.0e-4_PS,min(5.0_PS,((g%ms(i,j)%a_len**3)*mod_ratio(i))**0.33333333))
            g%ms(i,j)%a_len=a_new
            c_new=max(1.0e-4_PS,min(5.0_PS,((g%ms(i,j)%c_len**3)*mod_ratio(i))**0.33333333))
            g%ms(i,j)%c_len=c_new

            g%is(i,j)%d=max(0.0_PS,min(0.9_PS*a_new,((g%is(i,j)%d**3)*mod_ratio(i))**0.33333333))
            g%is(i,j)%ag=max(0.0_PS,min(a_new,((g%is(i,j)%ag**3)*mod_ratio(i))**0.33333333))
            g%is(i,j)%cg=max(0.0_PS,min(c_new,((g%is(i,j)%cg**3)*mod_ratio(i))**0.33333333))
            g%is(i,j)%n_exice=max(0.0_PS,min(max_exice,g%IS(i,j)%n_exice*mod_ratio(i)))

            if(g%IS(i,j)%V_cs<g%MS(i,j)%mass(1)/new_con/den_i) then
              g%IS(i,j)%V_cs=g%MS(i,j)%mass(1)/new_con/den_i
            endif
          end if
          g%ms(i,j)%mark = 2
        end if
      else if( g%ms(i,j)%con /= 0.0_ps .and. g%ms(i,j)%mass(1) == 0.0_ps) then
        iswitch(i)=2
        g%ms(i,j)%con = 0.0_ps
        if(g%token==2) then
          g%is(i,j)%v_cs=0.0_PS
          g%ms(i,j)%a_len=0.0_PS
          g%ms(i,j)%c_len=0.0_PS
          g%is(i,j)%d=0.0_PS
          g%is(i,j)%r=0.0_PS
          g%is(i,j)%e=0.0_PS
          g%is(i,j)%ag=0.0_PS
          g%is(i,j)%cg=0.0_PS
          g%is(i,j)%n_exice=0.0_PS
        end if
        g%ms(i,j)%mark = 3
      else if( g%ms(i,j)%con <= 0.0_ps .and. g%ms(i,j)%mass(1) > 0.0_ps) then
        iswitch(i)=3
        g%ms(i,j)%con=0.0_PS
        if(g%token==2) then
          g%is(i,j)%v_cs=0.0_PS
          g%ms(i,j)%a_len=0.0_PS
          g%ms(i,j)%c_len=0.0_PS
          g%is(i,j)%d=0.0_PS
          g%is(i,j)%r=0.0_PS
          g%is(i,j)%e=0.0_PS
          g%is(i,j)%ag=0.0_PS
          g%is(i,j)%cg=0.0_PS
          g%is(i,j)%n_exice=0.0_PS
        end if
        g%ms(i,j)%mark = 3
      else if( g%ms(i,j)%con < 0.0_ps ) then
        ierror(i)=1
      end if
    end do

!    do ik=1,g%N_vol*g%n_bin
!      k=(ik-1)/g%n_bin+1
!      i=ik-(k-1)*g%N_BIN
    do i = 1, g%N_BIN
       if(iswitch(i)==1) then
          do k = 1, g%N_vol
             g%ms(i,j)%vol(k)=g%ms(i,j)%vol(k)*mod_ratio(i)
          end do
       elseif(iswitch(i)==2.or.iswitch(i)==3) then
          do k = 1, g%N_vol
             g%ms(i,j)%vol(k)=0.0_PS
          end do
       endif
    enddo

!    do ik=1,(1+g%N_masscom)*g%n_bin
!      k=(ik-1)/g%n_bin+1
!      i=ik-(k-1)*g%N_BIN
    do i = 1, g%N_BIN
       if(iswitch(i)==2.or.iswitch(i)==3) then
          do k = 1, 1+g%N_masscom
             g%ms(i,j)%mass(k)=0.0_ps
          end do
       endif
    enddo

    if(any(ierror>0)) then
      do i=1,g%n_bin
        if(ierror(i)==1) then
           LOG_ERROR("check_con",*) "check_con > (4) something is wrong at ",i,j
           LOG_ERROR_CONT(*) "con,mass",g%ms(i,j)%con,g%ms(i,j)%mass(1)
        endif
      enddo
      call PRC_abort
    endif

  end subroutine check_con

!!$  subroutine check_con_scl( g, j)
!!$    implicit none
!!$    type (group), intent(inout)         :: g
!!$    integer,intent(in)            :: j
!!$    ! mean mass and middle mass
!!$    real(ps)       :: mean_mass, mid_mass, old_con,mod_ratio
!!$    integer        :: message
!!$    integer        :: i
!!$
!!$    ! minimum mass of ice crystal with 1 um radius
!!$    real(ps),parameter :: m_icmin=4.763209003e-12
!!$
!!$    ! realistic bounds for length predictions
!!$    real(ps),parameter :: max_exice=1.0
!!$
!!$    ! possible ratio of mean mass to bin boundaries
!!$    real(ps),parameter :: brat1=1.001,brat2=0.999
!!$
!!$    do i = 1, g%n_bin
!!$       message = 0
!!$       old_con = g%ms(i,j)%con
!!$       if( g%ms(i,j)%con > 0.0_ps .and. g%ms(i,j)%mass(1) > 0.0_ps) then
!!$          mean_mass = g%ms(i,j)%mass(1)/g%ms(i,j)%con
!!$          if( mean_mass < brat1*g%binb(i) .or. mean_mass > brat2*g%binb(i+1) ) then
!!$!!c             mid_mass =(g%binb(i)+g%binb(i+1))/2.0_ps
!!$             mid_mass =max(brat1*g%binb(i),min(mean_mass,brat2*g%binb(i+1)))
!!$             g%ms(i,j)%con = g%ms(i,j)%mass(1)/mid_mass
!!$             mod_ratio=old_con/max(g%ms(i,j)%con,1.0e-30_RP)
!!$             if(g%token==2) then
!!$                if(g%ms(i,j)%mass(imc)<g%ms(i,j)%con*m_icmin) then
!!$                   g%ms(i,j)%mass(imc)=g%ms(i,j)%con*m_icmin
!!$!!c                   write(*,*) "check_con> ice_mass fixed"
!!$                end if
!!$                g%is(i,j)%v_cs=g%is(i,j)%v_cs*mod_ratio
!!$                g%ms(i,j)%a_len=max(1.0e-4_PS,min(5.0_PS,((g%ms(i,j)%a_len**3.0)*mod_ratio)**0.33333333))
!!$                g%ms(i,j)%c_len=max(1.0e-4_PS,min(5.0_PS,((g%ms(i,j)%c_len**3.0)*mod_ratio)**0.33333333))
!!$                g%is(i,j)%d=max(0.0_PS,min(0.9_PS*g%MS(i,j)%a_len,((g%is(i,j)%d**3.0)*mod_ratio)**0.33333333))
!!$                g%ms(i,j)%vol=g%ms(i,j)%vol*mod_ratio
!!$                g%is(i,j)%ag=max(0.0_PS,min(g%MS(i,j)%a_len,((g%is(i,j)%ag**3.0)*mod_ratio)**0.33333333))
!!$                g%is(i,j)%cg=max(0.0_PS,min(g%MS(i,j)%c_len,((g%is(i,j)%cg**3.0)*mod_ratio)**0.33333333))
!!$                g%is(i,j)%n_exice=max(0.0_PS,min(max_exice,g%IS(i,j)%n_exice*mod_ratio))
!!$
!!$!!c                write(*,'("check_con> vol and len fixed",2I5,10ES15.6)') i,j,mod_ratio,&
!!$!!c                   g%IS(i,j)%V_cs,g%MS(i,j)%mass(1)/g%MS(i,j)%con,&
!!$!!c                   g%ms(i,j)%con,g%ms(i,j)%mass(1),g%ms(i,j)%a_len,g%ms(i,j)%c_len
!!$
!!$
!!$                if(g%IS(i,j)%V_cs<g%MS(i,j)%mass(1)/g%MS(i,j)%con/den_i) then
!!$!!c                   write(*,'("check_con> vol is wrong: bin,grid,vcs,vcsmin,modrat", 2I5,3ES15.6)') i,j,g%IS(i,j)%V_cs,g%MS(i,j)%mass(1)/g%MS(i,j)%con/den_i,mod_ratio
!!$                   g%IS(i,j)%V_cs=g%MS(i,j)%mass(1)/g%MS(i,j)%con/den_i
!!$                endif
!!$!!c                if(g%MS(i,j)%mass(imc)/g%MS(i,j)%con/(4.0*PI/3.0*(g%MS(i,j)%a_len**2+g%MS(i,j)%c_len**2)**1.5)&
!!$!!c                         <1.0e-3) then
!!$!!c                   write(*,'("check_con> small den",2I5,10ES15.6)') i,j,g%ms(i,j)%mass(imc),g%ms(i,j)%con,g%ms(i,j)%a_len,g%ms(i,j)%c_len,&
!!$!!c                      g%MS(i,j)%mass(imc)/g%MS(i,j)%con,&
!!$!!c                      g%MS(i,j)%mass(imc)/g%MS(i,j)%con/(4.0*PI/3.0*(g%MS(i,j)%a_len**2+g%MS(i,j)%c_len**2)**1.5)
!!$!!c                endif
!!$             end if
!!$             message = 2
!!$!!c             write(*,*) " check_con > concentration was fixed (1) at", i,j
!!$!!c             write(*,'(a18,e15.6,a4, e15.6)') &
!!$!!c                  "           > from ", old_con, " to ", g%ms(i,j)%con
!!$!!c             write(*,*) "           > the token is ", g%token
!!$          end if
!!$       else if( g%ms(i,j)%con /= 0.0_ps .and. g%ms(i,j)%mass(1) == 0.0_ps) then
!!$          g%ms(i,j)%con = 0.0_ps
!!$          g%ms(i,j)%mass=0.0_ps
!!$          if(g%token==2) then
!!$             g%is(i,j)%v_cs=0.0_PS
!!$             g%ms(i,j)%a_len=0.0_PS
!!$             g%ms(i,j)%c_len=0.0_PS
!!$             g%is(i,j)%d=0.0_PS
!!$             g%is(i,j)%r=0.0_PS
!!$             g%is(i,j)%e=0.0_PS
!!$             g%ms(i,j)%vol=0.0_PS
!!$             g%is(i,j)%ag=0.0_PS
!!$             g%is(i,j)%cg=0.0_PS
!!$             g%is(i,j)%n_exice=0.0_PS
!!$          end if
!!$          message = 3
!!$!!c          write(*,*) " check_con > concentration was fixed (2) at", i,j
!!$!!c          write(*,'(a18,e15.6,a4, e15.6)') &
!!$!!c               "           > from ", old_con, " to ", g%ms(i,j)%con
!!$!!c          write(*,*) "           > the token is ", g%token
!!$       else if( g%ms(i,j)%con <= 0.0_ps .and. g%ms(i,j)%mass(1) > 0.0_ps) then
!!$!!c          mid_mass = (g%binb(i)+g%binb(i+1))/2.0_ps
!!$!!c          g%ms(i)%con = g%ms(i)%mass(1)/mid_mass
!!$!!c          message = 2
!!$!!c          write(*,*) " check_con > concentration was fixed (3) at", i
!!$!!c          write(*,'(a18, e15.6,a4, e15.6)') &
!!$!!c               "           > from ", old_con, " to ", g%ms(i)%con
!!$!!c          write(*,*) "           > the token is ", g%token
!!$!!c          write(*,*) " check_con > con and mass were fixed to 0.0 (3) at", i,j
!!$!!c          write(*,'(a18, e15.6,a4, e15.6)') &
!!$!!c               "           > from con", old_con, " to ", 0.0
!!$!!c          write(*,'(a18, e15.6,a4, e15.6)') &
!!$!!c               "           > from mass", g%ms(i,j)%mass(1), " to ", 0.0
!!$          g%ms(i,j)%con=0.0_PS
!!$          g%ms(i,j)%mass=0.0_PS
!!$          if(g%token==2) then
!!$             g%is(i,j)%v_cs=0.0_PS
!!$             g%ms(i,j)%a_len=0.0_PS
!!$             g%ms(i,j)%c_len=0.0_PS
!!$             g%is(i,j)%d=0.0_PS
!!$             g%is(i,j)%r=0.0_PS
!!$             g%is(i,j)%e=0.0_PS
!!$             g%ms(i,j)%vol=0.0_PS
!!$             g%is(i,j)%ag=0.0_PS
!!$             g%is(i,j)%cg=0.0_PS
!!$             g%is(i,j)%n_exice=0.0_PS
!!$          end if
!!$          message = 3
!!$       else if( g%ms(i,j)%con < 0.0_ps ) then
!!$          write(*,*) "check_con > (4) something is wrong at ",i,j
!!$          write(*,*) "con,mass",g%ms(i,j)%con,g%ms(i,j)%mass(1)
!!$          stop
!!$       end if
!!$       g%ms(i,j)%mark = message
!!$    end do
!!$!!c    g%mark_cm(j) = message
!!$  end subroutine check_con_scl

!!$  subroutine check_length( g,n)
!!$    type (group), intent(inout)         :: g
!!$    integer, intent(in) :: n
!!$    ! axis ratio
!!$    real(ps)                           :: phi
!!$    ! possible axis ratio of solid hydrometeor
!!$    real(ps), parameter                :: max_phi = 1.0e+3
!!$    real(ps), parameter                :: min_phi = 1.0e-3
!!$    ! mean mass and middle mass
!!$    real(ps)                           :: mean_mass
!!$    ! possible minimum length by considering the density
!!$    real(ps)                           :: min_alen
!!$    ! possible maximum length by considering the density
!!$    real(ps)                           :: max_alen
!!$    ! possible maximum density of ice
!!$    real(ps), parameter                :: den_max = 0.91668
!!$    ! possible minimum density of ice
!!$    real(ps), parameter                :: den_min = 1.0e-8
!!$
!!$    integer    :: i
!!$    do i = 1, g%n_bin
!!$       if( g%MS(i,n)%con > 0.0_ps .and. g%MS(i,n)%mass(1) > 0.0_ps) then
!!$          mean_mass = g%MS(i,n)%mass(1)/g%MS(i,n)%con
!!$          if( g%MS(i,n)%a_len == 0.0_ps ) then
!!$             phi = max_phi
!!$          else
!!$             phi = min( max( g%MS(i,n)%c_len/g%MS(i,n)%a_len, min_phi), max_phi)
!!$          end if
!!$          min_alen = (3.0_ps*mean_mass/(4.0_ps*pi*phi*den_max))**(1.0/3.0)
!!$          max_alen = (3.0_ps*mean_mass/(4.0_ps*pi*phi*den_min))**(1.0/3.0)
!!$
!!$          if( g%MS(i,n)%a_len < min_alen ) then
!!$!!c             write(*,*) " check_length > minimum alen was used in ", i
!!$             g%MS(i,n)%a_len = min_alen
!!$          else if( g%MS(i,n)%a_len > max_alen ) then
!!$!!c             write(*,*) " check_length > maximum alen was used in ", i
!!$             g%MS(i,n)%a_len = max_alen
!!$          end if
!!$
!!$          g%MS(i,n)%c_len = g%MS(i,n)%a_len*phi
!!$       else if( g%MS(i,n)%con == 0.0_ps .and. g%MS(i,n)%mass(1) == 0.0_ps) then
!!$          g%MS(i,n)%a_len = 0.0_ps
!!$          g%MS(i,n)%c_len = 0.0_ps
!!$       else
!!$          write(*,*) "something wrong in initialiation!"
!!$          stop
!!$       end if
!!$    end do
!!$  end subroutine check_length

  function get_totalmass( g,n ) result(m_t)
    ! mass (liquid or solid) group
    type (group), intent(in)   :: g
    real(ps)   :: m_t
    integer,intent(in)  :: n
    integer  :: i

    m_t = 0.0_ps
    do i = 1, g%n_bin
       m_t = m_t + g%MS(i,n)%mass(1)
    end do
  end function get_totalmass


  subroutine print_dcondt( g, cur_time, ofname,id,jd,kd,istrt,ifunit)
    type (group), intent(in) :: g
    real(PS), intent(in)                      :: cur_time
    character (len = *)      :: ofname
    real(ps)                 :: total_tend
    integer,intent(in) :: istrt,ifunit
    integer :: id(g%l),jd(g%l),kd(g%l)
    integer                  :: i, j,n

200 format( f10.1, i4, 20( 2x, es17.8e3))
    if(istrt==1) then
       open( unit=ifunit, file=ofname,status='replace')
    else
       open( unit=ifunit, file=ofname,status='unknown',position='append')
    end if
!!c    write(ifunit,'(100a)') "# 1 time, 2 bin no, 3 total, 4 (1), 5 (2), 6 (3), 7 (4), &
!!c         8 (5), 9 (6), 10 (7), 11 (8), 12 (9) 13 (10) 14 (11) 15 (12)"

250 format('# ',3i6)
!!c    if(isect==1) then
!!c       write(ifunit,'(5i8)') g%n_bin,15,nmic
!!c    end if
    do n=1,g%l
       ! loop over grid points
!!c       if(g%mark_cm(n)==3) cycle
       write(ifunit,250) id(n),jd(n),kd(n)
       do i = 1, g%n_bin
          total_tend = 0.0_ps
          do j = 1, g%n_tendpros
             total_tend = total_tend + g%MS(i,n)%dcondt(j)
          end do
          write(ifunit,200) cur_time, i, total_tend, ( g%MS(i,n)%dcondt(j), j=1,g%n_tendpros)
       end do
!!c       write(ifunit,*) " "
    end do
    close(ifunit)
  end subroutine print_dcondt

  subroutine print_bdcondt( g, ofname,id,jd,kd,istrt,ifunit)
    type (group), intent(in) :: g
    !real(PS), intent(in)                      :: cur_time
    character (len = *)      :: ofname
    !real(ps)                 :: total_tend
    integer,intent(in) :: istrt,ifunit
    integer :: id(g%l),jd(g%l),kd(g%l)
    integer                  :: i, j,n

    if(istrt==1) then
       open( unit=ifunit, file=ofname,status='replace',form='unformatted')
    else
       open( unit=ifunit, file=ofname,status='unknown',form='unformatted',position='append')
    end if

!!c    # 1 id, 2 jd, 3 kd, &
!!c      4 (1), 5 (2), 6 (3), 7 (4), &
!!c      8 (5), 9 (6), 10 (7), 11 (8), 12 (9) 13 (10) 14 (11) 15 (12)"
!!c
    write(ifunit) (real(id(n),PS_KIND),real(jd(n),PS_KIND),real(kd(n),PS_KIND),&
          ((real(g%ms(i,n)%dcondt(j),PS_KIND), j=1,g%n_tendpros),i=1,g%N_BIN),n=1,g%L)
    close(ifunit)
  end subroutine print_bdcondt


  subroutine print_dmassdt( g, cur_time, ofname,id,jd,kd,istrt,ifunit)
    type (group), intent(in) :: g
    real(PS), intent(in)                      :: cur_time
    character (len = *)      :: ofname
    real(ps)                 :: total_tend
    integer,intent(in) :: istrt,ifunit
    integer :: id(g%l),jd(g%l),kd(g%l)
    integer                  :: i, j,n

200 format( f10.1, i4, 20( 2x, es17.8e3))
250 format('# ',3i6)

    if(istrt==1) then
       open( unit=ifunit, file=ofname,status='replace')
    else
       open( unit=ifunit, file=ofname,status='unknown',position='append')
    end if
!!c    write(ifunit,'(100a)') "# 1 time, 2 bin no, 3 total, 4 (1), 5 (2), 6 (3), 7 (4), &
!!c         8 (5), 9 (6), 10 (7), 11 (8), 12 (9) 13 (10) 14 (11) 15 (12)"
!!c    if(isect==1) then
!!c       write(ifunit,'(5i8)') g%n_bin,15,nmic
!!c    end if
    do n=1,g%l
       ! loop over grid points
!!c       if(g%mark_cm(n)==3) cycle
!!c       write(ifunit,*) "# grid number",n
       write(ifunit,250) id(n),jd(n),kd(n)
       do i = 1, g%n_bin
          total_tend = 0.0_ps
          do j = 1, g%n_tendpros
             total_tend = total_tend + g%MS(i,n)%dmassdt(1,j)
          end do
          write(ifunit,200) cur_time, i, total_tend, ( g%MS(i,n)%dmassdt(1,j), j=1,g%n_tendpros)
       end do
!!c       write(ifunit,*) " "
    end do
    close(ifunit)
  end subroutine print_dmassdt

  subroutine print_bdmassdt( g, ofname,id,jd,kd,istrt,ifunit)
    type (group), intent(in) :: g
    !real(PS), intent(in)                      :: cur_time
    character (len = *)      :: ofname
    !real(ps)                 :: total_tend
    integer,intent(in) :: istrt,ifunit
    integer :: id(g%l),jd(g%l),kd(g%l)
    integer                  :: i, j,n

    if(istrt==1) then
       open( unit=ifunit, file=ofname,status='replace',form='unformatted')
    else
       open( unit=ifunit, file=ofname,status='unknown',form='unformatted',position='append')
    end if
!!c    # 1 id, 2 jd, 3 kd, &
!!c      4 (1), 5 (2), 6 (3), 7 (4), &
!!c      8 (5), 9 (6), 10 (7), 11 (8), 12 (9) 13 (10) 14 (11) 15 (12)"
!!c
    write(ifunit) (real(id(n),PS_KIND),real(jd(n),PS_KIND),real(kd(n),PS_KIND),&
          ((real(g%ms(i,n)%dmassdt(1,j),PS_KIND), j=1,g%n_tendpros),i=1,g%N_BIN),n=1,g%L)
    close(ifunit)
  end subroutine print_bdmassdt

!!$  subroutine print_dvoldt( g, cur_time, ofname,id,jd,kd,istrt,ifunit)
!!$    type (group), intent(in) :: g
!!$    real(PS), intent(in)                      :: cur_time
!!$    character (len = *)      :: ofname
!!$    real(ps)                 :: total_tend_1,total_tend_2,total_tend_3
!!$    integer,intent(in) :: istrt,ifunit
!!$    integer :: id(g%l),jd(g%l),kd(g%l)
!!$    integer                  :: i, j,n
!!$
!!$200 format( f10.1, i4, 50( 2x, es17.8e3))
!!$250 format('# ',3i6)
!!$    if(istrt==1) then
!!$       open( unit=ifunit, file=ofname,status='replace')
!!$    else
!!$       open( unit=ifunit, file=ofname,status='unknown',position='append')
!!$    end if
!!$!!c    write(ifunit,'(100a)') "# 1 time, 2 bin no, 3 total_1, 4 (1), 5 (2), 6 (3), 7 (4), &
!!$!!c         8 (5), 9 (6), 10 (7), 11 (8), 12 (9), 13 (10), 14 (11), 15 (12)&
!!$!!c         16 total_2, 17 (1), 18 (2), 19 (3), 20 (4), &
!!$!!c         21 (5), 22 (6), 23 (7), 24 (8), 25 (9), 26 (10), 27 (11), 28 (12), &
!!$!!c         29 total_3, 30 (1), 31 (2), 32 (3), 33 (4), &
!!$!!c         34 (5), 35 (6), 36 (7), 37 (8), 38 (9), 39 (10), 40 (11), 41 (12)"
!!$
!!$!!c    if(isect==1) then
!!$!!c       write(ifunit,'(5i8)') g%n_bin,41,nmic
!!$!!c    end if
!!$    do n=1,g%l
!!$       ! +++ loop over grids +++
!!$!!c       if(g%mark_cm(n)==3) cycle
!!$!!c       write(ifunit,*) "# grid number",n
!!$       write(ifunit,250) id(n),jd(n),kd(n)
!!$       do i = 1, g%n_bin
!!$          total_tend_1 = 0.0_ps
!!$          total_tend_2 = 0.0_ps
!!$          total_tend_3 = 0.0_ps
!!$          do j = 1, g%n_tendpros
!!$             total_tend_1 = total_tend_1 + g%MS(i,n)%dvoldt(ivcs,j)
!!$             total_tend_2 = total_tend_2 + g%MS(i,n)%dvoldt(iacr,j)
!!$             total_tend_3 = total_tend_3 + g%MS(i,n)%dvoldt(iccr,j)
!!$          end do
!!$          write(ifunit,200) cur_time, i, &
!!$               total_tend_1, ( g%MS(i,n)%dvoldt(ivcs,j), j=1,g%n_tendpros),&
!!$               total_tend_2, ( g%MS(i,n)%dvoldt(iacr,j), j=1,g%n_tendpros),&
!!$               total_tend_3, ( g%MS(i,n)%dvoldt(iccr,j), j=1,g%n_tendpros)
!!$       end do
!!$!!c       write(ifunit,*) " "
!!$    end do
!!$    close(ifunit)
!!$
!!$  end subroutine print_dvoldt

!!$  subroutine print_rmatlab( g,ag,cur_time,ofname,id,jd,kd,istrt,ifunit)
!!$    type (group), intent(in) :: g
!!$    ! thermo variable object
!!$    type (airgroup) :: ag
!!$    real(PS), intent(in)                      :: cur_time
!!$    character (len = *)      :: ofname
!!$    integer :: id(g%l),jd(g%l),kd(g%l)
!!$    integer,intent(in) :: istrt,ifunit
!!$    integer                  :: i,n
!!$
!!$300 format( f10.1, 50( 2x, es17.8e3))
!!$250 format('# ',3i6)
!!$    if(istrt==1) then
!!$       open( unit=ifunit, file=ofname,status='replace')
!!$    else
!!$       open( unit=ifunit, file=ofname,status='unknown',position='append')
!!$    end if
!!$!!c    open( unit=2, file=ofname2,form='unformatted',position='append')
!!$    if( g%token == 2 ) then
!!$!!c       if(isect==1) then
!!$!!c          write(ifunit,'(5i8)') g%n_bin,44,nmic
!!$!!c          write(2) g%n_bin,28,nmic
!!$!!c       end if
!!$       do n=1,g%l
!!$          ! loop over grid points
!!$!!c          if(g%mark_cm(n)==3) cycle
!!$
!!$          write(ifunit,250) id(n),jd(n),kd(n)
!!$!!c          write(2) id(n),jd(n),kd(n)
!!$!!c       write(ifunit,'(100a)') "% 1 time, &
!!$!!c            2 con, 3 mass_t, 4 mass_r, 5 mass_a, 6 mass_ice,&
!!$!!c            7 mass_apt, 8 mass_aps, 9 mass_api,
!!$!!c            10 len, 11 den, 12 vtm, 13 tmp, 14 fv, 15 fac, 16 fh, &
!!$!!c            17 n_re, 18 cap, 19 binb,&
!!$!!c            20 a len, 21 c len, 22 d len, 23 r len, 24 e len,&
!!$!!c            25 v_cs, 26 semi_aip, 27 semi_cip, 28 v_r, 29 v_a, &
!!$!!c            30 t, 31 sw, 32 si, 33 dsw, 34 dsi, 35 mass_melt
!!$!!c            36 ag, 37 cg, 38 n_exice, 39 v_csw, 40 semi_a, 41 semi_c, 42 den_ip, 43 den_ic
!!$!!c            44 fkn, 45 mass_frez, 46 porocity
!!$
!!$          do i = 1, g%n_bin
!!$             write(ifunit,300) cur_time, &
!!$                  g%MS(i,n)%con,g%MS(i,n)%mass(imt),g%MS(i,n)%mass(imr),g%MS(i,n)%mass(ima),&
!!$                  g%MS(i,n)%mass(imc),g%MS(i,n)%mass(imat),g%MS(i,n)%mass(imas),g%MS(i,n)%mass(imai),&
!!$                  g%MS(i,n)%len,g%MS(i,n)%den,g%MS(i,n)%vtm,g%MS(i,n)%tmp,&
!!$                  g%MS(i,n)%fv,g%MS(i,n)%fac,g%MS(i,n)%fh,g%MS(i,n)%nre,&
!!$                  g%MS(i,n)%cap,g%binb(i),&
!!$                  g%MS(i,n)%a_len,g%MS(i,n)%c_len,g%IS(i,n)%d,g%IS(i,n)%r,g%IS(i,n)%e,&
!!$                  g%IS(i,n)%v_cs,g%IS(i,n)%semi_aip,g%IS(i,n)%semi_cip,&
!!$                  g%MS(i,n)%vol(1),g%MS(i,n)%vol(2),&
!!$                  ag%tv(n)%t,ag%tv(n)%s_v(1),ag%tv(n)%s_v(2),ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2),&
!!$                  g%MS(i,n)%mass(imw),&
!!$                  g%IS(i,n)%ag,g%IS(i,n)%cg,g%IS(i,n)%n_exice,&
!!$                  g%IS(i,n)%v_csw,g%MS(i,n)%semi_a,g%MS(i,n)%semi_c,g%IS(i,n)%den_ip,g%IS(i,n)%den_ic,&
!!$                  g%MS(i,n)%fkn,g%MS(i,n)%mass(imf),g%IS(i,n)%q_e
!!$
!!$          end do
!!$       end do
!!$    else if( g%token == 1 ) then
!!$!!c       if(isect==1) then
!!$!!c          write(ifunit,'(4i8)') g%n_bin,23,nmic
!!$!!c!!c          write(2) g%n_bin,24,nmic
!!$!!c       end if
!!$       do n=1,g%l
!!$          ! loop over grid points
!!$!!c          if(g%mark_cm(n)==3) cycle
!!$          write(ifunit,250) id(n),jd(n),kd(n)
!!$!!c          write(2) id(n),jd(n),kd(n)
!!$!!c       write(ifunit,'(100a)') "% 1 time, &
!!$!!c            2 con, 3 mass_t, 4 mass_apt, 5 mass_aps, 6 mass_api,
!!$!!c            7 len, 8 den, 9 vtm, 10 tmp, 11 fv, 12 fac, 13 fh, &
!!$!!c            14 n_re, 15 cap, 16 binb,&
!!$!!c            17 a len, 18 c len, 19 t, 20 sw, 21 si, 22 dsw, 23 dsi
!!$!!c            24 fkn
!!$
!!$          do i = 1, g%n_bin
!!$             write(ifunit,300) cur_time, &
!!$                  g%MS(i,n)%con,g%MS(i,n)%mass(1),g%MS(i,n)%mass(2),g%MS(i,n)%mass(3),&
!!$                  g%MS(i,n)%mass(4),&
!!$                  g%MS(i,n)%len,g%MS(i,n)%den,g%MS(i,n)%vtm,g%MS(i,n)%tmp,&
!!$                  g%MS(i,n)%fv,g%MS(i,n)%fac,g%MS(i,n)%fh,g%MS(i,n)%nre,&
!!$                  g%MS(i,n)%cap,g%binb(i),&
!!$                  g%MS(i,n)%a_len,g%MS(i,n)%c_len,&
!!$                  ag%tv(n)%t,ag%tv(n)%s_v(1),ag%tv(n)%s_v(2),ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2),&
!!$                  g%MS(i,n)%fkn
!!$          end do
!!$       end do
!!$    else if( g%token == 3 ) then
!!$!!c       if(isect==1) then
!!$!!c          write(ifunit,'(4i8)') g%n_bin,10,nmic
!!$!!c!!c          write(2) g%n_bin,21,nmic
!!$!!c       end if
!!$       do n=1,g%l
!!$          ! loop over grid points
!!$!!c          if(g%mark_cm(n)==3) cycle
!!$          write(ifunit,250) id(n),jd(n),kd(n)
!!$!!c          write(2) id(n),jd(n),kd(n)
!!$
!!$!!c       write(ifunit,'(100a)') "% 1 time, &
!!$!!c            2 con, 3 mass_t, 4 mass_aps, 5 mass_api,
!!$!!c            6 len, 7 den, 8 t, 9 sw, 10 si, 11 dsw, 12 dsi
!!$          do i = 1, g%n_bin
!!$             write(ifunit,300) cur_time, &
!!$                  g%MS(i,n)%con,g%MS(i,n)%mass(1),g%MS(i,n)%mass(2),g%MS(i,n)%mass(3),&
!!$                  g%MS(i,n)%len,g%MS(i,n)%den,&
!!$                  ag%tv(n)%t,ag%tv(n)%s_v(1),ag%tv(n)%s_v(2),ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2)
!!$          end do
!!$       end do
!!$    end if
!!$    close(ifunit)
!!$!!c    close(2)
!!$  end subroutine print_rmatlab

!!$  subroutine print_brmatlab( g,ag,cur_time,ofname,id,jd,kd,istrt,ifunit)
!!$    type (group), intent(in) :: g
!!$    ! thermo variable object
!!$    type (airgroup) :: ag
!!$    real(PS), intent(in)                      :: cur_time
!!$    character (len = *)      :: ofname
!!$    integer :: id(g%l),jd(g%l),kd(g%l)
!!$    integer,intent(in) :: istrt,ifunit
!!$    integer                  :: i,n
!!$
!!$    if(istrt==1) then
!!$       open( unit=ifunit, file=ofname,status='replace',form='unformatted')
!!$    else
!!$       open( unit=ifunit, file=ofname,status='unknown',form='unformatted',position='append')
!!$    end if
!!$    if( g%token == 2 ) then
!!$!!c            1 id, 2 jd, 3 kd, 4 time, 5 t, &
!!$!!c            6 sw, 7 si, 8 dsw, 9 dsi, &
!!$!!c                            +
!!$!!c            1 con, 2 mass_t, 3 mass_r, 4 mass_a, 5 mass_ice,&
!!$!!c            6 mass_apt, 7 mass_aps, 8 mass_api, 9 len, 10 den, &
!!$!!c           11 vtm, 12 tmp, 13 fv, 14 fac, 15 fh, &
!!$!!c           16 n_re, 17 cap, l8 binb, 19 a len, 20 c len, &
!!$!!c           21 d len, 22 r len, 23 e len, 24 v_cs, 25 semi_aip,&
!!$!!c           26 semi_cip, 27 v_r, 28 v_a, 29 mass_melt, 30 ag, &
!!$!!c           31 cg, 32 n_exice, 33 v_csw, 34 semi_a, 35 semi_c, &
!!$!!c           36 den_ip, 37 den_ic, 38 fkn, 39 mass_frez, 40 porocity
!!$
!!$          write(ifunit) (real(id(n),PS_KIND),real(jd(n),PS_KIND),real(kd(n),PS_KIND),cur_time,&
!!$                  ag%tv(n)%t,ag%tv(n)%s_v(1),ag%tv(n)%s_v(2),ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2),&
!!$                  !
!!$                 (g%ms(i,n)%con,g%ms(i,n)%mass(imt),g%ms(i,n)%mass(imr),g%ms(i,n)%mass(ima),&
!!$                  g%ms(i,n)%mass(imc),g%ms(i,n)%mass(imat),g%ms(i,n)%mass(imas),g%ms(i,n)%mass(imai),&
!!$                  g%ms(i,n)%len,g%ms(i,n)%den,g%ms(i,n)%vtm,g%ms(i,n)%tmp,&
!!$                  g%ms(i,n)%fv,g%ms(i,n)%fac,g%ms(i,n)%fh,g%ms(i,n)%nre,&
!!$                  g%ms(i,n)%cap,g%binb(i),&
!!$                  g%ms(i,n)%a_len,g%ms(i,n)%c_len,g%is(i,n)%d,g%is(i,n)%r,g%is(i,n)%e,&
!!$                  g%is(i,n)%v_cs,g%IS(i,n)%semi_aip,g%IS(i,n)%semi_cip,&
!!$                  g%ms(i,n)%vol(1),g%ms(i,n)%vol(2),&
!!$                  g%MS(i,n)%mass(imw),&
!!$                  g%IS(i,n)%ag,g%IS(i,n)%cg,g%IS(i,n)%n_exice,&
!!$                  g%IS(i,n)%v_csw,g%MS(i,n)%semi_a,g%MS(i,n)%semi_c,g%IS(i,n)%den_ip,g%IS(i,n)%den_ic,&
!!$                  g%MS(i,n)%fkn,g%MS(i,n)%mass(imf),g%IS(i,n)%q_e,i=1,g%N_BIN),n=1,g%L)
!!$
!!$
!!$    else if( g%token == 1 ) then
!!$!!c            1 id, 2 jd, 3 kd, 4 time, 5 t,&
!!$!!c            6 sw, 7 si, 8 dsw, 9 dsi,&
!!$!!c                           +
!!$!!c            1 con, 2 mass_t, 3 mass_apt, 4 mass_aps, 5 mass_api, &
!!$!!c            6 len, 7 den, 8 vtm, 9 tmp, 10 fv, &
!!$!!c           11 fac, 12 fh, 13 n_re, 14 cap, 15 binb,&
!!$!!c           16 a len, 17 c len, 18 fkn
!!$
!!$          write(ifunit) (real(id(n),PS_KIND),real(jd(n),PS_KIND),real(kd(n),PS_KIND),cur_time,&
!!$                  ag%tv(n)%t,ag%tv(n)%s_v(1),ag%tv(n)%s_v(2),ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2),&
!!$                  !
!!$                  (g%ms(i,n)%con,g%ms(i,n)%mass(rmt),g%ms(i,n)%mass(rmat),g%ms(i,n)%mass(rmas),&
!!$                  g%ms(i,n)%mass(rmai),&
!!$                  g%ms(i,n)%len,g%ms(i,n)%den,g%ms(i,n)%vtm,g%ms(i,n)%tmp,&
!!$                  g%ms(i,n)%fv,g%ms(i,n)%fac,g%ms(i,n)%fh,g%ms(i,n)%nre,&
!!$                  g%ms(i,n)%cap,g%binb(i),&
!!$                  g%ms(i,n)%a_len,g%ms(i,n)%c_len,&
!!$                  g%MS(i,n)%fkn,i=1,g%N_BIN),n=1,g%L)
!!$    else if( g%token == 3 ) then
!!c            1 id, 2 jd, 3 kd, 4 time, 5 time,&
!!c            6 sw, 7 si, 8 dsw, 9 dsi,&
!!c                           +
!!c            1 con, 2 mass_t, 3 mass_aps, 4 mass_api,
!!c            5 len, 6 den, 7 mean geo, 8 sd geo

!!$          write(ifunit) (real(id(n),PS_KIND),real(jd(n),PS_KIND),real(kd(n),PS_KIND),cur_time,&
!!$                  ag%tv(n)%t,ag%tv(n)%s_v(1),ag%tv(n)%s_v(2),ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2),&
!!$                  !
!!$                  (g%ms(i,n)%con,g%ms(i,n)%mass(amt),g%ms(i,n)%mass(ams),g%ms(i,n)%mass(ami),&
!!$                  g%ms(i,n)%len,g%ms(i,n)%den,&
!!$                  real(g%MS(i,n)%p(1),PS_KIND),real(exp(min(10.0_DS,g%MS(i,n)%p(2))),PS_KIND),&
!!$                  i=1,g%N_BIN),n=1,g%L)
!!$    end if
!!$    close(ifunit)
!!$!!c    close(2)
!!$  end subroutine print_brmatlab

!!$  subroutine print_imatlab( g,ofname,id,jd,kd,istrt,ifunit)
!!$    type (group), intent(in) :: g
!!$    !real(PS), intent(in)                      :: cur_time
!!$    character (len = *)      :: ofname
!!$    integer :: id(g%l),jd(g%l),kd(g%l)
!!$    integer,intent(in) :: istrt,ifunit
!!$    integer                  :: i,n
!!$
!!$300 format(8i3)
!!$250 format('# ',3i6)
!!$    if(istrt==1) then
!!$       open( unit=ifunit, file=ofname,status='replace')
!!$    else
!!$       open( unit=ifunit, file=ofname,status='unknown',position='append')
!!$    end if
!!$!!c    open( unit=2, file=ofname2,form='unformatted',position='append')
!!$    if( g%token == 2 ) then
!!$!!c       if(isect==1) then
!!$!!c          write(ifunit,'(5i8)') g%n_bin,6,nmic
!!$!!c!!c          write(2) g%n_bin,6,nmic
!!$!!c       end if
!!$       do n=1,g%l
!!$          ! loop over grid points
!!$!!c          if(g%mark_cm(n)==3) cycle
!!$          write(ifunit,250) id(n),jd(n),kd(n)
!!$!!c          write(2) id(n),jd(n),kd(n)
!!$          do i = 1, g%n_bin
!!$             write(ifunit,300) g%IS(i,n)%habit,g%IS(i,n)%sh_type,g%IS(i,n)%init_growth,g%IS(i,n)%growth_mode,g%IS(i,n)%is_mod
!!$!!c             write(2) g%IS(i,n)%habit,g%IS(i,n)%sh_type
!!$          end do
!!$       end do
!!$    else
!!$    end if
!!$    close(ifunit)
!!$!!c    close(2)
!!$  end subroutine print_imatlab

!!$  subroutine print_bimatlab( g,ofname,id,jd,kd,istrt,ifunit)
!!$    type (group), intent(in) :: g
!!$    !real(PS), intent(in)                      :: cur_time
!!$    character (len = *)      :: ofname
!!$    integer :: id(g%l),jd(g%l),kd(g%l)
!!$    integer,intent(in) :: istrt,ifunit
!!$    integer                  :: i,n
!!$
!!$    if(istrt==1) then
!!$       open( unit=ifunit, file=ofname,status='replace',form='unformatted')
!!$    else
!!$       open( unit=ifunit, file=ofname,status='unknown',form='unformatted',position='append')
!!$    end if
!!$    if( g%token == 2 ) then
!!$        write(ifunit) (id(n),jd(n),kd(n),&
!!$           (g%is(i,n)%habit,g%is(i,n)%sh_type,g%IS(i,n)%init_growth,g%IS(i,n)%growth_mode,g%IS(i,n)%is_mod,&
!!$            i=1,g%N_BIN),n=1,g%L)
!!$    else
!!$    end if
!!$    close(ifunit)
!!$  end subroutine print_bimatlab

  subroutine print_tendency_ap( g, cur_time,output_format,id,jd,kd,iproc,istrt,nc)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! print out the tendency of each process
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! number of categores for the group
    integer, intent(in) :: nc
    type (group), dimension(nc)       :: g
    real(PS), intent(in)     :: cur_time
    character(len=16),intent(in) :: output_format
    ! output file names
    character(len=20)  :: string_prefix1, string_prefix2,string_prefix3
    character(len=20)  :: string_postfix
    character(len=6)   :: string_time
    character(len=46)  :: ofname,ofname2
    !character(len=5)  :: ccat
    real(PS)               :: digit
    integer :: id(g(nc)%l),jd(g(nc)%l),kd(g(nc)%l)
    integer,intent(in) :: iproc,istrt

    integer            :: inum, jnum,ic,i,ifunit
    character(len=2) :: cproc

    ! +++ make names of output files +++
    !!        write(string_time, *) int(cur_time/60.0_ps)

    ! add processor number
    write(cproc,'(I2)') iproc
    do i=1,2
      if(cproc(i:i).eq.' ') cproc(i:i)='0'
    enddo
    ifunit=801+iproc

    do ic=1,nc
       if( g(nc)%token == 2 ) then
          string_prefix1 = "dcdt_s"
          string_prefix2 = "dmdt_s"
          string_prefix3 = "dvdt_s"
       else if( g(nc)%token == 11 ) then
          string_prefix1 = "dcdt_c"
          string_prefix2 = "dmdt_c"
       else if( g(nc)%token == 1 .or.  g(nc)%token == 12) then
          string_prefix1 = "dcdt_r"
          string_prefix2 = "dmdt_r"
       else if( g(nc)%token == 3) then
          string_prefix1 = "dcdt_a"
          string_prefix2 = "dmdt_a"
       end if
       string_prefix1=trim(string_prefix1)//achar(ic+48)//'_'
       string_prefix2=trim(string_prefix2)//achar(ic+48)//'_'
       string_prefix3=trim(string_prefix3)//achar(ic+48)//'_'

       string_time = "000000"
!!c    string_postfix="m.dat"          ! in minites
       string_postfix="s.dat"          ! in seconds


!!c       inum = int(cur_time/60.0_ps)
       inum = floor(cur_time+1.0e-5)
       digit = 1.0
       do
          inum = inum/10
          if( inum /= 0 ) then; digit = digit + 1.0
          else if( inum == 0 ) then; exit
          end if
       end do
   !!c    inum = int(cur_time/60.0_ps)
       inum = floor(cur_time+1.0e-5)
       do
          jnum = inum/int(10.0**(floor(digit+1.0e-5)-1))
          inum = inum - jnum*int(10.0**(floor(digit+1.0e-5)-1))
          string_time(6-int(digit)+1:6-int(digit)+1) = achar(jnum+48)
          digit = digit - 1.0
          if( digit <= 0.0 ) exit
       end do

       ofname = trim(string_prefix1) // trim(string_time)
       ofname = trim(ofname) // trim(string_postfix)//cproc

       select case(output_format)
       case('binary')
         ofname2='b'//ofname
         call print_bdcondt( g(ic), ofname2,id,jd,kd,istrt,ifunit)
       case default
         call print_dcondt( g(ic), cur_time, ofname,id,jd,kd,istrt,ifunit)
       end select

       ofname = trim(string_prefix2) // trim(string_time)
       ofname = trim(ofname) // trim(string_postfix)//cproc
       select case(output_format)
       case('binary')
         ofname2='b'//ofname
         call print_bdmassdt( g(ic), ofname2,id,jd,kd,istrt,ifunit)
       case default
         call print_dmassdt( g(ic), cur_time, ofname,id,jd,kd,istrt,ifunit)
       end select

!!c       if( g(ic)%token == 2 ) then
!!c          ofname = trim(string_prefix3) // trim(string_time)
!!c          ofname = trim(ofname) // trim(string_postfix)//cproc
!!c          call print_dvoldt( g(ic), cur_time, ofname,id,jd,kd,istrt,ifunit)
!!c       end if
    end do
  end subroutine print_tendency_ap

  subroutine print_tendency(g,cur_time,output_format,id,jd,kd,iproc,istrt)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! print out the tendency of each process
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (group), intent(in)       :: g
    real(PS), intent(in)     :: cur_time
    character(len=16),intent(in) :: output_format
    ! output file names
    character(len=20)  :: string_prefix1, string_prefix2,string_prefix3
    character(len=20)  :: string_postfix
    character(len=6)   :: string_time
    character(len=46)  :: ofname,ofname2
    !character(len=1)  :: ccat
    real(PS)               :: digit
    integer,intent(in) :: iproc,istrt
    integer :: id(g%l),jd(g%l),kd(g%l)

    integer            :: inum, jnum,i,ifunit!,ic
    character(len=2) :: cproc

    ! +++ make names of output files +++
    !!        write(string_time, *) int(cur_time/60.0_ps)

    if( g%token == 2 ) then
       string_prefix1 = "dcdt_s_"
       string_prefix2 = "dmdt_s_"
       string_prefix3 = "dvdt_s_"
    else if( g%token == 11 ) then
       string_prefix1 = "dcdt_c_"
       string_prefix2 = "dmdt_c_"
    else if( g%token == 1 .or.  g%token == 12) then
       string_prefix1 = "dcdt_r_"
       string_prefix2 = "dmdt_r_"
    end if

    string_time = "000000"
!!c    string_postfix="m.dat"          ! in minites
    string_postfix="s.dat"          ! in seconds


!!c     inum = int(cur_time/60.0_ps)
    inum = floor(cur_time+1.0e-5)
    digit = 1.0
    do
       inum = inum/10
       if( inum /= 0 ) then; digit = digit + 1.0
       else if( inum == 0 ) then; exit
       end if
    end do
   !!c    inum = int(cur_time/60.0_ps)
    inum = floor(cur_time+1.0e-5)
    do
       jnum = inum/int(10.0**(floor(digit+1.0e-5)-1))
       inum = inum - jnum*int(10.0**(floor(digit+1.0e-5)-1))
       string_time(6-int(digit)+1:6-int(digit)+1) = achar(jnum+48)
       digit = digit - 1.0
       if( digit <= 0.0 ) exit
    end do
    ! add processor number
    write(cproc,'(I2)') iproc
    do i=1,2
      if(cproc(i:i).eq.' ') cproc(i:i)='0'
    enddo
    ifunit=801+iproc

    ofname = trim(string_prefix1) // trim(string_time)
    ofname = trim(ofname) // trim(string_postfix)//cproc

    select case(output_format)
    case('binary')
      ofname2 = 'b'//ofname
      call print_bdcondt( g, ofname2,id,jd,kd,istrt,ifunit)
    case default
      call print_dcondt( g, cur_time, ofname,id,jd,kd,istrt,ifunit)
    end select

    ofname = trim(string_prefix2) // trim(string_time)
    ofname = trim(ofname) // trim(string_postfix)//cproc
    select case(output_format)
    case('binary')
      ofname2 = 'b'//ofname
      call print_bdmassdt( g, ofname2,id,jd,kd,istrt,ifunit)
    case default
      call print_dmassdt( g, cur_time, ofname,id,jd,kd,istrt,ifunit)
    end select

!!c    if( g%token == 2 ) then
!!c       ofname = trim(string_prefix3) // trim(string_time)
!!c       ofname = trim(ofname) // trim(string_postfix)//cproc
!!c       call print_dvoldt( g, cur_time, ofname,id,jd,kd,istrt,ifunit)
!!c    end if
  end subroutine print_tendency


!!$  function get_vapdepact(th_var,phase,alpha,A,B,Sd_max,dt,min_dt,r0,m0) result(m1)
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! use Twomey lower bound estimate for integral of supersaturation.
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! radius of aerosol
!!$    type (Thermo_Var), intent(in) :: th_var
!!$    integer,intent(in) :: phase
!!$    real(PS),intent(in) :: r0,dt,min_dt
!!$    real(PS),intent(in) :: alpha,A,B
!!$    real(PS) :: fkn
!!$    real(PS) :: m0,m1
!!$    real(PS) :: S_c2,Twomey,Sd_max
!!$
!!$    ! +++ calculate kinetic effect +++
!!$    fkn=get_fkn(th_var,phase,r0)
!!$
!!$
!!$    if(phase==1) then
!!$       ! +++ calculate critical saturation +++
!!$       S_c2=(2.0_PS/sqrt(B))*(A/(3.0_PS*r0))**3.0
!!$
!!$       Twomey=max(0.0_PS,(Sd_max**2.0-S_c2))/2.0/alpha/th_var%W+Sd_max*max(0.0_PS,dt-min_dt)
!!$
!!$       m1=max(m0,&
!!$            (coef4p*coef3i4p1i3*0.66666666666667_PS*th_var%GTP(1)*fkn*Twomey+m0**(2.0/3.0))**1.5)
!!$    elseif(phase==2) then
!!$       ! +++ calculate critical saturation +++
!!$       S_c2=((th_var%svw0+1.0_PS)*th_var%e_sat(1)/th_var%e_sat(2)-1.0_PS)**2.0
!!$
!!$       Twomey=max(0.0_PS,(Sd_max**2.0-S_c2))/2.0/alpha/th_var%W+Sd_max*max(0.0_PS,dt-min_dt)
!!$
!!$       m1=max(m0,&
!!$            (coef4p*coef3i4p1i3*0.66666666666667_PS*th_var%GTP(2)*fkn*Twomey+m0**(2.0/3.0))**1.5)
!!$    end if
!!$  end function get_vapdepact

!!$  function get_vapdepact_v2(th_var,phase,Sd_max,dt,r0,m0) result(m1)
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! use twomey lower bound estimate for integral of supersaturation.
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! radius of aerosol
!!$    type (thermo_var), intent(in) :: th_var
!!$    integer,intent(in) :: phase
!!$    real(ps),intent(in) :: r0
!!$    real(PS),intent(in) :: dt
!!$    real(ps) :: fkn
!!$    real(ps) :: m0,m1
!!$    real(ps) :: twomey,sd_max!,s_c2
!!$
!!$    ! +++ calculate kinetic effect +++
!!$    fkn=get_fkn(th_var,phase,r0)
!!$
!!$    ! +++ calculate critical saturation +++
!!$    twomey=Sd_max*dt
!!$    m1=max(m0,&
!!$         (coef4p*coef3i4p1i3*0.66666666666667_ps*th_var%gtp(phase)*fkn*twomey+m0**(2.0/3.0))**1.5)
!!$  end function get_vapdepact_v2

!!$  function get_inirad_ccn(W,AA,BB,S,rn_d) result(m1)
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! use Ivanova et al (1977)'s parameterization of wet radius at S=0.0
!!$    ! and see Kogan (1991)
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    real(PS),intent(in) :: W,AA,BB,S,rn_d
!!$    real(PS) :: m1
!!$    ! radius of the dry salt nucleus, ak is the factor
!!$    real(PS) :: rn,ak
!!$    real(PS) :: r_star
!!$    real(PS) :: a0,a2,p,q,dum,w1,w2,w3,w4,req1,req2,req3,req4,req0
!!$
!!$    rn=rn_d*1.0e+4_PS
!!$
!!$    ! +++ calculate threshold radius (in micron) +++
!!$    r_star=9.0e-2_PS*(W*0.01_PS)**(-0.16)
!!$
!!$    if( rn<r_star) then
!!$       ! --- use the equilibrium radius from Kohler curve ---
!!$       a2=-AA/S
!!$       a0=BB/S*rn_d**3.0
!!$       p=-a2*a2/3.0_PS
!!$       q=-a0-a2*a2*a2*2.0_PS/27.0_PS
!!$       dum=0.25_PS*q*q+p*p*p/27.0_PS
!!$       if(dum>0.0) then
!!$          w1=(0.5_PS*q+sqrt(dum))**(1.0/3.0)
!!$          w2=(0.5_PS*q-sqrt(dum))**(1.0/3.0)
!!$          req1=w1-p/3.0_PS/w1-a2/3.0_PS
!!$          req2=w2-p/3.0_PS/w2-a2/3.0_PS
!!$
!!$          w3=-w1
!!$          w4=-w2
!!$          req3=w3-p/3.0_PS/w3-a2/3.0_PS
!!$          req4=w4-p/3.0_PS/w4-a2/3.0_PS
!!$
!!$          req0=min(max(rn_d,req1,req2,req3,req4),5.0e-4_PS)
!!$       else
!!$          ! set it to the minimum bin size of liquid spectrum
!!$          req0=max(1.0e-5_PS,rn_d)
!!$       end if
!!$!!c       write(*,'("1: W,S,rn_d,req0",5ES15.6)') W,S,rn_d,req0
!!$       m1=coef4pi3*(req0)**3.0
!!$    else
!!$       ak=max(1.0_PS,5.8_PS*(W*0.01_PS)**(-0.12)*rn**(-0.214))
!!$       m1=coef4pi3*(ak*rn_d)**3.0
!!$!!c       write(*,'("2: W,S,rn_d,req0,ak",7ES15.6)') W,S,rn_d,ak*rn_d,ak
!!$    end if
!!$  end function get_inirad_ccn

!!$  subroutine cal_vapdepact_ice2(iswitch,coef_vap,Sd_max,dt,m0,m1,Qp)
!!$    ! switch for prognosing mass or assigning mass
!!$    integer, intent(in) :: iswitch
!!$    ! level of complexity
!!$    !integer, intent(in)           :: level
!!$    real(ps),intent(in) :: coef_vap,m0,dt
!!$    real(ps) :: m1
!!$    real(ps) :: sd_max
!!$    !integer :: growth_mode
!!$    ! non-mass variables of a representative particle in the shifted bin
!!$    ! argument 1 : volume of circumscribing sphere
!!$    !          2 : a-axis length
!!$    !          3 : c-axis length
!!$    !          4 : d-axis length
!!$    !          5 : r-axis length
!!$    real(PS), dimension(*)                     :: Qp
!!$
!!$    real(PS) :: d_mean_mass
!!$
!!$
!!$    if(iswitch==1) then
!!$       d_mean_mass=coef_vap*Sd_max*dt
!!$       m1=m0+d_mean_mass
!!$    else
!!$       d_mean_mass=m1-m0
!!$    endif
!!$
!!$    Qp(iacr)=(m1/den_i/coef3s)**(1.0/3.0)
!!$    Qp(iccr)=Qp(iacr)
!!$    Qp(ivcs)=coef4pi3*(Qp(iacr)**2+Qp(iccr)**2)**1.5
!!$
!!$  end subroutine cal_vapdepact_ice2

!!$  subroutine check_tendency_ap(g,ga,loc,KD,ID,JD)
!!$    ! assumptions
!!$    type (group)        :: g
!!$    !
!!$    !   for aerosol groups
!!$    type (group), dimension(*)        :: ga
!!$    !type (airgroup), intent(in) :: ag
!!$    integer,intent(in) :: loc
!!$    integer :: KD(*),ID(*),JD(*)
!!$    integer                 :: n
!!$
!!$    ! --- case of linear distribution ---
!!$    do n=1,ga(1)%l
!!$       if(sum(g%ms(:,n)%dmassdt(rmat,1))<0.0_PS) then
!!$          write(*,'("ck_tend_ap>",4I5,10ES15.6)') loc,KD(n),ID(n),JD(n),&
!!$           sum(g%ms(:,n)%dmassdt(rmat,1)),&
!!$           sum(ga(1)%ms(:,n)%dmassdt(amt,1)),sum(ga(2)%ms(:,n)%dmassdt(amt,1)),sum(ga(3)%ms(:,n)%dmassdt(amt,1))
!!$       end if
!!$    end do
!!$!    do n=1,ga(2)%l
!!$!       if(ga(2)%ms(1,n)%dcondt(1)>0.5_ps) then
!!$!          write(*,*) "loc is",loc
!!$!          write(*,11) n,ga(2)%ms(1,n)%dcondt(1),ag%tv(n)%s_v(1),ag%tv(n)%s_v(2)
!!$!          write(*,12) g%is(1,n)%sh_type,g%is(1,n)%habit,g%ms(1,n)%dcondt(1),g%ms(1,n)%dcondt(2)
!!$11     format("apc:n,dcondt2,sv1,sv2",i5,3es15.6)
!!$!12     format("apc:shtype,habit,dcondt_vp,col",2i5,2es15.6)
!!$!       end if
!!$!    end do
!!$  end subroutine check_tendency_ap

!!$  subroutine mcom_check(g,from)
!!$    ! assumptions
!!$    type (group)        :: g
!!$    character (len=*) :: from
!!$    integer                 :: i,n
!!$
!!$    ! --- case of linear distribution ---
!!$    do n=1,g%l
!!$       do i=1,g%n_bin
!!$          if(g%MS(i,n)%mass(imc)>g%MS(i,n)%mass(1)) then
!!$             write(*,*) "from",from
!!$             write(*,11) i,n,g%MS(i,n)%mass(1),g%MS(i,n)%mass(imc)
!!$             write(*,12) g%MS(i,n)%con,g%MS(i,n)%a_len,g%MS(i,n)%c_len
!!$11     format("mcomcheck:i,n,mass1,4",2i5,3es15.6)
!!$12     format("mcomcheck:con,a,c",3es15.6)
!!$          end if
!!$       end do
!!$    end do
!!$
!!$  end subroutine mcom_check

!!$  subroutine ini_group_midmass(g)
!!$    type (group)        :: g
!!$    integer                 :: i!,n
!!$    do i=1,g%n_bin
!!$       g%ms(i,1)%mass(1)=0.5_ps*(g%binb(i)+g%binb(i+1))
!!$       g%ms(i,1)%con=1.0_ps
!!$    end do
!!$  end subroutine ini_group_midmass

  subroutine ini_group_mp(g)
    implicit none
    type (group)        :: g
    integer                      :: i!, j,k,l,ick
    ! liquid water content (g/m^-3)
    real(ds), parameter              :: w_l = 0.31
    ! intercept (m^-3 mm^-1)
    real(ds)                         :: n_0!, x0
    ! slope (mm^-1)
    real(ds)                         :: lambda
    real(ds)                        :: d_1, d_2, t_1, t_2
    ! scaling diameter, normalization factor for truncated distribution
    real(ds)                         :: d_n!, s

!!c    ! apparent density
!!c    real(ds)                         :: den_a
!!c    ! bulk ice density
!!c    ! axis ratio
!!c    real(ds)                         :: phi,psi
!!c    ! mean mass
!!c    real(ds)                        :: mean_mass
!!c    real(ds)        :: vol

    ! initial mass of cloud
    !real(ds)      :: ini_mcloud

    ! mass of an ice crystal
    !real(ds)      :: m_ice,m_ini
    !real(ds)      :: dum, dum1, dum2

!!c    ! a-axis length and c-axis length
!!c    real(ds) :: a_len,c_len,d_len


    real(PS) :: rr!,m,left_bd,right_bd,rho,dia1,dia2,slp,m_m,r_m

    ! --- case of liquid phase ---
    ! rain rate (mm/hr)
    rr=54.0
    ! lambda, mm^(-1)
!!c     lambda = 3.033386824
    lambda=4.1*rr**(-0.21)
    ! mean diameter (mm)
    d_n = 1.0/lambda
    ! n_0, m^(-3) mm^(-1)
!!c     n_0 = 12697.46258
    n_0=8.0e+3

    do i = 1,g%n_bin
       ! +++ calculate the melt diameteor at bin boundaries +++
       d_1 = (6.0*g%binb(i)/(pi*1.0))**(1.0/3.0)
       d_2 = (6.0*g%binb(i+1)/(pi*1.0))**(1.0/3.0)
       ! +++ change the unit into mm +++
       d_1 = d_1*10.0
       d_2 = d_2*10.0

       t_1 = lambda*d_1
       t_2 = lambda*d_2

       g%ms(i,1)%con = (n_0/lambda) * ( exp(-lambda*d_1) - exp(-lambda*d_2))
       g%ms(i,1)%mass(1) = ((t_1**3.0+3.0*(t_1**2.0)+6.0*t_1+6.0)*exp(-t_1)&
            - (t_2**3.0+3.0*(t_2**2.0)+6.0*t_2+6.0)*exp(-t_2))*&
            (pi/6.0)*0.001*(n_0/(lambda**4.0))

       ! +++ change the unit into cm^-3 from m^-3 +++
       g%ms(i,1)%con=g%ms(i,1)%con/1000000.0
       g%ms(i,1)%mass(1)=g%ms(i,1)%mass(1)/1000000.0
       g%ms(i,1)%mass(rmat)=g%ms(i,rmt)%mass(rmt)*0.1
       g%ms(i,1)%mass(rmas)=g%ms(i,rmt)%mass(rmat)*0.9
       g%ms(i,1)%mass(rmai)=g%ms(i,rmt)%mass(rmat)*0.1

    end do
  end subroutine ini_group_mp

  subroutine gen_length(mass,alen,clen)
    real(PS),intent(in) :: mass!,tmp
    real(PS),intent(inout) :: alen,clen
    real(PS),parameter :: den_max=0.758088258

    ! assume hexagonal column
    alen=(mass/(2.0*pi*den_max))**0.33333333
    clen=alen

!!c    write(*,*) "gen_length > neg length for a or c, modified to",alen,clen

  end subroutine gen_length

!!$  subroutine cal_ts_ndrr(cur_time,g)
!!$    ! **********************************************************************
!!$    ! Calculate the time series of number density and rain rate.
!!$    ! **********************************************************************
!!$    type (Group), intent(inout)   :: g
!!$    real(PS),intent(in) :: cur_time
!!$    real(PS),dimension(3) :: total
!!$    integer :: n,i
!!$
!!$    open( unit=1, file="ts_ndrr.dat",POSITION='APPEND')
!!$    do n=1,g%L
!!$       total = 0.0_PS
!!$       do i = 1, g%N_BIN
!!$          total(1)=total(1)+g%MS(i,n)%con
!!$          total(2)=total(2)+g%MS(i,n)%mass(1)
!!$          total(3)=total(3)+g%MS(i,n)%mass(1)*g%MS(i,n)%vtm
!!$       end do
!!$       ! in m^-3
!!$       total(1)=total(1)*1.0e+6
!!$       ! in g m^-3
!!$       total(2)=total(2)*1.0e+6
!!$       ! in mm/hr
!!$       total(3)=total(3)*36000.0
!!$       write(1,100) n,cur_time,total(1:3)
!!$    end do
!!$100 format(I5,4ES15.6)
!!$    close(1)
!!$  end subroutine cal_ts_ndrr

!!$  subroutine check_rosrel(from,i,n,mass,alen,clen,n_exice,ag,cg,mark)
!!$    integer :: i,n,mark
!!$    real(PS) :: mass,alen,clen,n_exice,ag,cg,rlen_ana,maxd,maxd_ana
!!$    ! factor
!!$    real(PS),parameter :: fct1=0.1,fct2=0.9,N_ros=4.0
!!$    character (len=*) :: from
!!$    mark=0
!!$    if(n_exice>0.5_PS) then
!!$       if(ag<alen*fct1.and.cg>=clen*fct2) then
!!$          rlen_ana=(mass/(den_i*0.3257*N_ros**0.5206))**(1.0/3.0)
!!$          maxd=2.0*sqrt(1.0+0.25**2.0)*2.0*clen
!!$          maxd_ana=2.0*rlen_ana
!!$          if(maxd>2.0*maxd_ana) then
!!$             write(*,'("check_rosrel",A10,2I5,15ES15.6)') from,i,n,maxd/maxd_ana,maxd,maxd_ana,alen,clen,n_exice,ag,cg
!!$             mark=1
!!$          end if
!!$       end if
!!$    end if
!!$
!!$  end subroutine check_rosrel

!!$  subroutine ck_mlttend(gr,gs,ID,JD,KD)
!!$    ! mass (liquid or solid) group
!!$    type (Group), intent(inout)   :: gr,gs
!!$    integer,dimension(*) :: ID,JD,KD
!!$    real(PS) :: sum_dmdt_r,sum_dmdt_s
!!$    integer pro_type,n,i
!!$    pro_type=10
!!$    do n=1,gs%L
!!$       sum_dmdt_r=0.0_PS
!!$       do i=1,gr%N_BIN
!!$          sum_dmdt_r=sum_dmdt_r+gr%MS(i,n)%dmassdt(rmt,pro_type)
!!$       end do
!!$       sum_dmdt_s=0.0_PS
!!$       do i=1,gs%N_BIN
!!$          sum_dmdt_s=sum_dmdt_s+gs%MS(i,n)%dmassdt(imt,pro_type)
!!$       end do
!!$       if(sum_dmdt_s==0.0_PS.and.sum_dmdt_r==0.0_PS) cycle
!!$       write(*,'("ck_mlttend>k,i,j,sumdr,sumds",3I5,2ES15.6)') &
!!$            KD(n),ID(n),JD(n),sum_dmdt_r*gr%dt,sum_dmdt_s*gs%dt
!!$    enddo
!!$  end subroutine ck_mlttend

  subroutine update_mesrc(ag,gr,gs,ga,mes_rc,nacat)
    ! **********************************************************************
    ! diagnose hydrometeor flags
    ! **********************************************************************
    ! rain group, ice group
    type (Group), intent(inout) :: gr,gs
    type (group), dimension(*),intent(inout) :: ga
    ! thermo variable object
    type (AirGroup), intent(inout)  :: ag
    ! message from reality-check
    integer,dimension(*),intent(inout)   :: mes_rc
    integer,intent(in)  :: nacat
    ! 3D coordinates
    !integer :: ID(*),JD(*),KD(*)
    ! local vars
    integer :: n,i,ic
    real(PS),dimension(LMAX)         :: M_tot, M_v, M_tr, M_ts

    ! +++ initialization +++

    do n=1,ag%L
      mes_rc(n) = 1
      M_v(n) = ag%tv(n)%rv*ag%tv(n)%den
      M_tr(n)=0.0_PS
      M_ts(n)=0.0_PS
      gr%mark_cm(n)=0
      gs%mark_cm(n)=0
    enddo
    do ic=1,nacat
      do n=1,ag%L
        ga(ic)%mark_cm(n)=0
      enddo
    enddo

    do i=1,gr%n_bin
      do n=1,gr%L
        M_tr(n)=M_tr(n)+gr%MS(i,n)%mass(1)
      enddo
    enddo
    do i=1,gs%n_bin
      do n=1,gs%L
        M_ts(n)=M_ts(n)+gs%MS(i,n)%mass(1)
      enddo
    enddo
    do n=1,gs%L
      M_tot(n)=M_v(n)+M_tr(n)+M_ts(n)
    enddo

    do n=1,ag%L
      if( M_tot(n) <= 0.0 ) then
        mes_rc(n) = 0
      else
        if(M_tr(n)>0.0_PS) then
          if(M_ts(n)>0.0_PS) then
            mes_rc(n)=4
          else
            gs%mark_cm(n)=3
            mes_rc(n)=2
          end if
        else
          gr%mark_cm(n)=3
          if(M_ts(n)>0.0_PS) then
            mes_rc(n)=3
          else
            gs%mark_cm(n)=3
            mes_rc(n)=1
          end if
        end if
      endif
    enddo
  end subroutine update_mesrc

  subroutine cal_den_aclen_vec(g,icond1)
    implicit none
    type (Group), intent(inout)    :: g
    integer,dimension(g%N_BIN,g%L),intent(in)    ::  icond1
    real(PS) :: map,den_ap,r_n
    ! axis ratio, c-axis/a-axis
    real(PS)                           :: alpha
    real(PS)                           :: dum
    integer :: i,n

!    do in=1,g%n_bin*g%L
!      n=(in-1)/g%N_BIN+1
!      i=in-(n-1)*g%N_BIN
    do n = 1, g%L
    do i = 1, g%N_BIN
      if(icond1(i,n)==0) then

        ! calculate density of mixed aerosols
        den_ap=g%MS(i,n)%den_ai/(1.0_PS-g%MS(i,n)%eps_map* &
               (1.0_PS-g%MS(i,n)%den_ai/g%MS(i,n)%den_as))


        ! assume that v_d=v_w+v_n where v_w is volume of water and v_n is volume of dry aerosol.
        map=g%MS(i,n)%mass(rmat)/g%MS(i,n)%con
        g%MS(i,n)%den=g%MS(i,n)%mean_mass/((g%MS(i,n)%mean_mass-map)/den_w+&
                            map/den_ap)
        !tmp       den=1.0_PS

        r_n=(map/coef4pi3/den_ap)**(1.0/3.0)

        ! +++ calculate the equivalent diameter +++
        g%MS(i,n)%len = (6.0_PS*(g%MS(i,n)%mean_mass/(PI*g%MS(i,n)%den)))**(1.0/3.0)

        ! the above diagnosis of hydrometeor density occasionally restuls in a smaller density
        ! than the aerosol density. This causes a problem in diagnosing the equibrium saturation
        ! ratio. So, arbitrary factor 1.05 is set as a minimum.

        g%MS(i,n)%len = max(r_n*1.05_PS,g%MS(i,n)%len)
        g%MS(i,n)%den=g%MS(i,n)%mean_mass/(PI/6.0_PS*g%MS(i,n)%len**3)

        ! +++ calculate spheroidal shape
        !     based on Pruppacher and Klett (1997, p.397) +++
        if( g%MS(i,n)%len <= 280.0e-4 ) then
          alpha = 1.0_PS
          g%MS(i,n)%a_len = g%MS(i,n)%len/2.0_PS
          g%MS(i,n)%c_len = g%MS(i,n)%len/2.0_PS
        else if( 280.0e-4 < g%MS(i,n)%len .and. g%MS(i,n)%len <= 1.0e-1_PS ) then
          ! NOTE: This formula is assumed to be applicable for this range.
          ! It is supposed to be applicable to the range: from 1 mm to 9 mm.
          !
          dum = g%MS(i,n)%len

          alpha = max( 1.001668_PS - 0.098055_PS*dum &
               - 2.52686_PS*dum**2 + 3.75061_PS*dum**3 &
               - 1.68692_PS*dum**4, 1.0e-5_RP)

          dum=g%MS(i,n)%len/((8.0_PS*alpha)**(1.0/3.0))
          g%MS(i,n)%a_len = dum
          g%MS(i,n)%c_len = dum * alpha
        else if( 1.0e-1_PS < g%MS(i,n)%len ) then
          ! NOTE: This formula is assumed to be applicable for this range.
          ! It is supposed to be applicable to the range: from 1 mm to 9 mm.
          ! But it is not right for more than 5mm; alpha becomes negative.
          !
          dum = min(0.9_PS,g%MS(i,n)%len)
          alpha = max( 1.001668_PS - 0.098055_PS*dum &
               - 2.52686_PS*dum**2 + 3.75061_PS*dum**3 &
               - 1.68692_PS*dum**4, 1.0e-5_RP)
          dum = g%MS(i,n)%len/((8.0_PS*alpha)**(1.0_RP/3.0_RP))
          g%MS(i,n)%a_len = dum
          g%MS(i,n)%c_len = dum * alpha
        end if
        g%MS(i,n)%semi_a = g%MS(i,n)%a_len
        g%MS(i,n)%semi_c = g%MS(i,n)%c_len
      endif

    enddo
    enddo
  end subroutine cal_den_aclen_vec

!!$  subroutine denchk1(gs,i,n,id,jd,kd,from)
!!$    type (Group),intent(in)   :: gs
!!$    integer,dimension(*) :: id,jd,kd
!!$    integer,intent(in) :: i,n
!!$    character*(*) :: from
!!$
!!$    write(*,*) "denchk1,",trim(from),i,kd(n),id(n),jd(n),gs%MS(i,n)%mass(1)/&
!!$       max(1.0e-30_RP,gs%IS(i,n)%v_cs*gs%MS(i,n)%con),gs%MS(i,n)%con,gs%MS(i,n)%mass(1),gs%IS(i,n)%v_cs
!!$
!!$  end subroutine denchk1

!  subroutine lenchk1(gs,id,jd,kd,from)
!    type (Group),intent(in)   :: gs
!    integer,dimension(*) :: id,jd,kd
!    integer :: i,n,in
!    character*(*) :: from
!
!    do in=1,gs%n_bin*gs%L
!      n=(in-1)/gs%N_BIN+1
!      i=in-(n-1)*gs%N_BIN
!      if(gs%IS(i,n)%d>1.0e-20) then
!        write(*,*) "lenchk1:",trim(from),i,kd(n),id(n),jd(n),gs%MS(i,n)%mass(1) &
!         ,gs%MS(i,n)%con,gs%MS(i,n)%a_len,gs%MS(i,n)%c_len,gs%IS(i,n)%d
!      endif
!      if(isnan(gs%MS(i,n)%a_len)) then
!        write(*,*) "lenchk2:",trim(from),i,kd(n),id(n),jd(n),gs%MS(i,n)%mass(1) &
!         ,gs%MS(i,n)%con,gs%MS(i,n)%a_len,gs%MS(i,n)%c_len,gs%IS(i,n)%d
!      endif
!    enddo
!  end subroutine lenchk1


END MODULE class_Group
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!<   REFERENCE LIST  >
! Chen and Lamb (1994a)
!
! Bohm, J. P. , A general hydrodynamic theory for mixed-phase microphysics.
! Part I: drag and fall speed of hydrometeors.
!
! Hall, W. D. and H. R. Pruppacher, The Survival of Ice Particles Falling from Cirrus
! Clouds in Subsaturated Air, Journal of the Atmospheric Sciences, 1976.
!
! Low and List (1982)
! Manton and Cotton (1977)
