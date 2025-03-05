#include "scalelib.h"
!OCL SERIAL
MODULE class_Mass_Bin
  use scale_io
  use maxdims
  use mod_amps_const
  use class_Thermo_Var, only: &
     Thermo_Var
  use class_Ice_Shape, only: &
     Ice_Shape
  use class_AirGroup, only: &
     AirGroup
  use acc_amps
  use par_amps
  use com_amps, only: &
     debug
  implicit none  

  private

  public :: make_Mass_Bin
  public :: get_osm
  public :: get_critrad_itr
  public :: get_critrad_anal
  public :: get_hazerad_itr
  public :: get_hazerad_anal
  public :: get_aspect_ratio_sol
  public :: get_area_ratio_sol
  public :: interp_data1d_lut_big
  public :: cal_coef_Ts3
  public :: cal_best_number
  public :: cal_meanmass_vec
  public :: cal_mass_comp_ap
  public :: cal_melt_index
  public :: fix_mass_comp
  public :: diag_sh_type_v6
  public :: diag_pardis_ap

  public :: Mass_Bin
  public :: data1d_lut
  public :: data1d_lut_big

  TYPE Mass_Bin
     ! memory is aligned as declared below for use in common block
     sequence 

     ! parameters necessary to define concentration distribution of a bin
     ! - linear distribution for bin i (#/g/cm^3). 
     !      n(x) = n_0 + k ( x - x_0 )
     !           = p(1)   + p(3) ( x - p(2))
     !
     ! - cubic distribution (#/g/cm^3). 
     !      n(x) = p(1)*x^3+p(2)*x^2+p(3)*x+p(4)
     !
     ! - constant-size distribution
     !      p(1) : fixed total concentration
     !
     ! - Marshall-Palmer distribution
     !      p(1) : total concentration, #/cm^3
     !      p(2) : mean diameter, cm
     !      n(D) = (p(1)/p(2)) exp(-D/p(2))
     ! 
     ! - lognormal distribution
     !      p(1) : mean_ln
     !      p(2) : sigma_ln
     !      n(a) = N/sqrt(2PI)/a/sigma_ln*exp(-0.5*((ln(a)-mean_ln)/sigma_ln)**2.0)
     real(DS), dimension(3)            :: p      ! 8*3

     ! Tendency by each process for the bin, not for one particle
     !
     ! 1st argument for concentration variable, and
     ! 2nd argument for mass and volume variables are the process number 
     !
     !          1: vapor deposition
     !             - case of aerosols, source by evaporation
     !
     !          2: collision between the same phases
     !
     !          3: CCN activation for aerosol and liquid
     !             DHF for ice
     !
     !          4: collision between rain and ice hydrometeors
     !
     !          5: collision-breakup
     !
     !          6: hydrodynamic breakup
     !             DHF for aerosol
     !
     !          7: ice nucleation (sorption/deposition)
     !
     !          8: ice nucleation (contact)
     !             - case of cloud group (the 1st bin of rain),
     !               this is a transfer by
     !               contact nucleation
     !          
     !          9: secondary nucleation
     !             - case of liq group
     !               autoconversion process of the first bin
     !
     !         10: melting-shedding
     !
     !         11: ice nucleation (immersion freezing nucleation)
     !
     !
     ! 1st argument for mass tendency
     !          1: total mass
     !          2: mass by riming
     !          3: mass by aggregation
     !          4: mass of a mean ice crystal
     !          5: total mass of aerosols inside
     !          6: mass of a soluble material of aerosols
     !          7: mass of a insoluble material of aerosols
     !          8: mass of melt water
     ! 1st argument for volume tendency
     !          1: V_cs
     !          2: V_a
     !          3: V_c
     !          4: V_d
     !          5: V_r
     !          6: V_rime 
     !          7: V_agg  
     real(8), dimension(1+mxnmasscomp,mxntend)   :: dmassdt   ! 4*9*12
     real(8), dimension(mxnnonmc,mxntend)   :: dvoldt       ! 4*7*12
     real(8), dimension(mxntend)     :: dcondt        ! 4*12
!tmp     real(PS), pointer, dimension(:)     :: dcondt
!tmp     real(PS), pointer, dimension(:,:)   :: dmassdt
!tmp     real(PS), pointer, dimension(:,:)   :: dvoldt

     ! mass components
     ! for aerosol particles
     ! argument 1: total mass
     !          2: mass of a soluble material of aerosols
     !          3: mass of a insoluble material of aerosols
     ! for liquid particles
     ! argument 1: total mass
     !          2: total mass of aerosols inside
     !          3: mass of a soluble material of aerosols
     !          4: mass of a insoluble material of aerosols
     ! for solid particles
     ! argument 1: total mass
     !          2: total mass by riming in the bin
     !          3: total mass by aggregation in the bin
     !          4: mass of a mean hexagonal ice crystal + aerosol inside
     !          5: total mass of aerosols inside
     !          6: mass of a soluble material of aerosols
     !          7: mass of a insoluble material of aerosols
     !          8: mass of melt water
     !          m_1=m_2+m_3+m_4+m_8
     !          m_5=m_6+m_7
     real(PS), dimension(1+mxnmasscomp)   :: mass  ! 4*9

     ! concentration variable
     ! argument 1: activated IN concentration for contact parameter diagnosis
     !          2: not used 
     real(PS), dimension(mxnvol)   :: vol   ! 4*2

     ! Lagragian mass tendency for the bin, not for a particle
     !   - It is used for calculating sfc temperature of a hydrometeor 
     !     and melting-shedding process, so especially important for
     !     solid hydrometeors.
     ! 
     ! argument
     !          1: vapor deposition
     !          2: riming process
     real(PS), dimension(2)     :: Ldmassdt   ! 4*2
!tmp     real(PS), pointer, dimension(:)     :: Ldmassdt
!tmp     real(PS), allocatable, dimension(:)     :: Ldmassdt


     ! coeficients to reduce redundunt calculation
     ! 1: coefficient for vapor deposition calculation
     real(PS), dimension(2)            :: coef      ! 4*2
!tmp     real(PS), pointer, dimension(:)            :: coef
!tmp     real(PS), allocatable, dimension(:)            :: coef

     ! --- information for the bin
     ! lower and upper limits of the bin (g)
     real(PS)    :: lmt
     real(PS)    :: umt
     
     ! --- optionally provided ---
     ! concentration of the bin (#/cm3)
     real(PS)    :: con

     ! mean (concentration-weighted) mass defined as
     !           mean_mass = (total mass in the bin)/(total concentration in the bin)
     real(PS)    :: mean_mass

     ! characteristic length
     ! (total surface area)/
     ! (perimeter projected normal to the flow direction)
     ! 
     ! in case of sphere, it is the diameter
     real(PS)    :: len
     
     ! apparent density of the meteor (g/cm^3)
     real(PS)    :: den

     ! density of soluble and insoluble materials in an aerosol
     real(PS)    :: den_as,den_ai

     ! length along a-axis (cm) on mean mass in the bin
     real(PS)    :: a_len

     ! length along c-axis (cm) on mean mass in the bin
     real(PS)    :: c_len

     ! mass fraction of soluble material in the aerosol
     real(PS) :: eps_map
     
     ! --- diagnostic variables ---
     ! terminal velocity (cm/s)
     real(PS)    :: vtm
     
     ! particle temperature (K)
     real(PS)    :: tmp
     
     ! mean ventilation coefficient for mass
     real(PS)    :: fv

     ! kinetic effect
     real(PS)    :: fkn

     ! parameter of local ventilation coefficients for the c and a axes
     real(PS)    :: fac

     ! mean ventilation coefficient for heat
     real(PS)    :: fh

     ! Reynolds number
     real(PS)    :: Nre

     ! capacitance for vapor field (cm)
     real(PS)    :: CAP

     ! capacitance for vapor field (cm) for a hexagonal monocrystal
     real(PS)    :: CAP_hex

     ! lengths of semi-axis of a spheroid circumscribing all ice crystals
     real(PS)    :: semi_a, semi_c

     ! saturation vapor pressure at the surface of the hydrometeors
     real(PS)    :: e_sat

     ! de-activation radius for cloud droplets (used for evaporation)
     real(PS)    :: r_act
     ! critical radius
     real(PS)    :: r_crt

     ! lower bound of contact parameter [deg]
     real(PS) :: th00_cp

!!c     ! saturation mixing ratio at hydrometeor surface
!!c     real(PS)    :: rvs

     ! type of distribution
     ! 0 : constant-size disbtribution
     ! 1 : linear distribution w.r.t. mass
     ! 2 : Marshall-Palmer distribution w.r.t. diameter
     ! 3 : Lognormal distribution w.r.t. radius
     ! 4 : uniform distribution w.r.t. radius
     integer     :: dis_type


     ! marker
     ! this is used to reduce unnecessary computation
     integer          :: mark

     ! this is used to indicate evaporation of whole hydrometeors occured
     integer          :: inevp

     ! indication of melt water on ice surface
     ! 0: dry growth mode
     ! 1: wet growth mode
     ! 2: totally wet mode
     integer :: inmlt


  end type Mass_Bin

  integer,parameter :: nx_lut_max=100
  type data1d_lut
    sequence 
    integer :: n
    integer :: align_data1d_lut ! CHIARUI
    real(PS) :: xs,dx
    real(PS),dimension(nx_lut_max) :: y
  end type data1d_lut
  type data1d_lut_big
    sequence 
    integer :: n
    integer :: align_data1d_lut_big ! CHIARUI
    real(PS) :: xs,dx
    real(PS),dimension(501) :: y
  end type data1d_lut_big

  !
  ! Dust fraction (tuning parameter)
  !   This fraction of aerosols and cloud droplets 
  !   incldue ice nuclei.
  !
  !real(PS),parameter :: frac_dust=0.01_PS 
  
CONTAINS

  ! ************* Optional Constructor for all types **********************
!tmp  function make_Mass_Bin (dlmt, dumt, ddtype, dnaxis, dnrime, dnagg,dnap,dnmelt, &
!tmp       dnvol, dntp, &
!tmp       dcon, dmass_t, dmass_r, dmass_agg, &
!tmp       dlen, dden, dden_as, dden_ai,&
!tmp       dalen, dclen, dvol_r, dvol_a, &
!tmp       deps_map, dvtm, dfv, dfkn, dfac, dfh, dNre, dtmp, dcap, dmrk) &
!tmp       result (ms)
  subroutine make_Mass_Bin (ms, dlmt, dumt, ddtype, dnaxis, dnrime, dnagg,dnap,dnmelt, &
       dnvol, dntp, &
       dcon, dmass_t, dmass_r, dmass_agg, &
       dlen, dden, dden_as, dden_ai,&
       dalen, dclen, dvol_r, dvol_a, &
       deps_map,&
       dvtm, dfv, dfkn, dfac, dfh, dNre, dtmp, dcap, dmrk) 

!tmp    real(PS), optional, intent(in)    :: dlmt, dumt, dcon, dmass_t, &
!tmp         dmass_r, dmass_agg, dlen, dden, dden_as,dden_ai,dvol_r, dvol_a,deps_map
!tmp    integer, optional, intent(in)     :: dnaxis, dnrime, dnagg,dnap,dnmelt,dnvol, &
!tmp         ddtype, dmrk, dntp
!tmp    real(PS), optional, intent(in)    :: dalen, dclen, dvtm, dtmp
!tmp    real(PS), optional, intent(in)    :: dfv, dfkn, dfac, dfh, dNre, dcap

    real(PS), intent(in)    :: dlmt, dumt, dcon, dmass_t, &
         dmass_r, dmass_agg, dlen, dden, dden_as,dden_ai,dvol_r, dvol_a,deps_map
    integer, intent(in)     :: dnaxis, dnrime, dnagg,dnap,dnmelt,dnvol, &
         ddtype, dmrk, dntp
    real(PS), intent(in)    :: dalen, dclen, dvtm, dtmp
    real(PS), intent(in)    :: dfv, dfkn, dfac, dfh, dNre, dcap
    type (Mass_Bin)  :: ms

    integer          :: i, j, var_Status, naxis, nrime, nagg, nvol,nap,nmelt

    ms%lmt = dlmt
    ms%umt = dumt
    ms%con = dcon
    ms%len = dlen
    ms%den = dden
    ms%den_as = dden_as
    ms%den_ai = dden_ai

    ms%a_len = dalen
    ms%c_len = dclen

    ms%semi_a = 0.0_PS
    ms%semi_c = 0.0_PS

    ms%eps_map = deps_map
    ms%vtm = dvtm

    ms%tmp = dtmp
    ms%fv = dfv

    ms%fkn = dfkn

    ms%fac = dfac
    ms%fh = dfh

    ms%Nre = dNre
    ms%CAP = dcap

    ms%CAP_hex=0.0_PS

    ms%e_sat=0.0_PS

    ms%dis_type = ddtype

    ms%mean_mass=0.0_PS

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! If you want to add more mass components, you need to modify the following.
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! +++ number of axis +++
    naxis = dnaxis
    ! +++ riming component +++
    nrime = dnrime
    ! +++ aggregation component +++
    nagg = dnagg
    
    ! +++ aerosol component +++
    nap = dnap

    ! +++ melt component +++
    nmelt = dnmelt

    ! +++ volume component +++
    nvol = dnvol



    ! +++ allocate memory to value terms +++
!tmp    call allocate_make_mass_bin



    ms%mass=0.0_PS
    ms%mass(1) = dmass_t
    if( nrime == 1 ) then
       ms%mass(imr) = dmass_r
    end if
    if( nagg == 1 ) then
       ms%mass(ima) = dmass_agg
    end if

    ms%r_act=0.0_PS
    ms%r_crt=0.0_PS

    if( nvol == 2 ) then
       ms%vol(1) = dvol_r
       ms%vol(2) = dvol_a
    end if


    ms%mark = dmrk

    ms%inmlt=0


    ! ++++++ initialize it ++++++
    ms%dcondt = 0.0_PS
    ms%dmassdt = 0.0_PS
    ms%Ldmassdt = 0.0_PS
    ms%dvoldt = 0.0_PS

    ms%coef=0.0_PS

    ms%inevp=0

!tmp  contains
!tmp    subroutine allocate_make_mass_bin
!tmp      integer :: item(8)
!tmp      item(1:8)=0
!tmp      nullify(ms%mass)
!tmp      allocate( ms%mass(1+nrime+nagg+1+nap+nmelt), stat = item(1))
!tmp      nullify(ms%vol)
!tmp      allocate( ms%vol(nvol), stat = item(2))
!tmp      nullify(ms%p)
!tmp      if( ms%dis_type == 0 ) then
!tmp         ! --- constant-size distribution ---
!tmp         allocate( ms%p(1), stat = item(3))
!tmp      else if( ms%dis_type == 1) then
!tmp         ! --- linear distribution ---
!tmp         allocate( ms%p(3), stat = item(3))
!tmp      else if( ms%dis_type == 2 ) then
!tmp         ! --- Marshall-Palmear distribution ---
!tmp         allocate( ms%p(2), stat = item(3))
!tmp      else if( ms%dis_type == 3 ) then
!tmp         ! --- Log normal distribution ---
!tmp         allocate( ms%p(2), stat = item(3))
!tmp      else if( ms%dis_type == 4 ) then
!tmp         ! --- uniform distribution --
!tmp         allocate( ms%p(2), stat = item(3))
!tmp      end if
!tmp      ! +++ allocate the tendency matrix +++
!tmp      nullify(ms%dcondt)
!tmp      allocate( ms%dcondt(dntp), stat = item(4))
!tmp      nullify(ms%dmassdt)
!tmp      allocate( ms%dmassdt(1+nrime+nagg+1+nap+nmelt,dntp), stat = item(5))
!tmp      nullify(ms%Ldmassdt)
!tmp      allocate( ms%Ldmassdt(2), stat = item(6))
!tmp      nullify(ms%dvoldt)
!tmp      allocate( ms%dvoldt(1+naxis+1,dntp), stat = item(7))
!tmp      nullify(ms%coef)
!tmp      allocate( ms%coef(2), stat = item(8))
!tmp      if( any(item/=0) ) Stop 'Allocation failed at make_mass_bin'
!tmp    end subroutine allocate_make_mass_bin

!tmp  end function make_Mass_Bin
  end subroutine make_Mass_Bin

!tmp  function make_Mass_Bin2 () result (ms)
!tmp    type (Mass_Bin)  :: ms
!tmp
!tmp    integer          :: i, var_Status
!tmp
!tmp    allocate( ms%mass(3), stat = var_Status)
!tmp    allocate( ms%vol(2), stat = var_Status)
!tmp    allocate( ms%p(3), stat = var_Status)
!tmp
!tmp    ! +++ allocate the tendency matrix +++
!tmp    allocate( ms%dcondt(10), stat = var_Status)
!tmp    if(var_Status /= 0 ) stop "Memory not available for dcondt &
!tmp         in make_Mass_Bin"
!tmp    allocate( ms%dmassdt(3,10), stat = var_Status)
!tmp    if(var_Status /= 0 ) stop "Memory not available for dmassdt &
!tmp         in make_Mass_Bin"
!tmp    allocate( ms%Ldmassdt(2), stat = var_Status)
!tmp    if(var_Status /= 0 ) stop "Memory not available for Ldmassdt &
!tmp         in make_Mass_Bin"
!tmp    allocate( ms%dvoldt(7,10), stat = var_Status)
!tmp    if(var_Status /= 0 ) stop "Memory not available for dvoldt &
!tmp         in make_Mass_Bin"
!tmp  end function make_Mass_Bin2


  ! ***************** Public Deconstructor for a Mass_Bin type *****************
!tmp  subroutine delete_Mass_Bin (name)
!tmp    type (Mass_Bin), intent(inout) :: name
!tmp    integer :: item(8)
!tmp    item=0
!tmp    deallocate( name%coef, stat = item(8))
!tmp    nullify(name%coef)
!tmp    deallocate( name%dvoldt, stat = item(7))
!tmp    nullify(name%dvoldt)
!tmp    deallocate( name%Ldmassdt, stat = item(6))
!tmp    nullify(name%Ldmassdt)
!tmp    deallocate( name%dmassdt, stat = item(5))
!tmp    nullify(name%dmassdt)
!tmp    deallocate( name%dcondt, stat = item(4))
!tmp    nullify(name%dcondt)
!tmp    deallocate( name%p, stat = item(3))
!tmp    nullify(name%p)
!tmp    deallocate( name%vol, stat = item(2))
!tmp    nullify(name%vol)
!tmp    deallocate( name%mass, stat = item(1))
!tmp    nullify(name%mass)
!tmp    if( any(item/=0) ) Stop 'Deallocation failed at make_mass_bin'
!tmp
!tmp  end subroutine delete_Mass_Bin

!!c  subroutine assign_con_mass( MS, nN, nM)
!!c    type (Mass_Bin), intent(inout)    :: MS
!!c    real(PS), optional, intent(in)      :: nN, nM
!!c    if( nN < 0.0_PS .OR. nM < 0.0_PS ) then
!!c       write(*,*) "assign_con_mass> nN and nM are negative!"
!!c       stop
!!c    end if
!!c    MS%con = nN
!!c    MS%mass(1) = nM
!!c  end subroutine assign_con_mass

!!c  ! ***************** Public Constructor for a Mass_Bin type *****************
!!c  function Mass_Bin_ (dlmt, dumt, &
!!c       dmxr, dcon, dmass, ddia, dden, &
!!c       dche, dvtm, dtmp, dfv, drvs) result (ms)
!!c    real(PS), optional, intent(in)    :: dlmt, dumt, dmxr, dcon, dmass, ddia, dden
!!c    integer, optional, intent(in) :: dche
!!c    real(PS), optional, intent(in)    :: dvtm, dtmp, dfv, drvs
!!c    type (Mass_Bin)  :: ms
!!c    ms = Mass_Bin (dlmt, dumt, dmxr, dcon, dmass, ddia, dden, dche, &
!!c         dvtm, dtmp, dfv, drvs)
!!c  end function Mass_Bin_
!!c

!  FUNCTION addliquid(p, q)
!  TYPE (liquid), INTENT(IN) :: p,q
!  TYPE (liquid), INTENT(IN) :: addliquid
!  addliquid%dia = p%dia + q%dia
!  addliquid%mxr = p%mxr + q%mxr
!  END FUNCTION addliquid
!!$  subroutine undef_Mass_Bin(phase,ms)
!!$    integer, intent(in)  :: phase
!!$    type (Mass_Bin), intent(inout) :: ms
!!$    real(PS),parameter :: undef=999.9e+30
!!$    if(phase==1.or.phase==2) then
!!$      ms%p=undef
!!$    endif
!!$    ms%dmassdt=undef
!!$    ms%dvoldt=undef
!!$    ms%dcondt=undef
!!$    ms%mass=undef
!!$    ms%vol=undef
!!$    ms%Ldmassdt=undef
!!$    ms%coef=undef
!!$    ms%con=undef
!!$    ms%mean_mass=undef
!!$    ms%len=undef
!!$    ms%den=undef
!!$    ms%den_as=undef
!!$    ms%den_ai=undef
!!$    ms%a_len=undef
!!$    ms%c_len=undef
!!$    ms%eps_map=undef
!!$    ms%vtm=undef
!!$    ms%tmp=undef
!!$    ms%fv=undef
!!$    ms%fkn=undef
!!$    ms%fac=undef
!!$    ms%fh=undef
!!$    ms%Nre=undef
!!$    ms%CAP=undef
!!$    ms%CAP_hex=undef
!!$    ms%semi_a=undef
!!$    ms%semi_c=undef
!!$    ms%e_sat=undef
!!$    ms%r_act=undef
!!$    ms%r_crt=undef
!!$    ms%mark=0
!!$    ms%inevp=0
!!$    ms%inmlt=0
!!$
!!$  end subroutine undef_Mass_Bin

!!$  subroutine cal_tsfc_iter(T_inf,A1,B1,phase2,estbar,esitbar,out,em)
!!$    real(PS),intent(in) :: T_inf,A1,B1
!!$    integer,intent(in) :: phase2
!!$    ! saturation vapor lookup table
!!$    real(DS)  :: estbar(150),esitbar(111)
!!$    real(PS),intent(inout) :: out
!!$    integer,intent(inout) :: em
!!$
!!$    REAL(PS) :: tol=1.0e-6
!!$    integer :: ITMAX
!!$    real(PS) :: eps
!!$!    PARAMETER (ITMAX=50,EPS=epsilon(out))
!!$    PARAMETER (ITMAX=50,EPS=3.e-8)
!!$
!!$
!!$    !  Using Brent's method, find the root of a function func known to lie between x1 and x2.
!!$    !  The root, returned as zbrent, will be refined until its accuracy is tol.
!!$      
!!$    !  Parameters: Maximum allowed number of iterations, and machine floating-point precision.
!!$    INTEGER :: iter
!!$    REAL(PS) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
!!$
!!$!!c             tmp=th_var%T
!!$!!c             tmp=200.0
!!$          
!!$!!c             tmp=T_0
!!$!!c             write(*,*) phase
!!$!!c             do iter=1,10
!!$!!c                write(*,*) "tmp,esat",tmp,get_sat_vapor_pres_lk(phase,tmp,estbar,esitbar)
!!$!!c                tmp=tmp+5.0 
!!$!!c             end do
!!$!!c             stop
!!$    ! ++++++++++++++++++++++++++++
!!$    em=0
!!$
!!$    a=170.16_PS
!!$!!c    a=T_inf
!!$    a=273.16
!!$    b=300.16_PS
!!$!!c             b=273.16_PS
!!$    iter=1
!!$2221  continue
!!$    fa=func(a)
!!$2222  continue
!!$    fb=func(b)
!!$      
!!$    if(fa.gt.0.0_PS.and.fb.gt.0.0_PS) then
!!$       iter=iter+1
!!$       if(iter>50.or.a<123.0) then
!!$!!c          write(*,'("iter at Ts large 1",2I5,10ES15.6)') phase2,iter,T_inf,A1,B1,fa,fb,a,b
!!$          em=1
!!$          out=T_inf
!!$          return
!!$       end if
!!$!!c       a=a-50.0_PS
!!$!!c       a=a/1.2_PS
!!$       a=a-10.0_PS
!!$!!c       stop
!!$       goto 2221
!!$    elseif(fa.lt.0.0_PS.and.fb.lt.0.0_PS) then
!!$!!c       if((fa.gt.0.0_PS.and.fb.gt.0.0_PS).or.(fa.lt.0.0_PS.and.fb.lt.0.0_PS)) then
!!$       iter=iter+1
!!$       if(iter>50) then
!!$          write(*,'("iter at Ts large 2",2I5,10ES15.6)') phase2,iter,T_inf,A1,B1,fa,fb,a,b
!!$          em=2
!!$          out=T_inf
!!$          return
!!$       end if
!!$!!c       stop  
!!$!!c                b=b*100.0
!!$       b=b*2.0_PS
!!$!!c       b=b+10.0_PS
!!$!!c       b=b+50.0_PS
!!$       goto 2222
!!$    end if
!!$      
!!$    c=b
!!$    fc=fb
!!$    do iter=1,ITMAX
!!$       if((fb.gt.0.0_PS.and.fc.gt.0.0_PS).or.(fb.lt.0.0_PS.and.fc.lt.0.0_PS))then
!!$          c=a ! Rename a, b, c and adjust bounding interval d.
!!$          fc=fa
!!$          d=b-a
!!$          e=d
!!$       endif
!!$       if(abs(fc).lt.abs(fb)) then
!!$          a=b
!!$          b=c
!!$          c=a
!!$          fa=fb
!!$          fb=fc
!!$          fc=fa
!!$       endif
!!$       tol1=2.0_PS*EPS*abs(b)+0.5_PS*tol ! Convergence check.
!!$       xm=0.5_PS*(c-b)
!!$       if(abs(xm).le.tol1 .or. fb.eq.0.)then
!!$          out=b
!!$          goto 1111
!!$       endif
!!$       if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
!!$          s=fb/fa ! Attempt inverse quadratic interpolation.
!!$          if(a.eq.c) then
!!$             p=2.0_PS*xm*s
!!$             q=1.0_PS-s
!!$           else
!!$             q=fa/fc
!!$             r=fb/fc
!!$             p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.0_PS))
!!$             q=(q-1.0_PS)*(r-1.0_PS)*(s-1.0_PS)
!!$          endif
!!$          if(p.gt.0.) q=-q ! Check whether in bounds.
!!$          p=abs(p)
!!$          if(2.0_PS*p .lt. min(3.0_PS*xm*q-abs(tol1*q),abs(e*q))) then
!!$             e=d ! Accept interpolation.
!!$             d=p/q
!!$          else
!!$             d=xm ! Interpolation failed, use bisection.
!!$             e=d
!!$          endif
!!$       else ! Bounds decreasing too slowly, use bisection.
!!$          d=xm
!!$          e=d
!!$       endif
!!$       a=b ! Move last best guess to a.
!!$       fa=fb
!!$       if(abs(d) .gt. tol1) then ! Evaluate new trial root.
!!$          b=b+d
!!$       else
!!$          b=b+sign(tol1,xm)
!!$       endif
!!$       fb=func(b)
!!$    enddo
!!$    write(*,'("brent exceeding maximum iterations. n,iter",2I5,12ES15.6)') &
!!$         iter,a,b,fa,fb
!!$!!c             write(*,'("W,a_c,N_act,M_act,T,r_e",10ES15.6)') ag%TV(n)%W, a_c, N_act,M_act, ag%TV(n)%T,r_e
!!$      
!!$    out=T_inf
!!$      
!!$      
!!$1111  continue
!!$
!!$  contains      
!!$    function func(x) result(out)
!!$    use class_Thermo_Var, only: &
!!$       get_sat_vapor_pres_lk
!!$      real(PS),intent(in) :: x
!!$      real(PS) :: out
!!$!!c      out=x-B1+A1*get_sat_vapor_pres(phase2,x)/x
!!$      out=x-B1+A1*get_sat_vapor_pres_lk(phase2,x,estbar,esitbar)/x
!!$    end function func
!!$  end subroutine cal_tsfc_iter

!!$  subroutine cal_some_vel(ms,ishape,mode,alpha,OMEGA,A_e,A_c,P,q,char_len)
!!$    use class_Ice_Shape, only: &
!!$       get_total_sfc_area, &
!!$       get_effect_area, &
!!$       get_circum_area, &
!!$       get_perimeter, &
!!$       get_vip, &
!!$       den_i
!!$    type (Mass_Bin), intent(in)   :: ms
!!$    type (Ice_Shape), intent(in)  :: ishape
!!$    integer              :: mode
!!$    ! axial ratio, c/(2*(a+d))
!!$    real(PS)                    :: alpha
!!$    ! total surface area of the body
!!$    real(PS)                    :: OMEGA
!!$    ! effective cross-sectional area
!!$    real(PS)                    :: A_e
!!$    ! circumscribed cross-sectional area
!!$    real(PS)                    :: A_c
!!$    ! Perimeter 
!!$    real(PS)                    :: P
!!$    ! the ratio of the effective over the circumscribed cross-sectional area, Ae/A
!!$    real(PS)                    :: q,q_ic,q_cs
!!$    ! charasteristic length defined for Bohm (1992)
!!$    real(PS)                    :: char_len
!!$    ! bulk density of aggregates (Bohm 1989)
!!$    real(PS)                    :: b_den
!!$    ! aspect ratio assumed for rosettes and irregular crystals.
!!$    real(PS),parameter :: rho_re=0.25
!!$
!!$    ! this should correspond with aspect ratio diagnosis
!!$    real(PS) :: mass_q_agg,mass_q_grp
!!$
!!$    real(PS) :: mass1
!!$    real(PS),parameter :: mrat_agg=10.0_PS,mrat_grp=10.0_PS ! sense test grp=10
!!$
!!$    real(PS)  :: dum
!!$    if( ishape%sh_type == 1 .or. ishape%sh_type == 2 ) then
!!$       ! +++ get aspect ratio +++ 
!!$       alpha=ishape%phi_cs
!!$       
!!$       ! +++ calculate the total surface area +++
!!$       OMEGA = get_total_sfc_area( ishape, ms%semi_a, ms%semi_c, mode)
!!$       
!!$       ! +++ calculate the area ratio +++
!!$       ! ++++++ calculate the effective cross section ++++++
!!$       A_e = get_effect_area( ishape,ishape%sh_type,ms%semi_a, ms%semi_c, mode)
!!$       
!!$       ! ++++++ calculate the circumscribed cross-sectional area ++++++
!!$       A_c = get_circum_area( ishape,ishape%sh_type, ms%semi_a, ms%semi_c, mode)
!!$       
!!$!!c       if( ishape%sh_type == 2 ) then
!!$!!c          ! --- case of lightly rimed ice crystal ---
!!$!!c          ! NOTE: the definition by volume is necessary to be studied.
!!$!!c          dum=(4.0_PS*PI/3.0_PS)*ms%semi_c*ms%semi_a**2.0
!!$!!c          q = min( ( ishape%V_ic + ms%vol(1) )/dum, 1.0_PS)
!!$!!c       else
!!$       if( ishape%habit>=5 .and. ishape%is_mod(2)==1) then
!!$          ! use Heymsfield et al (2002)'s side planes
!!$          q=max(0.18_PS,min(0.88_PS,(ms%den/0.35_PS)**(1.0/2.34_PS)))
!!$       else
!!$          q=A_e/A_c
!!$       end if
!!$!!c       end if
!!$       ! +++ calculate the characteristic length based on
!!$       !     perimeter and total surface area +++
!!$       P = get_perimeter( ishape, ms%semi_a, ms%semi_c, mode)
!!$       
!!$       ! +++ calculate characteristic length +++
!!$       char_len = 2.0*sqrt(A_c/PI)
!!$    else if( ishape%sh_type <= 4 ) then
!!$!!c       if(mode==2) then
!!$!!c          alpha=0.25_PS
!!$!!c       else
!!$       alpha = ms%semi_c/ms%semi_a
!!$!!c       end if     
!!$   
!!$       ! +++ calculate the total surface area +++
!!$       OMEGA = get_total_sfc_area( ishape, ms%semi_a, ms%semi_c, mode)
!!$       
!!$       ! +++ calculate the area ratio +++
!!$       ! ++++++ calculate the effective cross section ++++++
!!$       A_e = get_effect_area( ishape, ishape%sh_type,ms%semi_a, ms%semi_c, mode)
!!$       
!!$       ! ++++++ calculate the circumscribed cross-sectional area ++++++
!!$       A_c = get_circum_area( ishape,ishape%sh_type, ms%semi_a, ms%semi_c, mode)
!!$
!!$       ! +++ calculate the characteristic length based on
!!$       !     perimeter and total surface area +++
!!$       P = get_perimeter( ishape, ms%semi_a, ms%semi_c, mode)
!!$
!!$       ! +++ calculate characteristic length +++
!!$       char_len = 2.0*sqrt(A_c/PI)
!!$
!!$       if(mode==2) then
!!$!!c          select case(ishape%habit)
!!$!!c          case (1)
!!$!!c             ! 1 : hexagonal plate
!!$!!c             b_den=0.9_PS
!!$!!c          case (2)
!!$!!c             ! 2 : broad-branch crystal
!!$!!c             b_den=0.5_PS
!!$!!c          case (3)
!!$!!c             ! 3 : columnar crystal
!!$!!c             b_den=0.8_PS
!!$!!c          case (4)
!!$!!c             ! 4 : rosette bullet
!!$!!c             b_den=0.8_PS
!!$!!c          case (5)
!!$!!c             ! 5 : irregular ice crystals
!!$!!c             b_den=0.9_PS
!!$!!c          end select
!!$!!c          q = 24.0_PS*ms%mean_mass/(PI*b_den*char_len**3.0)
!!$
!!$       
!!$          ! +++ calculate the area ratio of ice crystals +++
!!$          ! ++++++ calculate the effective cross section ++++++
!!$          if( ishape%habit>=5 ) then
!!$             ! use Heymsfield et al (2002)'s side planes
!!$             q_ic=max(0.18_PS,min(0.88_PS,(ms%mass(imc)/ms%con/ishape%V_ic/0.35_PS)&
!!$                  **(1.0/2.34_PS)))
!!$          else
!!$             q_ic=get_effect_area( ishape,1,ms%a_len, ms%c_len, mode)/&
!!$                  get_circum_area( ishape,1,ms%a_len, ms%c_len, mode)
!!$          end if
!!$!!c          q_cs=(ms%mean_mass/den_i/ishape%V_cs)**(2.0/3.0)
!!$!!c          q_cs=ms%mean_mass/den_i/ishape%V_cs
!!$!!c          q_cs=(ms%mass(imc)+ms%mass(ima)+ms%mass(imr))/ms%con/den_i/ishape%V_cs
!!$          q_cs=(ms%mass(imc)+ms%mass(ima)+ms%mass(imr))/ms%con/den_i/&
!!$                get_vip(ishape%is_mod(2),ishape%phi_cs,ishape%semi_aip)
!!$          
!!$          if(q_cs>1.0e-30_PS) then
!!$             mass_q_agg=ms%mass(imc)/ms%con*mrat_agg
!!$             mass1=(ms%mass(imc)+ms%mass(ima))/ms%con
!!$
!!$             if(mass1<mass_q_agg) then
!!$!!c                q=10.0**((log10(q_ic)-log10(q_cs))/(log10(ms%mass(imc)/ms%con)-log10(mass_q_agg))*&
!!$!!c                     (log10(mass1)-log10(mass_q_agg))+log10(q_cs))
!!$                q=q_cs*(mass1/mass_q_agg)**&
!!$                 (log10(q_ic/q_cs)/log10(ms%mass(imc)/ms%con/mass_q_agg))
!!$
!!$
!!$
!!$                if(q<0.0.or.q>2.0.or.(q<0.0.and.q>2.0)) then
!!$                   write(*,*) "q is not right for agg",ishape%habit,q,q_ic,q_cs,ms%mass(imc),ms%con,ms%mean_mass
!!$                   write(*,*) ms%a_len,ms%c_len,ishape%phi_ic,ishape%psi_ic,ishape%a,ishape%c
!!$                   write(*,*) ms%mass
!!$                   stop
!!$                end if
!!$             else
!!$                q=q_cs
!!$             end if
!!$           else
!!$              q=q_ic
!!$           endif
!!$       else
!!$          q = A_e/A_c
!!$       end if 
!!$
!!$    else if( ishape%sh_type == 5 ) then
!!$!!c       if(mode==2) then
!!$!!c          alpha=0.7_PS
!!$!!c       else
!!$          alpha = ms%semi_c/ms%semi_a
!!$!!c       end if     
!!$   
!!$       ! +++ calculate the total surface area +++
!!$       OMEGA = get_total_sfc_area( ishape, ms%semi_a, ms%semi_c, mode)
!!$       
!!$       ! +++ calculate the area ratio +++
!!$       ! ++++++ calculate the effective cross section ++++++
!!$       A_e = get_effect_area( ishape, ishape%sh_type,ms%semi_a, ms%semi_c, mode)
!!$       
!!$       ! ++++++ calculate the circumscribed cross-sectional area ++++++
!!$       A_c = get_circum_area( ishape, ishape%sh_type,ms%semi_a, ms%semi_c, mode)
!!$
!!$       ! +++ calculate the characteristic length based on
!!$       !     perimeter and total surface area +++
!!$       P = get_perimeter( ishape, ms%semi_a, ms%semi_c, mode)
!!$
!!$       ! +++ calculate characteristic length +++
!!$       char_len = 2.0*sqrt(A_c/PI)
!!$
!!$!!c       if(mode==2) then
!!$!!c          ! from Bohm (1991), which assume that 10% of diameter wide is
!!$!!c          ! covered by 50%.
!!$!!c          q = 0.82_PS
!!$       
!!$       if(mode==2) then
!!$          ! +++ calculate the area ratio of ice crystals +++
!!$          ! ++++++ calculate the effective cross section ++++++
!!$          if( ishape%habit>=5 ) then
!!$             ! use Heymsfield et al (2002)'s side planes
!!$             q_ic=max(0.18_PS,min(0.88_PS,(ms%mass(imc)/ms%con/ishape%V_ic/0.35_PS)&
!!$                  **(1.0/2.34_PS)))
!!$          else
!!$             q_ic=get_effect_area( ishape,1,ms%a_len, ms%c_len, mode)/&
!!$                  get_circum_area( ishape,1,ms%a_len, ms%c_len, mode)
!!$          end if
!!$!!c          q_cs=ms%mean_mass/den_i/ishape%V_cs
!!$!!c          q_cs=(ms%mass(imc)+ms%mass(ima)+ms%mass(imr))/ms%con/den_i/ishape%V_cs
!!$          q_cs=(ms%mass(imc)+ms%mass(ima)+ms%mass(imr))/ms%con/den_i/&
!!$                get_vip(ishape%is_mod(2),ishape%phi_cs,ishape%semi_aip)
!!$
!!$          if(q_cs>1.0e-30_PS) then
!!$             mass_q_grp=(ms%mass(imc)+ms%mass(ima))/ms%con*mrat_grp
!!$             mass1=(ms%mass(imc)+ms%mass(ima)+ms%mass(imr))/ms%con
!!$
!!$             if(mass1<mass_q_grp) then
!!$!!c                q=10.0**((log10(q_ic)-log10(q_cs))/(log10(ms%mass(imc)/ms%con)-log10(mass_q_grp))*&
!!$!!c                     (log10(mass1)-log10(mass_q_grp))+log10(q_cs))
!!$                q=q_cs*(mass1/mass_q_grp)**&
!!$                 (log10(q_ic/q_cs)/log10(ms%mass(imc)/ms%con/mass_q_grp))
!!$
!!$                if(q<0.0.or.q>2.0.or.(q<0.0.and.q>2.0)) then
!!$                   write(*,*) "q is not right for grp",q,q_ic,q_cs,ms%mass(imc),ms%con,ms%mean_mass
!!$                   write(*,*) ms%mass
!!$                   stop
!!$                end if
!!$             else
!!$                q=q_cs
!!$             end if
!!$           else
!!$              q=q_ic
!!$           end if
!!$       else
!!$          q = A_e/A_c
!!$       end if
!!$       
!!$    end if
!!$  end subroutine cal_some_vel

!!$  subroutine cal_mass_comp(token,level,th_var,eps_ap0, ms)
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! calculate 
!!$    ! 1. mass of rime component and aggregate component
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    type (Mass_Bin), intent(inout)   :: ms
!!$!!c    type (Ice_Shape), intent(in)    :: is
!!$    type (Thermo_Var), intent(in)    :: th_var
!!$    integer, intent(in)    ::  level,token
!!$    real(PS),intent(in) ::  eps_ap0
!!$    ! minimum mass of pristene crystal (1.0um)
!!$!!c    real,parameter :: min_mass=sq_three*3.0e-12_PS*den_i
!!$!!c    real(PS) :: ratio
!!$    ! minimum mass for fragments to be aggregates or rimed aggrates
!!$    real(PS),parameter :: min_mass=1.0e-5_PS,den_max=0.758088258_PS
!!$    ! mass limit to become liquid coated.
!!$    !  0.05 corresponds  1.6% of total radius of ice sphere with liquid coating.
!!$    ! 0.155 corresponds  5.0% of total radius of ice sphere with liquid coating.
!!$    ! 0.288 corresponds 10.0% of total radius of ice sphere with liquid coating.
!!$    ! 0.5               19.5%
!!$    real(PS),parameter :: fmliq_lmt=0.05
!!$    ! limit to be dealt as ice particle
!!$!!c    real(PS),parameter :: fmice_lmt=0.999
!!$!!c    real(PS),parameter :: fmice_lmt=0.75
!!$    real(PS),parameter :: fmice_lmt=0.25
!!$
!!$    real(PS) :: mass_left,mass_stot,oldmassc(mxnmasscomp+1),rat
!!$    integer :: ifix
!!$
!!$    ! initialize 
!!$    ms%inmlt=0
!!$    ms%inevp=0
!!$    ms%eps_map=eps_ap0
!!$
!!$    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) then
!!$       ms%mass=0.0_PS
!!$       return
!!$    endif
!!$    if(token==2) then
!!$       if(level<=3) then
!!$!!c          ms%mass(ima)=max(ms%mass(imt)-ms%mass(imr)-ms%mass(imc),0.0_PS)
!!$          if( th_var%T>=T_0 ) ms%inmlt=1
!!$       elseif(level<=5) then
!!$!!c          ms%mass(ima)=max(ms%mass(imt)-ms%mass(imr)-ms%mass(imc),0.0_PS)
!!$!!c          if(ms%mass(imat)<1.0e-30_PS) then
!!$!!c             ms%mass(imat)=min(max(1.0e-30_PS,ms%mass(imat)),ms%mass(imt))
!!$!!c          end if
!!$          ms%mass(imai)=max(ms%mass(imat)-ms%mass(imas),0.0_PS)
!!$          ms%eps_map=max(0.0_PS,min(1.0_PS,ms%mass(imas)/ms%mass(imat)))
!!$          if( th_var%T>=T_0 ) ms%inmlt=1
!!$       else
!!$!!c          ms%mass(ima)=max(ms%mass(imt)-ms%mass(imr)-ms%mass(imc)-ms%mass(imw),0.0_PS)
!!$!!c          if(ms%mass(imat)<1.0e-30_PS) then
!!$!!c             ms%mass(imat)=min(max(1.0e-30_PS,ms%mass(imat)),ms%mass(imt))
!!$!!c          end if
!!$
!!$          ms%mass(imai)=max(ms%mass(imat)-ms%mass(imas),0.0_PS)
!!$          ms%eps_map=max(0.0_PS,min(1.0_PS,ms%mass(imas)/ms%mass(imat)))
!!$          if( ms%mass(imw)/ms%mass(imt)>=fmice_lmt ) then
!!$             ms%inmlt=2
!!$! <<< 2015/03 T. Hashino added for KiD test
!!$!tmp          elseif( th_var%T>=T_0+5.0) then
!!$!tmp             ms%inmlt=2
!!$! >>> 2015/03 T. Hashino added for KiD test
!!$          elseif( ms%mass(imw)/ms%mass(imt)>=fmliq_lmt ) then
!!$             ms%inmlt=1
!!$          end if
!!$       end if
!!$!!c       if(ms%eps_map>=1.0e-6) then
!!$!!c          write(*,*) "cal_mass_comp > sol frac larger than 1.0e-6 in ice"
!!$!!c          write(*,'(8ES15.6)') ms%mass(imt),ms%mass(imat),ms%mass(imas)
!!$!!c       end if
!!$
!!$       ! maintain the relative relation of predicted mass component
!!$       mass_left=ms%mass(imt)-ms%mass(imc)
!!$       if(mass_left<0.0_PS) then
!!$!!c          write(*,*) "cal_masscomp: m4>m1",ms%mass
!!$!!c          write(*,*) "con, meanmass",ms%con,ms%mass(imt)/ms%con
!!$          ms%mass(imt)=ms%mass(imc)
!!$          mass_left=0.0_PS
!!$       endif
!!$       mass_stot=ms%mass(imr)+ms%mass(ima)+ms%mass(imw)
!!$       oldmassc=ms%mass
!!$       ifix=0
!!$       if(mass_stot>1.0e-25_PS) then 
!!$          rat=mass_left/mass_stot
!!$          ms%mass(imr)=ms%mass(imr)*rat
!!$          ms%mass(ima)=ms%mass(ima)*rat
!!$          ms%mass(imw)=ms%mass(imw)*rat
!!$          ifix=1
!!$       else
!!$          ms%mass(imc)=ms%mass(imt)
!!$          ms%mass(imr)=0.0_PS
!!$          ms%mass(ima)=0.0_PS
!!$          ms%mass(imw)=0.0_PS
!!$          ifix=2
!!$       end if
!!$
!!$       if(ms%mass(imr)+ms%mass(ima)+ms%mass(imw)+ms%mass(imc)<0.99*ms%mass(imt)) then
!!$          write(*,*) "fixing failed:",ms%mass(1:8)
!!$          write(*,*) "old masscomp",oldmassc
!!$          write(*,*) "ifix,rat",ifix,rat
!!$       endif
!!$
!!$    elseif(token==1) then
!!$       if(level>=4) then
!!$!!c          if(ms%mass(rmat)<1.0e-30_PS) then
!!$!!c             ms%mass(rmat)=min(max(1.0e-30_PS,ms%mass(rmat)),ms%mass(rmt))
!!$!!c          end if
!!$          ms%mass(rmai)=max(ms%mass(rmat)-ms%mass(rmas),0.0_PS)
!!$          ms%eps_map=max(0.0_PS,min(1.0_PS,ms%mass(rmas)/ms%mass(rmat)))
!!$       end if
!!$!!c       if(ms%eps_map<1.0e-6) then
!!$!!c          write(*,*) "cal_mass_comp > sol frac less than 1.0e-6 in liq"
!!$!!c          write(*,'(8ES15.6)') ms%mass(rmt),ms%mass(rmat),ms%mass(rmas)
!!$!!c       end if
!!$    elseif(token==3) then
!!$       if(level>=4) then
!!$          ms%eps_map=max(0.0_PS,min(1.0_PS,ms%mass(ams)/ms%mass(amt)))
!!$          ms%mass(ami)=max(ms%mass(amt)-ms%mass(ams),0.0_PS)
!!$       end if
!!$    end if
!!$
!!$
!!$!!c    if( ms%mass(1) <= min_mass*ms%con) then
!!$!!c       ms%mass(4)=ms%mass(4)+ms%mass(3)
!!$!!c       ms%mass(3)=0.0_PS
!!$!!c       write(*,*) "cal_mass_comp > mass_c per perticle fixed for agg",ms%mass(1:4)/ms%con
!!$!!c       ms%mass(4)=min_mass*ms%con
!!$!!c       ms%mass(1)=max(ms%mass(1),ms%mass(4))
!!$!!c       if( ms%mass(2)+ms%mass(3) > 1.0e-20 ) then
!!$!!c          ratio=(ms%mass(1)-ms%mass(4))/(ms%mass(2)+ms%mass(3))
!!$!!c          ms%mass(2)=ratio*ms%mass(2)
!!$!!c          ms%mass(3)=ratio*ms%mass(3)
!!$!!c       else
!!$!!c          ms%mass(2)=0.0
!!$!!c          ms%mass(3)=0.0
!!$!!c       end if
!!$!!c    end if
!!$    
!!$  end subroutine cal_mass_comp

  subroutine cal_mass_comp_ap(token,level,eps_ap0, ms)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate aerosol mass component
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Mass_Bin), intent(inout)   :: ms
    integer, intent(in)    ::  level,token
    real(PS),intent(in) ::  eps_ap0

    if(token==2) then
       if(level<=3) then

       elseif(level<=5) then
          ms%mass(imai)=max(ms%mass(imat)-ms%mass(imas),0.0_PS)
          ms%eps_map=max(0.0_PS,min(1.0_PS,ms%mass(imas)/ms%mass(imat)))
       else
          ms%mass(imai)=max(ms%mass(imat)-ms%mass(imas),0.0_PS)
          ms%eps_map=max(0.0_PS,min(1.0_PS,ms%mass(imas)/ms%mass(imat)))
       end if
    elseif(token==1) then
       if(level>=4) then
          ms%mass(rmai)=max(ms%mass(rmat)-ms%mass(rmas),0.0_PS)
          ms%eps_map=max(0.0_PS,min(1.0_PS,ms%mass(rmas)/ms%mass(rmat)))
       end if
    elseif(token==3) then
       if(level>=4) then
          ms%eps_map=max(0.0_PS,min(1.0_PS,ms%mass(ams)/ms%mass(amt)))
          ms%mass(ami)=max(ms%mass(amt)-ms%mass(ams),0.0_PS)
       end if
    end if
    
  end subroutine cal_mass_comp_ap

  subroutine cal_melt_index(level,th_var,ms,em)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate 
    ! 1. mass of rime component and aggregate component
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Mass_Bin), intent(inout)   :: ms
    type (Thermo_Var), intent(in)    :: th_var
    integer, intent(in)    ::  level
    integer,intent(inout) :: em
    ! minimum mass of pristene crystal (1.0um)
!!c    real,parameter :: min_mass=sq_three*3.0e-12_PS*den_i
!!c    real(PS) :: ratio
    ! minimum mass for fragments to be aggregates or rimed aggrates
    real(PS),parameter :: min_mass=1.0e-5_PS,den_max=0.758088258_PS
    ! mass limit to become liquid coated.
    !  0.05 corresponds  1.6% of total radius of ice sphere with liquid coating.
    ! 0.155 corresponds  5.0% of total radius of ice sphere with liquid coating.
    ! 0.288 corresponds 10.0% of total radius of ice sphere with liquid coating.
    ! 0.5               19.5%
    real(PS),parameter :: fmliq_lmt=0.05
    ! limit to be dealt as ice particle
!!c    real(PS),parameter :: fmice_lmt=0.999
!!c    real(PS),parameter :: fmice_lmt=0.75
    real(PS),parameter :: fmice_lmt=0.25

    real(PS) :: mass_left,mass_stot,oldmassc(mxnmasscomp+1),rat
    integer :: ifix
    integer :: i_T_geT0,i_mlt_ge_fmice,i_mlt_ge_fmliq,i_mass_left_lt0

    i_T_geT0=0.5*(1.0+sign(1.0_PS,th_var%T-T_0))

    if(level<=3) then
      ms%inmlt=i_T_geT0
    elseif(level<=5) then
      ms%inmlt=i_T_geT0
    else
      i_mlt_ge_fmice=0.5*(1.0+sign(1.0_PS,ms%mass(imw)/ms%mass(imt)-fmice_lmt))
      i_mlt_ge_fmliq=0.5*(1.0+sign(1.0_PS,ms%mass(imw)/ms%mass(imt)-fmliq_lmt))

      ms%inmlt=i_mlt_ge_fmice*2 +   &
               (1-i_mlt_ge_fmice)*i_mlt_ge_fmliq*1

!tmp          if( ms%mass(imw)/ms%mass(imt)>=fmice_lmt ) then
!tmp             ms%inmlt=2
! <<< 2015/03 T. Hashino added for KiD test
!tmp          elseif( th_var%T>=T_0+5.0) then
!tmp             ms%inmlt=2
! >>> 2015/03 T. Hashino added for KiD test
!tmp          elseif( ms%mass(imw)/ms%mass(imt)>=fmliq_lmt ) then
!tmp             ms%inmlt=1
!tmp          end if
    end if
    
  end subroutine cal_melt_index

  subroutine fix_mass_comp(level,ms,em)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate 
    ! 1. mass of rime component and aggregate component
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Mass_Bin), intent(inout)   :: ms
    integer, intent(in)    ::  level
    integer,intent(inout) :: em

    real(PS) :: mass_left,mass_stot,oldmassc(mxnmasscomp+1),rat
    integer :: ifix
    integer :: i_T_geT0,i_mlt_ge_fmice,i_mlt_ge_fmliq,i_mass_left_lt0

    ! maintain the relative relation of predicted mass component
    i_mass_left_lt0=0.5*(1.0-sign(1.0_PS,ms%mass(imt)-ms%mass(imc)))

    mass_left=max(0.0_PS, ms%mass(imt)-ms%mass(imc))
    ms%mass(imt)=(1.0-real(i_mass_left_lt0,PS_KIND))*ms%mass(imt)+ &
                          real(i_mass_left_lt0,PS_KIND)*ms%mass(imc)
       
    mass_stot=ms%mass(imr)+ms%mass(ima)+ms%mass(imw)
!tmp       oldmassc=ms%mass
    ifix=0
    if(mass_stot>1.0e-25_PS) then 
      rat=mass_left/mass_stot
      ms%mass(imr)=ms%mass(imr)*rat
      ms%mass(ima)=ms%mass(ima)*rat
      ms%mass(imw)=ms%mass(imw)*rat
      ifix=1
    else
      ms%mass(imc)=ms%mass(imt)
      ms%mass(imr)=0.0_PS
      ms%mass(ima)=0.0_PS
      ms%mass(imw)=0.0_PS
      ifix=2
    end if

    if(ms%mass(imr)+ms%mass(ima)+ms%mass(imw)+ms%mass(imc)<0.99*ms%mass(imt)) then
      em=ifix
!tmp          write(*,*) "fixing failed:",ms%mass(1:8)
!tmp          write(*,*) "old masscomp",oldmassc
!tmp          write(*,*) "ifix,rat",ifix,rat
    endif
  end subroutine fix_mass_comp

!!$  function rime_density() result(den_r)
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! calculate
!!$    ! 1. density of rime deposits or mean of density
!!$    !
!!$    ! NOTE: The radius, and the Reynolds number of cloud droplet has to be 
!!$    !       diagnosed beforehand.
!!$    !       So is the graupel.
!!$    !       
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    real (PS)      :: den_r
!!$    den_r = 0.454_PS
!!$!!!!!!!!!!! under construction
!!$  end function rime_density

!!$  subroutine det_sh_type( ms, is, level)
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! calculate 
!!$    ! 1. volume of rime component and aggregate component
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    type (Mass_Bin), intent(inout)   :: ms
!!$    type (Ice_Shape), intent(inout)    :: is
!!$    integer, intent(in)    ::  level
!!$    real(PS)  :: m_ice
!!$
!!$    if( ms%mass(imr) > ms%mass(imt)/100.0_PS ) then
!!$       if( ms%mass(ima) <= ms%mass(imc) ) then
!!$          if( ms%mass(imr) <= ms%mass(imc) ) then
!!$             ! rimed crystal
!!$             is%sh_type = 2
!!$          else
!!$             ! graupel
!!$             is%sh_type = 5
!!$          end if
!!$       else
!!$          if( ms%mass(imr) <= ms%mass(ima) ) then
!!$             ! rimed aggregate
!!$             is%sh_type = 4
!!$          else
!!$             ! graupel
!!$             is%sh_type = 5
!!$          end if
!!$       end if
!!$    else if( ms%mass(ima) > ms%mass(imc) ) then
!!$       is%sh_type = 3             
!!$    else
!!$       ! pristine crystal
!!$       is%sh_type = 1
!!$    end if
!!$  end subroutine det_sh_type
!!$
!!$  subroutine det_sh_type2( ms, is, level)
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! calculate 
!!$    ! 1. volume of rime component and aggregate component
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    type (Mass_Bin), intent(inout)   :: ms
!!$    type (Ice_Shape), intent(inout)    :: is
!!$    integer, intent(in)    ::  level
!!$    real(PS)  :: m_ice
!!$
!!$    if( ms%mass(imr) > ms%mass(imt)/100.0 ) then
!!$       if( ms%mass(ima) <= ms%mass(imt)/100.0 ) then
!!$          if( ms%mass(imr) <= ms%mass(imc) ) then
!!$             ! rimed crystal
!!$             is%sh_type = 2
!!$          else
!!$             ! graupel
!!$             is%sh_type = 5
!!$          end if
!!$       else
!!$          if( ms%mass(imr) <= ms%mass(ima) ) then
!!$             ! rimed aggregate
!!$             is%sh_type = 4
!!$          else
!!$             ! graupel
!!$             is%sh_type = 5
!!$          end if
!!$       end if
!!$    else if( ms%mass(ima) > ms%mass(imt)/100.0 ) then
!!$       is%sh_type = 3             
!!$    else
!!$       ! pristine crystal
!!$       is%sh_type = 1
!!$    end if
!!$  end subroutine det_sh_type2
!!$
!!$  subroutine det_sh_type3( ms, is, level)
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! calculate 
!!$    ! 1. volume of rime component and aggregate component
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    type (Mass_Bin), intent(inout)   :: ms
!!$    type (Ice_Shape), intent(inout)    :: is
!!$    integer, intent(in)    ::  level
!!$    real(PS)  :: m_ice
!!$    real(PS),parameter :: frac=0.1_PS
!!$!!c    real(PS),parameter :: frac=0.01_PS
!!$
!!$    if( ms%mass(imr) > ms%mass(imt)*frac ) then
!!$       if( ms%mass(ima) <= ms%mass(imt)*frac ) then
!!$          if( ms%mass(imr) <= ms%mass(imc) ) then
!!$             ! rimed crystal
!!$             is%sh_type = 2
!!$          else
!!$             ! graupel
!!$             is%sh_type = 5
!!$          end if
!!$       else
!!$          if( ms%mass(imr) <= ms%mass(ima) ) then
!!$             ! rimed aggregate
!!$             is%sh_type = 4
!!$          else
!!$             ! graupel
!!$             is%sh_type = 5
!!$          end if
!!$       end if
!!$    else if( ms%mass(ima) > ms%mass(imt)*frac ) then
!!$       is%sh_type = 3             
!!$    else
!!$       ! pristine crystal
!!$       is%sh_type = 1
!!$    end if
!!$  end subroutine det_sh_type3
!!$
!!$  subroutine det_sh_type4( ms, is, level)
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! calculate 
!!$    ! 1. volume of rime component and aggregate component
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    type (Mass_Bin), intent(inout)   :: ms
!!$    type (Ice_Shape), intent(inout)    :: is
!!$    integer, intent(in)    ::  level
!!$    real(PS)  :: m_ice
!!$    real(PS),parameter :: frac=0.1_PS
!!$!!c    real(PS),parameter :: frac=0.01_PS
!!$
!!$    if( ms%mass(imr)+ms%mass(imw) > ms%mass(imt)*frac ) then
!!$       if( ms%mass(ima) <= ms%mass(imt)*frac ) then
!!$          if( ms%mass(imr)+ms%mass(imw) <= ms%mass(imc) ) then
!!$             ! rimed crystal
!!$             is%sh_type = 2
!!$          else
!!$             ! graupel
!!$             is%sh_type = 5
!!$          end if
!!$       else
!!$          if( ms%mass(imr)+ms%mass(imw) <= ms%mass(ima) ) then
!!$             ! rimed aggregate
!!$             is%sh_type = 4
!!$          else
!!$             ! graupel
!!$             is%sh_type = 5
!!$          end if
!!$       end if
!!$    else if( ms%mass(ima) > ms%mass(imt)*frac ) then
!!$       is%sh_type = 3             
!!$    else
!!$       ! pristine crystal
!!$       is%sh_type = 1
!!$    end if
!!$  end subroutine det_sh_type4

!!$  subroutine diag_sh_type_v4(level, ratio_Mp, sh_type)
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! diagnose the type of ice particle
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    integer, intent(in)    ::  level
!!$    integer, intent(inout)    ::  sh_type
!!$    ! ratio of each component mass to the total mass in a bin
!!$    ! ratio_mp(1) : rime
!!$    ! ratio_mp(2) : aggregation
!!$    ! ratio_mp(3) : ice crystal
!!$    ! ratio_Mp(7) : melt mass
!!$    real(ps), pointer, dimension(:)      :: ratio_Mp
!!$
!!$    real(PS)  :: m_ice
!!$    real(PS),parameter :: frac=0.1_PS
!!$!!c    real(PS),parameter :: frac=0.01_PS
!!$
!!$    if( ratio_Mp(1)+ratio_Mp(7) > frac ) then
!!$       if( ratio_Mp(2) <= frac ) then
!!$          if( ratio_Mp(1)+ratio_Mp(7) <= ratio_Mp(3) ) then
!!$             ! rimed crystal
!!$             sh_type = 2
!!$          else
!!$             ! graupel
!!$             sh_type = 5
!!$          end if
!!$       else
!!$          if(ratio_Mp(1)+ratio_Mp(7) <= ratio_Mp(2) ) then
!!$             ! rimed aggregate
!!$             sh_type = 4
!!$          else
!!$             ! graupel
!!$             sh_type = 5
!!$          end if
!!$       end if
!!$    else if( ratio_Mp(2) > frac ) then
!!$       sh_type = 3             
!!$    else
!!$       ! pristine crystal
!!$       sh_type = 1
!!$    end if
!!$  end subroutine diag_sh_type_v4
!!$  subroutine diag_sh_type_v5(level, ratio_Mp, sh_type, is_mod)
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! diagnose the type of ice particle
!!$    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    integer, intent(in)    ::  level
!!$    integer, intent(inout)    ::  sh_type
!!$!tmp    integer,pointer, dimension(:)    ::  is_mod
!!$    integer, dimension(*)    ::  is_mod
!!$    ! 1 : cylinder
!!$    ! 2 : spheroid
!!$
!!$    ! ratio of each component mass to the total mass in a bin
!!$    ! ratio_mp(1) : rime
!!$    ! ratio_mp(2) : aggregation
!!$    ! ratio_mp(3) : ice crystal
!!$    ! ratio_Mp(7) : melt mass
!!$!tmp    real(ps), pointer, dimension(:)      :: ratio_Mp
!!$    real(ps), dimension(*)      :: ratio_Mp
!!$
!!$    real(PS)  :: m_ice
!!$    real(PS),parameter :: frac=0.1_PS
!!$!!c    real(PS),parameter :: frac=0.01_PS
!!$
!!$    ! default is a cylinder
!!$    is_mod(1)=1
!!$    is_mod(2)=1
!!$
!!$    if( ratio_Mp(1)+ratio_Mp(7) > frac ) then
!!$       if( ratio_Mp(2) <= frac ) then
!!$          if( ratio_Mp(1)+ratio_Mp(7) <= ratio_Mp(3) ) then
!!$             ! rimed crystal
!!$             sh_type = 2
!!$          else
!!$             ! graupel
!!$             sh_type = 5
!!$             is_mod(2)=2
!!$          end if
!!$       else
!!$          if(ratio_Mp(1)+ratio_Mp(7) <= ratio_Mp(2) ) then
!!$             ! rimed aggregate
!!$             sh_type = 4
!!$          else
!!$             ! graupel
!!$             sh_type = 5
!!$             is_mod(2)=2
!!$          end if
!!$       end if
!!$    else if( ratio_Mp(2) > frac ) then
!!$       sh_type = 3             
!!$    else
!!$       ! pristine crystal
!!$       sh_type = 1
!!$    end if
!!$  end subroutine diag_sh_type_v5
  subroutine diag_sh_type_v6(level, ratio_Mp, m_frz, sh_type, is_mod)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! diagnose the type of ice particle
    ! melt mass is excluded from the diagnosis
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, intent(in)    ::  level
    integer, intent(inout)    ::  sh_type
!tmp    integer,pointer, dimension(:)    ::  is_mod
    integer, dimension(*)    ::  is_mod
    ! 1 : cylinder
    ! 2 : spheroid

    ! ratio of each component mass to the total mass in a bin
    ! ratio_mp(1) : rime
    ! ratio_mp(2) : aggregation
    ! ratio_mp(3) : ice crystal
    ! ratio_Mp(7) : melt mass
!tmp    real(ps), pointer, dimension(:)      :: ratio_Mp
    real(ps), dimension(*)      :: ratio_Mp
    ! mass of frozen nucleation
    real(PS),intent(in)  :: m_frz

    real(PS),parameter :: frac=0.1_PS
!!c    real(PS),parameter :: frac=0.01_PS
    real(PS) :: ratio_dry

    ! default is a cylinder
    is_mod(1)=1
    is_mod(2)=1

    ratio_dry=max(0.0_PS,1.0_PS-ratio_Mp(imw_m))
    if( ratio_Mp(imr_m) > frac*ratio_dry ) then
       if( ratio_Mp(ima_m) <= frac*ratio_dry ) then
          if( ratio_Mp(imr_m) <= ratio_Mp(imc_m) ) then
             ! rimed crystal
             sh_type = 2
          else
             ! graupel
             sh_type = 5
             is_mod(2)=2
          end if
       else
          if(ratio_Mp(imr_m) <= ratio_Mp(ima_m) ) then
             ! rimed aggregate
             sh_type = 4
          else
             ! graupel
             sh_type = 5
             is_mod(2)=2
          end if
       end if
    else if( ratio_Mp(ima_m) > frac*ratio_dry ) then
       sh_type = 3             
    else
       ! pristine crystal
       sh_type = 1
    end if
    if(sh_type<=2) then
       ! condition for ice sphere
       !  1: vapor mass has to become 1000 times of frozen mass. (Bacon et al. 2003, QJRMS)
       !     but it seems to overestimate the presence of frozen ice, so changed to 8.
       !  2: diameter of frozen mass is larger than 100 micron. (Korolev 2004, JAM)
! <<<< 2015/01 T. Hashino mod
       if(ratio_Mp(imf_m)*8.0>ratio_Mp(imc_m).or. &
!tmp       if(ratio_Mp(imf_m)*1000.0>ratio_Mp(imc_m).or. &
!tmp       if(ratio_Mp(imf_m)>ratio_Mp(imc_m)*0.5.or. &
          m_frz>3.83978e-6_PS &
       ) then
!tmp       if(ratio_Mp(imf_m)*10.0>ratio_Mp(imc_m)) then
! >>> 2015/01 T. Hashino mod
          is_mod(2)=2   
       endif
    endif

  end subroutine diag_sh_type_v6

!!$  subroutine cal_aclen_liq(ms)
!!$    ! this is for liquid hydrometeors
!!$    type (Mass_Bin), intent(inout)   :: ms
!!$    ! axis ratio, c-axis/a-axis
!!$    real(PS)                           :: alpha
!!$    real(PS)                           :: mean_mass
!!$    real(PS)                           :: dum
!!$
!!$       ! +++ calculate the equivalent diameter +++
!!$       ms%len = (6.0_PS*(ms%mean_mass/(PI*ms%den)))**(1.0/3.0)
!!$       
!!$       ! +++ calculate spheroidal shape 
!!$       !     based on Pruppacher and Klett (1997, p.397) +++
!!$       if( ms%len <= 280.0e-4 ) then
!!$          alpha = 1.0_PS
!!$          ms%a_len = ms%len/2.0_PS
!!$          ms%c_len = ms%len/2.0_PS
!!$       else if( 280.0e-4 < ms%len .and. ms%len <= 1.0e-1_PS ) then
!!$          ! NOTE: This formula is assumed to be applicable for this range.
!!$          ! It is supposed to be applicable to the range: from 1 mm to 9 mm.
!!$          dum = ms%len
!!$
!!$          alpha = max( 1.001668_PS - 0.098055_PS*dum &
!!$               - 2.52686_PS*dum**2 + 3.75061_PS*dum**3 &
!!$               - 1.68692_PS*dum**4, 1.0e-5_RP)
!!$
!!$          dum=ms%len/((8.0_PS*alpha)**(1.0_RP/3.0_RP))
!!$          ms%a_len = dum
!!$          ms%c_len = dum * alpha 
!!$       else if( 1.0e-1_PS < ms%len ) then
!!$          ! NOTE: This formula is assumed to be applicable for this range.
!!$          ! It is supposed to be applicable to the range: from 1 mm to 9 mm.
!!$          ! But it is not right for more than 5mm; alpha becomes negative.
!!$          dum = min(0.9_PS,ms%len)
!!$          alpha = max( 1.001668_PS - 0.098055_PS*dum &
!!$               - 2.52686_PS*dum**2 + 3.75061_PS*dum**3 &
!!$               - 1.68692_PS*dum**4, 1.0e-5_RP)
!!$          dum = ms%len/((8.0_PS*alpha)**(1.0_RP/3.0_RP))
!!$          ms%a_len = dum
!!$          ms%c_len = dum * alpha 
!!$       end if
!!$       ms%semi_a = ms%a_len
!!$       ms%semi_c = ms%c_len
!!$  end subroutine cal_aclen_liq

!!$  subroutine cal_aclen(phase, ms)
!!$    ! phase of hydrometeor
!!$    ! 1: liquid, 2: solid, 3: aerosol
!!$    integer, intent(in)  :: phase
!!$    type (Mass_Bin), intent(inout)   :: ms
!!$    ! axis ratio, c-axis/a-axis
!!$    real(PS)                           :: alpha
!!$    real(PS)                           :: mean_mass
!!$    real(PS)                           :: dum
!!$
!!$    if( phase == 1 ) then
!!$       ! initialize
!!$       ms%len=0.0_PS
!!$       ms%a_len=0.0_PS
!!$       ms%c_len=0.0_PS
!!$       ms%semi_a=0.0_PS
!!$       ms%semi_c=0.0_PS
!!$       if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return
!!$
!!$       ! +++ calculate the equivalent diameter +++
!!$       ms%len = (6.0_PS*(ms%mean_mass/(PI*ms%den)))**(1.0/3.0)
!!$       
!!$       ! +++ calculate spheroidal shape 
!!$       !     based on Pruppacher and Klett (1997, p.397) +++
!!$       if( ms%len <= 280.0e-4 ) then
!!$          alpha = 1.0_PS
!!$          ms%a_len = ms%len/2.0_PS
!!$          ms%c_len = ms%len/2.0_PS
!!$       else if( 280.0e-4 < ms%len .and. ms%len <= 1.0e-1_PS ) then
!!$          ! NOTE: This formula is assumed to be applicable for this range.
!!$          ! It is supposed to be applicable to the range: from 1 mm to 9 mm.
!!$          dum = ms%len
!!$
!!$          alpha = max( 1.001668_PS - 0.098055_PS*dum &
!!$               - 2.52686_PS*dum**2 + 3.75061_PS*dum**3 &
!!$               - 1.68692_PS*dum**4, 1.0e-5_RP)
!!$          ms%a_len = ms%len/((8.0_PS*alpha)**(1.0_RP/3.0_RP))
!!$          ms%c_len = ms%a_len * alpha 
!!$       else if( 1.0e-1_PS < ms%len ) then
!!$          ! NOTE: This formula is assumed to be applicable for this range.
!!$          ! It is supposed to be applicable to the range: from 1 mm to 9 mm.
!!$          ! But it is not right for more than 5mm; alpha becomes negative.
!!$          dum = min(0.9_PS,ms%len)
!!$          alpha = max( 1.001668_PS - 0.098055_PS*dum &
!!$               - 2.52686_PS*dum**2.0 + 3.75061_PS*dum**3.0 &
!!$               - 1.68692_PS*dum**4.0, 1.0e-5_RP)
!!$          ms%a_len = ms%len/((8.0_PS*alpha)**(1.0_RP/3.0_RP))
!!$          ms%c_len = ms%a_len * alpha 
!!$       end if
!!$       ms%semi_a = ms%a_len
!!$       ms%semi_c = ms%c_len
!!$    end if
!!$  end subroutine cal_aclen

!!$  subroutine diag_pardis_ap_v0(icat,flagp,ap_lnsig,ap_mean,ms)
!!$!!c    use com_amps
!!$    integer,intent(in) :: icat,flagp
!!$    real(PS) :: ap_lnsig(*),ap_mean(*)
!!$    type (Mass_Bin), intent(inout)   :: ms
!!$
!!$    ! initialize
!!$    ms%p=0.0_PS
!!$    ms%a_len=0.0_PS
!!$    ms%c_len=0.0_PS
!!$    ms%len=0.0_PS
!!$!!c    write(*,*) "in diag_pardis",ms%con,ms%mass(1)
!!$    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return
!!$
!!$!tmp    write(*,*) "diag_ap, ap_lnsig,ap_mean",ap_lnsig(1:2),ap_mean(1:2)
!!$!!c    stop
!!$    if(ms%dis_type==3) then
!!$       ! aerosols 
!!$       ! lognormal distribtuion
!!$       if(flagp==2.or.flagp==-1.or.flagp==-2.or.flagp==-3) then
!!$          ! fixed standard deviation
!!$          ! m_ln=log(ms%p(1))
!!$          ms%p(1)=exp((log(ms%mean_mass/coef4pi3/ms%den)-4.5_PS*ap_lnsig(icat)**2.0)/3.0_PS)
!!$          ! standard deviation of lognormal distribution, sig_ln
!!$          ms%p(2)=ap_lnsig(icat)
!!$!!c          ms%a_len=exp(-ap_lnsig(icat)**2.0)*(ms%mean_mass/ms%den/coef4pi3)**(1.0/3.0)
!!$       elseif(flagp==3.or.flagp==-4.or.flagp==-5.or.flagp==-6) then
!!$          ! fixed mean radius
!!$          ! m_ln=log(ms%p(1))
!!$          ms%p(1)=ap_mean(icat)
!!$          ! standard deviation of lognormal distribution, sig_ln
!!$          ms%p(2)=sqrt(2.0_PS/9.0_PS*(-3.0_PS*log(ap_mean(icat))+log(ms%mean_mass/ms%den/coef4pi3)))
!!$       elseif(flagp==1) then
!!$          ! fixed standard deviation
!!$          ! m_ln=log(ms%p(1))
!!$          ms%p(1)=exp((log(ms%mean_mass/coef4pi3/ms%den)-4.5_PS*ap_lnsig(icat)**2.0)/3.0_PS)
!!$          ! standard deviation of lognormal distribution, sig_ln
!!$          ms%p(2)=ap_lnsig(icat)
!!$!!c          ms%a_len=exp(-ap_lnsig(icat)**2.0)*(ms%mean_mass/ms%den/coef4pi3)**(1.0/3.0)
!!$       else
!!$          write(*,*) "This flag for aerosols not defined in diag_pardis_ap",flagp
!!$          stop
!!$       end if
!!$       ms%a_len=ms%p(1)*exp(0.5_PS*ms%p(2)**2.0)
!!$       ms%len=ms%a_len*2.0_PS
!!$       ms%c_len=ms%a_len
!!$
!!$
!!$    elseif(ms%dis_type==4) then
!!$       ! aerosols 
!!$       ! uniform distribution      
!!$       ms%len=(ms%mean_mass/ms%den/coefpi6)**(1.0/3.0)
!!$       ms%p(1)=ms%len
!!$       ms%a_len=ms%len*0.5_PS
!!$       ms%c_len=ms%a_len
!!$    else
!!$       write(*,*) "This distribution type is not defined in diag_pardis_ap",ms%dis_type
!!$       stop
!!$    end if
!!$  end subroutine diag_pardis_ap_v0

  subroutine diag_pardis_ap(icat,flagp,ap_lnsig,ap_mean,ms,em)
!!c    use com_amps
    integer,intent(in) :: icat,flagp
    integer,intent(inout) :: em
    real(PS) :: ap_lnsig(*),ap_mean(*)
    type (Mass_Bin), intent(inout)   :: ms


    if(ms%dis_type==3) then
       ! aerosols 
       ! lognormal distribtuion
       if(flagp==2.or.flagp==-1.or.flagp==-2.or.flagp==-3) then
          ! fixed standard deviation
          ! m_ln=log(ms%p(1))
          ms%p(1)=exp((log(ms%mean_mass/coef4pi3/ms%den)-4.5_PS*ap_lnsig(icat)**2.0)/3.0_PS)
          ! standard deviation of lognormal distribution, sig_ln
          ms%p(2)=ap_lnsig(icat)
!!c          ms%a_len=exp(-ap_lnsig(icat)**2.0)*(ms%mean_mass/ms%den/coef4pi3)**(1.0/3.0)
       elseif(flagp==3.or.flagp==-4.or.flagp==-5.or.flagp==-6) then
          ! fixed mean radius
          ! m_ln=log(ms%p(1))
          ms%p(1)=ap_mean(icat)
          ! standard deviation of lognormal distribution, sig_ln
          ms%p(2)=sqrt(2.0_PS/9.0_PS*(-3.0_PS*log(ap_mean(icat))+log(ms%mean_mass/ms%den/coef4pi3)))
       elseif(flagp==1) then
          ! fixed standard deviation
          ! m_ln=log(ms%p(1))
          ms%p(1)=exp((log(ms%mean_mass/coef4pi3/ms%den)-4.5_PS*ap_lnsig(icat)**2.0)/3.0_PS)
          ! standard deviation of lognormal distribution, sig_ln
          ms%p(2)=ap_lnsig(icat)
!!c          ms%a_len=exp(-ap_lnsig(icat)**2.0)*(ms%mean_mass/ms%den/coef4pi3)**(1.0/3.0)
       else
         em=1
!!c          write(*,*) "This flag for aerosols not defined in diag_pardis_ap",flagp
!!c          stop
       end if
       ms%a_len=ms%p(1)*exp(0.5_PS*ms%p(2)**2.0)
       ms%len=ms%a_len*2.0_PS
       ms%c_len=ms%a_len


    elseif(ms%dis_type==4) then
       ! aerosols 
       ! uniform distribution      
       ms%len=(ms%mean_mass/ms%den/coefpi6)**(1.0/3.0)
       ms%p(1)=ms%len
       ms%a_len=ms%len*0.5_PS
       ms%c_len=ms%a_len
    else
       em=2
!!c       write(*,*) "This distribution type is not defined in diag_pardis_ap",ms%dis_type
!!c       stop
    end if
  end subroutine diag_pardis_ap

!!$  subroutine mod_tendency( ms, n_mcom, n_vcom, process, mod_c)
!!$    type (Mass_Bin), intent(inout)   :: ms
!!$    integer, intent(in)   :: n_mcom, n_vcom, process
!!$    real(PS)  :: mod_c
!!$    integer   :: i
!!$    do i = 1, (n_mcom+1)
!!$       ms%dmassdt(i,process) = ms%dmassdt(i,process)*mod_c
!!$    end do
!!$    ms%dcondt(process) = ms%dcondt(process)*mod_c
!!$
!!$    do i = 1, n_vcom
!!$       ms%dvoldt(i, process) = ms%dvoldt(i, process)*mod_c
!!$    end do
!!$  end subroutine mod_tendency

!!$  subroutine cal_meanmass( ms,em)
!!$    type (Mass_Bin), intent(inout)   :: ms
!!$    integer,intent(inout) :: em
!!$    if(ms%con>1.0e-30_PS.and.ms%mass(1)>1.0e-30_PS) then
!!$       ms%mean_mass=ms%mass(1)/ms%con
!!$
!!$       if(ms%mean_mass==0.0_PS)then
!!$          ms%mass=0.0_PS
!!$          ms%con=0.0_PS
!!$          em=10
!!$       endif
!!$    else
!!$       ms%mass=0.0_PS
!!$       ms%con=0.0_PS
!!$       ms%mean_mass=0.0_PS
!!$    end if
!!$  end subroutine cal_meanmass

  subroutine cal_meanmass_vec( ms,em)
    type (Mass_Bin), intent(inout)   :: ms
    integer,intent(inout) :: em

    ms%mean_mass=ms%mass(1)/ms%con

    if(ms%mean_mass==0.0_PS)then
      em=10
    endif
  end subroutine cal_meanmass_vec

!!$  subroutine cal_wvt_numint(phase, ms, th_var, mode, level, &
!!$       a, m1, m2, NDIS, WVT,ishape)
!!$    use class_Ice_Shape, only: &
!!$       make_ice_shape_zero
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! calculate the bin terminal velocity
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! phase of hydrometeor
!!$    integer, intent(in)  :: phase
!!$    type (Mass_Bin), intent(in)   :: ms
!!$    type (Thermo_Var), intent(in)    :: th_var
!!$    ! number of discretization
!!$    integer,intent(in) :: NDIS
!!$
!!$    ! parameter of linear distribution
!!$    real(PS),pointer,dimension(:)     :: a
!!$    type (Ice_Shape), optional       :: ishape
!!$    ! dummy mass_bin object
!!$    type (Mass_Bin)                  :: msdum
!!$    type (Ice_Shape)                 :: isdum
!!$    real(PS), pointer, dimension(:)  :: VT
!!$    real(PS), pointer, dimension(:)  :: WVT
!!$    real(PS), pointer, dimension(:)  :: mp
!!$    real(PS) :: m1,m2, phi, dalen, dclen
!!$    real(PS) :: n1, n2, val1, val2
!!$    real(PS) :: dm
!!$    real(PS), dimension(2)  :: sum_x
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
!!$
!!$    allocate( VT(NDIS), stat = var_Status)
!!$    if(var_Status /= 0 ) stop "Memory not available for VT  in vapor_deposition"
!!$        
!!$    allocate( mp(NDIS), stat = var_Status)
!!$    if(var_Status /= 0 ) stop "Memory not available for mp in vapor_deposition"
!!$         
!!$
!!$    dm = (m2 - m1)/real(NDIS,PS_KIND)
!!$    VT=0.0_PS
!!$
!!$!    msdum = make_Mass_Bin2()
!!$    msdum = ms
!!$    msdum%con=1.0_PS
!!$
!!$    if( phase == 2 ) then
!!$!!c    is_dum = make_Ice_Shape_zero ()
!!$!tmp       isdum = make_Ice_Shape_zero () 
!!$       call make_Ice_Shape_zero (isdum) 
!!$       isdum = ishape
!!$    end if
!!$    do i = 1, NDIS
!!$       mp(i) = m1 + real(i-1,PS_KIND)*dm
!!$       msdum%mass(1)=mp(i)
!!$       msdum%mean_mass=mp(i)
!!$
!!$       if( phase == 1 ) then
!!$
!!$          call cal_aclen(phase, msdum)
!!$          call cal_terminal_vel( phase, msdum, th_var, mode, level,em)
!!$       else if( phase == 2 ) then
!!$
!!$!!c          isdum%V_ic=ishape%V_ic*(mp(i)/ms%mean_mass)
!!$
!!$          msdum%a_len = ( mp(i)/get_mci2(1,ishape%phi_ic,ms%den))**(1.0/3.0)   
!!$               
!!$          isdum%c = isdum%phi_ic*msdum%a_len 
!!$          isdum%d = isdum%psi_ic*msdum%a_len 
!!$          isdum%a = msdum%a_len - isdum%d
!!$          msdum%c_len = isdum%c
!!$
!!$          if(ishape%sh_type<=4) then
!!$             msdum%semi_a = (mp(i)/ms%den/(coef4pi3*(1.0_PS+isdum%phi_cs**2)**1.5))**(1.0/3.0)
!!$             msdum%semi_c = msdum%semi_a*isdum%phi_cs
!!$          else
!!$             if(ishape%phi_cs<1.0_PS) then
!!$                msdum%semi_a = (mp(i)/ms%den/coef4pi3)**(1.0/3.0)
!!$                msdum%semi_c = msdum%semi_a*isdum%phi_cs
!!$             else
!!$                msdum%semi_c = (mp(i)/ms%den/coef4pi3)**(1.0/3.0)
!!$                msdum%semi_a = msdum%semi_c/isdum%phi_cs
!!$             end if
!!$          end if 
!!$          call cal_terminal_vel( phase, msdum, th_var, mode, level,em, isdum)
!!$       end if
!!$
!!$       VT(i) = msdum%vtm
!!$
!!$    end do
!!$
!!$    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! numerically integrate terminal velocity
!!$    WVT = 0.0_PS
!!$    sum_x=0.0_PS
!!$
!!$    do i = 1, (NDIS-1)
!!$       n1 = max( a(1), 0.0_PS) + a(3)*( mp(i) - a(2))
!!$       n2 = max( a(1), 0.0_PS) + a(3)*( mp(i+1) - a(2))
!!$       val1 = n1*mp(i)*VT(i)
!!$       val2 = n2*mp(i+1)*VT(i+1)          
!!$       WVT(1) = WVT(1) + (val1+val2)*dm/2.0_PS
!!$
!!$       val1 = n1*VT(i)
!!$       val2 = n2*VT(i+1)          
!!$       WVT(2)=WVT(2) + (val1+val2)*dm/2.0_PS
!!$
!!$       sum_x(1)=sum_x(1)+(n1*mp(i)+n2*mp(i+1))*dm/2.0_PS
!!$       sum_x(2)=sum_x(2)+(n1+n2)*dm/2.0_PS
!!$    end do
!!$    WVT(1) = WVT(1)/sum_x(1)
!!$    WVT(2) = WVT(2)/sum_x(2)
!!$
!!$    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$!!c    call delete_Mass_Bin (msdum)
!!$    deallocate( mp,stat=var_Status)
!!$    if(var_Status /= 0 ) stop "Memory not deallocated for mp in vapor_deposition"
!!$         
!!$    nullify(mp)
!!$    deallocate( VT,stat=var_Status)
!!$    if(var_Status /= 0 ) stop "Memory not deallocated for VT in vapor_deposition"
!!$        
!!$    nullify(VT)
!!$  end subroutine cal_wvt_numint

!!$  subroutine cal_wvt_numint_v2(phase, ms, th_var, mode, level, &
!!$       a, m1, m2, NDIS, WVT,ishape)
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! calculate the bin terminal velocity
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! phase of hydrometeor
!!$    integer, intent(in)  :: phase
!!$    type (Mass_Bin), intent(in)   :: ms
!!$    type (Thermo_Var), intent(in)    :: th_var
!!$    ! number of discretization
!!$    integer,intent(in) :: NDIS
!!$
!!$    ! parameter of linear distribution
!!$    real(PS),pointer,dimension(:)     :: a
!!$    type (Ice_Shape), optional       :: ishape
!!$    real(PS), pointer, dimension(:)  :: VT
!!$    real(PS), pointer, dimension(:)  :: WVT
!!$    real(PS), pointer, dimension(:)  :: mp
!!$    real(PS),intent(in) :: m1,m2
!!$    real(PS) :: A1,B1,C1,D1,alpha,beta,gamma,N_re_0,ratio,char_len,C_DP,C_DO,X,&
!!$         A_e,A_c,P,q,OMEGA,rad,U_S,lambda_a,Y,NP,NBO,k
!!$    real(PS), parameter         :: X_0=2.8e+06
!!$    real(PS) :: n1, n2, val1, val2
!!$    real(PS) :: dm
!!$    real(PS), dimension(2)  :: sum_x
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
!!$
!!$    allocate( VT(NDIS), stat = var_Status)
!!$    if(var_Status /= 0 ) stop "Memory not available for VT in vapor_deposition"
!!$         
!!$    allocate( mp(NDIS), stat = var_Status)
!!$    if(var_Status /= 0 ) stop "Memory not available for mp in vapor_deposition" 
!!$         
!!$
!!$    dm = (m2 - m1)/real(NDIS,PS_KIND)
!!$    VT=0.0_PS
!!$
!!$    if( phase == 1 ) then
!!$       char_len=0.5_PS*(ms%mean_mass/coefpi6)**(1.0_RP/3.0_RP)
!!$       do i = 1, NDIS
!!$          mp(i) = m1 + real(i-1,PS_KIND)*dm
!!$
!!$          rad=char_len*(mp(i)/ms%mean_mass)**(1.0_RP/3.0_RP)
!!$          if(0.5e-4_PS<=rad.and.rad<10.0e-4_PS) then
!!$             ! calculate Stokes terminal velocity
!!$             U_S=rad**2*gg*(den_w-th_var%den_a)/4.5_PS/th_var%d_vis
!!$             lambda_a=6.6e-6_PS*(th_var%d_vis/1.818e-4_PS)*(1013250.0_PS/th_var%P)*(th_var%T/293.15_PS)
!!$             VT(i)=(1.0_PS+1.26_PS*lambda_a/rad)*U_S
!!$          else if(rad<535.0e-4_PS) then
!!$             X=log(32.0_PS*rad**3*(den_w-th_var%den_a)*th_var%den_a*gg/3.0_PS/th_var%d_vis**2)
!!$             Y=-0.318657e+1_PS+0.992696_PS*X-0.153193e-2_PS*X**2-0.987059e-3_PS*X**3 &
!!$                  -0.578878e-3_PS*X**4+0.855176e-4*X**5-0.327815e-5*X**6
!!$             A1=exp(Y)
!!$             ! +++ calculate the terminal velocity +++
!!$             VT(i)=A1*th_var%d_vis/(2.0_PS*rad*th_var%den_a)
!!$          else
!!$             rad=min(rad,3500.0e-4_RP)
!!$             
!!$             NBO=gg*(den_w-th_var%den_a)*rad**2/th_var%sig_wa
!!$             NP=th_var%sig_wa**3*th_var%den_a**2/th_var%d_vis**4/gg
!!$             
!!$             X=log(NBO*NP**(1.0/6.0)*16.0_PS/3.0_PS)
!!$             Y=-0.500015e+1_PS+0.523778e+1_PS*X-0.204914e+1_PS*X**2+0.475294_PS*X**3 &
!!$                  -0.542819e-1_PS*X**4+0.238449e-2_PS*X**5
!!$             
!!$             A1=NP**(1.0/6.0)*exp(Y)
!!$             ! +++ calculate the terminal velocity +++
!!$             VT(i)=A1*th_var%d_vis/(2.0_PS*rad*th_var%den_a)
!!$             
!!$          end if
!!$       end do
!!$    else if( phase == 2 ) then
!!$             
!!$       ! +++ calculate some important variables, depending on level +++
!!$       call cal_some_vel(ms,ishape,mode,alpha,OMEGA,A_e,A_c,P,q,char_len)
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
!!$       
!!$       ! +++ calculate the pressure drag C_DP +++
!!$       if( 0.3_PS < alpha .and. alpha < 1.0_PS ) then
!!$          gamma = 3.76_PS - 8.41*alpha + 9.18*alpha**2.0 - 3.53*alpha**3.0
!!$          C_DP = 0.292*k*gamma
!!$       else if( alpha <= 0.3_PS ) then 
!!$          gamma = 1.98
!!$          C_DP = 0.292*k*gamma
!!$       else if( 1.0_PS <= alpha ) then
!!$          if( 1.0_PS <= q ) then
!!$             C_DP = (1.46*q-0.46)*(0.492-0.200/sqrt(alpha))
!!$          else
!!$!!c             gamma = 3.76_PS - 8.41*alpha + 9.18*alpha**2.0 - 3.53*alpha**3.0
!!$             gamma = 1.0
!!$             C_DP = 0.292*k*gamma
!!$          end if
!!$       end if
!!$       
!!$       ! +++ calculate the Oseen-type drag coefficient, C_DO +++
!!$       C_DO = 4.5*(k**2.0)*max(alpha, 1.0_PS)
!!$       ! +++ matching theory +++
!!$       gamma = (C_DO-C_DP)/(4.0_PS*C_DP)
!!$
!!$
!!$       ! +++ calculate the weighted terminal velocity +++
!!$       ! coef. from Reynolds number
!!$       A1=6.0*k*th_var%d_vis/(C_DP*th_var%den)
!!$       ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$       ! coef. from X
!!$       if( q <= 1.0_PS ) then
!!$          D1 = 8.0_PS*gg*th_var%den/&
!!$               (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*(q**0.25))
!!$       else if( q > 1.0_PS ) then
!!$          D1 = 8.0_PS*gg*th_var%den/&
!!$               (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*q)
!!$       end if
!!$       C1=sqrt(C_DP*D1)/(6.0_PS*k)
!!$ 
!!$      
!!$
!!$       do i = 1, NDIS
!!$          mp(i) = m1 + real(i-1,PS_KIND)*dm
!!$
!!$          ! +++ calculate the Davies or Best number +++
!!$          X=D1*mp(i)
!!$
!!$          ! +++ calculate the Reynolds number +++
!!$          beta = sqrt(1.0_PS + (C_DP/(6.0_PS*k))*sqrt(X/C_DP)) - 1.0_PS
!!$          N_re_0 = 6.0*k*beta**2.0/C_DP
!!$
!!$          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$          ! coef. from matching theory
!!$          B1=(1.0_PS+2.0_PS*beta*exp(-beta*gamma)/&
!!$               ((2.0_PS+beta)*(1.0_PS+beta)))
!!$
!!$          ! +++ correction for transition to turbulent flow +++
!!$          if( 1000.0 <= N_re_0 ) then
!!$             ratio=(1.0_PS + 1.6_PS*(X/X_0)**2.0)/(1.0_PS + (X/X_0)**2.0)
!!$             A1=A1/(ratio)
!!$             C1=C1*sqrt(ratio)
!!$!!c          C1=C1*sqrt((1.0_PS + 1.6_PS*(X/X_0)**2.0)/(1.0_PS + (X/X_0)**2.0))
!!$          end if
!!$          
!!$          
!!$          VT(i)=A1*B1/(char_len*(mp(i)/ms%mean_mass)**(1.0/3.0))*(&
!!$               2.0_PS+C1*sqrt(mp(i))-2.0_PS*sqrt(1.0_PS+C1*sqrt(mp(i))))
!!$       end do
!!$!!c       if(NDIS<=5) then
!!$!!c          open( unit=21, file="check_coef_svtm.dat",POSITION='APPEND')
!!$!!c          do i = 1, NDIS
!!$!!c             mp(i) = m1 + real(i-1)*dm
!!$!!c             write(21,'(10ES15.6)') ms%mean_mass,mp(i),C1,C1*sqrt(mp(i)),1.0+C1*sqrt(mp(i)),VT(i)
!!$!!c          end do
!!$!!c          close(21)
!!$!!c       end if
!!$    end  if
!!$    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! numerically integrate terminal velocity
!!$    WVT = 0.0_PS
!!$    sum_x=0.0_PS
!!$
!!$    do i = 1, (NDIS-1)
!!$       n1 = max( a(1), 0.0_PS) + a(3)*( mp(i) - a(2))
!!$       n2 = max( a(1), 0.0_PS) + a(3)*( mp(i+1) - a(2))
!!$       WVT(1) = WVT(1) + (n1*mp(i)*VT(i)+n2*mp(i+1)*VT(i+1))*dm/2.0_PS
!!$       WVT(2)=WVT(2) + (n1*VT(i)+n2*VT(i+1))*dm/2.0_PS
!!$       sum_x(1)=sum_x(1)+(n1*mp(i)+n2*mp(i+1))*dm/2.0_PS
!!$       sum_x(2)=sum_x(2)+(n1+n2)*dm/2.0_PS
!!$    end do
!!$    WVT(1) = WVT(1)/sum_x(1)
!!$    WVT(2) = WVT(2)/sum_x(2)
!!$
!!$    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$!!c    call delete_Mass_Bin (msdum)
!!$    deallocate( mp,stat=var_Status)
!!$    if(var_Status /= 0 ) stop "Memory not deallocated for mp in vapor_deposition"
!!$        
!!$    nullify(mp)
!!$    deallocate( VT,stat=var_Status)
!!$    if(var_Status /= 0 ) stop "Memory not deallocated for VT in vapor_deposition"
!!$         
!!$    nullify(VT)
!!$  end subroutine cal_wvt_numint_v2

!!$  subroutine cal_wvt_intpol(phase, ms, th_var, mode, level, &
!!$       a, m1, m2, WVT,ishape)
!!$    use class_Ice_Shape, only: &
!!$       make_ice_shape_zero
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! calculate the bin terminal velocity by interpolation function
!!$    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! phase of hydrometeor
!!$    integer, intent(in)  :: phase
!!$    type (Mass_Bin), intent(in)   :: ms
!!$    type (Thermo_Var), intent(in)    :: th_var
!!$    ! parameter of linear distribution
!!$    real(PS),pointer,dimension(:)     :: a
!!$    type (Ice_Shape), optional       :: ishape
!!$    ! dummy mass_bin object
!!$    type (Mass_Bin)                  :: msdum
!!$    type (Ice_Shape)                 :: isdum
!!$    real(PS), pointer, dimension(:)  :: VT
!!$    real(PS), pointer, dimension(:)  :: WVT
!!$    real(PS), dimension(2)  :: mp
!!$    real(PS) :: m1,m2, phi, dalen, dclen,A1,Am,A2
!!$    real(PS) :: n1, n2, nm, val1, val2, valm
!!$    real(PS) :: dm
!!$    real(PS), dimension(2)  :: sum_x
!!$
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
!!$    ! number of discretization
!!$    integer  :: NDIS
!!$
!!$
!!$    NDIS=2
!!$    allocate( VT(NDIS), stat = var_Status)
!!$    if(var_Status /= 0 ) stop "Memory not available for VT in vapor_deposition"
!!$        
!!$
!!$    VT=0.0_PS
!!$
!!$    mp(1)=m1
!!$    mp(2)=m2
!!$
!!$!    msdum = make_Mass_Bin2()
!!$    msdum = ms
!!$    msdum%con=1.0_PS
!!$
!!$    if( phase == 2 ) then
!!$!!c    is_dum = make_Ice_Shape_zero ()
!!$!tmp       isdum = make_Ice_Shape_zero () 
!!$       call make_Ice_Shape_zero (isdum) 
!!$       isdum = ishape
!!$    end if
!!$    do i=1,NDIS
!!$       msdum%mass(1)=mp(i)
!!$       msdum%mean_mass=mp(i)
!!$
!!$       if( phase == 1 ) then
!!$
!!$          call cal_aclen(phase, msdum)
!!$          call cal_terminal_vel( phase, msdum, th_var, mode, level,em)
!!$       else if( phase == 2 ) then
!!$
!!$          if(ishape%habit<=3) then
!!$             msdum%a_len=(mp(i)/ms%mean_mass)**(1.0/3.0)*ms%a_len
!!$             msdum%c_len=isdum%phi_ic*msdum%a_len 
!!$             isdum%c=msdum%c_len
!!$             isdum%d=isdum%psi_ic*msdum%a_len 
!!$             isdum%a=msdum%a_len-isdum%d
!!$             isdum%r=isdum%gam_ic*msdum%a_len
!!$             isdum%e=isdum%eta_ic*msdum%a_len
!!$          elseif(ishape%habit==4) then
!!$             ! rosette
!!$             isdum%r=(mp(i)/ms%mean_mass)**(1.0/3.0)*ishape%r
!!$             msdum%a_len=isdum%r/isdum%gam_ic
!!$             msdum%c_len=isdum%phi_ic*msdum%a_len 
!!$             isdum%c=msdum%c_len
!!$             isdum%d=isdum%psi_ic*msdum%a_len 
!!$             isdum%a=msdum%a_len-isdum%d
!!$             isdum%e=isdum%eta_ic*msdum%a_len
!!$          elseif(ishape%habit>=5) then
!!$             ! irregular crystal
!!$             isdum%e=(mp(i)/ms%mean_mass)**(1.0/3.0)*ishape%e
!!$             msdum%a_len=isdum%e/isdum%eta_ic
!!$             msdum%c_len=isdum%phi_ic*msdum%a_len 
!!$             isdum%c=msdum%c_len
!!$             isdum%d=isdum%psi_ic*msdum%a_len 
!!$             isdum%a=msdum%a_len-isdum%d
!!$             isdum%r=isdum%gam_ic*msdum%a_len
!!$          end if
!!$          msdum%semi_a=(mp(i)/ms%mean_mass)**(1.0/3.0)*ms%semi_a
!!$          msdum%semi_c=msdum%semi_a*isdum%phi_ic
!!$
!!$          call cal_terminal_vel( phase, msdum, th_var, mode, level,em, isdum)
!!$       end if
!!$
!!$       VT(i) = msdum%vtm
!!$
!!$    end do
!!$
!!$    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! fit 2nd order polynomial and integrate over the region.
!!$    WVT = 0.0_PS
!!$    sum_x=0.0_PS
!!$
!!$    n1=max( a(1), 0.0_PS) + a(3)*( mp(1) - a(2))
!!$    n2=max( a(1), 0.0_PS) + a(3)*( mp(2) - a(2))
!!$    nm=max( a(1), 0.0_PS) + a(3)*( ms%mean_mass - a(2))
!!$
!!$    val1 = n1*mp(1)*VT(1)
!!$    val2 = n2*mp(2)*VT(2)          
!!$    valm = nm*ms%mean_mass*ms%vtm
!!$
!!$!!c    A1=val1/(mp(1)-ms%mean_mass)/(mp(1)-mp(2))
!!$!!c    Am=valm/(ms%mean_mass-mp(1))/(ms%mean_mass-mp(2))
!!$!!c    A2=val2/(mp(2)-ms%mean_mass)/(mp(2)-mp(1))
!!$!!c
!!$!!c    WVT(1)=((mp(2)**3.0-mp(1)**3.0)*(A1+Am+A2)/3.0_PS&
!!$!!c         -(mp(2)**2.0-mp(1)**2.0)*&
!!$!!c         (A1*(ms%mean_mass+mp(2))+Am*(mp(1)+mp(2))+A2*(mp(1)+ms%mean_mass))/2.0_PS&
!!$!!c         +(mp(2)-mp(1))*(A1*ms%mean_mass*mp(2)+Am*mp(1)*mp(2)+A2*mp(1)*ms%mean_mass))/ms%mass(1)
!!$    WVT(1)=(val1/(ms%mean_mass-mp(1))*(&
!!$            (mp(2)**2.0+mp(2)*mp(1)+mp(1)**2.0)/3.0-(ms%mean_mass+mp(2))*(mp(2)+mp(1))/2.0+ms%mean_mass*mp(2))+&
!!$            valm*(mp(2)-mp(1))/(ms%mean_mass-mp(1))/(ms%mean_mass-mp(2))*&
!!$              ((mp(2)**2.0+mp(2)*mp(1)+mp(1)**2.0)/3.0-(mp(2)+mp(1))*(mp(2)+mp(1))/2.0+mp(1)*mp(2))+&
!!$            val2/(mp(2)-ms%mean_mass)*&
!!$              ((mp(2)**2.0+mp(2)*mp(1)+mp(1)**2.0)/3.0-&
!!$                  (ms%mean_mass+mp(1))*(mp(2)+mp(1))/2.0+ms%mean_mass*mp(1)))/ms%mass(1)
!!$
!!$
!!$
!!$
!!$    val1 = n1*VT(1)
!!$    val2 = n2*VT(2)          
!!$    valm = nm*ms%vtm
!!$
!!$!!c    A1=val1/(mp(1)-ms%mean_mass)/(mp(1)-mp(2))
!!$!!c    Am=valm/(ms%mean_mass-mp(1))/(ms%mean_mass-mp(2))
!!$!!c    A2=val2/(mp(2)-ms%mean_mass)/(mp(2)-mp(1))
!!$!!c
!!$!!c    WVT(2)=((mp(2)**3.0-mp(1)**3.0)*(A1+Am+A2)/3.0_PS&
!!$!!c         -(mp(2)**2.0-mp(1)**2.0)*&
!!$!!c         (A1*(ms%mean_mass+mp(2))+Am*(mp(1)+mp(2))+A2*(mp(1)+ms%mean_mass))/2.0_PS&
!!$!!c         +(mp(2)-mp(1))*(A1*ms%mean_mass*mp(2)+Am*mp(1)*mp(2)+A2*mp(1)*ms%mean_mass))/ms%con
!!$
!!$    WVT(2)=(val1/(ms%mean_mass-mp(1))*((mp(2)**2.0+mp(2)*mp(1)+mp(1)**2.0)/3.0-&
!!$               (ms%mean_mass+mp(2))*(mp(2)+mp(1))/2.0+ms%mean_mass*mp(2))+&
!!$            valm*(mp(2)-mp(1))/(ms%mean_mass-mp(1))/(ms%mean_mass-mp(2))*&
!!$               ((mp(2)**2.0+mp(2)*mp(1)+mp(1)**2.0)/3.0-(mp(2)+mp(1))*(mp(2)+mp(1))/2.0+mp(1)*mp(2))+&
!!$            val2/(mp(2)-ms%mean_mass)*((mp(2)**2.0+mp(2)*mp(1)+mp(1)**2.0)/3.0-&
!!$               (ms%mean_mass+mp(1))*(mp(2)+mp(1))/2.0+ms%mean_mass*mp(1)))/ms%con
!!$
!!$
!!$
!!$    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$!!c    call delete_Mass_Bin (msdum)
!!$    deallocate( VT,stat=var_Status)
!!$    if(var_Status /= 0 ) stop "Memory not deallocated for VT in vapor_deposition"
!!$        
!!$    nullify(VT)
!!$  end subroutine cal_wvt_intpol

!!$  function get_mci1(sh_type,phi,den) result(mci)
!!$    real(PS),intent(in) :: phi,den
!!$    integer, intent(in) :: sh_type
!!$    ! coefficient of mass-dimension relation, m = mci * D^qq
!!$    real(PS) :: mci
!!$    select case(sh_type)
!!$    case (1)
!!$       ! pristine crystal
!!$       ! for cylinder volume
!!$       mci = phi*PI*den/4.0_PS
!!$    case (2)
!!$       ! rimed pristine crystal
!!$       ! for cylinder volume
!!$       mci = phi*PI*den/4.0_PS
!!$    case default
!!$       ! for spheroid
!!$       mci = (PI/6.0_PS)*phi*den
!!$    end select
!!$  end function get_mci1
!!$  function get_mci2(sh_type,phi,den) result(mci)
!!$    real(PS),intent(in) :: phi,den
!!$    integer, intent(in) :: sh_type
!!$    ! coefficient of mass-dimension relation, m = mci * r ^qq
!!$    real(PS) :: mci
!!$    mci = coef4pi3*phi*den
!!$  end function get_mci2

!!$  subroutine cal_cseff(ms, th_var,y)
!!$    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    ! Calculate curvature and solute contribution
!!$    ! using (13-27b) Pruppacher and Klett (1996).
!!$    !
!!$    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$    use com_amps
!!$    implicit none
!!$    type (Mass_Bin), intent(inout) :: ms
!!$    type (Thermo_Var), intent(in) :: th_var
!!$    real(PS) :: mean_mass_s,y
!!$
!!$    mean_mass_s=ms%mass(3)/ms%con
!!$
!!$!!c    if(mean_mass_s/ms%mean_mass>0.99999) then
!!$!!c       y=2.0_PS*th_var%sig_wa/(R_v*th_var%T*den_w*ms%len*0.5_PS)
!!$    y=2.0_PS*th_var%sig_wa/(R_v*th_var%T*den_w*ms%len*0.5_PS)-&
!!$         nu_aps(1)*phi_aps(1)*mean_mass_s*M_W/M_aps(1)/(ms%mean_mass-mean_mass_s)
!!$
!!$    ! for now the solute effect is neglected because  I do not have
!!$    ! the skill to predict the water right after activation.         
!!$
!!$  end subroutine cal_cseff

!!$  subroutine cal_cond_mlt_grap(ms,th_var,a_tot,a_ice,T_s,dmdt,ick)
!!$    use class_Thermo_Var, only: &
!!$       get_sat_vapor_pres
!!$    implicit none
!!$    type (Mass_Bin), intent(in) :: ms
!!$    type (Thermo_Var), intent(in) :: th_var
!!$    real(PS),intent(in) :: a_tot,a_ice
!!$    real(PS) :: T_s,T_sb,dmdt
!!$    integer :: iter,ick
!!$
!!$!    T_s=th_var%T
!!$    T_s=T_0
!!$    T_sb=T_s
!!$    iter=0
!!$    do 
!!$       iter=iter+1
!!$       T_s=T_0+(th_var%k_a*(th_var%T-T_s)*ms%fh+&
!!$            th_var%D_v*L_e*M_w/R_u*(th_var%e/th_var%T-get_sat_vapor_pres(1, T_s )/T_s)*ms%fv)*(a_tot-a_ice)/a_ice/k_w
!!$  
!!$       if(abs(T_s-th_var%T)>50.0.or.iter==20) then
!!$          ick=1
!!$          dmdt=0.0
!!$          exit          
!!$       elseif(abs(T_s-T_sb)<0.001) then
!!$          ick=0
!!$          dmdt=4.0_PS*PI*a_tot*a_ice/(a_tot-a_ice)*k_w*(T_s-T_0)/L_f
!!$          exit
!!$       end if
!!$       T_sb=T_s
!!$    end do
!!$  end subroutine cal_cond_mlt_grap

!!$  subroutine cal_coef_Ts(ms,th_var,phase,TS_A1,TS_B11,TS_B12,phase2)
!!$    use class_Thermo_Var, only: &
!!$       get_mod_thermal_cond
!!$    type (Mass_Bin), intent(inout) :: ms
!!$    type (thermo_var), intent(in) :: th_var
!!$    real(PS),intent(inout) :: TS_A1,TS_B11,TS_B12
!!$    integer,intent(in) :: phase
!!$    integer,intent(inout) :: phase2
!!$    
!!$    real(PS) :: r_e
!!$
!!$    ! modified conductivity
!!$    real(PS) :: mk_a
!!$
!!$
!!$    TS_A1=0.0_PS
!!$    TS_B11=0.0_PS
!!$    TS_B12=0.0_PS
!!$    phase2=phase
!!$    if(ms%con<1.0e-30_PS.or.ms%mass(1)<1.0e-30_PS) return
!!$
!!$    ! calculate modified thermal conductivity
!!$    mk_a = get_mod_thermal_cond(sqrt(ms%semi_a**2+ms%semi_c**2),&
!!$                                th_var%K_a,th_var%T,th_var%den)
!!$
!!$    if(T_0>th_var%T) then
!!$!!c       r_e=th_var%e_sat_n(1)/th_var%e_sat_n(2)
!!$       if(ms%inmlt==0) then
!!$          phase2=phase
!!$               
!!$          TS_A1=L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          TS_B11=L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$               th_var%e_sat_n(phase2)/th_var%T_n/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          TS_B12=(L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$               th_var%e_sat_n(phase2)/th_var%T_n+&
!!$               (L_f+C_w*th_var%T_n)*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con+&
!!$               4.0_PS*PI*ms%CAP*mk_a*ms%fh*th_var%T_n)/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          
!!$       else
!!$          !
!!$          ! This formulation ignores conduction between the surface of the hydrometeor
!!$          ! and the surface of the ice core.
!!$          !
!!$          phase2=1
!!$          
!!$          TS_A1=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          TS_B11=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$               th_var%e_sat_n(phase2)/th_var%T_n/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          TS_B12=(L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$               th_var%e_sat_n(phase2)/th_var%T_n+&
!!$               (L_f+C_w*th_var%T_n)*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con+&
!!$               4.0_PS*PI*ms%CAP*mk_a*ms%fh*th_var%T_n)/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$       end if
!!$       
!!$    else
!!$       ! 
!!$       ! assume that the riming process produce surface water 
!!$       ! due to the liquid layer over the hydrometeors:
!!$       ! no fusion heat release by riming process
!!$       !
!!$       ! The liquid layer is turbulent, so the heat received at the
!!$       ! surface is balanced with the one of ice core.
!!$       !
!!$       ! Therefore, the surface temperature is T_0.
!!$       phase2=1
!!$!!c
!!$!!c             A1=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
!!$!!c                  (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$!!c
!!$!!c             B1=(L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$!!c                  th_var%e_sat_n(phase2)/th_var%T_n*&
!!$!!c                  (th_var%s_v_n(phase2)+1.0_PS)+&
!!$!!c                  (C_w*th_var%T_n)*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con+&
!!$!!c                  4.0_PS*PI*ms%CAP*mk_a*ms%fh*th_var%T_n)/&
!!$!!c                  (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$    end if
!!$  end subroutine cal_coef_Ts

!!$  subroutine cal_coef_Ts2(ms,th_var,phase,TS_A1,TS_B11_0,TS_B12_0,TS_B12_1,TS_B12_2,phase2)
!!$    use class_Thermo_Var, only: &
!!$       get_mod_thermal_cond
!!$    type (Mass_Bin), intent(inout) :: ms
!!$    type (thermo_var), intent(in) :: th_var
!!$    real(PS),intent(inout) :: TS_A1,TS_B11_0,TS_B12_0,TS_B12_1,TS_B12_2
!!$    integer,intent(in) :: phase
!!$    integer,intent(inout) :: phase2
!!$    
!!$    real(PS) :: r_e
!!$
!!$    ! modified conductivity
!!$    real(PS) :: mk_a
!!$
!!$    TS_A1=0.0_PS
!!$    TS_B11_0=0.0_PS
!!$    TS_B12_0=0.0_PS
!!$    TS_B12_1=0.0_PS
!!$    TS_B12_2=0.0_PS
!!$    phase2=phase
!!$    if(ms%con<1.0e-30_PS.or.ms%mass(1)<1.0e-30_PS) return
!!$
!!$    ! calculate modified thermal conductivity
!!$    mk_a = get_mod_thermal_cond(sqrt(ms%semi_a**2+ms%semi_c**2),&
!!$                                th_var%K_a,th_var%T,th_var%den)
!!$
!!$    if(T_0>th_var%T) then
!!$       !!c       r_e=th_var%e_sat(1)/th_var%e_sat(2)
!!$       if(ms%inmlt==0) then
!!$          phase2=phase
!!$          
!!$          TS_A1=L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          !!c          TS_B11=L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$          !!c               th_var%e_sat(phase2)/Ta/&
!!$          !!c               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          ! TS_B11=TS_B11_0*e_sat(Ta)/Ta
!!$          TS_B11_0=L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          !!c          TS_B12=(L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$          !!c               th_var%e_sat(phase2)/Ta+&
!!$          !!c               (L_f+C_w*Ta)*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con+&
!!$          !!c               4.0_PS*PI*ms%CAP*mk_a*ms%fh*Ta)/&
!!$          !!c               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          
!!$          ! TS_B12=TS_B12_0+TS_B12_1*e_sat(Ta)/Ta+TS_B12_2*Ta
!!$          !
!!$          TS_B12_0=(L_f*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con)/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          !TS_B12_1*e_sat(Ta)/Ta
!!$          TS_B12_1=L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          !TS_B12_2*Ta
!!$          TS_B12_2=(C_w*ms%Ldmassdt(2)+4.0_PS*PI*ms%CAP*mk_a*ms%fh)/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$       
!!$       
!!$       
!!$       else
!!$          !
!!$          ! This formulation ignores conduction between the surface of the hydrometeor
!!$          ! and the surface of the ice core.
!!$          !
!!$          phase2=1
!!$          
!!$          TS_A1=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          
!!$          !!c          TS_B11=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$          !!c               th_var%e_sat(phase2)/Ta/&
!!$          !!c               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          ! TS_B11=TS_B11_0*e_sat(Ta)/Ta
!!$          TS_B11_0=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          !!c          TS_B12=(L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$          !!c               th_var%e_sat(phase2)/Ta+&
!!$          !!c               (L_f+C_w*Ta)*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con+&
!!$          !!c               4.0_PS*PI*ms%CAP*mk_a*ms%fh*Ta)/&
!!$          !!c               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$          ! TS_B12=TS_B12_0+TS_B12_1*e_sat(Ta)/Ta+TS_B12_2*Ta
!!$          TS_B12_0=(L_f*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con)/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          ! 
!!$          TS_B12_1=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          !
!!$          TS_B12_2=(C_w*ms%Ldmassdt(2)+4.0_PS*PI*ms%CAP*mk_a*ms%fh)/&
!!$               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$          
!!$       end if
!!$    else
!!$       ! 
!!$       ! assume that the riming process produce surface water 
!!$       ! due to the liquid layer over the hydrometeors:
!!$       ! no fusion heat release by riming process
!!$       !
!!$       ! The liquid layer is turbulent, so the heat received at the
!!$       ! surface is balanced with the one of ice core.
!!$       !
!!$       ! Therefore, the surface temperature is T_0.
!!$       phase2=1
!!$!!c
!!$!!c             A1=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
!!$!!c                  (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$!!c
!!$!!c             B1=(L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!$!!c                  th_var%e_sat(phase2)/th_var%T*&
!!$!!c                  (th_var%s_v_n(phase2)+1.0_PS)+&
!!$!!c                  (C_w*th_var%T)*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con+&
!!$!!c                  4.0_PS*PI*ms%CAP*mk_a*ms%fh*th_var%T)/&
!!$!!c                  (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!$    end if
!!$
!!$  end subroutine cal_coef_Ts2

  subroutine cal_coef_Ts3(ms,th_var,phase,TS_A1,TS_B11,TS_B12,TS_B13,TS_D1 &
                         ,phase2,Tmax,mk_a)
    type (Mass_Bin), intent(inout) :: ms
    type (thermo_var), intent(in) :: th_var
    real(8),intent(inout) :: TS_A1,TS_B11,TS_B12,TS_B13,TS_D1
    integer,intent(in) :: phase
    ! modified conductivity
    real(PS),intent(in) :: mk_a
    integer,intent(inout) :: phase2
    real(PS),intent(inout) :: Tmax
    
    real(PS) :: L_x


    if(MS%inmlt==0) then
      phase2=phase
      L_x=L_s
      Tmax=T_0
    else
      !
      ! This formulation ignores conduction between the surface of the hydrometeor
      ! and the surface of the ice core.
      !
      phase2=1
      L_x=L_e
      Tmax=T_0+40.0
    endif
    if(ms%con>=1.0e-30_PS.and.ms%mass(1)>=1.0e-30_PS) then


               
      TS_A1=4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2)
                   
      TS_B11=L_x*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn

      TS_B12=L_f*ms%Ldmassdt(2) - &
            !
            L_f*ms%dmassdt(1,10)/ms%con 

      TS_B13=C_w*ms%Ldmassdt(2)+4.0_PS*PI*ms%CAP*mk_a*ms%fh


      TS_D1=L_x*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn
    else
      TS_A1=0.0_PS
      TS_B11=0.0_PS
      TS_B12=0.0_PS
      TS_B13=0.0_PS
      TS_D1=0.0_PS
    endif

  end subroutine cal_coef_Ts3

!!$  subroutine update_growth_mode(level,mes_rc,th_var)
!!$!tmp    use mod_amps_utility, only: get_growth_mode
!!$    use mod_amps_utility, only: get_growth_mode_max
!!$    integer,intent(in) :: level,mes_rc
!!$    type (Thermo_Var), intent(inout)  :: th_var
!!$
!!$    if(level==3.or.level==5.or.level==7) then
!!$       if((mes_rc==2.or.mes_rc==4).and.th_var%T<253.16_PS) then
!!$          ! if liquid co exist
!!$!tmp          th_var%nuc_gmode=get_growth_mode(th_var%T,th_var%e_sat(1)/th_var%e_sat(2)-1.0_PS)
!!$          th_var%nuc_gmode=get_growth_mode_max(th_var%T,th_var%e_sat(1)/th_var%e_sat(2)-1.0_PS)
!!$
!!$!!c          th_var%nuc_gmode=1
!!$
!!$       end if
!!$    end if
!!$  end subroutine update_growth_mode

!!$  function get_dmdt0(imode,tmp) result(out)
!!$    integer,intent(in) :: imode
!!$    real(PS),intent(in) :: tmp
!!$    real(PS) :: tmpc,out
!!$    real(PS),dimension(6) :: tmpr,irr_dvdt,pla_dvdt,col_dvdt,ros_dvdt,s3_dvdt
!!$    integer :: i1,i2
!!$    tmpr = (/-20.0, -30.0, -40.0, -50.0, -60.0, -70.0/)
!!$    irr_dvdt = (/13000.0, 4600.0, 0.0, 0.0, 0.0, 0.0/)
!!$    pla_dvdt = (/19000.0, 2600.0, 1100.0, 430.0, 340.0, 22.0/)
!!$    col_dvdt = (/0.0, 2300.0, 2200.0, 1000.0, 1400.0, 75.0/)
!!$    ros_dvdt = (/0.0, 0.0, 0.0, 2800.0, 2600.0, 85.0/)
!!$    s3_dvdt = (/22000.0, 5300.0, 850.0, 0.0, 0.0, 0.0/)
!!$    tmpc=tmp-273.16
!!$    if(imode==1) then
!!$       ! case of irregular polycrystals         
!!$       out=irr_dvdt(1)+(tmpc-tmpr(1))*(irr_dvdt(2)-irr_dvdt(1))/(tmpr(2)-tmpr(1))
!!$    elseif(imode==5) then
!!$       ! case of planar polycrystals         
!!$       i1=max(min(int((-20.0_RP-tmpc)/10.0_RP)+1,2),1)
!!$       i2=i1+1
!!$       out=s3_dvdt(i1)+(tmpc-tmpr(i1))*(s3_dvdt(i2)-s3_dvdt(i1))/(tmpr(i2)-tmpr(i1))
!!$    else if(imode==2) then
!!$       ! case of plates
!!$       i1=max(min(int((-20.0_RP-tmpc)/10.0_RP)+1,5),1)
!!$       i2=i1+1
!!$       out=pla_dvdt(i1)+(tmpc-tmpr(i1))*(pla_dvdt(i2)-pla_dvdt(i1))/(tmpr(i2)-tmpr(i1))
!!$    else if(imode==3) then
!!$       ! case of columns
!!$       i1=max(min(int((-20.0_RP-tmpc)/10.0_RP)+1,5),2)
!!$       i2=i1+1
!!$       out=col_dvdt(i1)+(tmpc-tmpr(i1))*(col_dvdt(i2)-col_dvdt(i1))/(tmpr(i2)-tmpr(i1))
!!$    else if(imode==4) then
!!$       ! case of rosetta bullets
!!$       i1=max(min(int((-20.0_RP-tmpc)/10.0_RP)+1,5),4)
!!$       i2=i1+1
!!$       out=ros_dvdt(i1)+(tmpc-tmpr(i1))*(ros_dvdt(i2)-ros_dvdt(i1))/(tmpr(i2)-tmpr(i1))
!!$    end if
!!$    out=min(max(out,22.0_RP),22000.0_RP)*den_i*1.0e-12_RP
!!$  end function get_dmdt0

  function get_osm(ica,molality) result(out)
    use scale_prc, only: &
       PRC_abort
    use com_amps, only: &
       APSNAME
    use mod_amps_utility, only: &
       osm_ammsul, &
       osm_sodchl
    implicit none
    integer,intent(in) :: ica
    real(PS),intent(in) :: molality
    real(PS) :: out

    select case(trim(APSNAME(ica)))
    case('NH42SO4')
       ! ammonium sulfate
       out=osm_ammsul(molality)
    case('NACL2')
       ! sodium chloride
       out=osm_sodchl(molality)
    case('NH4HSO4')
       ! ammonium bi-sulfate
       out=osm_ammsul(molality)
    case default
       LOG_ERROR("get_osm",*) "get_osm: No such chemicals",trim(APSNAME(ica))
       call PRC_abort
    end select
  end function get_osm

  function get_critrad_itr(AA,BB,n_aps,r_n) result(out)
  use acc_amps
  implicit none
  ! coefficients (micro meter and dimensionless)
  real(PS),intent(in) :: AA,BB
  ! moles of soluble mass
  real(PS),intent(in) :: n_aps
  ! radius of dry nuclei (micro meter)
  real(PS),intent(in) :: r_n
  real(PS) :: out

  !  Using Brent's method, find the root of a function func known to lie between x1 and x2.
  !  The root, returned as zbrent, will be refined until its accuracy is tol.
  REAL(PS) :: tol=1.0e-6
  integer :: ITMAX
  real(PS) :: eps
!  PARAMETER (ITMAX=50,EPS=epsilon(out))
  PARAMETER (ITMAX=50,EPS=3.e-8)

  !  Parameters: Maximum allowed number of iterations, and machine floating-point precision.
  INTEGER :: iter
  REAL(PS) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

  a=r_n*1.01
  b=r_n*1000.0


  iter=1
2221 continue
  fa=func(a)
2222 continue
  fb=func(b)

  if(fa.gt.0.0_PS.and.fb.gt.0.0_PS) then
     iter=iter+1
     if(iter>50) then
        if(debug) write(*,*) "iter at Ts large 1",fa,fb,a,b
     end if
!!$       a=a-50.0_PS
     a=a/2.0_PS
!try     a=max(r_n*1.00001,a*0.9_PS)
!!$       stop
     goto 2221
  elseif(fa.lt.0.0_PS.and.fb.lt.0.0_PS) then
!!$       if((fa.gt.0.0_PS.and.fb.gt.0.0_PS).or.(fa.lt.0.0_PS.and.fb.lt.0.0_PS)) then
     iter=iter+1
     if(iter>50) then
        if(debug) write(*,*) "iter at Ts large 1",fa,fb,a,b
     end if
!!$       stop  
!!$                b=b*100.0
     b=b*2.0_PS
!!$       b=b+50.0_PS
     goto 2222
  end if

  c=b
  fc=fb
  do iter=1,ITMAX
     if((fb.gt.0.0_PS.and.fc.gt.0.0_PS).or.(fb.lt.0.0_PS.and.fc.lt.0.0_PS))then
        c=a ! Rename a, b, c and adjust bounding interval d.
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
     tol1=2.0_PS*EPS*abs(b)+0.5_PS*tol ! Convergence check.
     xm=0.5_PS*(c-b)
     if(abs(xm).le.tol1 .or. fb.eq.0.)then
        out=b
        goto 1111
     endif
     if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
        s=fb/fa ! Attempt inverse quadratic interpolation.
        if(a.eq.c) then
           p=2.0_PS*xm*s
           q=1.0_PS-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.0_PS))
           q=(q-1.0_PS)*(r-1.0_PS)*(s-1.0_PS)
        endif
        if(p.gt.0.) q=-q ! Check whether in bounds.
        p=abs(p)
        if(2.0_PS*p .lt. min(3.0_PS*xm*q-abs(tol1*q),abs(e*q))) then
           e=d ! Accept interpolation.
           d=p/q
        else
           d=xm ! Interpolation failed, use bisection.
           e=d
        endif
     else ! Bounds decreasing too slowly, use bisection.
        d=xm
        e=d
     endif
     a=b ! Move last best guess to a.
     fa=fb
     if(abs(d) .gt. tol1) then ! Evaluate new trial root.
        b=b+d
     else
        b=b+sign(tol1,xm)
     endif
     fb=func(b)
  enddo
  if(debug) write(*,'("brent exceeding maximum iterations. n,iter",2I5,12ES15.6)') &
       iter,a,b,fa,fb
!!$             write(*,'("W,a_c,N_act,M_act,T,r_e",10ES15.6)') ag%TV(n)%W, a_c, N_act,M_act, ag%TV(n)%T,r_e

  out=r_n*1.1_PS


1111 continue

contains
  function func(x) result(out)
    real(PS),intent(in) :: x
    real(PS) :: out
    ! n_aps has to be in mol um^3/cm^3
    out=x**3-sqrt(3.0_PS*BB*get_osm(1,1.0e+3*n_aps/(den_w*coef4pi3*(x**3-r_n**3) )) *&
         r_n**3.0/AA)*x**2-r_n**3
  end function func

end function get_critrad_itr

function get_critrad_anal(AA,sb,beta,r_n) result(r_cr)
  ! analytical expresion (6.3.7) by Khvorostyanov and Curry (2014)
  use acc_amps
  implicit none
  ! coefficients (cm and dimensionless)
  real(PS),intent(in) :: AA,sb
  real(PS),intent(in) :: beta
  ! radius of dry nuclei (cm)
  real(PS),intent(in) :: r_n
  real(PS) :: r_cr  ! cm
  real(PS) :: Zd
  real(8) :: Zdc,Pp,Pm

  Zd=r_n**beta*sqrt(sb/(3.0_PS*AA))
  Zdc=Zd
  if(Zd>1.0e-2) then
    Pp=(Zdc*Zdc*Zdc+sqrt(Zdc*Zdc*Zdc+0.25_PS)+0.5_PS)**(1.0/3.0)
    Pm=(Zdc*Zdc*Zdc-sqrt(Zdc*Zdc*Zdc+0.25_PS)+0.5_PS)**(1.0/3.0)
    r_cr=r_n*(Zdc+Pp+Pm)
!tmp    if(Zdc*Zdc*Zdc-sqrt(Zdc*Zdc*Zdc+0.25_PS)+0.5_PS<0.0_PS) then
!tmp      write(*,*) "ck r_cr",i,n,r_cr,r_n,Zdc,Pp,Pm,r_n,beta,sb,AA, &
!tmp               g%MS(i,n)%a_len,den_ap,g%MS(i,n)%eps_map,&
!tmp               g%MS(i,n)%den,g%MS(i,n)%den_ai,g%MS(i,n)%den_as
!tmp    endif
  else
    r_cr=r_n*(1.0d+0+Zdc+2.0d+0*Zdc*Zdc*Zdc/3.0d+0)
  endif

end function get_critrad_anal

function get_hazerad_itr(AA,BB,SS,n_aps,r_n,rd_c,iswitch) result(out)
  !
  ! this one considers variable osmotic coefficient
  !
  use acc_amps
  implicit none
  ! coefficients (micro meter and dimensionless)
  real(PS),intent(in) :: AA,BB
  ! saturation ratio
  real(PS),intent(in) :: SS
  ! moles of soluble mass
  real(PS),intent(in) :: n_aps
  ! radius of dry nuclei (micro meter), critical radius of droplet
  real(PS),intent(in) :: r_n, rd_c
  integer,intent(in) :: iswitch
  real(PS) :: out


  !  Using Brent's method, find the root of a function func known to lie between x1 and x2.
  !  The root, returned as zbrent, will be refined until its accuracy is tol.
  REAL(PS) :: tol=1.0e-6
  integer :: ITMAX
  real(PS) :: eps
!  PARAMETER (ITMAX=50,EPS=epsilon(rd_c))
  PARAMETER (ITMAX=50,EPS=3.e-8)
    
  !  Parameters: Maximum allowed number of iterations, and machine floating-point precision.
  INTEGER :: iter
  REAL(PS) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
  
  if(iswitch==0) then
     a=rd_c
     b=r_n*1.01
  else
     a=rd_c*50.0
     b=rd_c*1.01
  end if

  iter=1
2221 continue
  fa=func(a)
2222 continue
  fb=func(b)
  
  if(fa.gt.0.0_PS.and.fb.gt.0.0_PS) then
     iter=iter+1
     if(iter>50) then
        if(debug) write(*,'("iter at get_hazerad large 1",I5,10ES15.6)') iswitch,fa,fb,a,b,rd_c,r_n
     end if
     a=a*2.0_PS
!!$       stop
     goto 2221
  elseif(fa.lt.0.0_PS.and.fb.lt.0.0_PS) then
!!$       if((fa.gt.0.0_PS.and.fb.gt.0.0_PS).or.(fa.lt.0.0_PS.and.fb.lt.0.0_PS)) then
     iter=iter+1
     if(iter>50) then
        if(debug) write(*,'("iter at get_hazerad large 2",I5,10ES15.6)') iswitch,fa,fb,a,b,rd_c,r_n
     end if
!!$       stop  
!!$                b=b*100.0
     b=b/2.0_PS
!!$       b=b+50.0_PS
     goto 2222
  end if
  
  c=b
  fc=fb
  do iter=1,ITMAX
     if((fb.gt.0.0_PS.and.fc.gt.0.0_PS).or.(fb.lt.0.0_PS.and.fc.lt.0.0_PS))then
        c=a ! Rename a, b, c and adjust bounding interval d.
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
     tol1=2.0_PS*EPS*abs(b)+0.5_PS*tol ! Convergence check.
     xm=0.5_PS*(c-b)
     if(abs(xm).le.tol1 .or. fb.eq.0.)then
        out=b
        goto 1111
     endif
     if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
        s=fb/fa ! Attempt inverse quadratic interpolation.
        if(a.eq.c) then
           p=2.0_PS*xm*s
           q=1.0_PS-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.0_PS))
           q=(q-1.0_PS)*(r-1.0_PS)*(s-1.0_PS)
        endif
        if(p.gt.0.) q=-q ! Check whether in bounds.
        p=abs(p)
        if(2.0_PS*p .lt. min(3.0_PS*xm*q-abs(tol1*q),abs(e*q))) then
           e=d ! Accept interpolation.
           d=p/q
        else
           d=xm ! Interpolation failed, use bisection.
           e=d
        endif
     else ! Bounds decreasing too slowly, use bisection.
        d=xm
        e=d
     endif
     a=b ! Move last best guess to a.
     fa=fb
     if(abs(d) .gt. tol1) then ! Evaluate new trial root.
        b=b+d
     else
        b=b+sign(tol1,xm)
     endif
     fb=func(b)
  enddo
  if(debug) write(*,'("brent exceeding maximum iterations. n,iter",2I5,12ES15.6)') &
       iter,a,b,fa,fb
!!$             write(*,'("W,a_c,N_act,M_act,T,r_e",10ES15.6)') ag%TV(n)%W, a_c, N_act,M_act, ag%TV(n)%T,r_e
  
  out=r_n*1.1_PS
  
  
1111 continue
!  write(*,*) "get_hazerad_itr> iter, rd",iter,out
  
contains      
  function func(x) result(out)
    real(PS),intent(in) :: x
    real(PS) :: out
    out=log(SS)*x**4-AA*x**3+(BB*get_osm(1,1.0e+3*n_aps/&
          (den_w*coef4pi3*(max(x,r_n*1.01)**3-r_n**3)))*&
         r_n**3-log(SS)*r_n**3)*x+AA*r_n**3
  end function func
  
end function get_hazerad_itr

function get_hazerad_anal(AA,sb,beta,S_w,r_n) result(r_w)
  !
  ! The formula is (6.2.12) and (6.2.28) of Khvorostyanov and
  ! Curry (2014).
  !
  ! this one considers variable osmotic coefficient
  !
  use acc_amps
  implicit none
  ! coefficients [cm]
  real(PS),intent(in) :: AA
  ! saturation ratio
  real(PS),intent(in) :: S_w
  ! b parameter (dimensionless or cm)
  real(PS),intent(in) :: sb
  ! radius of dry nuclei (cm)
  real(PS),intent(in) :: r_n
  ! beta parameter (0.5: soluble mass distributed to volume, 0: to surface)
  real(PS),intent(in) :: beta
  ! radius of wet aerosol particle (cm)
  real(PS) :: r_w
  ! correction term
  real(PS) :: C_w
  ! 1/(2**(1/3))
  real(PS), parameter :: c2=0.7937005259841
  real :: Zd,Zd2
  complex :: Zdc,Aci,Bci
  real(PS),parameter :: S_w_lmt=0.97
  real(PS) :: r1,a1,b1,r0
  integer :: i,n
!tmp  real(8),dimension(0:4) :: x,d
!tmp  real(8) :: t,p

  Zd=r_n**beta*sqrt(sb/(3.0*AA))
  if(S_w<S_w_lmt.and.r_n>1.0e-6) then
    ! Note: this solution is only valid for Sw < 0.97

    C_w=AA/(3.0*sb**(1.0/3.0)*r_n**(2.0*(1.0+beta)/3.0))
    r_w=r_n*(1.0_PS+sb*r_n**(2.0*(1.0+beta)-3.0)/(-log(S_w))* &
          (1.0_PS+C_w*(-log(S_w))**(-2.0/3.0))**(-3.0))**(1.0/3.0)

!tmp    write(*,*) "ck gethazerad1:",r_w,sb,r_n,S_w,-log(S_w)
  elseif(Zd<=0.1) then

    Zd2=Zd*Zd
    r_w=r_n*(1.0+Zd2-1.0/3.0*Zd2*Zd2*Zd2)

!tmp    write(*,*) "ck gethazerad3:",r_w,sb,r_n,S_w,Zd
  else
    ! Note: this solution is only valid for Sw=1.0

    Zdc=cmplx(Zd,0.0)
    Aci=(1.0+sqrt(1.0-4.0*Zdc**6))**(1.0/3.0)
    Bci=(1.0-sqrt(1.0-4.0*Zdc**6))**(1.0/3.0)

    r_w=c2*r_n*(Aci+Bci)

!tmp    write(*,*) "ck gethazerad2:",r_w,sb,r_n,Aci,Bci,Zdc,S_w

  endif
end function get_hazerad_anal

  function get_aspect_ratio_sol(ms,ishape) result(alpha)
    ! +++ return aspect ratio necesarry for calculating 
    !     terminal velocity +++
    type (Mass_Bin), intent(in)   :: ms
    type (Ice_Shape), intent(in)    :: ishape
    real(PS)                           :: alpha
    integer :: i_sh_le2,i_sh_le4

    i_sh_le2=0.5*(1.0-real(isign(1,ishape%sh_type-3),PS_KIND))
    i_sh_le4=0.5*(1.0-real(isign(1,ishape%sh_type-5),PS_KIND))

    alpha = real(i_sh_le2,PS_KIND)*ishape%phi_cs +&
            (1.0-real(i_sh_le2,PS_KIND))*real(i_sh_le4,PS_KIND)*ms%semi_c/ms%semi_a +&
            (1.0-real(i_sh_le4,PS_KIND))*ms%semi_c/ms%semi_a

  end function get_aspect_ratio_sol

  function get_area_ratio_sol(ms,ishape,mode,A_e,A_c,A_e_ic,A_c_ic,q_cs) result(q)
    ! +++ return area ratio necesarry for calculating 
    !     terminal velocity +++
    type (Mass_Bin), intent(in)   :: ms
    type (Ice_Shape), intent(in)    :: ishape
    integer,intent(in) :: mode
    real(PS)                           :: q
    real(PS),intent(in) :: A_e,A_e_ic,A_c,A_c_ic,q_cs

    integer :: i_sh_le2,i_sh_le4,i_habit_ge5,i_ismod2_eq1 &
        ,i_mode_eq1,i_qcs_gt0,i_mass1_lt_massqagg,i_mass2_lt_massqgrp

    ! this should correspond with aspect ratio diagnosis
    real(PS) :: mass_q_agg,mass_q_grp
    real(PS) :: mass1,mass2
    real(PS),parameter :: mrat_agg=10.0_PS,mrat_grp=10.0_PS ! sense test grp=10
    real(PS) :: q_sidepl,q_0,q_ic_sidepl,q_ic_0,q_ic,q_agg,q_grp
    real(8) :: dum1,dum2

    i_sh_le2=0.5*(1.0-real(isign(1,ishape%sh_type-3),PS_KIND))
    i_sh_le4=0.5*(1.0-real(isign(1,ishape%sh_type-5),PS_KIND))
    i_habit_ge5=0.5*(1.0+real(isign(1,ishape%habit-5),PS_KIND))

    i_ismod2_eq1=max(0.0,(-isign(1,ishape%is_mod(2)-1)*max(0,ishape%is_mod(2)-1)+1.0))*&
                 max(0.0,(-isign(1,ishape%is_mod(2)-1)*min(0,ishape%is_mod(2)-1)+1.0))

    i_mode_eq1=max(0.0,(-isign(1,mode-1)*max(0,mode-1)+1.0))*&
               max(0.0,(-isign(1,mode-1)*min(0,mode-1)+1.0))

    ! use Heymsfield et al (2002)'s side planes
    q_sidepl=max(0.18_PS,min(0.88_PS,(ms%den/0.35_PS)**(1.0/2.34_PS)))
    q_0=A_e/A_c

    ! use Heymsfield et al (2002)'s side planes
    q_ic_sidepl=max(0.18_PS,min(0.88_PS,(ms%mass(imc)/ms%con/ishape%V_ic/0.35_PS)&
             **(1.0/2.34_PS)))
    q_ic_0=A_e_ic/A_c_ic
    q_ic=real(i_habit_ge5,PS_KIND)*q_ic_sidepl+ &
         (1.0-real(i_habit_ge5,PS_KIND))*q_ic_0

    mass_q_agg=ms%mass(imc)/ms%con*mrat_agg
    mass1=(ms%mass(imc)+ms%mass(ima))/ms%con
    dum1=log10(q_ic/max(1.0e-30_RP,q_cs))/log10(ms%mass(imc)/ms%con/mass_q_agg) 
    dum2=mass1/mass_q_agg
    dum1=dum2**dum1
    q_agg=max(0.1_DS,min(1.0_DS,q_cs*dum1))
!!!    q_agg=q_cs*(mass1/mass_q_agg)**dum1

    mass_q_grp=(ms%mass(imc)+ms%mass(ima))/ms%con*mrat_grp
    mass2=(ms%mass(imc)+ms%mass(ima)+ms%mass(imr))/ms%con
    dum1=log10(q_ic/max(1.0e-30_RP,q_cs))/log10(ms%mass(imc)/ms%con/mass_q_grp)
    dum2=mass1/mass_q_grp
    dum1=dum2**dum1
    q_grp=max(0.1_DS,min(1.0_DS,q_cs*dum1))
!!!    q_grp=q_cs*(mass1/mass_q_grp)**&
!!!        (log10(q_ic/max(1.0e-30,q_cs))/log10(ms%mass(imc)/ms%con/mass_q_grp))

    i_qcs_gt0=0.5*(1.0+sign(1.0_PS,q_cs-1.0e-30_PS))
    i_mass1_lt_massqagg=0.5_RP*(1.0_RP-sign(1.0_PS,mass1-mass_q_agg))
    i_mass2_lt_massqgrp=0.5_RP*(1.0_RP-sign(1.0_PS,mass2-mass_q_grp))

    q = &
        ! for prestine and rimed crystals
            real(i_sh_le2,PS_KIND)*&
              (real(i_habit_ge5*i_ismod2_eq1,PS_KIND)*q_sidepl+&
               (1.0_PS-real(i_habit_ge5*i_ismod2_eq1,PS_KIND))*q_0 ) +&
        ! for aggregates
            (1.0_PS-real(i_sh_le2,PS_KIND))*real(i_sh_le4,PS_KIND)*&
            ( real(i_mode_eq1,PS_KIND)*q_0+&
              (1.0_PS-real(i_mode_eq1,PS_KIND))*&
              (&
              real(i_qcs_gt0,PS_KIND)*&
                ( &
                real(i_mass1_lt_massqagg,PS_KIND)*(  &
                  q_agg &
                 )+ &
                (1.0_PS-real(i_mass1_lt_massqagg,PS_KIND))*( &
                  q_cs &
                 ) &
                )+ &
              (1.0_PS-real(i_qcs_gt0,PS_KIND))*q_ic &
              ) &
            ) +&
        ! for graupel
            (1.0_PS-real(i_sh_le4,PS_KIND))* &
            ( real(i_mode_eq1,PS_KIND)*q_0+&
              (1.0_PS-real(i_mode_eq1,PS_KIND))*&
              (real(i_qcs_gt0,PS_KIND)*&
                ( &
                real(i_mass2_lt_massqgrp,PS_KIND)*(  &
                  q_grp &
                 )+ &
                (1.0_PS-real(i_mass2_lt_massqgrp,PS_KIND))*( &
                  q_cs &
                 ) &
                )+ &
              (1.0_PS-real(i_qcs_gt0,PS_KIND))*q_ic &
              ) &
            )

  end function get_area_ratio_sol

  subroutine cal_best_number(D1,Q1,th_var,q,alpha)
    ! +++ return Davis number necesarry for calculating 
    !   terminal velocity

    type (Thermo_Var), intent(in)    :: th_var
    ! aspect ratio
    real(PS),intent(in)    :: alpha
    ! area ratio
    real(PS),intent(in) :: q
    real(DS),intent(inout) :: D1,Q1
    real(DS) :: D1_0,D1_1,Q1_0,Q1_1

    integer :: i_q_le1

    i_q_le1=0.5_RP*(1.0_RP+sign(1.0_PS,1.0_PS-q))

    D1_0=8.0_PS*gg*th_var%den/&
          (PI*(th_var%d_vis**2)*max(alpha, 1.0_PS)*(q**0.25))
    Q1_0=q**0.25

    D1_1=8.0_PS*gg*th_var%den/&
          (PI*(th_var%d_vis**2)*max(alpha, 1.0_PS)*q)
    Q1_1=q


    D1=real(i_q_le1,PS_KIND)*D1_0 +&
       (1.0_RP-real(i_q_le1,PS_KIND))*D1_1

    Q1=real(i_q_le1,PS_KIND)*Q1_0 +&
       (1.0_RP-real(i_q_le1,PS_KIND))*Q1_1

  end subroutine cal_best_number

  function interp_data1d_lut_big(a,xi) result(yi)
    type(data1d_lut_big),intent(in) :: a
    real(PS) :: xi,yi,wx,x1
    integer :: i1

    i1=max(1,min(a%n-1 &
        ,int((xi-a%xs)/a%dx)+1))
    x1=real(i1-1,PS_KIND)*a%dx+a%xs
    wx=min(1.0_RP,max(0.0_RP,(xi-x1)/a%dx))
    yi=(1.0_RP-wx)*a%y(i1)+ &
             wx*a%y(i1+1) 
  
  end function interp_data1d_lut_big

!!$  function interp_data1d_lut_log10_big(a,xi) result(yi)
!!$    type(data1d_lut_big),intent(in) :: a
!!$    real(PS) :: xi,yi,wx,x1,xi_log
!!$    integer :: i1
!!$
!!$    xi_log=log10(max(1.0e-30_RP,xi))
!!$    
!!$    i1=max(1,min(a%n-1 &
!!$        ,int((xi_log-a%xs)/a%dx)+1))
!!$    x1=real(i1-1,PS_KIND)*a%dx+a%xs
!!$    wx=min(1.0_RP,max(0.0_RP,(xi_log-x1)/a%dx))
!!$    yi=(1.0_RP-wx)*a%y(i1)+ &
!!$             wx*a%y(i1+1) 
!!$  
!!$  end function interp_data1d_lut_log10_big

END MODULE class_Mass_Bin
