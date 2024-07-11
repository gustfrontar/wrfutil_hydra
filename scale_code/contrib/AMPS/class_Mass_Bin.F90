MODULE class_Mass_Bin
  use maxdims
  use class_Thermo_Var
  use class_Ice_Shape
  use acc_amps
  use par_amps
  implicit none  
  public
  public :: Mass_Bin
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
  subroutine undef_Mass_Bin(phase,ms)
    integer, intent(in)  :: phase
    type (Mass_Bin), intent(inout) :: ms
    real(PS),parameter :: undef=999.9e+30
    if(phase==1.or.phase==2) then
      ms%p=undef
    endif
    ms%dmassdt=undef
    ms%dvoldt=undef
    ms%dcondt=undef
    ms%mass=undef
    ms%vol=undef
    ms%Ldmassdt=undef
    ms%coef=undef
    ms%con=undef
    ms%mean_mass=undef
    ms%len=undef
    ms%den=undef
    ms%den_as=undef
    ms%den_ai=undef
    ms%a_len=undef
    ms%c_len=undef
    ms%eps_map=undef
    ms%vtm=undef
    ms%tmp=undef
    ms%fv=undef
    ms%fkn=undef
    ms%fac=undef
    ms%fh=undef
    ms%Nre=undef
    ms%CAP=undef
    ms%CAP_hex=undef
    ms%semi_a=undef
    ms%semi_c=undef
    ms%e_sat=undef
    ms%r_act=undef
    ms%r_crt=undef
    ms%mark=0
    ms%inevp=0
    ms%inmlt=0

  end subroutine undef_Mass_Bin

  subroutine cal_surface_temp( phase, ms, th_var,estbar,esitbar)
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
    type (Mass_Bin), intent(inout) :: ms
    type (Thermo_Var), intent(in) :: th_var
    integer :: em

    ! latent heat of condensation in ergs/g
    real(PS), parameter              :: L_e = 2.5e+10_PS
    ! latent heat of sublimation in ergs/g
    real(PS), parameter             :: L_s = 2.8337e+10_PS
    ! latent heat of freeze in ergs/g
    real(PS), parameter             :: L_f = 0.3337e+10_PS

!!c    real(PS) :: r_e
    ! maximum realisitc super saturation
    real(PS),parameter :: max_Sw=0.001

    integer             :: i,phase2
    real(PS) :: TS_A1,TS_B1,TS_B11,TS_B12,tmp,tmp_n,es_sfc
    real(PS), parameter              :: R_u = 8.3143e+07
    real(PS), parameter             :: M_w=18.0
    real(PS), parameter             :: MR=M_w/R_u
    real(PS), parameter         :: den_w = 1.0
    ! heat content of water
    real(PS),parameter :: c_w=4.187e+5

    ! saturation vapor lookup table
    real(DS)  :: estbar(150),esitbar(111)


    ! temperature at triple point 
    real(PS), parameter          :: T_0 = 273.16_PS

    ! modified conductivity
    real(PS) :: mk_a

! <<< 2014/12 T. Hashino bug found 
    phase2=phase
! >>>> 2014/12 T. Hashino bug found 
    if( ms%con > 1.0e-30_PS .and. ms%mass(1) > 1.0e-30_PS ) then
       if( phase == 1 ) then

!!c       ms%tmp = th_var%T + ( L_e * ms%Ldmassdt(1) )/&
!!c            (4.0_PS*PI*ms%CAP*th_var%k_a*ms%fh)

          mk_a = get_mod_thermal_cond(ms%a_len,th_var%K_a,th_var%T,th_var%den)

          ms%tmp=th_var%T +  L_e*(ms%coef(1)*th_var%s_v_n(phase)+ms%coef(2))/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh)            
!!c          write(*,'("surfT",10ES15.6)') th_var%s_v_n(phase),th_var%T,ms%tmp,ms%coef(1),ms%coef(2)
          
       else if( phase == 2 ) then

          call cal_coef_Ts(ms,th_var,phase,TS_A1,TS_B11,TS_B12,phase2)

          if(T_0>th_var%T) then
             if(phase2==1) then
!!c                TS_B1=TS_B11*min(max_Sw,th_var%s_v_n(phase2))+TS_B12
                TS_B1=TS_B11*th_var%s_v_n(phase2)+TS_B12
             else
                TS_B1=TS_B11*th_var%s_v_n(phase2)+TS_B12
             end if

             em=0
             call cal_tsfc_iter(th_var%T,TS_A1,TS_B1,phase2,estbar,esitbar,tmp,em)
             if(em/=0) then
                write(*,*) "sfc temp error"
                write(*,'(2I5,10Es15.6)') em,phase2,th_var%T,TS_A1,TS_B1,tmp
             endif
             ms%tmp=min(tmp,T_0)
          else
             ! 
             ! assume that the riming process produce surface water 
             ! due to the liquid layer over the hydrometeors:
             ! no fusion heat release by riming process
             !
             ! The liquid layer is turbulent, so the heat received at the
             ! surface is balanced with the one of ice core.
             !
             ! Therefore, the surface temperature is T_0.
! <<< 2015/03 T. Hashino mod for KiD
             ! above assumption seems not working. So set it to the ambient 
             ! temperature
             ms%tmp=T_0
!mod             ms%tmp=th_var%T
! >>> 2015/03 T. Hashino
!
          end if

!!c          ms%e_sat=get_sat_vapor_pres(phase2,ms%tmp)
          ms%e_sat=get_sat_vapor_pres_lk(phase2,ms%tmp,estbar,esitbar)
       end if
    else
       ms%tmp = th_var%T
       ms%e_sat=get_sat_vapor_pres_lk(phase2,ms%tmp,estbar,esitbar)
    end if
  end subroutine cal_surface_temp


  subroutine cal_tsfc_iter(T_inf,A1,B1,phase2,estbar,esitbar,out,em)
    real(PS),intent(in) :: T_inf,A1,B1
    integer,intent(in) :: phase2
    ! saturation vapor lookup table
    real(DS)  :: estbar(150),esitbar(111)
    real(PS),intent(inout) :: out
    integer,intent(inout) :: em

    REAL(PS) :: tol=1.0e-6
    integer :: ITMAX
    real(PS) :: eps
!    PARAMETER (ITMAX=50,EPS=epsilon(out))
    PARAMETER (ITMAX=50,EPS=3.e-8)


    !  Using Brent's method, find the root of a function func known to lie between x1 and x2.
    !  The root, returned as zbrent, will be refined until its accuracy is tol.
      
    !  Parameters: Maximum allowed number of iterations, and machine floating-point precision.
    INTEGER :: iter
    REAL(PS) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

!!c             tmp=th_var%T
!!c             tmp=200.0
          
!!c             tmp=T_0
!!c             write(*,*) phase
!!c             do iter=1,10
!!c                write(*,*) "tmp,esat",tmp,get_sat_vapor_pres_lk(phase,tmp,estbar,esitbar)
!!c                tmp=tmp+5.0 
!!c             end do
!!c             stop
    ! ++++++++++++++++++++++++++++
    em=0

    a=170.16_PS
!!c    a=T_inf
    a=273.16
    b=300.16_PS
!!c             b=273.16_PS
    iter=1
2221  continue
    fa=func(a)
2222  continue
    fb=func(b)
      
    if(fa.gt.0.0_PS.and.fb.gt.0.0_PS) then
       iter=iter+1
       if(iter>50.or.a<123.0) then
!!c          write(*,'("iter at Ts large 1",2I5,10ES15.6)') phase2,iter,T_inf,A1,B1,fa,fb,a,b
          em=1
          out=T_inf
          return
       end if
!!c       a=a-50.0_PS
!!c       a=a/1.2_PS
       a=a-10.0_PS
!!c       stop
       goto 2221
    elseif(fa.lt.0.0_PS.and.fb.lt.0.0_PS) then
!!c       if((fa.gt.0.0_PS.and.fb.gt.0.0_PS).or.(fa.lt.0.0_PS.and.fb.lt.0.0_PS)) then
       iter=iter+1
       if(iter>50) then
          write(*,'("iter at Ts large 2",2I5,10ES15.6)') phase2,iter,T_inf,A1,B1,fa,fb,a,b
          em=2
          out=T_inf
          return
       end if
!!c       stop  
!!c                b=b*100.0
       b=b*2.0_PS
!!c       b=b+10.0_PS
!!c       b=b+50.0_PS
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
    write(*,'("brent exceeding maximum iterations. n,iter",2I5,12ES15.6)') &
         iter,a,b,fa,fb
!!c             write(*,'("W,a_c,N_act,M_act,T,r_e",10ES15.6)') ag%TV(n)%W, a_c, N_act,M_act, ag%TV(n)%T,r_e
      
    out=T_inf
      
      
1111  continue

  contains      
    function func(x) result(out)
      real(PS),intent(in) :: x
      real(PS) :: out
!!c      out=x-B1+A1*get_sat_vapor_pres(phase2,x)/x
      out=x-B1+A1*get_sat_vapor_pres_lk(phase2,x,estbar,esitbar)/x
    end function func
  end subroutine cal_tsfc_iter


  subroutine cal_capacitance(level,phase,th_var, ms, is)
    ! +++ Calculate capacitance based on the spheroid assumption,
    !     for example, see Chen and Lamb (1994a).                     +++
    ! phase of hydrometeor
    ! 1: liquid, 2: solid
    integer, intent(in)  :: phase
    integer, intent(in)  :: level
    type (Thermo_Var), intent(in) :: th_var
    type (Mass_Bin), intent(inout) :: ms
    type (Ice_Shape), optional, intent(in)    :: is
    real(PS)                    :: eps, d, out, phi
    integer                     :: i
    ! coefficient for capacitance by Chiruta and Wang (2003)
!!c    real(PS),parameter   :: CAP_ros=0.434*N_lob**0.257
    real(PS),parameter   :: CAP_ros=0.619753727
    ! temperature of freezing
    real(PS), parameter           :: T_0 = 273.16

    ! initialize
    ms%CAP=0.0_PS
    ms%CAP_hex=0.0_PS
    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return

    if( phase == 1 ) then
       ! --- in case of liquid phase ---
       ms%CAP = ms%len/2.0_PS
    else if( phase == 2 ) then
       ! --- in case of solid phase ---
       ! +++ first caclulate for a hexagonal monocrystal +++
       if( is%phi_ic < 1.0_PS ) then
          ! --- in case of oblate spheroids ---
          d = sqrt( ms%a_len**2.0 - ms%c_len**2.0)
          eps = sqrt( 1.0_PS - is%phi_ic**2.0)
          ms%CAP_hex = d/asin(eps)
       else if( is%phi_ic> 1.0_PS ) then
          ! --- in case of prolate spheroids ---
          d = sqrt( ms%c_len**2.0 - ms%a_len**2.0)
          eps = sqrt( 1.0_PS - is%phi_ic**(-2.0))
          ms%CAP_hex = d/log((1.0_PS + eps)*is%phi_ic )
!!c             ms%CAP_hex = ( 0.708_PS + 0.615_PS * (is%phi_ic**(-0.76)))*ms%c_len
       else 
          ! --- in case of sphere ---
          ms%CAP_hex = ms%a_len
       end if


       ! ++++ second calculate the capacitance for polycrystals ++++
       if( is%sh_type <= 2 ) then
          ! --- case of ice crystals ---
          if(is%is_mod(2)==1) then
             if(is%habit<=3) then
                ms%CAP=ms%CAP_hex
             elseif(is%habit==4) then
                ms%CAP = CAP_ros*is%r
             elseif(is%habit==5) then
                ! assume phi=0.25 with a spheroid
                eps=0.968245837_PS
                d=is%e*eps
                ms%CAP=d/asin(eps)
             else
                ! assume phi=0.25 with a spheroid
                eps=0.968245837_PS
                d=max(is%e,is%r)*eps
                ms%CAP=d/asin(eps)
   !!c             ms%CAP = max(is%r,is%e)
             end if
          else
             ! --- in case of sphere ---
             ms%CAP = ms%semi_a
          endif
       elseif( is%sh_type <= 3 ) then
          ! --- case of aggregates and rimed aggregates ---
          ! use Capacitance/D_max = 0.25 derived by 
          !      Westbrook et al. (2008), JAS
          ms%CAP=(is%V_cs/coefpi6)**(1.0/3.0)*0.25_PS
       else
          ! --- case of graupels ---
          phi = is%phi_cs
          if( phi < 1.0_PS ) then
             ! --- in case of oblate spheroids ---
             d = sqrt( ms%semi_a**2.0 - ms%semi_c**2.0)
             eps = sqrt( 1.0_PS - phi**2.0)
             ms%CAP = d/asin(eps)
          else if( phi > 1.0_PS ) then
             ! --- in case of prolate spheroids ---
             d = sqrt( ms%semi_c**2.0 - ms%semi_a**2.0)
             eps = sqrt( 1.0_PS - phi**(-2.0))
             ms%CAP = d/log((1.0_PS + eps)*phi )
!!c             ms%CAP = ( 0.708_PS + 0.615_PS * (phi**(-0.76)))*ms%c_len
          else 
             ! --- in case of sphere ---
             ms%CAP = ms%semi_a
          end if
       end if
       if(level<=5.and.th_var%T>=T_0) then
          ms%CAP=max(ms%CAP,(ms%mean_mass/coef4pi3)**(1.0/3.0))
       elseif(th_var%T>=T_0) then
          ms%CAP=max(ms%CAP,(ms%mass(imw)/ms%con/coef4pi3)**(1.0/3.0))
       end if 
       
       ! +++ additional factor to represent the shape +++
       ! under construction
    end if
  end subroutine cal_capacitance

  subroutine cal_ventilation_coef( phase, ms, th_var,ishape)
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
    type (Mass_Bin), intent(inout) :: ms
    type (Thermo_Var), intent(in) :: th_var
    type (Ice_Shape), optional       :: ishape

    ! heat capacity of dry air
    real(PS), parameter              :: C_pd = 10.04e+06
    ! Schmidt number
    real(PS)                   :: N_sc
    ! Nusselt number
    real(PS)                   :: N_ns
    real(PS)                   :: Xv, Xh
    ! major semiaxis 
    real(PS)                   :: r_m
    ! gas constant for vapor in  ergs/deg/g
    real (PS), parameter             :: R_v = 4.615e+06_PS
    ! temperature of freezing
    real(PS), parameter           :: T_0 = 273.16

    ! initialize
    ms%Nre=0.0_PS
    ms%fv=1.0_PS
    ms%fh=1.0_PS
    ms%fkn=1.0_PS
    ms%fac=1.0_PS
    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return

    ! +++ calculate the Schmidt number +++
    N_sc = th_var%d_vis/th_var%den/th_var%D_v

    ! +++ calculate the Nusselt or Prandtl number +++
    ! +++ 
!!c  need to check
!!c    N_ns = th_var%d_vis/C_pd/th_var%k_a
    N_ns = th_var%d_vis/th_var%den/th_var%k_a
!!c    N_ns = 0.72_PS

    if( phase == 1 ) then
       ! --- in case of liquid phase ---
       ! assume the spherical water drops
       ! +++ calculate Reynolds number +++
       ms%Nre = ms%len*ms%vtm*th_var%den&
            /th_var%d_vis
         
       ! +++ calculate dimensionless number X +++
       Xv = (N_sc**(1.0/3.0))*sqrt(ms%Nre)
       Xh = (N_ns**(1.0/3.0))*sqrt(ms%Nre)

       if( Xv < 1.4_PS ) then
          ms%fv = 1.0_PS + 0.108_PS*(Xv**2.0)
       else if( 1.4_PS <= Xv .and. Xv <= 51.4_PS ) then
          ms%fv = 0.78_PS + 0.308_PS*Xv
       else
          ms%fv = 0.78_PS + 0.308_PS*51.4_PS
       end if

!!c       ms%fv=1.0_PS
       
       if( Xh < 1.4_PS ) then
          ms%fh = 1.0_PS + 0.108_PS*(Xh**2.0)
       else if( 1.4_PS <= Xh .and. Xh <= 51.4_PS ) then
          ms%fh = 0.78_PS + 0.308_PS*Xh
       else
          ms%fh = 0.78_PS + 0.308_PS*51.4_PS
       end if


       ! +++ calculate kinetic effect +++
       r_m=ms%len*0.5_PS
       ms%fkn = get_fkn(th_var,phase,r_m)

    else if( phase == 2 ) then
       ! --- in case of ice phase ---
       ! assume the ice crystal model
       if(.not.present(ishape)) then
          write(*,*) " cal_ventilation > Need ice shape object!"
          stop
       end if
       ! +++ calculate Reynolds number +++
       ms%Nre = ms%len*ms%vtm*th_var%den&
            /th_var%d_vis
       
       ! +++ calculate dimensionless number X +++
       Xv = (N_sc**(1.0/3.0))*sqrt(ms%Nre)
       Xh = (N_ns**(1.0/3.0))*sqrt(ms%Nre)

       if(ishape%sh_type<=2) then
          if( Xv < 1.0_PS ) then
             ms%fac = (1.0_PS + 0.14_PS*(Xv**2.0)*&
                  sqrt(ms%c_len/ms%len/2.0_PS))/&
                  (1.0_PS + 0.14_PS*(Xv**2.0)*&
                  sqrt(ms%a_len/ms%len/2.0_PS))
             
          else if( 1.0_PS <= Xv ) then
             ms%fac = (0.86_PS + 0.28_PS*Xv*&
                  sqrt(ms%c_len/ms%len/2.0_PS))/&
                  (0.86_PS + 0.28_PS*Xv*&
                  sqrt(ms%a_len/ms%len/2.0_PS))
          end if
       else
          ! assume crystals within aggregates, rimed aggregates, graupel
          ! won't have ventilation effect on the shape.
          ms%fac = 1.0_PS
!tmp          if( Xv < 1.0_PS ) then
!tmp             ms%fac = (1.0_PS + 0.14_PS*(Xv**2.0)*&
!tmp                  sqrt(ms%semi_c/ms%len/2.0_PS))/&
!tmp                  (1.0_PS + 0.14_PS*(Xv**2.0)*&
!tmp                  sqrt(ms%semi_a/ms%len/2.0_PS))
!tmp             
!tmp          else if( 1.0_PS <= Xv ) then
!tmp             ms%fac = (0.86_PS + 0.28_PS*Xv*&
!tmp                  sqrt(ms%semi_c/ms%len/2.0_PS))/&
!tmp                  (0.86_PS + 0.28_PS*Xv*&
!tmp                  sqrt(ms%semi_a/ms%len/2.0_PS))
!tmp          end if
       end if

       if(ishape%sh_type<=2.and.ishape%habit==4) then
          ! if the ice crystal is bullet rosettes
          ! use Lie et al (2003)'s parameterization.
          ms%fv = 1.0_PS + 0.3005_PS*Xv - 0.0022_PS*Xv*Xv
          ms%fh = 1.0_PS + 0.3005_PS*Xh - 0.0022_PS*Xh*Xh
       else
          ! if the ice crystal is hexagonal plates, dendrites, or columns,
          ! if it is polycrystal,
          ! or if it is aggregate,
          if( Xv < 1.0_PS ) then
             ms%fv = 1.0_PS + 0.14_PS*(Xv**2.0)
          else if( 1.0_PS <= Xv ) then
             ms%fv = 0.86_PS + 0.28_PS*Xv
          end if
          if( Xh < 1.0_PS ) then
             ms%fh = 1.0_PS + 0.14_PS*(Xh**2.0)
          else if( 1.0_PS <= Xh ) then
             ms%fh = 0.86_PS + 0.28_PS*Xh
          end if
       end if

       ! +++ calculate kinetic effect +++
       ms%fkn=1.0_PS

       r_m = sqrt(ms%semi_a**2+ms%semi_c**2)
       ms%fkn=get_fkn(th_var,phase,r_m)

       if(ms%inmlt==1.and.th_var%T>=T_0) then
          if( Xh < 1.4_PS ) then
             ms%fh = max(ms%fh,1.0_PS + 0.108_PS*(Xh**2.0))
          else if( 1.4_PS <= Xh .and. Xh <= 51.4_PS ) then
             ms%fh = max(ms%fh,0.78_PS + 0.308_PS*Xh)
          else
             ms%fh = max(ms%fh,0.78_PS + 0.308_PS*51.4_PS)
          end if
       end if
    end if

  END subroutine cal_ventilation_coef

  subroutine cal_terminal_vel( phase, ms, th_var, mode, level,em, ishape)
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
    ! phase of hydrometeor
    ! 1: liquid, 2: solid
    integer, intent(in)  :: phase
    type (Mass_Bin), intent(inout)   :: ms
    type (Ice_Shape), optional       :: ishape
    type (Thermo_Var), intent(in)    :: th_var
    integer                          :: mode
    integer, intent(in)              :: level

    ! effective cross-sectional area
    real(PS)                    :: A_e
    ! circumscribed cross-sectional area
    real(PS)                    :: A_c
    ! the ratio of the effective over the circumscribed cross-sectional area, Ae/A
    real(PS)                    :: q
    ! charasteristic length defined for Bohm (1992)
    real(PS)                    :: char_len

    ! axial ratio, c/(2*(a+d))
    real(PS)                    :: alpha
    ! shape factor
    real(PS)                    :: k
    ! the Davies or Best number
    real(PS)                    :: X, X_0, X_old
    ! pressure drag coefficient
    real(PS)                    :: C_DP
    ! 
    real(PS)                    :: beta, gamma, xx
    ! Oseen-type solution for the drag coefficient
    real(PS)                    :: C_DO
    ! Reynolds number
    real(PS)                    :: N_re_0

    ! total surface area of the body
    real(PS)                    :: OMEGA
    ! Perimeter 
    real(PS)                    :: P

    ! density of ice (bulk density)
    real(PS), parameter         :: den_i = 0.91668
    real(PS), parameter         :: gg = 980.0
    
    ! Stoke's terminal velocity
    real(PS) :: U_S
    ! radius of a drop
    real(PS) :: rad
    ! Bond number and physical property number
    real(PS) :: NBO,NP
    real(PS) :: Y,lambda_a
    real(PS), parameter         :: den_w = 1.0

    real(PS)  :: beta0
    integer                    :: i
    integer,intent(inout) :: em

    ! initialize
    ms%vtm=0.0_PS
    ms%Nre=0.0_PS
    if(phase==2) then
       ms%len=0.0_PS
    endif
    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return

    if( phase == 1 ) then
       ! +++ calculate the mass component for one particle +++
       rad=ms%len*0.5_PS
       
       
       if(rad<0.5e-4_PS) then  ! bug found 2016/05/30
          ms%vtm=0.0_PS
       elseif(0.5e-4_PS<=rad.and.rad<10.0e-4_PS) then
          ! calculate Stokes terminal velocity
          U_S=rad**2.0*gg*(den_w-th_var%den_a)/4.5_PS/th_var%d_vis
          lambda_a=6.6e-6_PS*(th_var%d_vis/1.818e-4_PS)*(1013250.0_PS/th_var%P)*(th_var%T/293.15_PS)
          ms%vtm=(1.0_PS+1.26_PS*lambda_a/rad)*U_S

       elseif(rad<535.0e-4_PS) then
          X=log(32.0_PS*rad**3.0*(den_w-th_var%den_a)*th_var%den_a*gg/3.0_PS/th_var%d_vis**2.0)

   
!          Y=-0.318657e+1_PS+0.992696_PS*X-0.153193e-2_PS*X**2.0-0.987059e-3_PS*X**3.0&
!               -0.578878e-3_PS*X**4.0+0.855176e-4*X**5.0-0.327815e-5*X**6.0

          Y=-0.318657e+1_PS+0.992696_PS*X-0.153193e-2_PS*X*X-0.987059e-3_PS*X*X*X &
               -0.578878e-3_PS*X*X*X*X+0.855176e-4*X*X*X*X*X &
               -0.327815e-5*X*X*X*X*X*X
          ms%Nre=exp(Y)
          ! +++ calculate the terminal velocity +++
          ms%vtm=ms%Nre*th_var%d_vis/(2.0_PS*rad*th_var%den_a)

       else
          rad=min(rad,3500.0e-4_RP)
          
          NBO=gg*(den_w-th_var%den_a)*rad**2/th_var%sig_wa
          NP=th_var%sig_wa**3*th_var%den_a**2/th_var%d_vis**4/gg
          
          X=log(NBO*NP**(1.0_RP/6.0_RP)*16.0_PS/3.0_PS)
!          Y=-0.500015e+1_PS+0.523778e+1_PS*X-0.204914e+1_PS*X**2.0+0.475294_PS*X**3.0&
!               -0.542819e-1_PS*X**4.0+0.238449e-2_PS*X**5.0
          Y=-0.500015e+1_PS+0.523778e+1_PS*X-0.204914e+1_PS*X*X+0.475294_PS*X*X*X &
               -0.542819e-1_PS*X*X*X*X+0.238449e-2_PS*X*X*X*X*X

          
          ms%Nre=NP**(1.0_RP/6.0_RP)*exp(Y)
          ! +++ calculate the terminal velocity +++
          ms%vtm=ms%Nre*th_var%d_vis/(2.0_PS*rad*th_var%den_a)

          
       end if

    else if( phase == 2 ) then
       ! --- in case of ice phase ---
       if(.not.present(ishape)) then
          write(*,*) " cal_terminal_vel > Need ice shape object!"
          stop
       end if
       
       ! +++ calculate some important variables, depending on level +++
       call cal_some_vel(ms,ishape,mode,alpha,OMEGA,A_e,A_c,P,q,char_len)
       ms%len = OMEGA/P
       ishape%q_e=q

!!c       if(ishape%sh_type>=3) then
!!c          open( unit=21, file="check_svtm.dat",POSITION='APPEND')
!!c          write(21,'(10ES15.6)') ms%mean_mass,ms%semi_a,ms%semi_c,alpha,q,char_len,ms%Nre
!!c          close(21)
!!c       end if


       
       ! +++ calculate the shape factor +++
       if( alpha <= 0.16_PS ) then
          k = 0.85_PS
       else if ( 0.16_PS < alpha .and. alpha <= 1.0_PS ) then
          k = 0.82_PS + 0.18_PS*alpha
       else if ( 1.0_PS < alpha .and. alpha <= 3.0_PS ) then
          k = 0.37_PS + 0.63_PS/alpha
       else if ( 3.0_PS < alpha ) then
          k = 1.33_PS/( log(alpha) + 1.19_PS)
       end if
       
       ! +++ calculate the Davies or Best number +++
       if( q <= 1.0_PS ) then
          X = 8.0_PS*ms%mean_mass*gg*th_var%den_a/&
               (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*(q**0.25))
       else if( q > 1.0_PS ) then
          X = 8.0_PS*ms%mean_mass*gg*th_var%den_a/&
               (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*q)
       end if
       
       ! +++ calculate the pressure drag C_DP +++
       if( 0.3_PS < alpha .and. alpha < 1.0_PS ) then
          gamma = 3.76_PS - 8.41*alpha + 9.18*alpha**2.0 - 3.53*alpha**3.0
          C_DP = 0.292*k*gamma
       else if( alpha <= 0.3_PS ) then 
          gamma = 1.98
          C_DP = 0.292*k*gamma
!!c          else if( 1.0_PS <= alpha .and. 1.0_PS < q ) then
       else if( 1.0_PS <= alpha ) then
          if( 1.0_PS <= q ) then
             C_DP = (1.46*q-0.46)*(0.492-0.200/sqrt(alpha))
          else
!!c             gamma = 3.76_PS - 8.41*alpha + 9.18*alpha**2.0 - 3.53*alpha**3.0
             gamma = 1.0
             C_DP = 0.292*k*gamma
          end if
       end if
!!c       write(*,'(I3,5ES14.5)') i,k,gamma,alpha,q,C_DP
       
       ! +++ calculate the Oseen-type drag coefficient, C_DO +++
       C_DO = 4.5*(k**2.0)*max(alpha, 1.0_PS)
       
       ! +++ calculate the Reynolds number +++
       beta = sqrt(1.0_PS + (C_DP/(6.0_PS*k))*sqrt(X/C_DP)) - 1.0_PS
       beta0=beta
       N_re_0 = 6.0*k*beta**2.0/C_DP
       
       ! +++ correction for transition to turbulent flow +++
       if( 1000.0 <= N_re_0 ) then
          X_0 = 2.8e+06
          C_DP=C_DP*(1.0_PS + 1.6_PS*(X/X_0)**2.0)/(1.0_PS + (X/X_0)**2.0)

          beta = sqrt(1.0_PS + (C_DP/(6.0_PS*k))*sqrt(X/C_DP)) - 1.0_PS             
          N_re_0 = 6.0*k*beta**2.0/C_DP
          
       end if
       
       ! +++ matching theory +++
       gamma = (C_DO-C_DP)/(4.0_PS*C_DP)
       ms%Nre = N_re_0*(1.0_PS+2.0_PS*beta*exp(-beta*gamma)/&
            ((2.0_PS+beta)*(1.0_PS+beta)))
       if( ms%Nre < 0.0_PS ) then
          write(*,*) "Reynolds number is negative!!"
          write(*,'(10ES14.6)') ms%mass,ms%con,ms%a_len,ms%c_len,q,k,gamma
          write(*,'(10ES14.6)') N_re_0,exp(-beta*gamma),beta,X,C_DP
          write(*,'(2I5,10ES14.6)') ishape%sh_type,ishape%habit,ms%semi_a,ms%semi_c,ishape%phi_cs
          write(*,'(10ES14.6)') ishape%semi_aip,ishape%semi_cip,ishape%V_cs
          em=2
          return
       end if
       
       ! +++ calculate the terminal velocity +++
       ms%vtm = ms%Nre*th_var%d_vis/(char_len*th_var%den_a)
       if(ms%vtm>5.0e+3.and.ms%den>0.8_PS) then
          write(*,*) "terminal_vel > The ice vel is set to 50.0m/s"
          write(*,*) "old vel",ms%vtm
          write(*,'("mass,con,alen,clen,q,k,gamma",10ES14.6)') ms%mass(1),ms%con,ms%a_len,ms%c_len,q,k,gamma
          write(*,'("N_re0,ex,beta,X,C_DP",10ES14.6)') N_re_0,exp(-beta*gamma),beta,X,C_DP
          write(*,'("semi_a,semi_c,mean_mass,beta0",10ES14.6)') ms%semi_a,ms%semi_c,ms%mean_mass,beta0
          write(*,'("char_len,d_vis,den,Nre,A_c",10ES14.6)') char_len,th_var%d_vis,th_var%den,ms%Nre,A_c
          write(*,'("pressure,T",10ES14.6)') th_var%P,th_var%T
          ms%vtm=5.0e+3_PS
       elseif(ms%vtm>2.0e+3) then
          ms%vtm=2.0e+3_PS
       elseif(ms%vtm<0.0) then
          write(*,*) "terminal_vel > The ice vel is set to 0.0m/s"
          write(*,*) "old vel",ms%vtm
          write(*,'("mass,con,alen,clen,q,k,gamma",10ES14.6)') ms%mass(1),ms%con,ms%a_len,ms%c_len,q,k,gamma
          write(*,'("N_re0,ex,beta,X,C_DP",10ES14.6)') N_re_0,exp(-beta*gamma),beta,X,C_DP
          write(*,'("semi_a,semi_c,mean_mass,beta0",10ES14.6)') ms%semi_a,ms%semi_c,ms%mean_mass,beta0
          write(*,'("char_len,d_vis,den,Nre,A_c",10ES14.6)') char_len,th_var%d_vis,th_var%den,ms%Nre,A_c
          write(*,'("pressure,T",10ES14.6)') th_var%P,th_var%T
          ms%vtm=0.0_PS
       endif
    end if

  end subroutine cal_terminal_vel

  subroutine cal_some_vel(ms,ishape,mode,alpha,OMEGA,A_e,A_c,P,q,char_len)
    type (Mass_Bin), intent(in)   :: ms
    type (Ice_Shape), intent(in)  :: ishape
    integer              :: mode
    ! axial ratio, c/(2*(a+d))
    real(PS)                    :: alpha
    ! total surface area of the body
    real(PS)                    :: OMEGA
    ! effective cross-sectional area
    real(PS)                    :: A_e
    ! circumscribed cross-sectional area
    real(PS)                    :: A_c
    ! Perimeter 
    real(PS)                    :: P
    ! the ratio of the effective over the circumscribed cross-sectional area, Ae/A
    real(PS)                    :: q,q_ic,q_cs
    ! charasteristic length defined for Bohm (1992)
    real(PS)                    :: char_len
    ! bulk density of aggregates (Bohm 1989)
    real(PS)                    :: b_den
    ! aspect ratio assumed for rosettes and irregular crystals.
    real(PS),parameter :: rho_re=0.25

    ! this should correspond with aspect ratio diagnosis
    real(PS) :: mass_q_agg,mass_q_grp

    real(PS) :: mass1
    real(PS),parameter :: mrat_agg=10.0_PS,mrat_grp=10.0_PS ! sense test grp=10

    real(PS)  :: dum
    if( ishape%sh_type == 1 .or. ishape%sh_type == 2 ) then
       ! +++ get aspect ratio +++ 
       alpha=ishape%phi_cs
       
       ! +++ calculate the total surface area +++
       OMEGA = get_total_sfc_area( ishape, ms%semi_a, ms%semi_c, mode)
       
       ! +++ calculate the area ratio +++
       ! ++++++ calculate the effective cross section ++++++
       A_e = get_effect_area( ishape,ishape%sh_type,ms%semi_a, ms%semi_c, mode)
       
       ! ++++++ calculate the circumscribed cross-sectional area ++++++
       A_c = get_circum_area( ishape,ishape%sh_type, ms%semi_a, ms%semi_c, mode)
       
!!c       if( ishape%sh_type == 2 ) then
!!c          ! --- case of lightly rimed ice crystal ---
!!c          ! NOTE: the definition by volume is necessary to be studied.
!!c          dum=(4.0_PS*PI/3.0_PS)*ms%semi_c*ms%semi_a**2.0
!!c          q = min( ( ishape%V_ic + ms%vol(1) )/dum, 1.0_PS)
!!c       else
       if( ishape%habit>=5 .and. ishape%is_mod(2)==1) then
          ! use Heymsfield et al (2002)'s side planes
          q=max(0.18_PS,min(0.88_PS,(ms%den/0.35_PS)**(1.0/2.34_PS)))
       else
          q=A_e/A_c
       end if
!!c       end if
       ! +++ calculate the characteristic length based on
       !     perimeter and total surface area +++
       P = get_perimeter( ishape, ms%semi_a, ms%semi_c, mode)
       
       ! +++ calculate characteristic length +++
       char_len = 2.0*sqrt(A_c/PI)
    else if( ishape%sh_type <= 4 ) then
!!c       if(mode==2) then
!!c          alpha=0.25_PS
!!c       else
       alpha = ms%semi_c/ms%semi_a
!!c       end if     
   
       ! +++ calculate the total surface area +++
       OMEGA = get_total_sfc_area( ishape, ms%semi_a, ms%semi_c, mode)
       
       ! +++ calculate the area ratio +++
       ! ++++++ calculate the effective cross section ++++++
       A_e = get_effect_area( ishape, ishape%sh_type,ms%semi_a, ms%semi_c, mode)
       
       ! ++++++ calculate the circumscribed cross-sectional area ++++++
       A_c = get_circum_area( ishape,ishape%sh_type, ms%semi_a, ms%semi_c, mode)

       ! +++ calculate the characteristic length based on
       !     perimeter and total surface area +++
       P = get_perimeter( ishape, ms%semi_a, ms%semi_c, mode)

       ! +++ calculate characteristic length +++
       char_len = 2.0*sqrt(A_c/PI)

       if(mode==2) then
!!c          select case(ishape%habit)
!!c          case (1)
!!c             ! 1 : hexagonal plate
!!c             b_den=0.9_PS
!!c          case (2)
!!c             ! 2 : broad-branch crystal
!!c             b_den=0.5_PS
!!c          case (3)
!!c             ! 3 : columnar crystal
!!c             b_den=0.8_PS
!!c          case (4)
!!c             ! 4 : rosette bullet
!!c             b_den=0.8_PS
!!c          case (5)
!!c             ! 5 : irregular ice crystals
!!c             b_den=0.9_PS
!!c          end select
!!c          q = 24.0_PS*ms%mean_mass/(PI*b_den*char_len**3.0)

       
          ! +++ calculate the area ratio of ice crystals +++
          ! ++++++ calculate the effective cross section ++++++
          if( ishape%habit>=5 ) then
             ! use Heymsfield et al (2002)'s side planes
             q_ic=max(0.18_PS,min(0.88_PS,(ms%mass(imc)/ms%con/ishape%V_ic/0.35_PS)&
                  **(1.0/2.34_PS)))
          else
             q_ic=get_effect_area( ishape,1,ms%a_len, ms%c_len, mode)/&
                  get_circum_area( ishape,1,ms%a_len, ms%c_len, mode)
          end if
!!c          q_cs=(ms%mean_mass/den_i/ishape%V_cs)**(2.0/3.0)
!!c          q_cs=ms%mean_mass/den_i/ishape%V_cs
!!c          q_cs=(ms%mass(imc)+ms%mass(ima)+ms%mass(imr))/ms%con/den_i/ishape%V_cs
          q_cs=(ms%mass(imc)+ms%mass(ima)+ms%mass(imr))/ms%con/den_i/&
                get_vip(ishape%is_mod(2),ishape%phi_cs,ishape%semi_aip)
          
          if(q_cs>1.0e-30_PS) then
             mass_q_agg=ms%mass(imc)/ms%con*mrat_agg
             mass1=(ms%mass(imc)+ms%mass(ima))/ms%con

             if(mass1<mass_q_agg) then
!!c                q=10.0**((log10(q_ic)-log10(q_cs))/(log10(ms%mass(imc)/ms%con)-log10(mass_q_agg))*&
!!c                     (log10(mass1)-log10(mass_q_agg))+log10(q_cs))
                q=q_cs*(mass1/mass_q_agg)**&
                 (log10(q_ic/q_cs)/log10(ms%mass(imc)/ms%con/mass_q_agg))



                if(q<0.0.or.q>2.0.or.(q<0.0.and.q>2.0)) then
                   write(*,*) "q is not right for agg",ishape%habit,q,q_ic,q_cs,ms%mass(imc),ms%con,ms%mean_mass
                   write(*,*) ms%a_len,ms%c_len,ishape%phi_ic,ishape%psi_ic,ishape%a,ishape%c
                   write(*,*) ms%mass
                   stop
                end if
             else
                q=q_cs
             end if
           else
              q=q_ic
           endif
       else
          q = A_e/A_c
       end if 

    else if( ishape%sh_type == 5 ) then
!!c       if(mode==2) then
!!c          alpha=0.7_PS
!!c       else
          alpha = ms%semi_c/ms%semi_a
!!c       end if     
   
       ! +++ calculate the total surface area +++
       OMEGA = get_total_sfc_area( ishape, ms%semi_a, ms%semi_c, mode)
       
       ! +++ calculate the area ratio +++
       ! ++++++ calculate the effective cross section ++++++
       A_e = get_effect_area( ishape, ishape%sh_type,ms%semi_a, ms%semi_c, mode)
       
       ! ++++++ calculate the circumscribed cross-sectional area ++++++
       A_c = get_circum_area( ishape, ishape%sh_type,ms%semi_a, ms%semi_c, mode)

       ! +++ calculate the characteristic length based on
       !     perimeter and total surface area +++
       P = get_perimeter( ishape, ms%semi_a, ms%semi_c, mode)

       ! +++ calculate characteristic length +++
       char_len = 2.0*sqrt(A_c/PI)

!!c       if(mode==2) then
!!c          ! from Bohm (1991), which assume that 10% of diameter wide is
!!c          ! covered by 50%.
!!c          q = 0.82_PS
       
       if(mode==2) then
          ! +++ calculate the area ratio of ice crystals +++
          ! ++++++ calculate the effective cross section ++++++
          if( ishape%habit>=5 ) then
             ! use Heymsfield et al (2002)'s side planes
             q_ic=max(0.18_PS,min(0.88_PS,(ms%mass(imc)/ms%con/ishape%V_ic/0.35_PS)&
                  **(1.0/2.34_PS)))
          else
             q_ic=get_effect_area( ishape,1,ms%a_len, ms%c_len, mode)/&
                  get_circum_area( ishape,1,ms%a_len, ms%c_len, mode)
          end if
!!c          q_cs=ms%mean_mass/den_i/ishape%V_cs
!!c          q_cs=(ms%mass(imc)+ms%mass(ima)+ms%mass(imr))/ms%con/den_i/ishape%V_cs
          q_cs=(ms%mass(imc)+ms%mass(ima)+ms%mass(imr))/ms%con/den_i/&
                get_vip(ishape%is_mod(2),ishape%phi_cs,ishape%semi_aip)

          if(q_cs>1.0e-30_PS) then
             mass_q_grp=(ms%mass(imc)+ms%mass(ima))/ms%con*mrat_grp
             mass1=(ms%mass(imc)+ms%mass(ima)+ms%mass(imr))/ms%con

             if(mass1<mass_q_grp) then
!!c                q=10.0**((log10(q_ic)-log10(q_cs))/(log10(ms%mass(imc)/ms%con)-log10(mass_q_grp))*&
!!c                     (log10(mass1)-log10(mass_q_grp))+log10(q_cs))
                q=q_cs*(mass1/mass_q_grp)**&
                 (log10(q_ic/q_cs)/log10(ms%mass(imc)/ms%con/mass_q_grp))

                if(q<0.0.or.q>2.0.or.(q<0.0.and.q>2.0)) then
                   write(*,*) "q is not right for grp",q,q_ic,q_cs,ms%mass(imc),ms%con,ms%mean_mass
                   write(*,*) ms%mass
                   stop
                end if
             else
                q=q_cs
             end if
           else
              q=q_ic
           end if
       else
          q = A_e/A_c
       end if
       
    end if
  end subroutine cal_some_vel


  subroutine cal_mass_comp(token,level,th_var,eps_ap0, ms)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate 
    ! 1. mass of rime component and aggregate component
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Mass_Bin), intent(inout)   :: ms
!!c    type (Ice_Shape), intent(in)    :: is
    type (Thermo_Var), intent(in)    :: th_var
    integer, intent(in)    ::  level,token
    real(PS),intent(in) ::  eps_ap0
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
    ! temperature at triple point 
    real (PS), parameter          :: T_0 = 273.16_PS

    real(PS) :: mass_left,mass_stot,oldmassc(mxnmasscomp+1),rat
    integer :: ifix

    ! initialize 
    ms%inmlt=0
    ms%inevp=0
    ms%eps_map=eps_ap0

    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) then
       ms%mass=0.0_PS
       return
    endif
    if(token==2) then
       if(level<=3) then
!!c          ms%mass(ima)=max(ms%mass(imt)-ms%mass(imr)-ms%mass(imc),0.0_PS)
          if( th_var%T>=T_0 ) ms%inmlt=1
       elseif(level<=5) then
!!c          ms%mass(ima)=max(ms%mass(imt)-ms%mass(imr)-ms%mass(imc),0.0_PS)
!!c          if(ms%mass(imat)<1.0e-30_PS) then
!!c             ms%mass(imat)=min(max(1.0e-30_PS,ms%mass(imat)),ms%mass(imt))
!!c          end if
          ms%mass(imai)=max(ms%mass(imat)-ms%mass(imas),0.0_PS)
          ms%eps_map=max(0.0_PS,min(1.0_PS,ms%mass(imas)/ms%mass(imat)))
          if( th_var%T>=T_0 ) ms%inmlt=1
       else
!!c          ms%mass(ima)=max(ms%mass(imt)-ms%mass(imr)-ms%mass(imc)-ms%mass(imw),0.0_PS)
!!c          if(ms%mass(imat)<1.0e-30_PS) then
!!c             ms%mass(imat)=min(max(1.0e-30_PS,ms%mass(imat)),ms%mass(imt))
!!c          end if

          ms%mass(imai)=max(ms%mass(imat)-ms%mass(imas),0.0_PS)
          ms%eps_map=max(0.0_PS,min(1.0_PS,ms%mass(imas)/ms%mass(imat)))
          if( ms%mass(imw)/ms%mass(imt)>=fmice_lmt ) then
             ms%inmlt=2
! <<< 2015/03 T. Hashino added for KiD test
!tmp          elseif( th_var%T>=T_0+5.0) then
!tmp             ms%inmlt=2
! >>> 2015/03 T. Hashino added for KiD test
          elseif( ms%mass(imw)/ms%mass(imt)>=fmliq_lmt ) then
             ms%inmlt=1
          end if
       end if
!!c       if(ms%eps_map>=1.0e-6) then
!!c          write(*,*) "cal_mass_comp > sol frac larger than 1.0e-6 in ice"
!!c          write(*,'(8ES15.6)') ms%mass(imt),ms%mass(imat),ms%mass(imas)
!!c       end if

       ! maintain the relative relation of predicted mass component
       mass_left=ms%mass(imt)-ms%mass(imc)
       if(mass_left<0.0_PS) then
!!c          write(*,*) "cal_masscomp: m4>m1",ms%mass
!!c          write(*,*) "con, meanmass",ms%con,ms%mass(imt)/ms%con
          ms%mass(imt)=ms%mass(imc)
          mass_left=0.0_PS
       endif
       mass_stot=ms%mass(imr)+ms%mass(ima)+ms%mass(imw)
       oldmassc=ms%mass
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
          write(*,*) "fixing failed:",ms%mass(1:8)
          write(*,*) "old masscomp",oldmassc
          write(*,*) "ifix,rat",ifix,rat
       endif

    elseif(token==1) then
       if(level>=4) then
!!c          if(ms%mass(rmat)<1.0e-30_PS) then
!!c             ms%mass(rmat)=min(max(1.0e-30_PS,ms%mass(rmat)),ms%mass(rmt))
!!c          end if
          ms%mass(rmai)=max(ms%mass(rmat)-ms%mass(rmas),0.0_PS)
          ms%eps_map=max(0.0_PS,min(1.0_PS,ms%mass(rmas)/ms%mass(rmat)))
       end if
!!c       if(ms%eps_map<1.0e-6) then
!!c          write(*,*) "cal_mass_comp > sol frac less than 1.0e-6 in liq"
!!c          write(*,'(8ES15.6)') ms%mass(rmt),ms%mass(rmat),ms%mass(rmas)
!!c       end if
    elseif(token==3) then
       if(level>=4) then
          ms%eps_map=max(0.0_PS,min(1.0_PS,ms%mass(ams)/ms%mass(amt)))
          ms%mass(ami)=max(ms%mass(amt)-ms%mass(ams),0.0_PS)
       end if
    end if


!!c    if( ms%mass(1) <= min_mass*ms%con) then
!!c       ms%mass(4)=ms%mass(4)+ms%mass(3)
!!c       ms%mass(3)=0.0_PS
!!c       write(*,*) "cal_mass_comp > mass_c per perticle fixed for agg",ms%mass(1:4)/ms%con
!!c       ms%mass(4)=min_mass*ms%con
!!c       ms%mass(1)=max(ms%mass(1),ms%mass(4))
!!c       if( ms%mass(2)+ms%mass(3) > 1.0e-20 ) then
!!c          ratio=(ms%mass(1)-ms%mass(4))/(ms%mass(2)+ms%mass(3))
!!c          ms%mass(2)=ratio*ms%mass(2)
!!c          ms%mass(3)=ratio*ms%mass(3)
!!c       else
!!c          ms%mass(2)=0.0
!!c          ms%mass(3)=0.0
!!c       end if
!!c    end if
    
  end subroutine cal_mass_comp

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
    ! temperature at triple point 
    real (PS), parameter          :: T_0 = 273.16_PS

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

  function rime_density() result(den_r)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate
    ! 1. density of rime deposits or mean of density
    !
    ! NOTE: The radius, and the Reynolds number of cloud droplet has to be 
    !       diagnosed beforehand.
    !       So is the graupel.
    !       
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real (PS)      :: den_r
    den_r = 0.454_PS
!!!!!!!!!!! under construction
  end function rime_density

  subroutine det_sh_type( ms, is, level)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate 
    ! 1. volume of rime component and aggregate component
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Mass_Bin), intent(inout)   :: ms
    type (Ice_Shape), intent(inout)    :: is
    integer, intent(in)    ::  level
    real(PS)  :: m_ice

    if( ms%mass(imr) > ms%mass(imt)/100.0_PS ) then
       if( ms%mass(ima) <= ms%mass(imc) ) then
          if( ms%mass(imr) <= ms%mass(imc) ) then
             ! rimed crystal
             is%sh_type = 2
          else
             ! graupel
             is%sh_type = 5
          end if
       else
          if( ms%mass(imr) <= ms%mass(ima) ) then
             ! rimed aggregate
             is%sh_type = 4
          else
             ! graupel
             is%sh_type = 5
          end if
       end if
    else if( ms%mass(ima) > ms%mass(imc) ) then
       is%sh_type = 3             
    else
       ! pristine crystal
       is%sh_type = 1
    end if
  end subroutine det_sh_type

  subroutine det_sh_type2( ms, is, level)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate 
    ! 1. volume of rime component and aggregate component
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Mass_Bin), intent(inout)   :: ms
    type (Ice_Shape), intent(inout)    :: is
    integer, intent(in)    ::  level
    real(PS)  :: m_ice

    if( ms%mass(imr) > ms%mass(imt)/100.0 ) then
       if( ms%mass(ima) <= ms%mass(imt)/100.0 ) then
          if( ms%mass(imr) <= ms%mass(imc) ) then
             ! rimed crystal
             is%sh_type = 2
          else
             ! graupel
             is%sh_type = 5
          end if
       else
          if( ms%mass(imr) <= ms%mass(ima) ) then
             ! rimed aggregate
             is%sh_type = 4
          else
             ! graupel
             is%sh_type = 5
          end if
       end if
    else if( ms%mass(ima) > ms%mass(imt)/100.0 ) then
       is%sh_type = 3             
    else
       ! pristine crystal
       is%sh_type = 1
    end if
  end subroutine det_sh_type2

  subroutine det_sh_type3( ms, is, level)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate 
    ! 1. volume of rime component and aggregate component
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Mass_Bin), intent(inout)   :: ms
    type (Ice_Shape), intent(inout)    :: is
    integer, intent(in)    ::  level
    real(PS)  :: m_ice
    real(PS),parameter :: frac=0.1_PS
!!c    real(PS),parameter :: frac=0.01_PS

    if( ms%mass(imr) > ms%mass(imt)*frac ) then
       if( ms%mass(ima) <= ms%mass(imt)*frac ) then
          if( ms%mass(imr) <= ms%mass(imc) ) then
             ! rimed crystal
             is%sh_type = 2
          else
             ! graupel
             is%sh_type = 5
          end if
       else
          if( ms%mass(imr) <= ms%mass(ima) ) then
             ! rimed aggregate
             is%sh_type = 4
          else
             ! graupel
             is%sh_type = 5
          end if
       end if
    else if( ms%mass(ima) > ms%mass(imt)*frac ) then
       is%sh_type = 3             
    else
       ! pristine crystal
       is%sh_type = 1
    end if
  end subroutine det_sh_type3

  subroutine det_sh_type4( ms, is, level)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate 
    ! 1. volume of rime component and aggregate component
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Mass_Bin), intent(inout)   :: ms
    type (Ice_Shape), intent(inout)    :: is
    integer, intent(in)    ::  level
    real(PS)  :: m_ice
    real(PS),parameter :: frac=0.1_PS
!!c    real(PS),parameter :: frac=0.01_PS

    if( ms%mass(imr)+ms%mass(imw) > ms%mass(imt)*frac ) then
       if( ms%mass(ima) <= ms%mass(imt)*frac ) then
          if( ms%mass(imr)+ms%mass(imw) <= ms%mass(imc) ) then
             ! rimed crystal
             is%sh_type = 2
          else
             ! graupel
             is%sh_type = 5
          end if
       else
          if( ms%mass(imr)+ms%mass(imw) <= ms%mass(ima) ) then
             ! rimed aggregate
             is%sh_type = 4
          else
             ! graupel
             is%sh_type = 5
          end if
       end if
    else if( ms%mass(ima) > ms%mass(imt)*frac ) then
       is%sh_type = 3             
    else
       ! pristine crystal
       is%sh_type = 1
    end if
  end subroutine det_sh_type4

  subroutine diag_sh_type_v4(level, ratio_Mp, sh_type)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! diagnose the type of ice particle
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, intent(in)    ::  level
    integer, intent(inout)    ::  sh_type
    ! ratio of each component mass to the total mass in a bin
    ! ratio_mp(1) : rime
    ! ratio_mp(2) : aggregation
    ! ratio_mp(3) : ice crystal
    ! ratio_Mp(7) : melt mass
    real(ps), pointer, dimension(:)      :: ratio_Mp

    real(PS)  :: m_ice
    real(PS),parameter :: frac=0.1_PS
!!c    real(PS),parameter :: frac=0.01_PS

    if( ratio_Mp(1)+ratio_Mp(7) > frac ) then
       if( ratio_Mp(2) <= frac ) then
          if( ratio_Mp(1)+ratio_Mp(7) <= ratio_Mp(3) ) then
             ! rimed crystal
             sh_type = 2
          else
             ! graupel
             sh_type = 5
          end if
       else
          if(ratio_Mp(1)+ratio_Mp(7) <= ratio_Mp(2) ) then
             ! rimed aggregate
             sh_type = 4
          else
             ! graupel
             sh_type = 5
          end if
       end if
    else if( ratio_Mp(2) > frac ) then
       sh_type = 3             
    else
       ! pristine crystal
       sh_type = 1
    end if
  end subroutine diag_sh_type_v4
  subroutine diag_sh_type_v5(level, ratio_Mp, sh_type, is_mod)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! diagnose the type of ice particle
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

    real(PS)  :: m_ice
    real(PS),parameter :: frac=0.1_PS
!!c    real(PS),parameter :: frac=0.01_PS

    ! default is a cylinder
    is_mod(1)=1
    is_mod(2)=1

    if( ratio_Mp(1)+ratio_Mp(7) > frac ) then
       if( ratio_Mp(2) <= frac ) then
          if( ratio_Mp(1)+ratio_Mp(7) <= ratio_Mp(3) ) then
             ! rimed crystal
             sh_type = 2
          else
             ! graupel
             sh_type = 5
             is_mod(2)=2
          end if
       else
          if(ratio_Mp(1)+ratio_Mp(7) <= ratio_Mp(2) ) then
             ! rimed aggregate
             sh_type = 4
          else
             ! graupel
             sh_type = 5
             is_mod(2)=2
          end if
       end if
    else if( ratio_Mp(2) > frac ) then
       sh_type = 3             
    else
       ! pristine crystal
       sh_type = 1
    end if
  end subroutine diag_sh_type_v5
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

  subroutine cal_cs_spheroid(level,alen,clen,rlen,elen,V_ic,semi_a_i,semi_c_i,habit,&
       sh_type,semi_a,semi_c,phi_cs,V_cs,&
       mean_mass,m_ic,m_ag)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate 
    ! 1. semi-axis lengths of circumscribing sphere
    ! 2. axis ratio of it.
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, intent(in)    ::  level,habit,sh_type
    real(PS),intent(in) :: alen,clen,rlen,elen,V_ic,semi_a_i,semi_c_i,m_ic,m_ag,mean_mass
    real(PS),intent(inout) :: semi_a,semi_c,phi_cs,V_cs
    ! 4 pi /3, pi/6
    real(PS),parameter :: coef4pi3=4.188790205_PS,coefpi6=0.523598776

    ! assumption on aspect ratio
!!c    real(PS),parameter :: phi_agg=0.25_PS,mass_phi_agg=1.0e-4_PS
!!c    real(PS),parameter :: phi_agg=0.33333333333_PS,mass_phi_agg=1.0e-4_PS
    real(PS),parameter :: phi_agg=0.35_PS, mrat_agg=10.0_PS
!!c    real(PS),parameter :: phi_grp=1.0_PS,mrat_grp=10.0_PS
    real(PS),parameter :: phi_grp=0.8_PS,mrat_grp=10.0_PS ! sense test grp=10
    real(PS) :: mass_phi_agg,mass_phi_grp,mass1
    ! assumption for planar and irregular polycrystals
    real(PS),parameter :: spx_p=0.25


    real(PS) :: phi_ic,phi,dia

    select case(sh_type)
    case (1)
       ! pristine crystal
       ! assume the axis ratio is kept the same as pristine crystal inside.
       semi_a=semi_a_i
       semi_c=semi_c_i
       phi_cs=semi_c/semi_a

       V_cs=V_ic


    case (2)
       ! rimed crystal
       ! assume the axis ratio is kept the same as pristine crystal inside.
       semi_a=semi_a_i
       semi_c=semi_c_i
       phi_cs=semi_c/semi_a

       V_cs=V_ic



    case (3)
       ! 3. unrimed aggregates
       !
       !  assume cylinder

       dia=10.0_PS*(V_cs/coefpi6)**(1.0/3.0)
       mass_phi_agg=m_ic*mrat_agg

       mass1=m_ic+m_ag

       if(mass1<mass_phi_agg) then
          if(habit<=3) then
             phi_ic=clen/alen
          else
             phi_ic=spx_p
          end if
          if(phi_agg<=0.0.or.phi_ic<=0.0.or.m_ic<=0.0.or.&
               mass1<=0.0) then
             write(*,15) phi_ic,m_ic,mass1
15           format("Something is wrong in cs_spheroid for agg:",4ES15.6)
             write(*,'("mass1,m_ic",8ES15.6)') mass1,m_ic
             write(*,'("alen,clen,v_cs",8ES15.6)') alen,clen
             stop
          end if


          phi=(log10(phi_agg)-log10(phi_ic))/(log10(mass_phi_agg)-log10(m_ic))*&
               (log10(mass1)-log10(mass_phi_agg))+log10(phi_agg)
          phi=10.0**phi
       else
          phi=phi_agg
       end if
       phi_cs=phi

       semi_a=max(dia/20.0_PS/sqrt(1.0_PS+phi_cs**2.0),semi_a_i)
       semi_c=max(phi_cs*semi_a,semi_c_i)
       V_cs=coef4pi3*(semi_a**2.0+semi_c**2.0)**1.5
       phi_cs=semi_c/semi_a

    case (4)
       ! 4. rimed aggregates
       !
       !  assume cylinder
       !
       !  assume cylinder

       dia=10.0_PS*(V_cs/coefpi6)**(1.0/3.0)

       mass_phi_agg=m_ic*mrat_agg

       mass1=m_ic+m_ag

       if(mass1<mass_phi_agg) then
          if(habit<=3) then
             phi_ic=clen/alen
          else
             phi_ic=spx_p
          end if
          if(phi_agg<=0.0.or.phi_ic<=0.0.or.m_ic<=0.0.or.&
               mass1<=0.0) then
             write(*,16) phi_ic,m_ic,mass1
16           format("Something is wrong in cs_spheroid for rimed agg:",4ES15.6)
             write(*,'("mass1,m_ic",8ES15.6)') mass1,m_ic
             write(*,'("alen,clen,v_cs",8ES15.6)') alen,clen
             stop
          end if


          phi=(log10(phi_agg)-log10(phi_ic))/(log10(mass_phi_agg)-log10(m_ic))*&
               (log10(mass1)-log10(mass_phi_agg))+log10(phi_agg)
          phi=10.0**phi
       else
          phi=phi_agg
       end if
       phi_cs=phi

       semi_a=max(dia/20.0_PS/sqrt(1.0_PS+phi_cs**2.0),semi_a_i)
       semi_c=max(phi_cs*semi_a,semi_c_i)
       V_cs=coef4pi3*(semi_a**2.0+semi_c**2.0)**1.5
       phi_cs=semi_c/semi_a

    case (5)
       ! graupel
       !
       !    assume that a spheroid is the shape of rimed aggregates.
       !
       dia=10.0_PS*(V_cs/coefpi6)**(1.0/3.0) ! in mm
!!c       write(*,*) "bef get phi 5"

       mass_phi_grp=(m_ic+m_ag)*mrat_grp

       if(mean_mass<mass_phi_grp) then
          if(habit<=3) then
             phi_ic=clen/alen
          else
             phi_ic=spx_p
          end if

          if(phi_grp<=0.0.or.phi_ic<=0.0.or.m_ic<=0.0.or.&
               mean_mass<=0.0) then
             write(*,17) phi_ic,m_ic,mean_mass
17           format("Something is wrong in cs_spheroid at graupels:",4ES15.6)
             write(*,'("mean_mass,m_ic",8ES15.6)') mean_mass,m_ic
             write(*,'("alen,clen,v_cs",8ES15.6)') alen,clen
             stop
          end if

          phi=(log10(phi_grp)-log10(phi_ic))/(log10(mass_phi_grp)-log10(m_ic))*&
               (log10(mean_mass)-log10(mass_phi_grp))+log10(phi_grp)
          phi=10.0**phi
       else
          phi=phi_grp
       end if
!!c       write(*,*) "aft get phi 5"

       phi_cs=phi
       if(phi_cs<=1.0) then
          ! maximum dimension corresponds to the semi-a length
          semi_a=max(dia/20.0_PS,semi_a_i)
          semi_c=max(phi_cs*semi_a,semi_c_i)
       else
          ! maximum dimension corresponds to the semi-c length
          semi_c=max(dia/20.0_PS,semi_c_i)
          semi_a=max(semi_c/phi_cs,semi_a_i)
       end if
       V_cs=coef4pi3*max(semi_a,semi_c)**3.0
       phi_cs=semi_c/semi_a

    case default
       write(*,*) "ERROR! cal_cs_sheroid."
       stop
    end select
  end subroutine cal_cs_spheroid

  subroutine cal_cs_spheroid2(level,alen,clen,rlen,elen,V_ic,semi_a_i,semi_c_i,habit,&
       sh_type,semi_aip,semi_cip,phi_cs,V_cs,&
       semi_a,semi_c,V_csw,&
       mean_mass,m_ic,m_ag,m_mlt,&
       is_mod,iswitch)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate 
    ! 1. semi-axis lengths of circumscribing sphere
    ! 2. axis ratio of it.
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, intent(in)    ::  level,habit,sh_type,iswitch
    real(PS),intent(in) :: alen,clen,rlen,elen,V_ic,semi_a_i,semi_c_i,m_ic,m_ag,m_mlt,mean_mass
    real(PS),intent(inout) :: semi_aip,semi_cip,phi_cs,V_cs,semi_a,semi_c,V_csw
    integer,dimension(2) :: is_mod

    ! assumption on aspect ratio
!!c    real(PS),parameter :: phi_agg=0.25_PS,mass_phi_agg=1.0e-4_PS
!!c    real(PS),parameter :: phi_agg=0.33333333333_PS,mass_phi_agg=1.0e-4_PS
    real(PS),parameter :: phi_agg=0.35_PS, mrat_agg=10.0_PS
!!c    real(PS),parameter :: phi_grp=1.0_PS,mrat_grp=10.0_PS
    real(PS),parameter :: phi_grp=0.8_PS,mrat_grp=10.0_PS ! sense test grp=10
    real(PS) :: mass_phi_agg,mass_phi_grp,mass1
    ! assumption for planar and irregular polycrystals
    real(PS),parameter :: spx_p=0.25

    real(PS), parameter         :: den_w = 1.0

    real(PS) :: phi_ic,phi,dia,V_ip,V_space,ak3

    select case(sh_type)
    case (1)
       ! pristine crystal
       ! assume the axis ratio is kept the same as pristine crystal inside.
       semi_aip=semi_a_i
       semi_cip=semi_c_i
       phi_cs=semi_cip/semi_aip
       V_cs=V_ic

!!c       V_ip=coef2p*phi_cs*semi_aip**3.0
       V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
       V_space=max(0.0_PS,V_ip-(mean_mass-m_mlt)/den_i)
       ak3=1.0_PS+max(m_mlt/den_w-V_space,0.0_PS)/V_ip
       V_csw=ak3*V_cs
       semi_a=ak3**(1.0/3.0)*semi_aip
       semi_c=semi_a*phi_cs
!!c
!!c       semi_a=semi_aip
!!c       semi_c=semi_cip
!!c       V_csw=V_cs

    case (2)
       ! rimed crystal
       !
       V_cs=max(V_cs,V_ic)
       if(iswitch==0) then
          ! assume the cylinder volume
!!c          semi_aip=(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0_PS+phi_cs**2.0)
!!c          semi_cip=semi_aip*phi_cs
          call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
       elseif(iswitch==1) then
          ! assume the same growth as in riming scheme.
          if(habit<=3) then
             phi_ic=clen/alen
          else
             phi_ic=spx_p
          end if

          if(semi_a_i>semi_c_i) then
             if(is_mod(2)==1) then
                phi_cs=sqrt((V_cs/V_ic)**(2.0/3.0)*(1.0_PS+phi_ic**2)-1.0_PS)
             elseif(is_mod(2)==2) then
                phi_cs=phi_ic*(V_cs/V_ic)
             end if
             semi_aip=semi_a_i
             semi_cip=semi_aip*phi_cs
          elseif(semi_a_i<semi_c_i) then
             if(is_mod(2)==1) then
                phi_cs=1.0/sqrt((V_cs/V_ic)**(2.0/3.0)*(1.0_PS+phi_ic**(-2))-1.0_PS)
             elseif(is_mod(2)==2) then
                phi_cs=phi_ic*sqrt(V_ic/V_cs)
             end if
             semi_cip=semi_c_i
             semi_aip=semi_cip/phi_cs
          else
             phi_cs=phi_ic
             call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
          end if

!!c          if(semi_a_i>semi_c_i) then
!!c             phi_cs=sqrt((V_cs/V_ic)**(2.0/3.0)*(1.0_PS+phi_ic**2)-1.0_PS)
!!c             semi_aip=semi_a_i
!!c             semi_cip=semi_aip*phi_cs
!!c          elseif(semi_a_i<semi_c_i) then
!!c             phi_cs=1.0/sqrt((V_cs/V_ic)**(2.0/3.0)*(1.0_PS+phi_ic**(-2))-1.0_PS)
!!c             semi_cip=semi_c_i
!!c             semi_aip=semi_cip/phi_cs
!!c          else
!!c             phi_cs=phi_ic
!!c             ! assume the cylinder volume
!!c             semi_aip=(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0_PS+phi_cs**2.0)
!!c             semi_cip=semi_aip*phi_cs
!!c          end if

!!c
!!c       ! assume the axis ratio is kept the same as pristine crystal inside.
!!c       semi_aip=semi_a_i
!!c       semi_cip=semi_c_i
!!c       phi_cs=semi_cip/semi_aip
!!c       V_cs=V_ic
       end if
!!c       semi_aip=max(semi_aip,semi_a_i)
!!c       semi_cip=max(semi_cip,semi_c_i)
       phi_cs=semi_cip/semi_aip

       ! assume the cylinder volume
!!c       v_cs=coef4pi3*semi_aip**3.0*(1.0_PS+phi_cs**2.0)**1.5
       V_cs=get_vcs(is_mod(2),phi_cs,semi_aip)


!!c       V_ip=coef2p*phi_cs*semi_aip**3.0
       V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
       V_space=max(0.0_PS,V_ip-(mean_mass-m_mlt)/den_i)
       ak3=1.0_PS+max(m_mlt/den_w-V_space,0.0_PS)/V_ip
       V_csw=ak3*V_cs
       semi_a=ak3**(1.0/3.0)*semi_aip
       semi_c=semi_a*phi_cs
!!c
!!c       semi_a=semi_aip
!!c       semi_c=semi_cip
!!c       V_csw=V_cs

    case (3)
       ! 3. unrimed aggregates
       !
       V_cs=max(V_cs,V_ic)
       if(iswitch==1) then
          mass_phi_agg=m_ic*mrat_agg
          
          mass1=m_ic+m_ag
          
          if(mass1<mass_phi_agg) then
             if(habit<=3) then
                phi_ic=clen/alen
             else
                phi_ic=spx_p
             end if
             if(phi_agg<=0.0.or.phi_ic<=0.0.or.m_ic<=0.0.or.&
                  mass1<=0.0) then
                write(*,15) phi_ic,m_ic,mass1
15              format("Something is wrong in cs_spheroid for agg:",4ES15.6)
                write(*,'("mass1,m_ic",8ES15.6)') mass1,m_ic
                write(*,'("alen,clen,v_cs",8ES15.6)') alen,clen
                stop
             end if
             
             
             phi=(log10(phi_agg)-log10(phi_ic))/(log10(mass_phi_agg)-log10(m_ic))*&
                  (log10(mass1)-log10(mass_phi_agg))+log10(phi_agg)
             phi=10.0**phi
          else
             phi=phi_agg
          end if
          phi_cs=phi
       end if

       ! assume the cylinder volume
!!c       semi_aip=max(semi_a_i,(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0_PS+phi_cs**2.0))
!!c       semi_cip=max(semi_c_i,semi_aip*phi_cs)
       call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
!!c       semi_aip=max(semi_aip,semi_a_i)
!!c       semi_cip=max(semi_cip,semi_c_i)
       phi_cs=semi_cip/semi_aip

       ! assume the cylinder volume
!!c       v_cs=coef4pi3*semi_aip**3.0*(1.0_PS+phi_cs**2.0)**1.5
       V_cs=get_vcs(is_mod(2),phi_cs,semi_aip)


!!c       V_ip=coef2p*phi_cs*semi_aip**3.0
       V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
       V_space=max(0.0_PS,V_ip-(mean_mass-m_mlt)/den_i)
       ak3=1.0_PS+max(m_mlt/den_w-V_space,0.0_PS)/V_ip
       V_csw=ak3*V_cs
       semi_a=ak3**(1.0/3.0)*semi_aip
       semi_c=semi_a*phi_cs
!!c
!!c       semi_a=semi_aip
!!c       semi_c=semi_cip
!!c       V_csw=V_cs

    case (4)
       ! 4. rimed aggregates
       !
       !  assume cylinder
       !
       !  assume cylinder

       V_cs=max(V_cs,V_ic)
       if(iswitch==1) then
          mass_phi_agg=m_ic*mrat_agg
          
          mass1=m_ic+m_ag
          
          if(mass1<mass_phi_agg) then
             if(habit<=3) then
                phi_ic=clen/alen
             else
                phi_ic=spx_p
             end if
             if(phi_agg<=0.0.or.phi_ic<=0.0.or.m_ic<=0.0.or.&
                  mass1<=0.0) then
                write(*,16) phi_ic,m_ic,mass1
16              format("Something is wrong in cs_spheroid for rimed agg:",4ES15.6)
                write(*,'("mass1,m_ic",8ES15.6)') mass1,m_ic
                write(*,'("alen,clen,v_cs",8ES15.6)') alen,clen
                stop
             end if
             
             
             phi=(log10(phi_agg)-log10(phi_ic))/(log10(mass_phi_agg)-log10(m_ic))*&
                  (log10(mass1)-log10(mass_phi_agg))+log10(phi_agg)
             phi=10.0**phi
          else
             phi=phi_agg
          end if
          phi_cs=phi
       end if

       ! assume the cylinder volume
!!c       semi_aip=max(semi_a_i,(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0_PS+phi_cs**2.0))
!!c       semi_cip=max(semi_c_i,semi_aip*phi_cs)

       call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
!!c       semi_aip=max(semi_aip,semi_a_i)
!!c       semi_cip=max(semi_cip,semi_c_i)

       phi_cs=semi_cip/semi_aip

       ! assume the cylinder volume
!!c       v_cs=coef4pi3*semi_aip**3.0*(1.0_PS+phi_cs**2.0)**1.5
       V_cs=get_vcs(is_mod(2),phi_cs,semi_aip)



!!c       V_ip=coef2p*phi_cs*semi_aip**3.0
       V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
       V_space=max(0.0_PS,V_ip-(mean_mass-m_mlt)/den_i)
       ak3=1.0_PS+max(m_mlt/den_w-V_space,0.0_PS)/V_ip
       V_csw=ak3*V_cs
       semi_a=ak3**(1.0/3.0)*semi_aip
       semi_c=semi_a*phi_cs
!!c
!!c       semi_a=semi_aip
!!c       semi_c=semi_cip
!!c       V_csw=V_cs

    case (5)
       ! graupel
       !
       !    assume that a spheroid is the shape of rimed aggregates.
       !

       V_cs=max(V_cs,V_ic)
!!c       write(*,*) "bef get phi 5"

       if(iswitch==1) then          
          mass_phi_grp=(m_ic+m_ag)*mrat_grp
          
          if(mean_mass<mass_phi_grp) then
             if(habit<=3) then
                phi_ic=clen/alen
             else
                phi_ic=spx_p
             end if
             
             if(phi_grp<=0.0.or.phi_ic<=0.0.or.m_ic<=0.0.or.&
                  mean_mass<=0.0) then
                write(*,17) phi_ic,m_ic,mean_mass
17              format("Something is wrong in cs_spheroid at graupels:",4ES15.6)
                write(*,'("mean_mass,m_ic",8ES15.6)') mean_mass,m_ic
                write(*,'("alen,clen,v_cs",8ES15.6)') alen,clen
                stop
             end if
             
             phi=(log10(phi_grp)-log10(phi_ic))/(log10(mass_phi_grp)-log10(m_ic))*&
                  (log10(mean_mass)-log10(mass_phi_grp))+log10(phi_grp)
             phi=10.0**phi
          else
             phi=phi_grp
          end if
!!c       write(*,*) "aft get phi 5"

          phi_cs=phi

       end if

       ! assume the spheroidal volume
!!c       if(phi_cs<1.0) then
!!c          semi_aip=max(semi_a_i,(V_cs/coef4pi3)**(1.0/3.0))
!!c          semi_cip=max(semi_c_i,semi_aip*phi_cs)
!!c       else
!!c          semi_cip=max(semi_c_i,(V_cs/coef4pi3)**(1.0/3.0))
!!c          semi_aip=max(semi_a_i,semi_cip/phi_cs)
!!c       end if

       call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
!!c       semi_aip=max(semi_aip,semi_a_i)
!!c       semi_cip=max(semi_cip,semi_c_i)

       phi_cs=semi_cip/semi_aip

!!c       ! assume the spheroidal volume
!!c       if(phi_cs<1.0_PS) then
!!c          v_cs=coef4pi3*semi_aip**3.0
!!c       else
!!c          v_cs=coef4pi3*(semi_aip*phi_cs)**3.0
!!c       end if
       V_cs=get_vcs(is_mod(2),phi_cs,semi_aip)


!!c       V_ip=coef4pi3*phi_cs*semi_aip**3.0
       V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
       V_space=max(0.0_PS,V_ip-(mean_mass-m_mlt)/den_i)
       ak3=1.0_PS+max(m_mlt/den_w-V_space,0.0_PS)/V_ip
       V_csw=ak3*V_cs
       semi_a=ak3**(1.0/3.0)*semi_aip
       semi_c=semi_a*phi_cs
!!c
!!c       semi_a=semi_aip
!!c       semi_c=semi_cip
!!c       V_csw=V_cs



    case default
       write(*,*) "ERROR! cal_cs_sheroid."
       stop
    end select
  end subroutine cal_cs_spheroid2

  subroutine cal_cs_spheroid3(level,alen,clen,rlen,elen,V_ic,&
       semi_a_i,semi_c_i,habit,&
       sh_type,semi_aip,semi_cip,phi_cs,V_cs,&
       semi_a,semi_c,V_csw,&
       mean_mass,m_ic,m_ag,&
       m_rm,m_mlt,&
       is_mod,iswitch)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate 
    ! 1. semi-axis lengths of circumscribing sphere
    ! 2. axis ratio of it.
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, intent(in)    ::  level,habit,sh_type,iswitch
    real(PS),intent(in) :: alen,clen,rlen,elen,V_ic,semi_a_i,semi_c_i,m_ic,m_ag,m_rm,m_mlt,mean_mass
    real(PS),intent(inout) :: semi_aip,semi_cip,phi_cs,V_cs,semi_a,semi_c,V_csw
    integer,dimension(2) :: is_mod

    ! assumption on aspect ratio
!!c    real(PS),parameter :: phi_agg=0.25_PS,mass_phi_agg=1.0e-4_PS
!!c    real(PS),parameter :: phi_agg=0.33333333333_PS,mass_phi_agg=1.0e-4_PS
    real(PS),parameter :: phi_agg=1.0/3.0, mrat_agg=10.0_PS
!!c    real(PS),parameter :: phi_grp=1.0_PS,mrat_grp=10.0_PS
    real(PS),parameter :: phi_grp=0.8_PS,mrat_grp=10.0_PS ! sense test grp=10
    real(PS) :: mass_phi_agg,mass_phi_grp,mass1
    ! assumption for planar and irregular polycrystals
    real(PS),parameter :: spx_p=0.25

    real(PS), parameter         :: den_w = 1.0

    ! minimum (maximum) volume of circumscribing sphere corresponds to 1 um (10cm) radius
    real(PS),parameter :: V_csmin=4.18879020478639e-12,V_csmax=4.18879020478639e+3

    real(PS) :: phi_ic,phi,dia,V_ip,V_space,ak3

    select case(sh_type)
    case (1)
       ! pristine crystal
       if(is_mod(2)==1) then
          ! assume the axis ratio is kept the same as pristine crystal inside.
          semi_aip=semi_a_i
          semi_cip=semi_c_i
          phi_cs=semi_cip/semi_aip
          V_cs=V_ic
       else
          ! ice sphere
          phi_cs=1.0_PS
          V_cs=max(V_csmin,(mean_mass-m_mlt)/den_i)
          call cal_semiac_ip(is_mod(2),phi_cs,V_cs,semi_aip,semi_cip)
       endif

!!c          V_ip=coef2p*phi_cs*semi_aip**3.0
       V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
       V_space=max(0.0_PS,V_ip-(mean_mass-m_mlt)/den_i)
       ak3=1.0_PS+max(m_mlt/den_w-V_space,0.0_PS)/V_ip
       V_csw=ak3*V_cs
       semi_a=ak3**(1.0/3.0)*semi_aip
       semi_c=semi_a*phi_cs
!!c
!!c       semi_a=semi_aip
!!c       semi_c=semi_cip
!!c       V_csw=V_cs

    case (2)
       ! rimed crystal
       !
       V_cs=max(V_cs,V_ic)
       if(iswitch==0) then
          ! assume the cylinder volume
!!c          semi_aip=(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0_PS+phi_cs**2.0)
!!c          semi_cip=semi_aip*phi_cs
          call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
       elseif(iswitch==1) then
          ! assume the same growth as in riming scheme.
          if(is_mod(2)==1) then
             ! cyclinder ice particle model
             if(habit<=3) then
                phi_ic=clen/alen
             else
                phi_ic=spx_p
             end if
             if(semi_a_i>semi_c_i) then
                phi_cs=sqrt((V_cs/V_ic)**(2.0/3.0)*(1.0_PS+phi_ic**2)-1.0_PS)
                semi_aip=semi_a_i
                semi_cip=semi_aip*phi_cs
             elseif(semi_a_i<semi_c_i) then
                phi_cs=1.0/sqrt((V_cs/V_ic)**(2.0/3.0)*(1.0_PS+phi_ic**(-2))-1.0_PS)
                semi_cip=semi_c_i
                semi_aip=semi_cip/phi_cs
             else
                phi_cs=phi_ic
                call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
             end if
          else
             ! ice sphere
             phi_cs=1.0_PS
             call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
          endif
!!c          if(semi_a_i>semi_c_i) then
!!c             phi_cs=sqrt((V_cs/V_ic)**(2.0/3.0)*(1.0_PS+phi_ic**2)-1.0_PS)
!!c             semi_aip=semi_a_i
!!c             semi_cip=semi_aip*phi_cs
!!c          elseif(semi_a_i<semi_c_i) then
!!c             phi_cs=1.0/sqrt((V_cs/V_ic)**(2.0/3.0)*(1.0_PS+phi_ic**(-2))-1.0_PS)
!!c             semi_cip=semi_c_i
!!c             semi_aip=semi_cip/phi_cs
!!c          else
!!c             phi_cs=phi_ic
!!c             ! assume the cylinder volume
!!c             semi_aip=(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0_PS+phi_cs**2.0)
!!c             semi_cip=semi_aip*phi_cs
!!c          end if

!!c
!!c       ! assume the axis ratio is kept the same as pristine crystal inside.
!!c       semi_aip=semi_a_i
!!c       semi_cip=semi_c_i
!!c       phi_cs=semi_cip/semi_aip
!!c       V_cs=V_ic
       end if
!!c       semi_aip=max(semi_aip,semi_a_i)
!!c       semi_cip=max(semi_cip,semi_c_i)
       phi_cs=semi_cip/semi_aip

       ! assume the cylinder volume
!!c       v_cs=coef4pi3*semi_aip**3.0*(1.0_PS+phi_cs**2.0)**1.5
       V_cs=get_vcs(is_mod(2),phi_cs,semi_aip)


!!c       V_ip=coef2p*phi_cs*semi_aip**3.0
       V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
       V_space=max(0.0_PS,V_ip-(mean_mass-m_mlt)/den_i)
       ak3=1.0_PS+max(m_mlt/den_w-V_space,0.0_PS)/V_ip
       V_csw=ak3*V_cs
       semi_a=ak3**(1.0/3.0)*semi_aip
       semi_c=semi_a*phi_cs
!!c
!!c       semi_a=semi_aip
!!c       semi_c=semi_cip
!!c       V_csw=V_cs

    case (3)
       ! 3. unrimed aggregates
       !
       V_cs=max(V_cs,V_ic)
       if(iswitch==1) then
          mass_phi_agg=m_ic*mrat_agg
          
          mass1=m_ic+m_ag
          
          if(mass1<mass_phi_agg) then
             if(habit<=3) then
                phi_ic=clen/alen
             else
                phi_ic=spx_p
             end if
             if(phi_agg<=0.0.or.phi_ic<=0.0.or.m_ic<=0.0.or.&
                  mass1<=0.0) then
                write(*,15) phi_ic,m_ic,mass1
15              format("Something is wrong in cs_spheroid for agg:",4ES15.6)
                write(*,'("mass1,m_ic",8ES15.6)') mass1,m_ic
                write(*,'("alen,clen,v_cs",8ES15.6)') alen,clen
                stop
             end if
             
             
!!c             phi=(log10(phi_agg)-log10(phi_ic))/(log10(mass_phi_agg)-log10(m_ic))*&
!!c                  (log10(mass1)-log10(mass_phi_agg))+log10(phi_agg)
!!c             phi=10.0**phi

             phi=phi_agg*(mass1/mass_phi_agg)**(log10(phi_agg/phi_ic)/log10(mass_phi_agg/m_ic))

          else
             phi=phi_agg
          end if
          phi_cs=phi
       end if

       ! assume the cylinder volume
!!c       semi_aip=max(semi_a_i,(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0_PS+phi_cs**2.0))
!!c       semi_cip=max(semi_c_i,semi_aip*phi_cs)
       call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
!!c       semi_aip=max(semi_aip,semi_a_i)
!!c       semi_cip=max(semi_cip,semi_c_i)
       phi_cs=semi_cip/semi_aip

       ! assume the cylinder volume
!!c       v_cs=coef4pi3*semi_aip**3.0*(1.0_PS+phi_cs**2.0)**1.5
       V_cs=get_vcs(is_mod(2),phi_cs,semi_aip)


!!c       V_ip=coef2p*phi_cs*semi_aip**3.0
       V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
       V_space=max(0.0_PS,V_ip-(mean_mass-m_mlt)/den_i)
       ak3=1.0_PS+max(m_mlt/den_w-V_space,0.0_PS)/V_ip
       V_csw=ak3*V_cs
       semi_a=ak3**(1.0/3.0)*semi_aip
       semi_c=semi_a*phi_cs
!!c
!!c       semi_a=semi_aip
!!c       semi_c=semi_cip
!!c       V_csw=V_cs

    case (4)
       ! 4. rimed aggregates
       !
       !  assume cylinder
       !
       !  assume cylinder

       V_cs=max(V_cs,V_ic)
       if(iswitch==1) then
          mass_phi_agg=m_ic*mrat_agg
          
          mass1=m_ic+m_ag
          
          if(mass1<mass_phi_agg) then
             if(habit<=3) then
                phi_ic=clen/alen
             else
                phi_ic=spx_p
             end if
             if(phi_agg<=0.0.or.phi_ic<=0.0.or.m_ic<=0.0.or.&
                  mass1<=0.0) then
                write(*,16) phi_ic,m_ic,mass1
16              format("Something is wrong in cs_spheroid for rimed agg:",4ES15.6)
                write(*,'("mass1,m_ic",8ES15.6)') mass1,m_ic
                write(*,'("alen,clen,v_cs",8ES15.6)') alen,clen
                stop
             end if
             
             
!!c             phi=(log10(phi_agg)-log10(phi_ic))/(log10(mass_phi_agg)-log10(m_ic))*&
!!c                  (log10(mass1)-log10(mass_phi_agg))+log10(phi_agg)
!!c             phi=10.0**phi

             phi=phi_agg*(mass1/mass_phi_agg)**(log10(phi_agg/phi_ic)/log10(mass_phi_agg/m_ic))
          else
             phi=phi_agg
          end if
          phi_cs=phi
       end if

       ! assume the cylinder volume
!!c       semi_aip=max(semi_a_i,(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0_PS+phi_cs**2.0))
!!c       semi_cip=max(semi_c_i,semi_aip*phi_cs)

       call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
!!c       semi_aip=max(semi_aip,semi_a_i)
!!c       semi_cip=max(semi_cip,semi_c_i)

       phi_cs=semi_cip/semi_aip

       ! assume the cylinder volume
!!c       v_cs=coef4pi3*semi_aip**3.0*(1.0_PS+phi_cs**2.0)**1.5
       V_cs=get_vcs(is_mod(2),phi_cs,semi_aip)



!!c       V_ip=coef2p*phi_cs*semi_aip**3.0
       V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
       V_space=max(0.0_PS,V_ip-(mean_mass-m_mlt)/den_i)
       ak3=1.0_PS+max(m_mlt/den_w-V_space,0.0_PS)/V_ip
       V_csw=ak3*V_cs
       semi_a=ak3**(1.0/3.0)*semi_aip
       semi_c=semi_a*phi_cs
!!c
!!c       semi_a=semi_aip
!!c       semi_c=semi_cip
!!c       V_csw=V_cs

    case (5)
       ! graupel
       !
       !    assume that a spheroid is the shape of rimed aggregates.
       !

       V_cs=max(V_cs,V_ic)
!!c       write(*,*) "bef get phi 5"

       if(iswitch==1) then          
          mass_phi_grp=(m_ic+m_ag)*mrat_grp
          mass1=m_ic+m_ag+m_rm
          
          if(mass1<mass_phi_grp) then
             if(habit<=3) then
                phi_ic=clen/alen
             else
                phi_ic=spx_p
             end if
             
             if(phi_grp<=0.0.or.phi_ic<=0.0.or.m_ic<=0.0.or.&
                  mass1<=0.0) then
                write(*,17) phi_ic,m_ic,mass1
17              format("Something is wrong in cs_spheroid at graupels:",4ES15.6)
                write(*,'("mass1,m_ic",8ES15.6)') mass1,m_ic,m_ag,m_rm
                write(*,'("alen,clen,v_cs",8ES15.6)') alen,clen
                stop
             end if
             
!!c             phi=(log10(phi_grp)-log10(phi_ic))/(log10(mass_phi_grp)-log10(m_ic))*&
!!c                  (log10(mass1)-log10(mass_phi_grp))+log10(phi_grp)
!!c             phi=10.0**phi
             phi=phi_grp*(mass1/mass_phi_grp)**(log10(phi_grp/phi_ic)/log10(mass_phi_grp/m_ic))
          else
             phi=phi_grp
          end if
!!c       write(*,*) "aft get phi 5"

          phi_cs=phi

       end if

       ! assume the spheroidal volume
!!c       if(phi_cs<1.0) then
!!c          semi_aip=max(semi_a_i,(V_cs/coef4pi3)**(1.0/3.0))
!!c          semi_cip=max(semi_c_i,semi_aip*phi_cs)
!!c       else
!!c          semi_cip=max(semi_c_i,(V_cs/coef4pi3)**(1.0/3.0))
!!c          semi_aip=max(semi_a_i,semi_cip/phi_cs)
!!c       end if

       call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
!!c       semi_aip=max(semi_aip,semi_a_i)
!!c       semi_cip=max(semi_cip,semi_c_i)

       phi_cs=semi_cip/semi_aip

!!c       ! assume the spheroidal volume
!!c       if(phi_cs<1.0_PS) then
!!c          v_cs=coef4pi3*semi_aip**3.0
!!c       else
!!c          v_cs=coef4pi3*(semi_aip*phi_cs)**3.0
!!c       end if
       V_cs=get_vcs(is_mod(2),phi_cs,semi_aip)


!!c       V_ip=coef4pi3*phi_cs*semi_aip**3.0
       V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
       V_space=max(0.0_PS,V_ip-(mean_mass-m_mlt)/den_i)
       ak3=1.0_PS+max(m_mlt/den_w-V_space,0.0_PS)/V_ip
       V_csw=ak3*V_cs
       semi_a=ak3**(1.0/3.0)*semi_aip
       semi_c=semi_a*phi_cs
!!c
!!c       semi_a=semi_aip
!!c       semi_c=semi_cip
!!c       V_csw=V_cs



    case default
       write(*,*) "ERROR! cal_cs_sheroid."
       stop
    end select
  end subroutine cal_cs_spheroid3


  subroutine cal_cs_sphere(ms, is, level)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate 
    ! 1. semi-axis lengths of circumscribing sphere
    ! 2. axis ratio of it.
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Mass_Bin), intent(inout)   :: ms
    type (Ice_Shape), intent(inout)    :: is
    integer, intent(in)    ::  level
    real(PS)               :: dia
    !     m = a * dia^b
    ! mass in mg, diameter in mm
    real(PS)               :: a, b
    real(PS)               :: phi,max_ice,min_ice,phi_ic
    ! minimum mass for fragments to be aggregates or rimed aggrates
    real(PS),parameter :: min_mass=1.0e-5_PS,den_max=0.758088258_PS

    real(PS),parameter :: phi_agg=0.25_PS,mass_phi_agg=1.0e-4_PS
    real(PS),parameter :: phi_grp=0.7_PS,mass_phi_grp=1.0e-3_PS
    ! Matson and Huggins (1980)'s meadian of hail stones
!!c    real(PS),parameter :: phi_grp=0.77_PS,mass_phi_grp=0.479972525615449


    real(PS),parameter :: phi_re=0.25_PS
    real(PS),parameter :: den_a=0.3_PS


    select case(is%sh_type)
    case (1)
       ! pristine crystal
       ! assume the axis ratio is kept the same as pristine crystal inside.
       if(is%habit<=3) then
          ms%semi_a=ms%a_len
          ms%semi_c=ms%c_len
          is%phi_cs = ms%semi_c/ms%semi_a
       elseif(is%habit>=4) then
          ms%semi_a=max(is%r,is%e)
          is%phi_cs=0.25_PS
          ms%semi_c=ms%semi_a*is%phi_cs
       end if
    case (2)
       ! rimed crystal
       ! assume the axis ratio is kept the same as pristine crystal inside.
       if(is%habit<=3) then
          ms%semi_a=ms%a_len
          ms%semi_c=ms%c_len
          is%phi_cs = ms%semi_c/ms%semi_a
       elseif(is%habit>=4) then
          ms%semi_a=max(is%r,is%e)
          is%phi_cs=0.25_PS
          ms%semi_c=ms%semi_a*is%phi_cs
       end if

    case (3)
       ! 3. unrimed aggregates
       !
       !    assume that a cylinder is the shape of aggregates.
       !
       dia=(is%V_cs/coefpi6)**(1.0/3.0)
       is%phi_cs=1.0_PS

       ms%semi_a=dia/2.0_PS
       ms%semi_c=ms%semi_a
    case (4)
       ! 4. rimed aggregates
       !
       !    assume that a spheroid is the shape of rimed aggregates.
       !
       dia=(is%V_cs/coefpi6)**(1.0/3.0)
       is%phi_cs=1.0_PS

       ms%semi_a=dia/2.0_PS
       ms%semi_c=ms%semi_a
    case (5)
       ! graupel
       !
       !    assume that a spheroid is the shape of rimed aggregates.
       !
       dia=(is%V_cs/coefpi6)**(1.0/3.0)
       is%phi_cs=1.0_PS

       ms%semi_a=dia/2.0_PS
       ms%semi_c=ms%semi_a
    case (6)
       ! hail
       ! assume that the axis ratio is 1.
       dia=(is%V_cs/coefpi6)**(1.0/3.0)
       is%phi_cs=1.0_PS

       ms%semi_a=dia/2.0_PS
       ms%semi_c=ms%semi_a
    case default
       write(*,*) "ERROR! cal_cs_sheroid."
       stop
    end select

    ! +++ bound the semi_a, and semi_c estimates +++
!!c    if(is%habit<=3) then
!!c       min_ice=min(ms%a_len,ms%c_len)
!!c       max_ice=max(ms%a_len,ms%c_len)
!!c    else
!!c       min_ice=max(is%r,is%e)
!!c       max_ice=max(is%r,is%e)
!!c    end if
!!c    if(ms%semi_a>=ms%semi_c) then
!!c       ms%semi_a=max(max_ice,ms%semi_a)
!!c       ms%semi_c=max(min_ice,ms%semi_c)
!!c    else
!!c       ms%semi_c=max(max_ice,ms%semi_c)
!!c       ms%semi_a=max(min_ice,ms%semi_a)
!!c    end if
!!c    is%phi_cs=ms%semi_c/ms%semi_a
!!c
!!c    is%V_cs=4.188790205_PS*((1.0+is%phi_cs**2.0)**1.5)*ms%semi_a**3.0


  end subroutine cal_cs_sphere

  subroutine cal_bulk_density(level,token,mean_mass,den,&
       m_ic,den_ic,&
       habit,alen,clen,dlen,rlen,elen,phi_ic,psi_ic,v_ic,&
       sh_type,semi_a,semi_c,phi_cs,v_cs,&
       m_mlt, &     
       from)

    integer,intent(in) :: level,token,habit,sh_type
    real(PS),intent(in) :: mean_mass
    real(PS),intent(inout) :: den
    real(PS),optional,intent(in) :: m_ic
    real(PS),optional,intent(inout) :: den_ic,alen,clen,dlen,rlen,elen,phi_ic,psi_ic,v_ic,&
         semi_a,semi_c,phi_cs,v_cs
    real(PS),optional,intent(in) :: m_mlt
    character (len=*) :: from
    
    real(PS) :: den_max,v_ic_hex
    integer :: ierror

    ! 3 sqrt(3), 4pi/3
    real(PS),parameter :: coef3s=5.196152423, coef4pi3=4.18879020478639

    ! assume that the possible minimum density is 1.0e-5
    real(PS),parameter :: den_min=1.0e-5

    if(token==1) then
       ! liquid hydrometeors
       ! assume den=1.0
       den=1.0_PS
    elseif(token==2) then

       ! 1. calculate bulk crystal density
       ! pristine crystal and rimed crystals
       v_ic_hex=coef4pi3*(alen**2+clen**2)**1.5
       den_ic=m_ic/v_ic_hex

       ! maximum possible bulk sphere density of pristine hexagonal crystal
       den_max=coef3s*phi_ic*den_i*(1.0-psi_ic)/&
            (coef4pi3*(1.0+phi_ic**2.0)**1.5)

       ierror=0
       if(den_ic/den_max-1.0>0.1) then
          write(*,*) from
          write(*,'("cal_bulk_density:error 1, bf, alen,clen,dlen,m_ic,denic",5ES15.6)') alen,clen,dlen,m_ic,den_ic

          den_ic=den_max
          V_ic_hex=m_ic/den_max



          alen=(V_ic_hex/coef4pi3)**(1.0/3.0)/sqrt(1.0+phi_ic**2.0)
          clen=alen*phi_ic
          dlen=alen*psi_ic

          write(*,'("                       1, af alen,clen,dlen,m_ic,denx",5ES15.6)') alen,clen,dlen,m_ic,den_max
          ierror=1
       end if

       ! 2. calculate bulk sphere density       
       if(sh_type==1) then
          if(habit<=3) then
             den=den_ic
             if(ierror==1) then
                semi_a=alen
                semi_c=clen
                phi_cs=semi_c/semi_a
                V_ic=V_ic_hex
                V_cs=V_ic_hex
             end if
          elseif(habit>=4) then
             den=mean_mass/v_cs
             if(den>den_i) then
                ! assume it is sphere
                write(*,*) from

                den=den_i
                V_ic=mean_mass/den
                V_cs=V_ic

                semi_a=(V_cs/coef4pi3)**(1.0/3.0)
                semi_c=semi_a
                phi_cs=1.0_PS

                if(habit==4) then
                   rlen=semi_a
                elseif(habit==5) then
                   elen=semi_a
                else
                   elen=semi_a
                   rlen=elen
                end if

                write(*,*) "cal_bulk_density:error 2, alen,clen,dlen,m_ic",alen,clen,dlen,m_ic
                write(*,*) "cal_bulk_density:error 2, mm,semia,semic",mean_mass,semi_a,semi_c

             end if
          end if

       elseif(sh_type==2) then
          den=mean_mass/V_cs
!!c          if(den>den_i) then
          if((mean_mass-m_mlt)/V_cs>den_i) then
             ! assume it is sphere


!!c             den=den_i
!!c             V_cs=mean_mass/den_i
             V_cs=(mean_mass-m_mlt)/den_i
             den=mean_mass/V_cs
             semi_a=(V_cs/coef4pi3)**(1.0/3.0)
             semi_c=semi_a
             phi_cs=1.0_PS


             alen=(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0+phi_ic**2.0)
             clen=alen*phi_ic
             dlen=alen*psi_ic

             if(habit==4) then
                rlen=semi_a
             elseif(habit==5) then
                elen=semi_a
             else
                rlen=semi_a
                elen=semi_a
             end if
             V_ic=V_cs

             write(*,*) from
             write(*,*) "cal_bulk_density:error 3, alen,clen,dlen,m_ic",alen,clen,dlen,m_ic
             write(*,*) "cal_bulk_density:error 3, mm,semia,semic",mean_mass,semi_a,semi_c

          end if

       elseif(sh_type<=4) then
          den=mean_mass/V_cs
          if((mean_mass-m_mlt)/V_cs>den_i) then

             write(*,*) from
             write(*,'("cal_bulk_density:error 4, bf, semi_a,semi_c,den",3ES15.6)') semi_a,semi_c,den
             V_cs=(mean_mass-m_mlt)/den_i
             den=mean_mass/V_cs
             semi_a=(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0_PS+phi_cs**2.0)
             semi_c=semi_a*phi_cs

             write(*,'("cal_bulk_density:error 4, af, semi_a,semi_c,den",3ES15.6)') semi_a,semi_c,den
          elseif(den<den_min) then
             write(*,*) from
             write(*,'("cal_bulk_density:error 6 bf, den,semi_a,semi_c",3ES15.6)') den,semi_a,semi_c
             den=den_min
             V_cs=mean_mass/den_min
             semi_a=(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0_PS+phi_cs**2.0)
             semi_c=semi_a*phi_cs

             alen=min(semi_a,alen)
             clen=min(semi_c,clen)
             phi_ic=clen/alen
             dlen=alen*psi_ic
             V_ic=coef4pi3*alen**3.0*(1.0_PS+phi_ic**2.0)**1.5

             write(*,'("cal_bulk_density:error 6 af, den,semi_a,semi_c",3ES15.6)') den,semi_a,semi_c
          end if
       else
          den=mean_mass/V_cs

          if((mean_mass-m_mlt)/V_cs>den_i) then
             write(*,*) from
             write(*,'("cal_bulk_density:error 5, bf, semi_a,semi_c,den",3ES15.6)') semi_a,semi_c,den

             V_cs=(mean_mass-m_mlt)/den_i
             den=mean_mass/V_cs


             if(phi_cs<1.0) then
                semi_a=(V_cs/coef4pi3)**(1.0/3.0)
                semi_c=semi_a*phi_cs
             else
                semi_c=(V_cs/coef4pi3)**(1.0/3.0)
                semi_a=semi_c/phi_cs
             end if



             write(*,'("cal_bulk_density:error 5, af, semi_a,semi_c,den",3ES15.6)') semi_a,semi_c,den
          elseif(den<den_min) then
             write(*,*) from
             write(*,'("cal_bulk_density:error 7 bf, den,semi_a,semi_c",3ES15.6)') den,semi_a,semi_c
             den=den_min
             V_cs=mean_mass/den_min

             if(phi_cs<1.0) then
                semi_a=(V_cs/coef4pi3)**(1.0/3.0)
                semi_c=semi_a*phi_cs
             else
                semi_c=(V_cs/coef4pi3)**(1.0/3.0)
                semi_a=semi_c/phi_cs
             end if


             alen=min(semi_a,alen)
             clen=min(semi_c,clen)
             phi_ic=clen/alen
             dlen=alen*psi_ic
             V_ic=coef4pi3*alen**3.0*(1.0_PS+phi_ic**2.0)**1.5

             write(*,'("cal_bulk_density:error 7 af, den,semi_a,semi_c",3ES15.6)') den,semi_a,semi_c

          end if
       end if
    end if
  end subroutine cal_bulk_density
  subroutine cal_bulk_density2(level,token,mean_mass,den,&
       m_ic,den_ic,&
       habit,alen,clen,dlen,rlen,elen,phi_ic,psi_ic,v_ic,&
       sh_type,semi_aip,semi_cip,phi_cs,V_cs,den_ip,&
       semi_a,semi_c,V_csw,&
       m_rim,m_agg,m_mlt, &     
       is_mod,iswitch,iswitch_grp,&
       from)
    real(PS),intent(inout) :: den,den_ip
    integer,intent(in) :: level,token,habit,sh_type,iswitch,iswitch_grp
    real(PS),intent(in) :: mean_mass
    real(PS),optional,intent(in) :: m_ic
    real(PS),optional,intent(inout) :: den_ic,alen,clen,dlen,rlen,elen,phi_ic,psi_ic,v_ic,&
         semi_aip,semi_cip,phi_cs,v_cs,semi_a,semi_c,v_csw
    real(PS),optional,intent(in) :: m_mlt,m_rim,m_agg
    integer,dimension(2) :: is_mod
    character (len=*) :: from
    
    real(PS) :: den_max,v_ic_hex,v_space,v_ip,ak3, phi_cs_max,m_dry
    integer :: ierror,ierror2

    ! assume that the possible minimum density is 1.0e-5
    real(PS),parameter :: den_min=1.0e-5
    ! assume that the possible minimum density of hex ice crystal is 1.0e-4
    ! The hex column with dendritic arms (phi=20, psi=0.9) can reach 
    ! 0.0002832.
    real(PS),parameter :: den_min_hex=1.0e-4

    ! bulk water density
    real(PS),parameter :: den_w=1.0
    ! fraction of maximum possible ice-water mixture
    real(PS),parameter :: fdenmx=0.9999

    if(token==1) then
       ! liquid hydrometeors
       ! assume den=1.0
       den=1.0_PS
    elseif(token==2) then

       ! 1. calculate bulk crystal density
       ! pristine crystal and rimed crystals
       v_ic_hex=coef4pi3*(alen**2+clen**2)**1.5
       den_ic=m_ic/v_ic_hex

       ! maximum possible bulk sphere density of pristine hexagonal crystal
       den_max=coef3sq3*phi_ic*den_i*(1.0-psi_ic)/&
            (coef4pi3*(1.0+phi_ic**2.0)**1.5)

       ierror=0
       if(den_ic/den_max-1.0>0.1) then
!!c          write(*,*) from
!!c          write(*,'("cal_bulk_density:error 1, bf, alen,clen,dlen,m_ic,denic",5ES15.6)') alen,clen,dlen,m_ic,den_ic

          den_ic=den_max
          V_ic_hex=m_ic/den_max

          alen=(V_ic_hex/coef4pi3)**(1.0/3.0)/sqrt(1.0+phi_ic**2.0)
          clen=alen*phi_ic
          dlen=alen*psi_ic

!!c          write(*,'("                       1, af alen,clen,dlen,m_ic,denx",5ES15.6)') alen,clen,dlen,m_ic,den_max
          ierror=1
       end if
       if(den_ic<den_min_hex) then
!!c          write(*,*) from
!!c          write(*,'("cal_bulk_density:error 1, bf, alen,clen,dlen,m_ic,denic",5ES15.6)') alen,clen,dlen,m_ic,den_ic

          den_ic=den_min_hex
          V_ic_hex=m_ic/den_min_hex

          alen=(V_ic_hex/coef4pi3)**(1.0/3.0)/sqrt(1.0+phi_ic**2.0)
          clen=alen*phi_ic
          dlen=alen*psi_ic

!!c          write(*,'("                       1, af alen,clen,dlen,m_ic,denx",5ES15.6)') alen,clen,dlen,m_ic,den_max
          ierror=1
       end if

       ! 2. calculate bulk sphere density for dry ice particle and wet ice particle
       ierror2=0
       m_dry=m_ic+m_agg+m_rim
       if(sh_type<=2) then

          den_ip=m_dry/v_cs
          if(habit<=3) then
             if(ierror==1) then
                semi_aip=alen
                semi_cip=clen
                phi_cs=semi_cip/semi_aip
                V_ic=V_ic_hex
                V_cs=V_ic_hex
             end if
          elseif(habit>=4) then
!!c             if(den_ip>=den_i.or.semi_aip**2+semi_cip**2<4.0e-6) then
             if(den_ip>=den_i) then
                !
                ! assume it is sphere
                ! this may happen because the empirical equation used for these habits
                ! can produce density more than den_i in small range.
                !
                write(*,*) from,habit

                V_ic=m_dry/den_i
                V_cs=V_ic
                den_ip=den_i

                semi_aip=(V_cs/coef4pi3)**(1.0/3.0)
                semi_cip=semi_aip
                phi_cs=1.0_PS

                if(habit==4) then
                   rlen=semi_aip
                elseif(habit==5) then
                   elen=semi_aip
                elseif(habit==6) then
                   elen=semi_aip
                   rlen=elen
                end if

                write(*,'("cal_bulk_density:error 1_2, alen,clen,dlen,m_ic,m_mlt",10ES15.6)') &
                     alen,clen,dlen,m_ic,m_mlt,mean_mass,semi_aip,semi_cip
                is_mod(1)=2
                is_mod(2)=2

             end if
          end if

          if(iswitch==1) then
             if(is_mod(2)==1) then
                den_max=1.5_PS*phi_cs/(1.0_PS+phi_cs**2)**1.5*den_i
             elseif(is_mod(2)==2) then
                den_max=min(phi_cs,phi_cs**(-2))*den_i
             end if
             
             call den_ck1
             if(ierror2==1) then
                alen=(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0+phi_ic**2.0)
                clen=alen*phi_ic
                dlen=alen*psi_ic
                
                if(habit==4) then
                   rlen=semi_aip
                elseif(habit==5) then
                   elen=semi_aip
                elseif(habit==6) then
                   rlen=semi_aip
                   elen=semi_aip
                end if
                V_ic=V_cs
             end if
          end if

          den=mean_mass/V_csw
          if(iswitch==1) then
             if(is_mod(2)==1) then
                den_max=1.5_PS*phi_cs/(1.0_PS+phi_cs**2)**1.5*den_w*fdenmx
             elseif(is_mod(2)==2) then
                den_max=min(phi_cs,phi_cs**(-2))*den_w*fdenmx
             end if
             call den_ck3
             call den_ck4
          end if
       elseif(sh_type<=4) then

          den_ip=m_dry/v_cs
          if(iswitch==1) then
             if(is_mod(2)==1) then
                den_max=1.5_PS*phi_cs/(1.0_PS+phi_cs**2)**1.5*den_i
             elseif(is_mod(2)==2) then
                den_max=min(phi_cs,phi_cs**(-2))*den_i
             end if
             call den_ck1
             call den_ck2
          end if

          den=mean_mass/V_csw
          if(iswitch==1) then
             if(is_mod(2)==1) then
                den_max=1.5_PS*phi_cs/(1.0_PS+phi_cs**2)**1.5*den_w*fdenmx
             elseif(is_mod(2)==2) then
                den_max=min(phi_cs,phi_cs**(-2))*den_w*fdenmx
             end if
             call den_ck3
             call den_ck4
          end if
       else

          den_ip=m_dry/v_cs
          if(iswitch==1) then

             if(iswitch_grp==1) then
                ! This is because the aspect ratio is accurately not predicted by _rim5, so that
                ! den_max is not really correct.
!!c          phi_cs_max=den_ip/den_i
                if(is_mod(2)==1) then
                   den_max=1.5_PS*phi_cs/(1.0_PS+phi_cs**2)**1.5*den_i
                elseif(is_mod(2)==2) then
                   den_max=min(phi_cs,phi_cs**(-2))*den_i
                end if
                call den_ck1
             end if
             call den_ck2
          end if

          den=mean_mass/V_csw
          if(iswitch==1) then
             call den_ck3
             if(iswitch_grp==1) then
                if(is_mod(2)==1) then
                   den_max=1.5_PS*phi_cs/(1.0_PS+phi_cs**2)**1.5*den_w*fdenmx
                elseif(is_mod(2)==2) then
                   den_max=min(phi_cs,phi_cs**(-2))*den_w*fdenmx
                end if
                call den_ck4
             end if
          end if
       end if
    end if
   contains
    subroutine den_ck1
      if(den_ip-den_max>1.0e-4*den_max) then
!!c      if(den_ip>den_max) then
!!c         write(*,*) from,habit,sh_type
!!c         write(*,'("cal_bulk_density:error 1, bf, semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw

         V_cs=m_dry/den_max
         call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
         den_ip=den_max
         ierror2=1
!!c         write(*,'("cal_bulk_density:error 1, af, semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw
      end if
    end subroutine den_ck1
    subroutine den_ck2
!!c      if(den_ip<den_min) then
      if(den_ip-den_min<-1.0e-4*den_min) then
!!c         write(*,*) from,habit,sh_type
!!c         write(*,'("cal_bulk_density:error 2, bf, semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw
         V_cs=m_dry/den_min
         call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
         den_ip=den_min

         alen=min(semi_aip,alen)
         clen=min(semi_cip,clen)
         phi_ic=clen/alen
         dlen=alen*psi_ic
         V_ic=coef4pi3*alen**3.0*(1.0_PS+phi_ic**2.0)**1.5

         ierror2=1
!!c         write(*,'("cal_bulk_density:error 2, af, semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw
      end if
    end subroutine den_ck2
    subroutine den_ck3
      if(ierror2==1) then
!!c      if(den>den_w.or.ierror2==1) then
!!c         write(*,*) from,habit,sh_type
!!c         write(*,'("cal_bulk_density:error 3, bf, semi_a,semi_c,den,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,den,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw
!!c         end if

         V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
         V_space=max(0.0_PS,V_ip-m_dry/den_i)
         ak3=1.0_PS+max(m_mlt-V_space*den_w,0.0_PS)/V_ip
         V_csw=ak3*V_cs
         semi_a=ak3**(1.0/3.0)*semi_aip
         semi_c=semi_a*phi_cs
         den=mean_mass/V_csw
!!c         write(*,'("cal_bulk_density:error 3, af, semi_a,semi_c,den,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,den,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw
      end if
    end subroutine den_ck3
    subroutine den_ck4
      if(den-den_max>1.0e-4*den_max) then
         write(*,*) from,habit,sh_type
         write(*,'("cal_bulk_density:error 4-1, bf, semi_a,semi_c,den,semi_aip,semi_cip",20ES15.6)') &
             semi_a,semi_c,den,semi_aip,semi_cip
         write(*,'("cal_bulk_density:error 4-2, bf, den_ip,mean_mass,m_ic,m_mlt,m_rim",20ES15.6)') &
             den_ip,mean_mass,m_ic,m_mlt,m_rim
         write(*,'("cal_bulk_density:error 4-3, bf, V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') &
             V_ic,V_ip,V_space,V_cs,V_csw

         V_csw=mean_mass/den_max
         call cal_semiac_ip(is_mod(2),phi_cs,v_csw,semi_a,semi_c)
         den=den_max

         write(*,'("cal_bulk_density:error 4-1, af, semi_a,semi_c,den,semi_aip,semi_cip",20ES15.6)') &
             semi_a,semi_c,den,semi_aip,semi_cip
         write(*,'("cal_bulk_density:error 4-2, af, den_ip,mean_mass,m_ic,m_mlt,m_rim",20ES15.6)') &
             den_ip,mean_mass,m_ic,m_mlt,m_rim
         write(*,'("cal_bulk_density:error 4-3, af, V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') &
             V_ic,V_ip,V_space,V_cs,V_csw

!!c         if(V_space>1.0e-6) then
!!c         ! increase V_ip to make V_space to hold the water, then some of the melt water has to be decreased.
!!c         den_ism=m_dry/V_ip
!!c         if(mean_mass<=m_mlt.or.den_ism-den_i>den_i*1.0e-4) then
!!c            write(*,*) "den_ck4>woops",mean_mass,m_mlt,den_ism,V_ip
!!c         end if
!!c         
!!c         dV_ip=(mean_mass/den_i+(1.0_PS/den_w-1.0_PS/den_i)*m_mlt-V_ip)/&
!!c                     ((1.0_PS/den_w-1.0_PS/den_i)*den_ism+1.0_PS)
!!c         if(dV_ip<=0.0_PS) then
!!c            write(*,*) "den_ck4>something not right:",V_ip,dV_ip,mean_mass,m_mlt,den_ism,semi_aip,semi_cip,phi_cs
!!c!!c            stop
!!c            return
!!c         end if
!!c         m_mlt=m_mlt-den_ism*dV_ip
!!c         m_rim=m_rim+den_ism*dV_ip
!!c
!!c         V_ip=V_ip+dV_ip
!!c         semi_aip=(V_ip/get_coef_ip(is_mod(2))/phi_cs)**(1.0/3.0)
!!c         semi_cip=semi_aip*phi_cs
!!c         V_cs=get_vcs(is_mod(2),phi_cs,semi_aip)
!!c         V_csw=V_cs
!!c         semi_a=semi_aip
!!c         semi_c=semi_cip
!!c         end if  
!!c         den=mean_mass/V_csw
!!c         if(habit/=6) then
!!c         end if
!!c         if(den>den_w) then
!!c            write(*,*) "cal_bulk_density>something is wrong:",den,semi_a,semi_c,V_csw,V_ip,V_space,ak3,m_mlt
!!c         end if
      end if
    end subroutine den_ck4

  end subroutine cal_bulk_density2

  subroutine cal_bulk_density3(level,token,mean_mass,den,&
       m_ic,den_ic,&
       habit,alen,clen,dlen,rlen,elen,phi_ic,psi_ic,v_ic,&
       sh_type,semi_aip,semi_cip,phi_cs,V_cs,den_ip,&
       semi_a,semi_c,V_csw,&
       m_rim,m_agg,m_mlt, &     
       is_mod,iswitch,iswitch_grp,&
       id,jd,kd,from)
    real(PS),intent(inout) :: den,den_ip
    integer,intent(in) :: level,token,habit,sh_type,iswitch,iswitch_grp
    integer,intent(in) :: id,jd,kd
    real(PS),intent(in) :: mean_mass
    real(PS),optional,intent(in) :: m_ic
    real(PS),optional,intent(inout) :: den_ic,alen,clen,dlen,rlen,elen,phi_ic,psi_ic,v_ic,&
         semi_aip,semi_cip,phi_cs,v_cs,semi_a,semi_c,v_csw
    real(PS),optional,intent(in) :: m_mlt,m_rim,m_agg
    integer,dimension(2) :: is_mod
    character (len=*) :: from
    
    real(PS) :: den_max,v_ic_hex,v_space,v_ip,ak3, phi_cs_max,m_dry
    integer :: ierror,ierror2

    ! assume that the possible minimum density is 1.0e-4
    real(PS),parameter :: den_min=1.0e-4
    ! assume that the possible minimum density of hex ice crystal is 1.0e-4
    ! The hex column with dendritic arms (phi=20, psi=0.9) can reach 
    ! 0.0002832.
    real(PS),parameter :: den_min_hex=1.0e-4

    ! bulk water density
    real(PS),parameter :: den_w=1.0
    ! fraction of maximum possible ice-water mixture
    real(PS),parameter :: fdenmx=0.9999

    if(token==1) then
       ! liquid hydrometeors
       ! assume den=1.0
       den=1.0_PS
    elseif(token==2) then

       ! 1. calculate bulk crystal density
       ! pristine crystal and rimed crystals
       v_ic_hex=coef4pi3*(alen**2+clen**2)**1.5
       den_ic=m_ic/v_ic_hex

       ! maximum possible bulk sphere density of pristine hexagonal crystal
       den_max=coef3sq3*phi_ic*den_i*(1.0-psi_ic)/&
            (coef4pi3*(1.0+phi_ic**2.0)**1.5)

       ierror=0
       if(den_ic/den_max-1.0>0.1) then
!!c          write(*,*) from
!!c          write(*,'("cal_bulk_density:error 1a, bf, id,jd,kd,alen,clen,dlen,m_ic,denic,den_max,phi",3I5,10ES15.6)') &
!!c                     id,jd,kd,alen,clen,dlen,m_ic,den_ic,den_max,clen/alen

          den_ic=den_max
          V_ic_hex=m_ic/den_max

          alen=(V_ic_hex/coef4pi3)**(1.0/3.0)/sqrt(1.0+phi_ic**2.0)
          clen=alen*phi_ic
          dlen=alen*psi_ic

!!c          write(*,'("                       1a, af id,jd,kd,alen,clen,dlen,m_ic,denic,den_mxd,phi",3I5,10ES15.6)') &
!!c                     id,jd,kd,alen,clen,dlen,m_ic,den_ic,den_max,clen/alen
          ierror=1
       end if
       if(den_ic<den_min_hex) then
!!c          write(*,*) from
!!c          write(*,'("cal_bulk_density:error 1b, bf, id,jd,kd,alen,clen,dlen,m_ic,denic,phi",3I5,10ES15.6)') &
!!c                   id,jd,kd,alen,clen,dlen,m_ic,den_ic,clen/alen

          den_ic=den_min_hex
          V_ic_hex=m_ic/den_min_hex

          alen=(V_ic_hex/coef4pi3)**(1.0/3.0)/sqrt(1.0+phi_ic**2.0)
          clen=alen*phi_ic
          dlen=alen*psi_ic

!!c          write(*,'("                       1b, af id,jd,kd,alen,clen,dlen,m_ic,denic,phi",3I5,10ES15.6)') &
!!c                    id,jd,kd,alen,clen,dlen,m_ic,den_ic,clen/alen
          ierror=1
       end if

       ! 2. calculate bulk sphere density for dry ice particle and wet ice particle
       ierror2=0
       m_dry=m_ic+m_agg+m_rim
       if(sh_type<=2) then

          den_ip=m_dry/v_cs

          if(is_mod(2)==1) then
             if(habit<=3) then
                if(ierror==1) then
                   semi_aip=alen
                   semi_cip=clen
                   phi_cs=semi_cip/semi_aip
                   V_ic=V_ic_hex
                   V_cs=V_ic_hex

                   den_ip=m_dry/v_cs

                   V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
                   V_space=max(0.0_PS,V_ip-(mean_mass-m_mlt)/den_i)
                   ak3=1.0_PS+max(m_mlt/den_w-V_space,0.0_PS)/V_ip
                   V_csw=ak3*V_cs

                end if
             elseif(habit>=4) then
!!c             if(den_ip>=den_i.or.semi_aip**2+semi_cip**2<4.0e-6) then
                if(den_ip>=den_i) then
                   !
                   ! assume it is sphere
                   ! this may happen because the empirical equation used for these habits
                   ! can produce density more than den_i in small range.
                   !
!!c                   write(*,*) from,habit

                   V_ic=m_dry/den_i
                   V_cs=V_ic
                   den_ip=den_i

                   semi_aip=(V_cs/coef4pi3)**(1.0/3.0)
                   semi_cip=semi_aip
                   phi_cs=1.0_PS
   
                   if(habit==4) then
                      rlen=semi_aip
                   elseif(habit==5) then
                      elen=semi_aip
                   elseif(habit==6) then
                      elen=semi_aip
                      rlen=elen
                   end if

!!c                   write(*,*) "cal_bulk_density:error 1_2, alen,clen,dlen,m_ic,m_mlt",&
!!c                             alen,clen,dlen,m_ic,m_mlt,mean_mass,semi_aip,semi_cip

                   is_mod(1)=2
                   is_mod(2)=2
                endif
             endif
          elseif(is_mod(2)==2) then
             if(den_ip*0.9999>den_i) then
                !
                ! assume it is sphere
                !
!!c                write(*,*) from,habit,sh_type
!!c                write(*,'("cal_bulk_density:error 1_2_3, den_ip,a,c,d,m_ic,m_mlt,mm,saip,scip",10ES15.6)') &
!!c                        den_ip,alen,clen,dlen,m_ic,m_mlt,mean_mass,semi_aip,semi_cip

                V_cs=m_dry/den_i
                den_ip=den_i

                semi_aip=(V_cs/coef4pi3)**(1.0/3.0)
                semi_cip=semi_aip
                phi_cs=1.0_PS
   
             endif
          endif

          if(iswitch==1) then
             if(is_mod(2)==1) then
                den_max=1.5_PS*phi_cs/(1.0_PS+phi_cs**2)**1.5*den_i
             elseif(is_mod(2)==2) then
                den_max=min(phi_cs,phi_cs**(-2))*den_i
             end if
             
             call den_ck1
             if(sh_type==2) call den_ck2
             if(ierror2==1) then
                alen=(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0+phi_ic**2.0)
                clen=alen*phi_ic
                dlen=alen*psi_ic
                
                if(habit==4) then
                   rlen=semi_aip
                elseif(habit==5) then
                   elen=semi_aip
                elseif(habit==6) then
                   rlen=semi_aip
                   elen=semi_aip
                end if
                V_ic=V_cs
             end if
          end if

          den=mean_mass/V_csw
          if(iswitch==1) then
             if(is_mod(2)==1) then
                den_max=1.5_PS*phi_cs/(1.0_PS+phi_cs**2)**1.5*den_w*fdenmx
             elseif(is_mod(2)==2) then
                den_max=min(phi_cs,phi_cs**(-2))*den_w*fdenmx
             end if
             call den_ck3
             call den_ck4
          end if
       elseif(sh_type<=4) then

          den_ip=m_dry/v_cs
          if(iswitch==1) then
             if(is_mod(2)==1) then
                den_max=1.5_PS*phi_cs/(1.0_PS+phi_cs**2)**1.5*den_i
             elseif(is_mod(2)==2) then
                den_max=min(phi_cs,phi_cs**(-2))*den_i
             end if
             call den_ck1
             call den_ck2
          end if

          den=mean_mass/V_csw
          if(iswitch==1) then
             if(is_mod(2)==1) then
                den_max=1.5_PS*phi_cs/(1.0_PS+phi_cs**2)**1.5*den_w*fdenmx
             elseif(is_mod(2)==2) then
                den_max=min(phi_cs,phi_cs**(-2))*den_w*fdenmx
             end if
             call den_ck3
             call den_ck4
          end if
       else

          den_ip=m_dry/v_cs
          if(iswitch==1) then

             if(iswitch_grp==1) then
                ! This is because the aspect ratio is accurately not predicted by _rim5, so that
                ! den_max is not really correct.
!!c          phi_cs_max=den_ip/den_i
                if(is_mod(2)==1) then
                   den_max=1.5_PS*phi_cs/(1.0_PS+phi_cs**2)**1.5*den_i
                elseif(is_mod(2)==2) then
                   den_max=min(phi_cs,phi_cs**(-2))*den_i
                end if
                call den_ck1
             end if
             call den_ck2
          end if

          den=mean_mass/V_csw
          if(iswitch==1) then
             call den_ck3
             if(iswitch_grp==1) then
                if(is_mod(2)==1) then
                   den_max=1.5_PS*phi_cs/(1.0_PS+phi_cs**2)**1.5*den_w*fdenmx
                elseif(is_mod(2)==2) then
                   den_max=min(phi_cs,phi_cs**(-2))*den_w*fdenmx
                end if
                call den_ck4
             end if
          end if
       end if
    end if
   contains
    subroutine den_ck1
      if(den_ip-den_max>1.0e-4*den_max) then
!!c      if(den_ip>den_max) then
!!c         write(*,*) from,habit,sh_type
!!c         write(*,'("cal_bulk_density:error 1, bf, semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw

         V_cs=m_dry/den_max
         call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
         den_ip=den_max
         ierror2=1
!!c         write(*,'("cal_bulk_density:error 1, af, semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw
      end if
    end subroutine den_ck1
    subroutine den_ck2
!!c      if(den_ip<den_min) then
      if(den_ip-den_min<-1.0e-4*den_min) then
!!c         write(*,*) from,habit,sh_type
!!c         write(*,'("cal_bulk_density:error 2, bf, semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw
         V_cs=m_dry/den_min
         call cal_semiac_ip(is_mod(2),phi_cs,v_cs,semi_aip,semi_cip)
         den_ip=den_min

         alen=min(semi_aip,alen)
         clen=min(semi_cip,clen)
         phi_ic=clen/alen
         dlen=alen*psi_ic
         V_ic=coef4pi3*alen**3.0*(1.0_PS+phi_ic**2.0)**1.5

         ierror2=1
!!c         write(*,'("cal_bulk_density:error 2, af, semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw
      end if
    end subroutine den_ck2
    subroutine den_ck3
      if(ierror2==1) then
!!c      if(den>den_w.or.ierror2==1) then
!!c         write(*,*) from,habit,sh_type
!!c         write(*,'("cal_bulk_density:error 3, bf, semi_a,semi_c,den,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,den,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw
!!c         end if

         V_ip=get_vip(is_mod(2),phi_cs,semi_aip)
         V_space=max(0.0_PS,V_ip-m_dry/den_i)
         ak3=1.0_PS+max(m_mlt-V_space*den_w,0.0_PS)/V_ip
         V_csw=ak3*V_cs
         semi_a=ak3**(1.0/3.0)*semi_aip
         semi_c=semi_a*phi_cs
         den=mean_mass/V_csw
!!c         write(*,'("cal_bulk_density:error 3, af, semi_a,semi_c,den,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,den,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw
      end if
    end subroutine den_ck3
    subroutine den_ck4
      if(den-den_max>1.0e-4*den_max) then
!!c         write(*,*) from,habit,sh_type
!!c         write(*,'("cal_bulk_density:error 4, bf, semi_a,semi_c,den,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,den,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw

         V_csw=mean_mass/den_max
         call cal_semiac_ip(is_mod(2),phi_cs,v_csw,semi_a,semi_c)
         den=den_max

!!c         write(*,'("cal_bulk_density:error 4, af, semi_a,semi_c,den,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw",20ES15.6)') semi_a,semi_c,den,semi_aip,semi_cip,den_ip,mean_mass,m_ic,m_mlt,m_rim,V_ic,V_ip,V_space,V_cs,V_csw

!!c         if(V_space>1.0e-6) then
!!c         ! increase V_ip to make V_space to hold the water, then some of the melt water has to be decreased.
!!c         den_ism=m_dry/V_ip
!!c         if(mean_mass<=m_mlt.or.den_ism-den_i>den_i*1.0e-4) then
!!c            write(*,*) "den_ck4>woops",mean_mass,m_mlt,den_ism,V_ip
!!c         end if
!!c         
!!c         dV_ip=(mean_mass/den_i+(1.0_PS/den_w-1.0_PS/den_i)*m_mlt-V_ip)/&
!!c                     ((1.0_PS/den_w-1.0_PS/den_i)*den_ism+1.0_PS)
!!c         if(dV_ip<=0.0_PS) then
!!c            write(*,*) "den_ck4>something not right:",V_ip,dV_ip,mean_mass,m_mlt,den_ism,semi_aip,semi_cip,phi_cs
!!c!!c            stop
!!c            return
!!c         end if
!!c         m_mlt=m_mlt-den_ism*dV_ip
!!c         m_rim=m_rim+den_ism*dV_ip
!!c
!!c         V_ip=V_ip+dV_ip
!!c         semi_aip=(V_ip/get_coef_ip(is_mod(2))/phi_cs)**(1.0/3.0)
!!c         semi_cip=semi_aip*phi_cs
!!c         V_cs=get_vcs(is_mod(2),phi_cs,semi_aip)
!!c         V_csw=V_cs
!!c         semi_a=semi_aip
!!c         semi_c=semi_cip
!!c         end if  
!!c         den=mean_mass/V_csw
!!c         if(habit/=6) then
!!c         end if
!!c         if(den>den_w) then
!!c            write(*,*) "cal_bulk_density>something is wrong:",den,semi_a,semi_c,V_csw,V_ip,V_space,ak3,m_mlt
!!c         end if
      end if
    end subroutine den_ck4

  end subroutine cal_bulk_density3



  subroutine cal_bulk_density_ap(level,token,den,&
       deps_ap,den_ai,den_as)
    integer,intent(in) :: level,token
    real(PS),intent(in) :: deps_ap,den_ai,den_as
    real(PS),intent(inout) :: den
    
    ! aerosols
    den=den_ai/(1.0_PS-deps_ap*(1.0_PS-den_ai/den_as))
  end subroutine cal_bulk_density_ap

  subroutine cal_bulk_density_liq(level,mm,map,den_ai,den_as,eps_map,den)
    real(PS),intent(inout) :: den
    integer,intent(in) :: level
    real(PS),intent(in) :: mm,map,den_ai,den_as,eps_map

    ! bulk water density
    real(PS),parameter :: den_w=1.0
    real(PS) :: den_ap

    ! liquid hydrometeors
    if(mm<1.0e-30_PS) then
       den=1.0_PS
    else
       ! calculate density of mixed aerosols
       den_ap=den_ai/(1.0_PS-eps_map*(1.0_PS-den_ai/den_as))
       ! assume that v_d=v_w+v_n where v_w is volume of water and v_n is volume of dry aerosol.
       den=min(den_ap*0.95,mm/((mm-map)/den_w+map/den_ap))
!tmp       den=1.0_PS
    endif
  end subroutine cal_bulk_density_liq

  subroutine cal_aclen_liq(ms)
    ! this is for liquid hydrometeors
    type (Mass_Bin), intent(inout)   :: ms
    ! axis ratio, c-axis/a-axis
    real(PS)                           :: alpha
    ! density of water (g/cm^3)
    real (PS), parameter               :: den_w = 1.0_PS
    real(PS)                           :: mean_mass
    real(PS)                           :: dum

       ! +++ calculate the equivalent diameter +++
       ms%len = (6.0_PS*(ms%mean_mass/(PI*ms%den)))**(1.0/3.0)
       
       ! +++ calculate spheroidal shape 
       !     based on Pruppacher and Klett (1997, p.397) +++
       if( ms%len <= 280.0e-4 ) then
          alpha = 1.0_PS
          ms%a_len = ms%len/2.0_PS
          ms%c_len = ms%len/2.0_PS
       else if( 280.0e-4 < ms%len .and. ms%len <= 1.0e-1_PS ) then
          ! NOTE: This formula is assumed to be applicable for this range.
          ! It is supposed to be applicable to the range: from 1 mm to 9 mm.
          dum = ms%len

          alpha = max( 1.001668_PS - 0.098055_PS*dum &
               - 2.52686_PS*dum**2 + 3.75061_PS*dum**3 &
               - 1.68692_PS*dum**4, 1.0e-5_RP)

          dum=ms%len/((8.0_PS*alpha)**(1.0_RP/3.0_RP))
          ms%a_len = dum
          ms%c_len = dum * alpha 
       else if( 1.0e-1_PS < ms%len ) then
          ! NOTE: This formula is assumed to be applicable for this range.
          ! It is supposed to be applicable to the range: from 1 mm to 9 mm.
          ! But it is not right for more than 5mm; alpha becomes negative.
          dum = min(0.9_PS,ms%len)
          alpha = max( 1.001668_PS - 0.098055_PS*dum &
               - 2.52686_PS*dum**2 + 3.75061_PS*dum**3 &
               - 1.68692_PS*dum**4, 1.0e-5_RP)
          dum = ms%len/((8.0_PS*alpha)**(1.0_RP/3.0_RP))
          ms%a_len = dum
          ms%c_len = dum * alpha 
       end if
       ms%semi_a = ms%a_len
       ms%semi_c = ms%c_len
  end subroutine cal_aclen_liq

  subroutine cal_aclen(phase, ms)
    ! phase of hydrometeor
    ! 1: liquid, 2: solid, 3: aerosol
    integer, intent(in)  :: phase
    type (Mass_Bin), intent(inout)   :: ms
    ! axis ratio, c-axis/a-axis
    real(PS)                           :: alpha
    ! density of water (g/cm^3)
    real (PS), parameter               :: den_w = 1.0_PS
    real(PS)                           :: mean_mass
    real(PS)                           :: dum

    if( phase == 1 ) then
       ! initialize
       ms%len=0.0_PS
       ms%a_len=0.0_PS
       ms%c_len=0.0_PS
       ms%semi_a=0.0_PS
       ms%semi_c=0.0_PS
       if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return

       ! +++ calculate the equivalent diameter +++
       ms%len = (6.0_PS*(ms%mean_mass/(PI*ms%den)))**(1.0/3.0)
       
       ! +++ calculate spheroidal shape 
       !     based on Pruppacher and Klett (1997, p.397) +++
       if( ms%len <= 280.0e-4 ) then
          alpha = 1.0_PS
          ms%a_len = ms%len/2.0_PS
          ms%c_len = ms%len/2.0_PS
       else if( 280.0e-4 < ms%len .and. ms%len <= 1.0e-1_PS ) then
          ! NOTE: This formula is assumed to be applicable for this range.
          ! It is supposed to be applicable to the range: from 1 mm to 9 mm.
          dum = ms%len

          alpha = max( 1.001668_PS - 0.098055_PS*dum &
               - 2.52686_PS*dum**2 + 3.75061_PS*dum**3 &
               - 1.68692_PS*dum**4, 1.0e-5_RP)
          ms%a_len = ms%len/((8.0_PS*alpha)**(1.0_RP/3.0_RP))
          ms%c_len = ms%a_len * alpha 
       else if( 1.0e-1_PS < ms%len ) then
          ! NOTE: This formula is assumed to be applicable for this range.
          ! It is supposed to be applicable to the range: from 1 mm to 9 mm.
          ! But it is not right for more than 5mm; alpha becomes negative.
          dum = min(0.9_PS,ms%len)
          alpha = max( 1.001668_PS - 0.098055_PS*dum &
               - 2.52686_PS*dum**2.0 + 3.75061_PS*dum**3.0 &
               - 1.68692_PS*dum**4.0, 1.0e-5_RP)
          ms%a_len = ms%len/((8.0_PS*alpha)**(1.0_RP/3.0_RP))
          ms%c_len = ms%a_len * alpha 
       end if
       ms%semi_a = ms%a_len
       ms%semi_c = ms%c_len
    end if
  end subroutine cal_aclen

  subroutine diag_pardis_ap_v0(icat,flagp,ap_lnsig,ap_mean,ms)
!!c    use com_amps
    integer,intent(in) :: icat,flagp
    real(PS) :: ap_lnsig(*),ap_mean(*)
    type (Mass_Bin), intent(inout)   :: ms

    ! initialize
    ms%p=0.0_PS
    ms%a_len=0.0_PS
    ms%c_len=0.0_PS
    ms%len=0.0_PS
!!c    write(*,*) "in diag_pardis",ms%con,ms%mass(1)
    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return

!tmp    write(*,*) "diag_ap, ap_lnsig,ap_mean",ap_lnsig(1:2),ap_mean(1:2)
!!c    stop
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
          write(*,*) "This flag for aerosols not defined in diag_pardis_ap",flagp
          stop
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
       write(*,*) "This distribution type is not defined in diag_pardis_ap",ms%dis_type
       stop
    end if
  end subroutine diag_pardis_ap_v0

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

  subroutine mod_tendency( ms, n_mcom, n_vcom, process, mod_c)
    type (Mass_Bin), intent(inout)   :: ms
    integer, intent(in)   :: n_mcom, n_vcom, process
    real(PS)  :: mod_c
    integer   :: i
    do i = 1, (n_mcom+1)
       ms%dmassdt(i,process) = ms%dmassdt(i,process)*mod_c
    end do
    ms%dcondt(process) = ms%dcondt(process)*mod_c

    do i = 1, n_vcom
       ms%dvoldt(i, process) = ms%dvoldt(i, process)*mod_c
    end do
  end subroutine mod_tendency

  subroutine cal_meanmass( ms,em)
    type (Mass_Bin), intent(inout)   :: ms
    integer,intent(inout) :: em
    if(ms%con>1.0e-30_PS.and.ms%mass(1)>1.0e-30_PS) then
       ms%mean_mass=ms%mass(1)/ms%con

       if(ms%mean_mass==0.0_PS)then
          ms%mass=0.0_PS
          ms%con=0.0_PS
          em=10
       endif
    else
       ms%mass=0.0_PS
       ms%con=0.0_PS
       ms%mean_mass=0.0_PS
    end if
  end subroutine cal_meanmass

  subroutine cal_meanmass_vec( ms,em)
    type (Mass_Bin), intent(inout)   :: ms
    integer,intent(inout) :: em

    ms%mean_mass=ms%mass(1)/ms%con

    if(ms%mean_mass==0.0_PS)then
      em=10
    endif
  end subroutine cal_meanmass_vec

  subroutine cal_wvt_numint(phase, ms, th_var, mode, level, &
       a, m1, m2, NDIS, WVT,ishape)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate the bin terminal velocity
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! phase of hydrometeor
    integer, intent(in)  :: phase
    type (Mass_Bin), intent(in)   :: ms
    type (Thermo_Var), intent(in)    :: th_var
    ! number of discretization
    integer,intent(in) :: NDIS

    ! parameter of linear distribution
    real(PS),pointer,dimension(:)     :: a
    type (Ice_Shape), optional       :: ishape
    ! dummy mass_bin object
    type (Mass_Bin)                  :: msdum
    type (Ice_Shape)                 :: isdum
    real(PS), pointer, dimension(:)  :: VT
    real(PS), pointer, dimension(:)  :: WVT
    real(PS), pointer, dimension(:)  :: mp
    real(PS) :: m1,m2, phi, dalen, dclen
    real(PS) :: n1, n2, val1, val2
    real(PS) :: dm
    real(PS), dimension(2)  :: sum_x

    ! +++ mode of shape to calculate the ventilation coefficient 
    !     and terminal velocity +++
    ! 1 : spheroid with a_len and c_len
    ! 2 : ice crystal model
    integer,intent(in)   :: mode

    ! level of complexity
    integer,intent(in)   :: level
    integer  :: i, var_Status
    integer :: em


    allocate( VT(NDIS), stat = var_Status)
    if(var_Status /= 0 ) stop "Memory not available for VT  in vapor_deposition"
        
    allocate( mp(NDIS), stat = var_Status)
    if(var_Status /= 0 ) stop "Memory not available for mp in vapor_deposition"
         

    dm = (m2 - m1)/real(NDIS,PS_KIND)
    VT=0.0_PS

!    msdum = make_Mass_Bin2()
    msdum = ms
    msdum%con=1.0_PS

    if( phase == 2 ) then
!!c    is_dum = make_Ice_Shape_zero ()
!tmp       isdum = make_Ice_Shape_zero () 
       call make_Ice_Shape_zero (isdum) 
       isdum = ishape
    end if
    do i = 1, NDIS
       mp(i) = m1 + real(i-1,PS_KIND)*dm
       msdum%mass(1)=mp(i)
       msdum%mean_mass=mp(i)

       if( phase == 1 ) then

          call cal_aclen(phase, msdum)
          call cal_terminal_vel( phase, msdum, th_var, mode, level,em)
       else if( phase == 2 ) then

!!c          isdum%V_ic=ishape%V_ic*(mp(i)/ms%mean_mass)

          msdum%a_len = ( mp(i)/get_mci2(1,ishape%phi_ic,ms%den))**(1.0/3.0)   
               
          isdum%c = isdum%phi_ic*msdum%a_len 
          isdum%d = isdum%psi_ic*msdum%a_len 
          isdum%a = msdum%a_len - isdum%d
          msdum%c_len = isdum%c

          if(ishape%sh_type<=4) then
             msdum%semi_a = (mp(i)/ms%den/(coef4pi3*(1.0_PS+isdum%phi_cs**2)**1.5))**(1.0/3.0)
             msdum%semi_c = msdum%semi_a*isdum%phi_cs
          else
             if(ishape%phi_cs<1.0_PS) then
                msdum%semi_a = (mp(i)/ms%den/coef4pi3)**(1.0/3.0)
                msdum%semi_c = msdum%semi_a*isdum%phi_cs
             else
                msdum%semi_c = (mp(i)/ms%den/coef4pi3)**(1.0/3.0)
                msdum%semi_a = msdum%semi_c/isdum%phi_cs
             end if
          end if 
          call cal_terminal_vel( phase, msdum, th_var, mode, level,em, isdum)
       end if

       VT(i) = msdum%vtm

    end do

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! numerically integrate terminal velocity
    WVT = 0.0_PS
    sum_x=0.0_PS

    do i = 1, (NDIS-1)
       n1 = max( a(1), 0.0_PS) + a(3)*( mp(i) - a(2))
       n2 = max( a(1), 0.0_PS) + a(3)*( mp(i+1) - a(2))
       val1 = n1*mp(i)*VT(i)
       val2 = n2*mp(i+1)*VT(i+1)          
       WVT(1) = WVT(1) + (val1+val2)*dm/2.0_PS

       val1 = n1*VT(i)
       val2 = n2*VT(i+1)          
       WVT(2)=WVT(2) + (val1+val2)*dm/2.0_PS

       sum_x(1)=sum_x(1)+(n1*mp(i)+n2*mp(i+1))*dm/2.0_PS
       sum_x(2)=sum_x(2)+(n1+n2)*dm/2.0_PS
    end do
    WVT(1) = WVT(1)/sum_x(1)
    WVT(2) = WVT(2)/sum_x(2)

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!c    call delete_Mass_Bin (msdum)
    deallocate( mp,stat=var_Status)
    if(var_Status /= 0 ) stop "Memory not deallocated for mp in vapor_deposition"
         
    nullify(mp)
    deallocate( VT,stat=var_Status)
    if(var_Status /= 0 ) stop "Memory not deallocated for VT in vapor_deposition"
        
    nullify(VT)
  end subroutine cal_wvt_numint

  subroutine cal_wvt_numint_v2(phase, ms, th_var, mode, level, &
       a, m1, m2, NDIS, WVT,ishape)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate the bin terminal velocity
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! phase of hydrometeor
    integer, intent(in)  :: phase
    type (Mass_Bin), intent(in)   :: ms
    type (Thermo_Var), intent(in)    :: th_var
    ! number of discretization
    integer,intent(in) :: NDIS

    ! parameter of linear distribution
    real(PS),pointer,dimension(:)     :: a
    type (Ice_Shape), optional       :: ishape
    real(PS), pointer, dimension(:)  :: VT
    real(PS), pointer, dimension(:)  :: WVT
    real(PS), pointer, dimension(:)  :: mp
    real(PS),intent(in) :: m1,m2
    real(PS) :: A1,B1,C1,D1,alpha,beta,gamma,N_re_0,ratio,char_len,C_DP,C_DO,X,&
         A_e,A_c,P,q,OMEGA,rad,U_S,lambda_a,Y,NP,NBO,k
    real(PS), parameter         :: gg = 980.0,X_0=2.8e+06
    real(PS), parameter         :: den_w = 1.0
    real(PS) :: n1, n2, val1, val2
    real(PS) :: dm
    real(PS), dimension(2)  :: sum_x

    ! +++ mode of shape to calculate the ventilation coefficient 
    !     and terminal velocity +++
    ! 1 : spheroid with a_len and c_len
    ! 2 : ice crystal model
    integer,intent(in)   :: mode

    ! level of complexity
    integer,intent(in)   :: level
    integer  :: i, var_Status
    integer :: em


    allocate( VT(NDIS), stat = var_Status)
    if(var_Status /= 0 ) stop "Memory not available for VT in vapor_deposition"
         
    allocate( mp(NDIS), stat = var_Status)
    if(var_Status /= 0 ) stop "Memory not available for mp in vapor_deposition" 
         

    dm = (m2 - m1)/real(NDIS,PS_KIND)
    VT=0.0_PS

    if( phase == 1 ) then
       char_len=0.5_PS*(ms%mean_mass/coefpi6)**(1.0_RP/3.0_RP)
       do i = 1, NDIS
          mp(i) = m1 + real(i-1,PS_KIND)*dm

          rad=char_len*(mp(i)/ms%mean_mass)**(1.0_RP/3.0_RP)
          if(0.5e-4_PS<=rad.and.rad<10.0e-4_PS) then
             ! calculate Stokes terminal velocity
             U_S=rad**2*gg*(den_w-th_var%den_a)/4.5_PS/th_var%d_vis
             lambda_a=6.6e-6_PS*(th_var%d_vis/1.818e-4_PS)*(1013250.0_PS/th_var%P)*(th_var%T/293.15_PS)
             VT(i)=(1.0_PS+1.26_PS*lambda_a/rad)*U_S
          else if(rad<535.0e-4_PS) then
             X=log(32.0_PS*rad**3*(den_w-th_var%den_a)*th_var%den_a*gg/3.0_PS/th_var%d_vis**2)
             Y=-0.318657e+1_PS+0.992696_PS*X-0.153193e-2_PS*X**2-0.987059e-3_PS*X**3 &
                  -0.578878e-3_PS*X**4+0.855176e-4*X**5-0.327815e-5*X**6
             A1=exp(Y)
             ! +++ calculate the terminal velocity +++
             VT(i)=A1*th_var%d_vis/(2.0_PS*rad*th_var%den_a)
          else
             rad=min(rad,3500.0e-4_RP)
             
             NBO=gg*(den_w-th_var%den_a)*rad**2/th_var%sig_wa
             NP=th_var%sig_wa**3*th_var%den_a**2/th_var%d_vis**4/gg
             
             X=log(NBO*NP**(1.0/6.0)*16.0_PS/3.0_PS)
             Y=-0.500015e+1_PS+0.523778e+1_PS*X-0.204914e+1_PS*X**2+0.475294_PS*X**3 &
                  -0.542819e-1_PS*X**4+0.238449e-2_PS*X**5
             
             A1=NP**(1.0/6.0)*exp(Y)
             ! +++ calculate the terminal velocity +++
             VT(i)=A1*th_var%d_vis/(2.0_PS*rad*th_var%den_a)
             
          end if
       end do
    else if( phase == 2 ) then
             
       ! +++ calculate some important variables, depending on level +++
       call cal_some_vel(ms,ishape,mode,alpha,OMEGA,A_e,A_c,P,q,char_len)
       
       ! +++ calculate the shape factor +++
       if( alpha <= 0.16_PS ) then
          k = 0.85_PS
       else if ( 0.16_PS < alpha .and. alpha <= 1.0_PS ) then
          k = 0.82_PS + 0.18_PS*alpha
       else if ( 1.0_PS < alpha .and. alpha <= 3.0_PS ) then
          k = 0.37_PS + 0.63_PS/alpha
       else if ( 3.0_PS < alpha ) then
          k = 1.33_PS/( log(alpha) + 1.19_PS)
       end if
       
       
       ! +++ calculate the pressure drag C_DP +++
       if( 0.3_PS < alpha .and. alpha < 1.0_PS ) then
          gamma = 3.76_PS - 8.41*alpha + 9.18*alpha**2.0 - 3.53*alpha**3.0
          C_DP = 0.292*k*gamma
       else if( alpha <= 0.3_PS ) then 
          gamma = 1.98
          C_DP = 0.292*k*gamma
       else if( 1.0_PS <= alpha ) then
          if( 1.0_PS <= q ) then
             C_DP = (1.46*q-0.46)*(0.492-0.200/sqrt(alpha))
          else
!!c             gamma = 3.76_PS - 8.41*alpha + 9.18*alpha**2.0 - 3.53*alpha**3.0
             gamma = 1.0
             C_DP = 0.292*k*gamma
          end if
       end if
       
       ! +++ calculate the Oseen-type drag coefficient, C_DO +++
       C_DO = 4.5*(k**2.0)*max(alpha, 1.0_PS)
       ! +++ matching theory +++
       gamma = (C_DO-C_DP)/(4.0_PS*C_DP)


       ! +++ calculate the weighted terminal velocity +++
       ! coef. from Reynolds number
       A1=6.0*k*th_var%d_vis/(C_DP*th_var%den)
       ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! coef. from X
       if( q <= 1.0_PS ) then
          D1 = 8.0_PS*gg*th_var%den/&
               (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*(q**0.25))
       else if( q > 1.0_PS ) then
          D1 = 8.0_PS*gg*th_var%den/&
               (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*q)
       end if
       C1=sqrt(C_DP*D1)/(6.0_PS*k)
 
      

       do i = 1, NDIS
          mp(i) = m1 + real(i-1,PS_KIND)*dm

          ! +++ calculate the Davies or Best number +++
          X=D1*mp(i)

          ! +++ calculate the Reynolds number +++
          beta = sqrt(1.0_PS + (C_DP/(6.0_PS*k))*sqrt(X/C_DP)) - 1.0_PS
          N_re_0 = 6.0*k*beta**2.0/C_DP

          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! coef. from matching theory
          B1=(1.0_PS+2.0_PS*beta*exp(-beta*gamma)/&
               ((2.0_PS+beta)*(1.0_PS+beta)))

          ! +++ correction for transition to turbulent flow +++
          if( 1000.0 <= N_re_0 ) then
             ratio=(1.0_PS + 1.6_PS*(X/X_0)**2.0)/(1.0_PS + (X/X_0)**2.0)
             A1=A1/(ratio)
             C1=C1*sqrt(ratio)
!!c          C1=C1*sqrt((1.0_PS + 1.6_PS*(X/X_0)**2.0)/(1.0_PS + (X/X_0)**2.0))
          end if
          
          
          VT(i)=A1*B1/(char_len*(mp(i)/ms%mean_mass)**(1.0/3.0))*(&
               2.0_PS+C1*sqrt(mp(i))-2.0_PS*sqrt(1.0_PS+C1*sqrt(mp(i))))
       end do
!!c       if(NDIS<=5) then
!!c          open( unit=21, file="check_coef_svtm.dat",POSITION='APPEND')
!!c          do i = 1, NDIS
!!c             mp(i) = m1 + real(i-1)*dm
!!c             write(21,'(10ES15.6)') ms%mean_mass,mp(i),C1,C1*sqrt(mp(i)),1.0+C1*sqrt(mp(i)),VT(i)
!!c          end do
!!c          close(21)
!!c       end if
    end  if
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! numerically integrate terminal velocity
    WVT = 0.0_PS
    sum_x=0.0_PS

    do i = 1, (NDIS-1)
       n1 = max( a(1), 0.0_PS) + a(3)*( mp(i) - a(2))
       n2 = max( a(1), 0.0_PS) + a(3)*( mp(i+1) - a(2))
       WVT(1) = WVT(1) + (n1*mp(i)*VT(i)+n2*mp(i+1)*VT(i+1))*dm/2.0_PS
       WVT(2)=WVT(2) + (n1*VT(i)+n2*VT(i+1))*dm/2.0_PS
       sum_x(1)=sum_x(1)+(n1*mp(i)+n2*mp(i+1))*dm/2.0_PS
       sum_x(2)=sum_x(2)+(n1+n2)*dm/2.0_PS
    end do
    WVT(1) = WVT(1)/sum_x(1)
    WVT(2) = WVT(2)/sum_x(2)

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!c    call delete_Mass_Bin (msdum)
    deallocate( mp,stat=var_Status)
    if(var_Status /= 0 ) stop "Memory not deallocated for mp in vapor_deposition"
        
    nullify(mp)
    deallocate( VT,stat=var_Status)
    if(var_Status /= 0 ) stop "Memory not deallocated for VT in vapor_deposition"
         
    nullify(VT)
  end subroutine cal_wvt_numint_v2


  subroutine cal_wterm_vel_v3(phase, ms, th_var, mode, level, &
       a, m1, m2, WVT,ishape)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate the terminal velocity weighted by concentration and mass.
    ! use the approximation of NRe= a*X**b fitted to results by Bohm equation.
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! phase of hydrometeor
    integer, intent(in)  :: phase
    type (Mass_Bin), intent(in)   :: ms
    type (Thermo_Var), intent(in)    :: th_var

    ! parameter of linear distribution
    real(DS),dimension(*)     :: a
    real(DS), dimension(*)  :: WVT
    type (Ice_Shape), optional       :: ishape
    !
    real(DS),intent(in) :: m1,m2
    real(PS) :: alpha,char_len,OMEGA,P,&
         A_e,A_c,q,X1,X2
    real(DS) :: A1,B1,D1,Q1,M,N,MT,NT
    real(DS), parameter         :: gg = 980.0d+0
    real(DS), parameter         :: den_w = 1.0d+0
    real(PS), dimension(4) :: a_rex,b_rex

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

    ! level of complexity
    integer,intent(in)   :: level
    integer  :: i, var_Status
    integer :: em

    real(PS), parameter :: chiarui_1=0.0_PS, chiarui_2=1.0e+30


    WVT(1:2)=0.0_DS
    NT=0.0_DS
    MT=0.0_DS


    if( phase == 1 ) then
       ! --- in case of liquid phase ---
       a_rex=(/la1,la2,la3,la4/)
       b_rex=(/lb1,lb2,lb3,lb4/)
       X_b=(/chiarui_1,LXlim1,LXlim2,LXlim3,chiarui_2/)

       ! calculate Best number for the bin limits
       ! +++ calculate the axial ratio +++
       alpha=ms%c_len/ms%a_len

       ! +++ calculate the effective cross section +++
       A_e = PI * ms%a_len*ms%a_len
          
       ! +++ calculate the circumscribed cross-sectional area +++
       A_c = A_e
       
       q = 1.0_PS

       char_len=2.0_PS*ms%a_len

       ! +++ calculate the Davies or Best number +++
       D1=8.0_DS*gg*th_var%den/&
            (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*(q**0.25))

       Q1=q**0.25

       X1=D1*m1
       X2=D1*m2


       do i=1,4

          if(X_b(i+1)<X1) cycle
          if(X_b(i)>X2) exit

          m_sub1=max(m1,real(X_b(i)/D1,DS))
          m_sub2=min(m2,real(X_b(i+1)/D1,DS))

          D_1=(m_sub1/coedpi6)**(1.0_RP/3.0_RP)
          D_2=(m_sub2/coedpi6)**(1.0_RP/3.0_RP)

          N=(m_sub2-m_sub1)*(a(1)-a(3)*(a(2)-(m_sub2+m_sub1)/2.0))
          M=a(1)*(m_sub2**2.0_RP-m_sub1**2.0_RP)/2.0_RP &
               +a(3)*((m_sub2**3-m_sub1**3)/3.0_RP-a(2)*(m_sub2**2-m_sub1**2)/2.0_RP)


          A1=th_var%d_vis**(1.0_RP-2.0_RP*b_rex(i))*th_var%den_a**(b_rex(i)-1.0_RP)*&
               (8.0_PS*gg/PI)**b_rex(i)*(Q1*max(alpha,1.0_RP))**(-b_rex(i))*&
               ms%mean_mass**(1.0_RP/3.0_RP)/char_len*a_rex(i)

          B1=b_rex(i)-1.0_DS/3.0_DS

          if(D_2<D_max) then
             ! 1. mass-weighted terminal velocity
             WVT(1)=WVT(1)+max(0.0_DS,&
                  A1*( (a(1)-a(3)*a(2))/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS))+&
                  a(3)/(B1+3.0_DS)*(m_sub2**(B1+3.0_DS)-m_sub1**(B1+3.0_DS))))
             
             ! 2. con-weighted terminal velocity
             WVT(2)=WVT(2)+max(0.0_DS,&
                  A1*( (a(1)-a(3)*a(2))/(B1+1.0_DS)*(m_sub2**(B1+1.0_DS)-m_sub1**(B1+1.0_DS))+&
                  a(3)/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS))))
          elseif(D_1<D_max) then

             v_max=A1*m_max**B1
             
             ! 1. mass-weighted terminal velocity
             WVT(1)=WVT(1)+max(0.0_DS,&
                  A1*( (a(1)-a(3)*a(2))/(B1+2.0_DS)*(m_max**(B1+2.0_DS)-m_sub1**(B1+2.0_DS))+&
                  a(3)/(B1+3.0_DS)*(m_max**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)))+&
                  v_max*( (a(1)-a(3)*a(2))/2.0_DS*(m_sub2**2.0-m_max**2.0)+&
                  a(3)/3.0_DS*(m_sub2**3.0-m_max**3.0)))
                  
             ! 2. con-weighted terminal velocity
             WVT(2)=WVT(2)+max(0.0_DS,&
                  A1*( (a(1)-a(3)*a(2))/(B1+1.0_DS)*(m_max**(B1+1.0_DS)-m_sub1**(B1+1.0_DS))+&
                  a(3)/(B1+2.0_DS)*(m_max**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)))+&
                  v_max*( (a(1)-a(3)*a(2))*(m_sub2-m_max)+&
                  a(3)/2.0_DS*(m_sub2**2.0-m_max**2.0)))
          else
             v_max=A1*m_max**B1
             
             ! 1. mass-weighted terminal velocity
             WVT(1)=WVT(1)+max(0.0_DS,&
                  v_max*( (a(1)-a(3)*a(2))/2.0_DS*(m_sub2**2.0-m_sub1**2.0)+&
                  a(3)/3.0_DS*(m_sub2**3.0-m_sub1**3.0)))
                  
             ! 2. con-weighted terminal velocity
             WVT(2)=WVT(2)+max(0.0_DS,&
                  v_max*( (a(1)-a(3)*a(2))*(m_sub2-m_sub1)+&
                  a(3)/2.0_DS*(m_sub2**2.0-m_sub1**2.0)))

          end if
          NT=NT+N
          MT=MT+M
       end do
       WVT(1)=WVT(1)/MT
       WVT(2)=WVT(2)/NT
!!c       WVT(1)=max(min(WVT(1)/MT,ms%vtm*2.0),ms%vtm*0.5)
!!c       WVT(2)=max(min(WVT(2)/NT,ms%vtm*2.0),ms%vtm*0.5)

    else if( phase == 2 ) then
       ! --- in case of liquid phase ---
       a_rex=(/sa1,sa2,sa3,sa4/)
       b_rex=(/sb1,sb2,sb3,sb4/)
       X_b=(/chiarui_1,SXlim1,SXlim2,SXlim3,chiarui_2/)

       ! calculate Best number for the bin limits
       call cal_some_vel(ms,ishape,mode,alpha,OMEGA,A_e,A_c,P,q,char_len)
       
       ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! coef. from X
       if( q <= 1.0_PS ) then
          D1 = 8.0_PS*gg*th_var%den/&
               (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*(q**0.25))
          Q1=q**0.25
       else if( q > 1.0_PS ) then
          D1 = 8.0_PS*gg*th_var%den/&
               (PI*(th_var%d_vis**2.0)*max(alpha, 1.0_PS)*q)
          Q1=q
       end if

       X1=D1*m1
       X2=D1*m2

       do i=1,4

          if(X_b(i+1)<X1) cycle
          if(X_b(i)>X2) exit

          m_sub1=max(m1,real(X_b(i)/D1,DS))
          m_sub2=min(m2,real(X_b(i+1)/D1,DS))

          N=(m_sub2-m_sub1)*(a(1)-a(3)*(a(2)-(m_sub2+m_sub1)/2.0_RP))
          M=a(1)*(m_sub2**2-m_sub1**2.0)/2.0_RP &
               +a(3)*((m_sub2**3-m_sub1**3)/3.0_RP-a(2)*(m_sub2**2-m_sub1**2)/2.0_RP)


          A1=th_var%d_vis**(1.0_RP-2.0_RP*b_rex(i))*th_var%den_a**(b_rex(i)-1.0_RP)*&
               (8.0_DS*gg/PI)**b_rex(i)*(Q1*max(alpha,1.0_RP))**(-b_rex(i))*&
               ms%mean_mass**(1.0_RP/3.0_RP)/char_len*a_rex(i)

          B1=b_rex(i)-1.0_DS/3.0_DS

          ! 1. mass-weighted terminal velocity
          WVT(1)=WVT(1)+max(0.0_DS,&
               A1*( (a(1)-a(3)*a(2))/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS))+&
               a(3)/(B1+3.0_DS)*(m_sub2**(B1+3.0_DS)-m_sub1**(B1+3.0_DS))))
          
          ! 2. con-weighted terminal velocity
          WVT(2)=WVT(2)+max(0.0_DS,&
               A1*( (a(1)-a(3)*a(2))/(B1+1.0_DS)*(m_sub2**(B1+1.0_DS)-m_sub1**(B1+1.0_DS))+&
               a(3)/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS))))

          NT=NT+N
          MT=MT+M
       end do
       WVT(1)=WVT(1)/MT
       WVT(2)=WVT(2)/NT
    end  if


  end subroutine cal_wterm_vel_v3


  subroutine cal_wvt_intpol(phase, ms, th_var, mode, level, &
       a, m1, m2, WVT,ishape)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate the bin terminal velocity by interpolation function
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! phase of hydrometeor
    integer, intent(in)  :: phase
    type (Mass_Bin), intent(in)   :: ms
    type (Thermo_Var), intent(in)    :: th_var
    ! parameter of linear distribution
    real(PS),pointer,dimension(:)     :: a
    type (Ice_Shape), optional       :: ishape
    ! dummy mass_bin object
    type (Mass_Bin)                  :: msdum
    type (Ice_Shape)                 :: isdum
    real(PS), pointer, dimension(:)  :: VT
    real(PS), pointer, dimension(:)  :: WVT
    real(PS), dimension(2)  :: mp
    real(PS) :: m1,m2, phi, dalen, dclen,A1,Am,A2
    real(PS) :: n1, n2, nm, val1, val2, valm
    real(PS) :: dm
    real(PS), dimension(2)  :: sum_x


    ! +++ mode of shape to calculate the ventilation coefficient 
    !     and terminal velocity +++
    ! 1 : spheroid with a_len and c_len
    ! 2 : ice crystal model
    integer,intent(in)   :: mode

    ! level of complexity
    integer,intent(in)   :: level
    integer  :: i, var_Status
    integer :: em

    ! number of discretization
    integer  :: NDIS


    NDIS=2
    allocate( VT(NDIS), stat = var_Status)
    if(var_Status /= 0 ) stop "Memory not available for VT in vapor_deposition"
        

    VT=0.0_PS

    mp(1)=m1
    mp(2)=m2

!    msdum = make_Mass_Bin2()
    msdum = ms
    msdum%con=1.0_PS

    if( phase == 2 ) then
!!c    is_dum = make_Ice_Shape_zero ()
!tmp       isdum = make_Ice_Shape_zero () 
       call make_Ice_Shape_zero (isdum) 
       isdum = ishape
    end if
    do i=1,NDIS
       msdum%mass(1)=mp(i)
       msdum%mean_mass=mp(i)

       if( phase == 1 ) then

          call cal_aclen(phase, msdum)
          call cal_terminal_vel( phase, msdum, th_var, mode, level,em)
       else if( phase == 2 ) then

          if(ishape%habit<=3) then
             msdum%a_len=(mp(i)/ms%mean_mass)**(1.0/3.0)*ms%a_len
             msdum%c_len=isdum%phi_ic*msdum%a_len 
             isdum%c=msdum%c_len
             isdum%d=isdum%psi_ic*msdum%a_len 
             isdum%a=msdum%a_len-isdum%d
             isdum%r=isdum%gam_ic*msdum%a_len
             isdum%e=isdum%eta_ic*msdum%a_len
          elseif(ishape%habit==4) then
             ! rosette
             isdum%r=(mp(i)/ms%mean_mass)**(1.0/3.0)*ishape%r
             msdum%a_len=isdum%r/isdum%gam_ic
             msdum%c_len=isdum%phi_ic*msdum%a_len 
             isdum%c=msdum%c_len
             isdum%d=isdum%psi_ic*msdum%a_len 
             isdum%a=msdum%a_len-isdum%d
             isdum%e=isdum%eta_ic*msdum%a_len
          elseif(ishape%habit>=5) then
             ! irregular crystal
             isdum%e=(mp(i)/ms%mean_mass)**(1.0/3.0)*ishape%e
             msdum%a_len=isdum%e/isdum%eta_ic
             msdum%c_len=isdum%phi_ic*msdum%a_len 
             isdum%c=msdum%c_len
             isdum%d=isdum%psi_ic*msdum%a_len 
             isdum%a=msdum%a_len-isdum%d
             isdum%r=isdum%gam_ic*msdum%a_len
          end if
          msdum%semi_a=(mp(i)/ms%mean_mass)**(1.0/3.0)*ms%semi_a
          msdum%semi_c=msdum%semi_a*isdum%phi_ic

          call cal_terminal_vel( phase, msdum, th_var, mode, level,em, isdum)
       end if

       VT(i) = msdum%vtm

    end do

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! fit 2nd order polynomial and integrate over the region.
    WVT = 0.0_PS
    sum_x=0.0_PS

    n1=max( a(1), 0.0_PS) + a(3)*( mp(1) - a(2))
    n2=max( a(1), 0.0_PS) + a(3)*( mp(2) - a(2))
    nm=max( a(1), 0.0_PS) + a(3)*( ms%mean_mass - a(2))

    val1 = n1*mp(1)*VT(1)
    val2 = n2*mp(2)*VT(2)          
    valm = nm*ms%mean_mass*ms%vtm

!!c    A1=val1/(mp(1)-ms%mean_mass)/(mp(1)-mp(2))
!!c    Am=valm/(ms%mean_mass-mp(1))/(ms%mean_mass-mp(2))
!!c    A2=val2/(mp(2)-ms%mean_mass)/(mp(2)-mp(1))
!!c
!!c    WVT(1)=((mp(2)**3.0-mp(1)**3.0)*(A1+Am+A2)/3.0_PS&
!!c         -(mp(2)**2.0-mp(1)**2.0)*&
!!c         (A1*(ms%mean_mass+mp(2))+Am*(mp(1)+mp(2))+A2*(mp(1)+ms%mean_mass))/2.0_PS&
!!c         +(mp(2)-mp(1))*(A1*ms%mean_mass*mp(2)+Am*mp(1)*mp(2)+A2*mp(1)*ms%mean_mass))/ms%mass(1)
    WVT(1)=(val1/(ms%mean_mass-mp(1))*(&
            (mp(2)**2.0+mp(2)*mp(1)+mp(1)**2.0)/3.0-(ms%mean_mass+mp(2))*(mp(2)+mp(1))/2.0+ms%mean_mass*mp(2))+&
            valm*(mp(2)-mp(1))/(ms%mean_mass-mp(1))/(ms%mean_mass-mp(2))*&
              ((mp(2)**2.0+mp(2)*mp(1)+mp(1)**2.0)/3.0-(mp(2)+mp(1))*(mp(2)+mp(1))/2.0+mp(1)*mp(2))+&
            val2/(mp(2)-ms%mean_mass)*&
              ((mp(2)**2.0+mp(2)*mp(1)+mp(1)**2.0)/3.0-&
                  (ms%mean_mass+mp(1))*(mp(2)+mp(1))/2.0+ms%mean_mass*mp(1)))/ms%mass(1)




    val1 = n1*VT(1)
    val2 = n2*VT(2)          
    valm = nm*ms%vtm

!!c    A1=val1/(mp(1)-ms%mean_mass)/(mp(1)-mp(2))
!!c    Am=valm/(ms%mean_mass-mp(1))/(ms%mean_mass-mp(2))
!!c    A2=val2/(mp(2)-ms%mean_mass)/(mp(2)-mp(1))
!!c
!!c    WVT(2)=((mp(2)**3.0-mp(1)**3.0)*(A1+Am+A2)/3.0_PS&
!!c         -(mp(2)**2.0-mp(1)**2.0)*&
!!c         (A1*(ms%mean_mass+mp(2))+Am*(mp(1)+mp(2))+A2*(mp(1)+ms%mean_mass))/2.0_PS&
!!c         +(mp(2)-mp(1))*(A1*ms%mean_mass*mp(2)+Am*mp(1)*mp(2)+A2*mp(1)*ms%mean_mass))/ms%con

    WVT(2)=(val1/(ms%mean_mass-mp(1))*((mp(2)**2.0+mp(2)*mp(1)+mp(1)**2.0)/3.0-&
               (ms%mean_mass+mp(2))*(mp(2)+mp(1))/2.0+ms%mean_mass*mp(2))+&
            valm*(mp(2)-mp(1))/(ms%mean_mass-mp(1))/(ms%mean_mass-mp(2))*&
               ((mp(2)**2.0+mp(2)*mp(1)+mp(1)**2.0)/3.0-(mp(2)+mp(1))*(mp(2)+mp(1))/2.0+mp(1)*mp(2))+&
            val2/(mp(2)-ms%mean_mass)*((mp(2)**2.0+mp(2)*mp(1)+mp(1)**2.0)/3.0-&
               (ms%mean_mass+mp(1))*(mp(2)+mp(1))/2.0+ms%mean_mass*mp(1)))/ms%con



    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!c    call delete_Mass_Bin (msdum)
    deallocate( VT,stat=var_Status)
    if(var_Status /= 0 ) stop "Memory not deallocated for VT in vapor_deposition"
        
    nullify(VT)
  end subroutine cal_wvt_intpol

  function get_mci1(sh_type,phi,den) result(mci)
    real(PS),intent(in) :: phi,den
    integer, intent(in) :: sh_type
    ! coefficient of mass-dimension relation, m = mci * D^qq
    real(PS) :: mci
    select case(sh_type)
    case (1)
       ! pristine crystal
       ! for cylinder volume
       mci = phi*PI*den/4.0_PS
    case (2)
       ! rimed pristine crystal
       ! for cylinder volume
       mci = phi*PI*den/4.0_PS
    case default
       ! for spheroid
       mci = (PI/6.0_PS)*phi*den
    end select
  end function get_mci1
  function get_mci2(sh_type,phi,den) result(mci)
    real(PS),intent(in) :: phi,den
    integer, intent(in) :: sh_type
    ! coefficient of mass-dimension relation, m = mci * r ^qq
    real(PS) :: mci
    mci = coef4pi3*phi*den
  end function get_mci2

  subroutine cal_cseff(ms, th_var,y)
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Calculate curvature and solute contribution
    ! using (13-27b) Pruppacher and Klett (1996).
    !
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    use com_amps
    implicit none
    type (Mass_Bin), intent(inout) :: ms
    type (Thermo_Var), intent(in) :: th_var
    real(PS), parameter :: den_w=1.0_PS
    real(PS) :: mean_mass_s,y

    mean_mass_s=ms%mass(3)/ms%con

!!c    if(mean_mass_s/ms%mean_mass>0.99999) then
!!c       y=2.0_PS*th_var%sig_wa/(R_v*th_var%T*den_w*ms%len*0.5_PS)
    y=2.0_PS*th_var%sig_wa/(R_v*th_var%T*den_w*ms%len*0.5_PS)-&
         nu_aps(1)*phi_aps(1)*mean_mass_s*M_W/M_aps(1)/(ms%mean_mass-mean_mass_s)

    ! for now the solute effect is neglected because  I do not have
    ! the skill to predict the water right after activation.         


         
  end subroutine cal_cseff


  subroutine cal_cond_mlt_grap(ms,th_var,a_tot,a_ice,T_s,dmdt,ick)
    implicit none
    type (Mass_Bin), intent(in) :: ms
    type (Thermo_Var), intent(in) :: th_var
    real(PS),intent(in) :: a_tot,a_ice
    real(PS) :: T_s,T_sb,dmdt
    integer :: iter,ick

    ! latent heat of freeze in ergs/g
    real (PS), parameter          :: L_f = 0.3337e+10_PS
    ! latent heat of condensation in ergs/g
    real (PS), parameter              :: L_e = 2.5e+10_PS
    real(PS), parameter              :: R_u = 8.3143e+07
    real (PS), parameter             :: M_w=18.0
    ! thermal conductivity
    real(PS),parameter :: k_w=0.58e+5    

    ! temperature at triple point 
    real (PS), parameter          :: T_0 = 273.16_PS

!    T_s=th_var%T
    T_s=T_0
    T_sb=T_s
    iter=0
    do 
       iter=iter+1
       T_s=T_0+(th_var%k_a*(th_var%T-T_s)*ms%fh+&
            th_var%D_v*L_e*M_w/R_u*(th_var%e/th_var%T-get_sat_vapor_pres(1, T_s )/T_s)*ms%fv)*(a_tot-a_ice)/a_ice/k_w
  
       if(abs(T_s-th_var%T)>50.0.or.iter==20) then
          ick=1
          dmdt=0.0
          exit          
       elseif(abs(T_s-T_sb)<0.001) then
          ick=0
          dmdt=4.0_PS*PI*a_tot*a_ice/(a_tot-a_ice)*k_w*(T_s-T_0)/L_f
          exit
       end if
       T_sb=T_s
    end do
  end subroutine cal_cond_mlt_grap

  subroutine cal_growth_mode(level,mes_rc,ms,is, th_var)
!tmp    use mod_amps_utility, only: get_growth_mode
    use mod_amps_utility, only: get_growth_mode_max
    integer,intent(in) :: level,mes_rc
    type (Mass_Bin), intent(in) :: ms
    type (Ice_Shape), intent(inout)  :: is
    type (Thermo_Var), intent(in)  :: th_var
    real(PS) :: svi

    ! initialize
    IS%growth_mode=0
    IS%init_growth=0
    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return

    if(level==3.or.level==5.or.level==7) then
!!c       if(max(IS%r,IS%e,MS%a_len,MS%c_len)<10.0e-4_PS.and.th_var%T<253.16) then
       if(max(IS%r,IS%e,sqrt(MS%a_len**2.0+MS%c_len**2.0))<10.0e-4_PS.and.th_var%T<253.16_PS) then

          svi=th_var%s_v(2)
          if(mes_rc==2.or.mes_rc==4) then
             ! if liquid co exist
             svi=th_var%e_sat(1)/th_var%e_sat(2)-1.0_PS
          end if
!tmp          IS%growth_mode=get_growth_mode(th_var%T,svi)
          IS%growth_mode=get_growth_mode_max(th_var%T,svi)
          IS%init_growth=1

!!c          IS%growth_mode=1

       else
          IS%init_growth=0
          if(IS%habit==1.or.IS%habit==2) then
             IS%growth_mode=2 
          else if(IS%habit==3) then
             IS%growth_mode=3
          else if(IS%habit==4) then
             IS%growth_mode=4
          else if(IS%habit==5) then
             IS%growth_mode=5
          else if(IS%habit>=6) then
             IS%growth_mode=1
          end if
       end if
    else
       if(max(IS%r,IS%e,MS%a_len,MS%c_len)<10.0e-4_PS.and.th_var%T<253.16_PS) then
          IS%init_growth=1
       else
          IS%init_growth=0
       end if
       if(IS%habit==1.or.IS%habit==2) then
          IS%growth_mode=2 
       else if(IS%habit==3) then
          IS%growth_mode=3
       else
          write(*,*) "habit is not right:level,habit",level,IS%habit
          stop
       end if
    end if
  end subroutine cal_growth_mode

  subroutine cal_coef_vapdep(phase,level,ms,th_var,is)
    integer, intent(in)  :: phase,level
    type (Mass_Bin), intent(inout) :: ms
    type (Thermo_Var), intent(in) :: th_var
    type (Ice_Shape), optional       :: is
    real(PS), parameter              :: R_u = 8.3143e+07
    real(PS), parameter             :: M_w=18.0
    real(PS), parameter             :: MR=M_w/R_u
    real(PS), parameter         :: den_w = 1.0
    ! temperature at triple point 
    real (PS), parameter          :: T_0 = 273.16_PS

    ! initialize
    ms%coef=0.0_PS
    if( ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return
  

    if(phase==1) then
       ms%coef(1)=4.0_PS*PI*th_var%GTP(phase)*ms%CAP*ms%fv*ms%fkn
       ms%coef(2)=0.0_PS
    elseif(phase==2) then
       if(th_var%T<T_0.and.ms%inmlt==0) then
          ms%coef(1)=  4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
               th_var%e_sat(phase)/th_var%T
          ms%coef(2)= -4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
               ms%e_sat/ms%tmp
       else
          if(level>=6) then
             if(ms%mass(imw)/ms%con>1.0e-15_PS) then
                ms%coef(1)=  4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
                     th_var%e_sat(1)/th_var%T
                ms%coef(2)= -4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
                     ms%e_sat/ms%tmp
             else
                ms%coef=0.0_PS
             end if
          else
             ms%coef=0.0_PS
          end if
       end if
    end if
  end subroutine cal_coef_vapdep

  subroutine cal_coef_vapdep2(phase,level,ms,th_var,nu_aps,phi_aps,m_aps,is)
    integer, intent(in)  :: phase,level
    type (Mass_Bin), intent(inout) :: ms
    type (Thermo_Var), intent(in) :: th_var
    type (Ice_Shape), optional       :: is
    real(PS), parameter              :: R_u = 8.3143e+07
    real(PS), parameter             :: M_w=18.0
    real(PS), parameter             :: MR=M_w/R_u
    real(PS), parameter         :: den_w = 1.0
    real(PS),dimension(*),intent(in) :: phi_aps,nu_aps,m_aps

    ! latent heat of condensation in ergs/g
    real(PS), parameter              :: L_e = 2.5e+10_PS

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

    ! temperature at triple point 
    real (PS), parameter          :: T_0 = 273.16_PS

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
          rd_c=get_critrad_itr(AA*1.0e+4,BB,n_aps,(r_n3**(1.0/3.0))*1.0e+4)
          ms%r_crt=1.0e-4_PS*rd_c

          if(rd_c<(r_n3**(1.0/3.0))*1.0e+4) then
            write(*,*) "rd_c is smaller than dry radius",rd_c,(r_n3**(1.0/3.0))*1.0e+4
          endif

          ! calcluate haze size with SRW.   This is used for evaporation only.
          ms%r_act=1.0e-4_PS*get_hazerad_itr(AA*1.0e+4_PS,BB,SRW,n_aps &
                  ,(r_n3**(1.0/3.0))*1.0e+4,rd_c,0)

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

  subroutine cal_coef_vapdep_ap(level,ica,ms,th_var)
    integer, intent(in)  :: level,ica
    type (Mass_Bin), intent(inout) :: ms
    type (Thermo_Var), intent(in) :: th_var
    integer :: i,n
    ! coefficient for capacitance by Chiruta and Wang (2003)
!!c    real(PS),parameter   :: CAP_ros=0.434*N_lob**0.257
    real(PS),parameter   :: CAP_ros=0.619753727
    real(PS) :: fkn

    ! initialize
    ms%coef=0.0_PS
    if(ms%con<=1.0e-30_PS.or.ms%mass(1)<=1.0e-30_PS) return

    if(ica==2) then
       ! +++ calculate kinetic effect +++
       fkn=get_fkn(th_var,2,ms%a_len)
       if(level<=2.or.level==4.or.level==6) then
          ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*ms%a_len*fkn
       elseif(level==3.or.level==5.or.level==7) then
          if(th_var%T>-20.0+273.16) then
             ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*ms%a_len*fkn
          else
             if(th_var%nuc_gmode==2.or.th_var%nuc_gmode==3) then
                ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*ms%a_len*fkn
             elseif(th_var%nuc_gmode==4) then
                ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*ms%a_len*fkn
!!c                ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*CAP_ros*ms%a_len*fkn
!!c             if(th_var%nuc_gmode==4) then
!!c                ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*CAP_ros*ms%a_len*fkn
             else
                ms%coef(1)=4.0_PS*PI*th_var%gtp(2)*ms%a_len*fkn
!!c                ms%coef(1)=get_dmdt0(th_var%nuc_gmode,th_var%T)
             end if
          end if
       end if 
    end if
  end subroutine cal_coef_vapdep_ap


  subroutine cal_coef_Ts(ms,th_var,phase,TS_A1,TS_B11,TS_B12,phase2)
    type (Mass_Bin), intent(inout) :: ms
    type (thermo_var), intent(in) :: th_var
    real(PS),intent(inout) :: TS_A1,TS_B11,TS_B12
    integer,intent(in) :: phase
    integer,intent(inout) :: phase2
    real(PS), parameter              :: R_u = 8.3143e+07
    real(PS), parameter             :: M_w=18.0
    real(PS), parameter             :: MR=M_w/R_u
    ! heat content of water
    real(PS),parameter :: c_w=4.187e+5
    ! latent heat of condensation in ergs/g
    real(PS), parameter              :: L_e = 2.5e+10_PS
    ! latent heat of sublimation in ergs/g
    real(PS), parameter             :: L_s = 2.8337e+10_PS
    ! latent heat of freeze in ergs/g
    real(PS), parameter             :: L_f = 0.3337e+10_PS
    ! temperature at triple point 
    real(PS), parameter          :: T_0 = 273.16_PS
    
    real(PS) :: r_e

    ! modified conductivity
    real(PS) :: mk_a


    TS_A1=0.0_PS
    TS_B11=0.0_PS
    TS_B12=0.0_PS
    phase2=phase
    if(ms%con<1.0e-30_PS.or.ms%mass(1)<1.0e-30_PS) return

    ! calculate modified thermal conductivity
    mk_a = get_mod_thermal_cond(sqrt(ms%semi_a**2+ms%semi_c**2),&
                                th_var%K_a,th_var%T,th_var%den)

    if(T_0>th_var%T) then
!!c       r_e=th_var%e_sat_n(1)/th_var%e_sat_n(2)
       if(ms%inmlt==0) then
          phase2=phase
               
          TS_A1=L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          TS_B11=L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
               th_var%e_sat_n(phase2)/th_var%T_n/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          TS_B12=(L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
               th_var%e_sat_n(phase2)/th_var%T_n+&
               (L_f+C_w*th_var%T_n)*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con+&
               4.0_PS*PI*ms%CAP*mk_a*ms%fh*th_var%T_n)/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          
       else
          !
          ! This formulation ignores conduction between the surface of the hydrometeor
          ! and the surface of the ice core.
          !
          phase2=1
          
          TS_A1=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          TS_B11=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
               th_var%e_sat_n(phase2)/th_var%T_n/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          TS_B12=(L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
               th_var%e_sat_n(phase2)/th_var%T_n+&
               (L_f+C_w*th_var%T_n)*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con+&
               4.0_PS*PI*ms%CAP*mk_a*ms%fh*th_var%T_n)/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
       end if
       
    else
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
!!c
!!c             A1=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
!!c                  (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!c
!!c             B1=(L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!c                  th_var%e_sat_n(phase2)/th_var%T_n*&
!!c                  (th_var%s_v_n(phase2)+1.0_PS)+&
!!c                  (C_w*th_var%T_n)*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con+&
!!c                  4.0_PS*PI*ms%CAP*mk_a*ms%fh*th_var%T_n)/&
!!c                  (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
    end if
  end subroutine cal_coef_Ts

  subroutine cal_coef_Ts2(ms,th_var,phase,TS_A1,TS_B11_0,TS_B12_0,TS_B12_1,TS_B12_2,phase2)
    type (Mass_Bin), intent(inout) :: ms
    type (thermo_var), intent(in) :: th_var
    real(PS),intent(inout) :: TS_A1,TS_B11_0,TS_B12_0,TS_B12_1,TS_B12_2
    integer,intent(in) :: phase
    integer,intent(inout) :: phase2
    real(PS), parameter              :: R_u = 8.3143e+07
    real(PS), parameter             :: M_w=18.0
    real(PS), parameter             :: MR=M_w/R_u
    ! heat content of water
    real(PS),parameter :: c_w=4.187e+5
    ! latent heat of condensation in ergs/g
    real(PS), parameter              :: L_e = 2.5e+10_PS
    ! latent heat of sublimation in ergs/g
    real(PS), parameter             :: L_s = 2.8337e+10_PS
    ! latent heat of freeze in ergs/g
    real(PS), parameter             :: L_f = 0.3337e+10_PS
    ! temperature at triple point 
    real(PS), parameter          :: T_0 = 273.16_PS
    
    real(PS) :: r_e

    ! modified conductivity
    real(PS) :: mk_a

    TS_A1=0.0_PS
    TS_B11_0=0.0_PS
    TS_B12_0=0.0_PS
    TS_B12_1=0.0_PS
    TS_B12_2=0.0_PS
    phase2=phase
    if(ms%con<1.0e-30_PS.or.ms%mass(1)<1.0e-30_PS) return

    ! calculate modified thermal conductivity
    mk_a = get_mod_thermal_cond(sqrt(ms%semi_a**2+ms%semi_c**2),&
                                th_var%K_a,th_var%T,th_var%den)

    if(T_0>th_var%T) then
       !!c       r_e=th_var%e_sat(1)/th_var%e_sat(2)
       if(ms%inmlt==0) then
          phase2=phase
          
          TS_A1=L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          !!c          TS_B11=L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
          !!c               th_var%e_sat(phase2)/Ta/&
          !!c               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          ! TS_B11=TS_B11_0*e_sat(Ta)/Ta
          TS_B11_0=L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          !!c          TS_B12=(L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
          !!c               th_var%e_sat(phase2)/Ta+&
          !!c               (L_f+C_w*Ta)*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con+&
          !!c               4.0_PS*PI*ms%CAP*mk_a*ms%fh*Ta)/&
          !!c               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          
          ! TS_B12=TS_B12_0+TS_B12_1*e_sat(Ta)/Ta+TS_B12_2*Ta
          !
          TS_B12_0=(L_f*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con)/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          !TS_B12_1*e_sat(Ta)/Ta
          TS_B12_1=L_s*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          !TS_B12_2*Ta
          TS_B12_2=(C_w*ms%Ldmassdt(2)+4.0_PS*PI*ms%CAP*mk_a*ms%fh)/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
       
       
       
       else
          !
          ! This formulation ignores conduction between the surface of the hydrometeor
          ! and the surface of the ice core.
          !
          phase2=1
          
          TS_A1=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          
          !!c          TS_B11=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
          !!c               th_var%e_sat(phase2)/Ta/&
          !!c               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          ! TS_B11=TS_B11_0*e_sat(Ta)/Ta
          TS_B11_0=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          !!c          TS_B12=(L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
          !!c               th_var%e_sat(phase2)/Ta+&
          !!c               (L_f+C_w*Ta)*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con+&
          !!c               4.0_PS*PI*ms%CAP*mk_a*ms%fh*Ta)/&
          !!c               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
          ! TS_B12=TS_B12_0+TS_B12_1*e_sat(Ta)/Ta+TS_B12_2*Ta
          TS_B12_0=(L_f*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con)/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          ! 
          TS_B12_1=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          !
          TS_B12_2=(C_w*ms%Ldmassdt(2)+4.0_PS*PI*ms%CAP*mk_a*ms%fh)/&
               (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
          
       end if
    else
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
!!c
!!c             A1=L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn/&
!!c                  (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
!!c
!!c             B1=(L_e*4.0_PS*PI*th_var%D_v*ms%CAP*MR*ms%fv*ms%fkn*&
!!c                  th_var%e_sat(phase2)/th_var%T*&
!!c                  (th_var%s_v_n(phase2)+1.0_PS)+&
!!c                  (C_w*th_var%T)*ms%Ldmassdt(2)-L_f*ms%dmassdt(1,10)/ms%con+&
!!c                  4.0_PS*PI*ms%CAP*mk_a*ms%fh*th_var%T)/&
!!c                  (4.0_PS*PI*ms%CAP*mk_a*ms%fh+C_w*ms%Ldmassdt(2))
    end if

  end subroutine cal_coef_Ts2

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
    real(PS), parameter              :: R_u = 8.3143e+07
    real(PS), parameter             :: M_w=18.0
    real(PS), parameter             :: MR=M_w/R_u
    ! heat content of water
    real(PS),parameter :: c_w=4.187e+5
    ! latent heat of condensation in ergs/g
    real(PS), parameter              :: L_e = 2.5e+10_PS
    ! latent heat of sublimation in ergs/g
    real(PS), parameter             :: L_s = 2.8337e+10_PS
    ! latent heat of freeze in ergs/g
    real(PS), parameter             :: L_f = 0.3337e+10_PS
    ! temperature at triple point 
    real(PS), parameter          :: T_0 = 273.16_PS
    
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

  subroutine update_growth_mode(level,mes_rc,th_var)
!tmp    use mod_amps_utility, only: get_growth_mode
    use mod_amps_utility, only: get_growth_mode_max
    integer,intent(in) :: level,mes_rc
    type (Thermo_Var), intent(inout)  :: th_var

    if(level==3.or.level==5.or.level==7) then
       if((mes_rc==2.or.mes_rc==4).and.th_var%T<253.16_PS) then
          ! if liquid co exist
!tmp          th_var%nuc_gmode=get_growth_mode(th_var%T,th_var%e_sat(1)/th_var%e_sat(2)-1.0_PS)
          th_var%nuc_gmode=get_growth_mode_max(th_var%T,th_var%e_sat(1)/th_var%e_sat(2)-1.0_PS)

!!c          th_var%nuc_gmode=1

       end if
    end if
  end subroutine update_growth_mode
  function get_dmdt0(imode,tmp) result(out)
    integer,intent(in) :: imode
    real(PS),intent(in) :: tmp
    real(PS) :: tmpc,out
    real(PS),dimension(6) :: tmpr,irr_dvdt,pla_dvdt,col_dvdt,ros_dvdt,s3_dvdt
    integer :: i1,i2
    tmpr = (/-20.0, -30.0, -40.0, -50.0, -60.0, -70.0/)
    irr_dvdt = (/13000.0, 4600.0, 0.0, 0.0, 0.0, 0.0/)
    pla_dvdt = (/19000.0, 2600.0, 1100.0, 430.0, 340.0, 22.0/)
    col_dvdt = (/0.0, 2300.0, 2200.0, 1000.0, 1400.0, 75.0/)
    ros_dvdt = (/0.0, 0.0, 0.0, 2800.0, 2600.0, 85.0/)
    s3_dvdt = (/22000.0, 5300.0, 850.0, 0.0, 0.0, 0.0/)
    tmpc=tmp-273.16
    if(imode==1) then
       ! case of irregular polycrystals         
       out=irr_dvdt(1)+(tmpc-tmpr(1))*(irr_dvdt(2)-irr_dvdt(1))/(tmpr(2)-tmpr(1))
    elseif(imode==5) then
       ! case of planar polycrystals         
       i1=max(min(int((-20.0_RP-tmpc)/10.0_RP)+1,2),1)
       i2=i1+1
       out=s3_dvdt(i1)+(tmpc-tmpr(i1))*(s3_dvdt(i2)-s3_dvdt(i1))/(tmpr(i2)-tmpr(i1))
    else if(imode==2) then
       ! case of plates
       i1=max(min(int((-20.0_RP-tmpc)/10.0_RP)+1,5),1)
       i2=i1+1
       out=pla_dvdt(i1)+(tmpc-tmpr(i1))*(pla_dvdt(i2)-pla_dvdt(i1))/(tmpr(i2)-tmpr(i1))
    else if(imode==3) then
       ! case of columns
       i1=max(min(int((-20.0_RP-tmpc)/10.0_RP)+1,5),2)
       i2=i1+1
       out=col_dvdt(i1)+(tmpc-tmpr(i1))*(col_dvdt(i2)-col_dvdt(i1))/(tmpr(i2)-tmpr(i1))
    else if(imode==4) then
       ! case of rosetta bullets
       i1=max(min(int((-20.0_RP-tmpc)/10.0_RP)+1,5),4)
       i2=i1+1
       out=ros_dvdt(i1)+(tmpc-tmpr(i1))*(ros_dvdt(i2)-ros_dvdt(i1))/(tmpr(i2)-tmpr(i1))
    end if
    out=min(max(out,22.0_RP),22000.0_RP)*den_i*1.0e-12_RP
  end function get_dmdt0
  function get_osm(ica,molality) result(out)
    use com_amps, only: APSNAME
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
       write(*,*) "get_osm: No such chemicals",trim(APSNAME(ica))
       stop
    end select
  end function get_osm

  function osm_ammsul(molality) result(out)
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
    
    write(*,*) "ERROR! osm_ammsul: Molality is out of bound",molality

      

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
        write(*,*) "iter at Ts large 1",fa,fb,a,b
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
        write(*,*) "iter at Ts large 1",fa,fb,a,b
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
  write(*,'("brent exceeding maximum iterations. n,iter",2I5,12ES15.6)') &
       iter,a,b,fa,fb
!!$             write(*,'("W,a_c,N_act,M_act,T,r_e",10ES15.6)') ag%TV(n)%W, a_c, N_act,M_act, ag%TV(n)%T,r_e

  out=r_n*1.1_PS


1111 continue

contains
  function func(x) result(out)
    real(PS),intent(in) :: x
    real(PS) :: out
    ! g cm^{-3}
    real(PS),parameter :: den_w=1.0_PS
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
        write(*,'("iter at get_hazerad large 1",I5,10ES15.6)') iswitch,fa,fb,a,b,rd_c,r_n
     end if
     a=a*2.0_PS
!!$       stop
     goto 2221
  elseif(fa.lt.0.0_PS.and.fb.lt.0.0_PS) then
!!$       if((fa.gt.0.0_PS.and.fb.gt.0.0_PS).or.(fa.lt.0.0_PS.and.fb.lt.0.0_PS)) then
     iter=iter+1
     if(iter>50) then
        write(*,'("iter at get_hazerad large 2",I5,10ES15.6)') iswitch,fa,fb,a,b,rd_c,r_n
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
  write(*,'("brent exceeding maximum iterations. n,iter",2I5,12ES15.6)') &
       iter,a,b,fa,fb
!!$             write(*,'("W,a_c,N_act,M_act,T,r_e",10ES15.6)') ag%TV(n)%W, a_c, N_act,M_act, ag%TV(n)%T,r_e
  
  out=r_n*1.1_PS
  
  
1111 continue
!  write(*,*) "get_hazerad_itr> iter, rd",iter,out
  
contains      
  function func(x) result(out)
    real(PS),intent(in) :: x
    real(PS) :: out
    real(PS),parameter :: den_w=1.0_PS
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
    real(PS), parameter         :: gg = 980.0

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

  function interp_data1d_lut_log10_big(a,xi) result(yi)
    type(data1d_lut_big),intent(in) :: a
    real(PS) :: xi,yi,wx,x1,xi_log
    integer :: i1

    xi_log=log10(max(1.0e-30_RP,xi))
    
    i1=max(1,min(a%n-1 &
        ,int((xi_log-a%xs)/a%dx)+1))
    x1=real(i1-1,PS_KIND)*a%dx+a%xs
    wx=min(1.0_RP,max(0.0_RP,(xi_log-x1)/a%dx))
    yi=(1.0_RP-wx)*a%y(i1)+ &
             wx*a%y(i1+1) 
  
  end function interp_data1d_lut_log10_big

END MODULE class_Mass_Bin
