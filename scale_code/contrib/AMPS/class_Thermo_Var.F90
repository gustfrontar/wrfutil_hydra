!OCL SERIAL
MODULE class_Thermo_Var
  use acc_amps
  use mod_amps_const
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! The object Thermo_Var includes the thermodynamical variables
!
! All units are CGS.
! 
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  implicit none  

  private

  public :: make_Thermo_Var3
  public :: make_Thermo_Var3_2
  public :: get_diffusivity
  public :: get_mod_diffusivity
  public :: get_mod_thermal_cond
  public :: get_thermal_conductivity
  public :: get_dynamic_viscosity
  public :: get_sfc_tension
  public :: get_GTP
  public :: get_sat_vapor_pres
  public :: get_sat_vapor_pres_lk
  public :: get_T_fesv_lk
  public :: get_fkn
  public :: renew_rv_var

  public :: Thermo_Var

  TYPE Thermo_Var
     ! memory is aligned as declared below.
     sequence

!     private
     ! (moist) ambient temperature in K
     real (PS)            :: T,T_n,T_m
     ! (moist) ambient pressure in g /s^2/cm = 0.1 Pa = 0.001 mb
     real (PS)            :: P

     ! (moist) ambient density in g /cm^3
     real (PS)            :: den

     ! density of ambient dry air in g/cm^3
     real(PS)             :: den_a


     ! diffusivity of water vapor in air
     real (PS)            :: D_v

     ! thermal conductivity of dry air
     real (PS)            :: k_a

     ! dynamic viscosity of air
     real (PS)            :: d_vis

     ! vapor pressure at infinity
     real (PS)            :: e
     ! saturation vapor pressure over water and ice at infinity
     ! argument 1 : liquid
     !          2 : ice
     real (PS), dimension(2)            :: e_sat
     real (PS), dimension(2)            :: e_sat_n

     ! super saturation of moist air with respect to a plane
     ! water surface, with respect to a plane ice surface
     ! argument 1 : liquid
     !          2 : ice
     real (PS), dimension(2)            :: s_v,s_v_n

     ! supersaturation over liquid at the beginning of the time step.
     real (PS)       :: svw0

     ! diffusivity function for
     ! argument 1 : liquid
     !          2 : ice
     real (PS), dimension(2)            :: GTP

     ! modified diffusivity function for
     ! argument 1 : liquid
     !          2 : ice
     real (PS), dimension(2)            :: MGTP

     ! mixing ratio of vapor at infinity
     real (PS)               :: rv
     ! saturation mixing ratio of vapor over water and ice at infinity
     ! argument 1 : liquid
     !          2 : ice
     real (PS), dimension(2)            :: rv_sat

     ! surface tension of water
     real (PS) :: sig_wa

     ! vertical velocity
     real (PS)                         :: W
     
     ! total mass tendency of vapor (g/s)
     real (PS) :: dmassdt
     ! tendency from evaporating hydrometeors
     ! 1: evaporation of liquid
     ! 2: evaporation of solid
     ! 3: immersion freezing
     real (PS),dimension(3) :: dmassdt_v

     ! growth mode of nucleated ice crystals
     integer :: nuc_gmode 

     integer :: align_Thermo_Var ! CHIARUI

  endtype Thermo_Var

CONTAINS

  ! ************* Optional Constructor for a Thermo_Var type *****************
  function make_Thermo_Var (temp, pres, dv, ka, dvis, &
       de, esw, esi, svw, svi, dw) &
       result (th_var)
    real (PS), optional, intent(in)    :: temp, pres, dv, ka
    real (PS), optional, intent(in)    :: dvis, de, esw, esi, svw, svi, dw
    type (Thermo_Var)  :: th_var
    integer            :: i

    if( temp /= -9.9_PS ) then; th_var%T = temp
    else ; th_var%T = 293.15_PS
    end if

    if( pres /= -9.9_PS ) then; th_var%P = pres
    else ; th_var%T = 1013250.0_PS
    end if

    if( dv   /= -9.9_PS ) then; th_var%D_v = dv
    else 
       th_var%D_v = get_diffusivity(th_var%P,th_var%T)
    end if

    if( ka  /= -9.9_PS ) then; th_var%k_a = ka
    else
       th_var%k_a = get_thermal_conductivity(th_var%T)
    end if

    if( dvis  /= -9.9_PS ) then; th_var%d_vis = dvis
    else
       th_var%d_vis = get_dynamic_viscosity( th_var%T)
    end if


    if( de /= -9.9_PS ) then; th_var%e = de
    else 
       !!!! under construction!!!!
       th_var%e = de
       !       th_var%e = get_vapor_pres(th_var%T)
    end if
    if( esw /= -9.9_PS ) then; th_var%e_sat(1) = esw
    else 
       th_var%e_sat(1) = get_sat_vapor_pres(1, th_var%T)
    end if
    if( esi /= -9.9_PS ) then; th_var%e_sat(2) = esi
    else 
       th_var%e_sat(2) = get_sat_vapor_pres(2, th_var%T)
    end if
    if( svw /= 999.9_PS ) then; th_var%s_v(1) = svw
    else 
       th_var%s_v(1) = th_var%e/th_var%e_sat(1) - 1.0_PS
    end if
    if( svi /= 999.9_PS ) then; th_var%s_v(2) = svi
    else 
       th_var%s_v(2) = th_var%e/th_var%e_sat(2) - 1.0_PS
    end if

    if( de /= -9.9_PS ) then
       th_var%den = M_a*(th_var%P - (1.0_PS - Rdvchiarui)*th_var%e )/(R_u*th_var%T)
    end if

    th_var%den_a=M_a*(th_var%P-th_var%e)/(R_u*th_var%T)

    do i=1, 2
       th_var%GTP(i) = get_GTP( th_var%T, th_var%P, th_var%D_v, th_var%k_a, &
         th_var%e_sat(i), i)
    end do

    th_var%rv = get_rv( th_var%T, th_var%e, th_var%den)
    do i=1, 2
       th_var%rv_sat(i)=Rdvchiarui*th_var%e_sat(i)/max((th_var%P-th_var%e_sat(i)),th_var%e_sat(i))
!!c       th_var%rv_sat(i) = get_rv( th_var%T, th_var%e_sat(i), th_var%den)
    end do

    if( dw  /= 999.9_PS ) then; th_var%W = dw
    else
       th_var%W = 0.0_PS
    end if

    th_var%sig_wa=get_sfc_tension(th_var%T)
    th_var%dmassdt=0.0_PS
    th_var%dmassdt_v=0.0_PS

    th_var%e_sat(2)=min(th_var%e_sat(1),th_var%e_sat(2))
    th_var%rv_sat(2)=min(th_var%rv_sat(1),th_var%rv_sat(2))

    th_var%e_sat_n=th_var%e_sat
    th_var%T_n=th_var%T
    th_var%T_m=th_var%T
  end function make_Thermo_Var

  function make_Thermo_Var2 (drv, dden, pres, temp, dw ) result (th_var)
    real(MP_KIND), intent(in)   :: drv, dden, pres, temp, dw
    type (Thermo_Var)  :: th_var
    integer            :: i
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    ! Input variables are in MKS unit.
    ! They has to be changed to CGS unit.
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    th_var%rv = drv
    th_var%den = dden*1.0e-3_PS
    th_var%P = pres*1.0e+1_PS
    th_var%T = temp
    th_var%W = dw*1.0e+2_PS

    th_var%D_v = get_diffusivity(th_var%P,th_var%T)
    th_var%k_a = get_thermal_conductivity(th_var%T)
    th_var%d_vis = get_dynamic_viscosity( th_var%T)
    th_var%sig_wa=get_sfc_tension(th_var%T)
!!c    th_var%e = (th_var%rv*th_var%den*R_v)*th_var%T
    th_var%e=th_var%P*th_var%rv/(Rdvchiarui+th_var%rv)

    th_var%den_a=M_a*(th_var%P-th_var%e)/(R_u*th_var%T)

    do i=1,2
       th_var%e_sat(i) = get_sat_vapor_pres(i, th_var%T)
       th_var%s_v(i) = th_var%e/th_var%e_sat(i) - 1.0_PS
       th_var%GTP(i) = get_GTP( th_var%T, th_var%P, th_var%D_v, th_var%k_a, &
         th_var%e_sat(i), i)
!!c       th_var%rv_sat(i) = get_rv( th_var%T, th_var%e_sat(i), th_var%den)
       th_var%rv_sat(i)=Rdvchiarui*th_var%e_sat(i)/max((th_var%P-th_var%e_sat(i)),th_var%e_sat(i))
       th_var%s_v_n(i)=th_var%s_v(i)
    end do

    th_var%dmassdt=0.0_PS
    th_var%dmassdt_v=0.0_PS

    th_var%e_sat(2)=min(th_var%e_sat(1),th_var%e_sat(2))
    th_var%rv_sat(2)=min(th_var%rv_sat(1),th_var%rv_sat(2))

    th_var%e_sat_n=th_var%e_sat
    th_var%T_n=th_var%T
    th_var%T_m=th_var%T
  end function make_Thermo_Var2
  function make_Thermo_Var3 (drv, dden, pres, temp, dw,destbar,desitbar ) result (th_var)
!tmp    use mod_amps_utility, only: get_growth_mode
    use mod_amps_utility, only: get_growth_mode_max
    real(MP_KIND), intent(in)   :: drv, dden, pres, temp, dw
    real(DS), intent(in)   :: destbar(150),desitbar(111)
    type (Thermo_Var)  :: th_var
    integer            :: i
    real(PS) :: check

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    ! Input variables are in MKS unit.
    ! They has to be changed to CGS unit.
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    th_var%rv = drv
    th_var%den = dden*1.0e-3_PS
    th_var%P = pres*1.0e+1_PS
    th_var%T = temp
    th_var%W = dw*1.0e+2_PS

    th_var%D_v = get_diffusivity(th_var%P,th_var%T)
    th_var%k_a = get_thermal_conductivity(th_var%T)
    th_var%d_vis = get_dynamic_viscosity( th_var%T)
    th_var%sig_wa=get_sfc_tension(th_var%T)
!!c    th_var%e = (th_var%rv*th_var%den*R_v)*th_var%T
    th_var%e=th_var%P*th_var%rv/(Rdvchiarui+th_var%rv)

    th_var%den_a=M_a*(th_var%P-th_var%e)/(R_u*th_var%T)

    do i=1,2
       th_var%e_sat(i) = get_sat_vapor_pres_lk(i, th_var%T, destbar, desitbar )
!!c       th_var%e_sat(i) = get_sat_vapor_pres(i, th_var%T)
       th_var%s_v(i) = th_var%e/th_var%e_sat(i) - 1.0_PS
       th_var%GTP(i) = get_GTP( th_var%T, th_var%P, th_var%D_v, th_var%k_a, &
         th_var%e_sat(i), i)
       th_var%rv_sat(i)=Rdvchiarui*th_var%e_sat(i)/max((th_var%P-th_var%e_sat(i)),th_var%e_sat(i))
!!c       th_var%rv_sat(i) = get_rv( th_var%T, th_var%e_sat(i), th_var%den)

       if(th_var%rv_sat(i)<=th_var%rv) then
          th_var%s_v(i)=max(0.0_PS,th_var%s_v(i))
       else
          th_var%s_v(i)=min(0.0_PS,th_var%s_v(i))
       end if

       th_var%s_v_n(i)=th_var%s_v(i)
    end do
    th_var%dmassdt=0.0_PS
    th_var%dmassdt_v=0.0_PS

    th_var%e_sat(2)=min(th_var%e_sat(1),th_var%e_sat(2))
    th_var%rv_sat(2)=min(th_var%rv_sat(1),th_var%rv_sat(2))
!!c    write(*,*) "rv1,rv2",th_var%rv,0.622*th_var%e/max((th_var%P-th_var%e),1.e-10)

!tmp    th_var%nuc_gmode=get_growth_mode(th_var%T,th_var%s_v(2))    
    th_var%nuc_gmode=get_growth_mode_max(th_var%T,th_var%s_v(2))    

    th_var%e_sat_n=th_var%e_sat
    th_var%T_n=th_var%T
    th_var%T_m=th_var%T
  end function make_Thermo_Var3
  function make_Thermo_Var3_2 (pres, temp, dw,destbar,desitbar,phase,RH ) result (th_var)
!tmp    use mod_amps_utility, only: get_growth_mode
    use mod_amps_utility, only: get_growth_mode_max
    real(PS), intent(in)   :: pres, temp, dw,RH
    real(DS), intent(in)   :: destbar(150),desitbar(111)
    integer, intent(in) :: phase
    type (Thermo_Var)  :: th_var
    integer            :: i
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    ! Input variables are in MKS unit.
    ! They has to be changed to CGS unit.
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    th_var%P = pres*1.0e+1_PS
    th_var%T = temp
    th_var%W = dw*1.0e+2_PS

    th_var%D_v = get_diffusivity(th_var%P,th_var%T)
    th_var%k_a = get_thermal_conductivity(th_var%T)
    th_var%d_vis = get_dynamic_viscosity( th_var%T)

    th_var%sig_wa=get_sfc_tension(th_var%T)


    do i=1,2
       th_var%e_sat(i) = get_sat_vapor_pres_lk(i, th_var%T, destbar, desitbar )
!!c       th_var%e_sat(i) = get_sat_vapor_pres(i, th_var%T)
       th_var%GTP(i) = get_GTP( th_var%T, th_var%P, th_var%D_v, th_var%k_a, &
         th_var%e_sat(i), i)
       th_var%rv_sat(i)=Rdvchiarui*th_var%e_sat(i)/max((th_var%P-th_var%e_sat(i)),th_var%e_sat(i))
!!c       th_var%rv_sat(i) = get_rv( th_var%T, th_var%e_sat(i), th_var%den)
    end do



    th_var%rv = RH*th_var%rv_sat(phase)/100.0

    th_var%e = th_var%rv*th_var%P/(Rdvchiarui+th_var%rv)

    do i=1,2
       th_var%s_v(i) = th_var%e/th_var%e_sat(i) - 1.0_PS
       th_var%s_v_n(i)=th_var%s_v(i)
    end do

    th_var%den = M_a*(th_var%P - (1.0_PS - Rdvchiarui)*th_var%e )/(R_u*th_var%T)
    th_var%den_a=M_a*(th_var%P-th_var%e)/(R_u*th_var%T)
    th_var%dmassdt=0.0_PS
    th_var%dmassdt_v=0.0_PS

    th_var%e_sat(2)=min(th_var%e_sat(1),th_var%e_sat(2))
    th_var%rv_sat(2)=min(th_var%rv_sat(1),th_var%rv_sat(2))

!tmp    th_var%nuc_gmode=get_growth_mode(th_var%T,th_var%s_v(2))    
    th_var%nuc_gmode=get_growth_mode_max(th_var%T,th_var%s_v(2))    

    th_var%e_sat_n=th_var%e_sat
    th_var%T_n=th_var%T
    th_var%T_m=th_var%T
  end function make_Thermo_Var3_2
!!c  ! ************* Optional Constructor for a Thermo_Var type *****************
!!c  function Thermo_Var_ (temp, pres, dv, ka, &
!!c       de, esw, esi, svw, svi) &
!!c       result (th_var)
!!c    real (PS), optional, intent(in)    :: temp, pres, dv, ka, de, esw, esi, svw, svi
!!c    type (Thermo_Var)  :: th_var
!!c
!!c    th_var = Thermo_Var( temp, pres, dv, ka, &
!!c         de, esw, esi, svw, svi)
!!c  end function Thermo_Var_

  function get_diffusivity(P,T) result (D_v)
    ! **********************************************************************
    ! Calculate the diffusivity of water vapor in air for temperatures
    ! between -40 to 40 C according to Hall and Pruppacher (1976)
    ! **********************************************************************
    ! ambient pressure and temperature in g/s^2/cm and Kelvin
    real (PS), intent(in)   :: P, T
    ! diffusivity of water vapor in cm^2/sec
    real (PS)               :: D_v
    real (PS), parameter:: P_0 = 1013250.0_PS ! g/s^2/cm
    D_v = 0.211_PS * ((T/T_0)**1.94_PS) * (P_0/P)
  end function get_diffusivity
  function get_mod_diffusivity(iphase,radius,OD_v,T) result (out)
    ! **********************************************************************
    ! Calculate the modified diffusivity of water vapor in air
    ! **********************************************************************
    integer,intent(in) :: iphase
    real (PS), intent(in) :: radius
    ! ambient pressure and temperature in g/s^2/cm and Kelvin
    real (PS), intent(in)   :: T
    ! original diffusivity of water vapor in cm^2/sec
    real (PS), intent(in)               :: OD_v
    ! modified diffusivity of water vapor in cm^2/sec
    real (PS)                :: out
    ! condensation coefficient
    real(PS) :: a_c
    ! condensation coefficient  for liquid
    real (PS), parameter             :: a_cliq=1.0
!test    real (PS), parameter             :: a_cliq=0.036
    ! deposition coefficient  for ice
    real (PS), parameter             :: a_cice1=0.5
    ! deposition coefficient  for ice for T<-40C
    real (PS), parameter             :: a_cice2=0.006
    ! vapor jump length
    real (PS), parameter             :: del_v=1.0e-5

    if(iphase==1) then
       ! liquid
       a_c=a_cliq
    else
       !if(iphase==2) then
       ! ice
       a_c=a_cice1
    endif
    out=radius/(radius+del_v)+sqrt(2.0*PI/(R_v*T))*OD_v/(radius*a_c)
    out=OD_v/out
  end function get_mod_diffusivity
  function get_mod_thermal_cond(radius,OK_a,T,den) result (out)
    ! **********************************************************************
    ! Calculate the modified diffusivity of water vapor in air
    ! **********************************************************************
    real (PS), intent(in) :: radius
    ! ambient pressure and temperature in g/s^2/cm and Kelvin
    real (PS), intent(in)   :: T,den
    ! original thermal diffusivity of water vapor in cm^2/sec
    real (PS), intent(in)               :: OK_a
    ! modified thermal diffusivity of water vapor in cm^2/sec
    real (PS)                :: out
    ! thermal accommodation coefficient
    real (PS), parameter             :: a_t=0.96
    ! vapor jump length
    real (PS), parameter             :: del_t=2.16e-5


    out=radius/(radius+del_t)+sqrt(2.0*PI*M_a/(R_u*T))*OK_a/(den*radius*a_t*c_pa)
!!c    out=radius/(radius+del_t)+sqrt(2.0*PI*M_a/(R_u*T))*OK_a/(radius*a_t*c_pa)
    out=OK_a/out
  end function get_mod_thermal_cond

  function get_thermal_conductivity(T) result (k_a)
    ! **********************************************************************
    ! Calculate the thermal conductivity of dry air
    ! according to Beard and Pruppacher (1971a)
    ! **********************************************************************
    ! ambient temperature in Kelvin
    real (PS), intent(in)   :: T
    ! thermal conductivity in 
    real (PS)               :: k_a
    k_a = (5.69 + 0.017*(T-273.16))*4.1868*1.0e+02 ! gcm^2/s^2 /cm/sec/deg
  end function get_thermal_conductivity

  function get_sat_vapor_pres(phase, d_T) result (e_sat)
    ! **********************************************************************
    ! Calculate the saturation vapor pressure over water and ice
    ! according to Lowe and Ficke (1974)
    ! **********************************************************************
    ! ambient temperature (K)
    real (PS), intent(in) :: d_T
    ! phase of water which the vapor pressure exist over
    !    water => 1   ice => 2 
    integer          :: phase
    ! saturation vapor pressure ( g/s^2/cm )
    real (PS)                   :: e_sat
    ! coefficient
    real (PS), dimension(7)     :: a
!    real (PS)                   :: T
!    T = d_T - 273.16

    if( phase == 1 ) then
       ! over liquid
!       a(:) = (/6.107799961, 4.436518521e-01, 1.428945805e-02, 2.650648471e-04,&
!            3.031240396e-06, 2.034080948e-08, 6.136820929e-11/)
       e_sat = 6.1070_PS*exp(17.15_PS*(d_T-273.16_PS)/(d_T-38.25_PS))
    else if( phase == 2) then
       ! over ice
!       a(:) = (/6.10690449, 5.02660639e-01, 1.87743264e-02, 4.13476180e-04,&
!            5.72333773e-06, 4.71651246e-08, 1.78086695e-10/)
       e_sat = 6.1064_PS*exp(21.88_PS*(d_T-273.16_PS)/(d_T-7.65_PS))
    end if
!    e_sat = a(1) + T*(a(2)+T*(a(3)+T*(a(4)+T*(a(5)+T*(a(6)+T*a(7)))))) ! hPa

    e_sat = 1000.0_PS*e_sat ! g/s^2/cm
    
  end function get_sat_vapor_pres

  function get_sat_vapor_pres_lk(phase, T, estbar, esitbar ) result(e_sat)
    !_______________________________________________________________________
    !     THIS ROUTINE COMPUTES SATURATION MIXING RATIO OVER ICE
    !     WATER BASED ON A TABLE LOOK UP PROCEEDURE DEVELOPED BY
    !     DERICKSON AND COTTON (1977) FOR THE LIQUID PHASE
    !_______________________________________________________________________
    implicit none
    integer,intent(in)  :: phase
    real(PS),intent(in) :: T
    real(DS),intent(in) :: estbar(150),esitbar(111)
    integer :: I,J
    real :: wt
    real(PS) :: e_sat
    ! g/s^2/cm
    if( phase==1 ) then
       I=MAX(1,MIN(int(T)-163,149))
       wt=MAX(MIN(T-real(I+163,PS_KIND),1.0_RP),0.0_RP)
       e_sat=(estbar(I)*(1.0_RP-wt)+estbar(I+1)*wt)*10.0_RP
    else
       ! phase == 2
       J=MAX(1,MIN(int(T)-163,110))
       wt=MAX(MIN(T-real(J+163,PS_KIND),1.0_RP),0.0_RP)
       e_sat=(esitbar(J)*(1.0_RP-wt)+esitbar(J+1)*wt)*10.0_RP
    end if
!tmp    write(*,*) "ck phase,est",phase,estbar(1:10),esitbar(1:10),e_sat
!!c    do i=1,150
!!c       write(*,*) "estbar2",i,estbar(i)
!!c    end do
!!c    do i=1,111
!!c       write(*,*) "esitbar2",i,esitbar(i)
!!c    end do
!!c    stop
  end function get_sat_vapor_pres_lk

  function get_dynamic_viscosity( d_T) result (mu)
    ! calculate the dynamic viscosity, g/cm/s
    ! ambient temperature (K)
    real (PS), intent(in) :: d_T
    real (PS)             :: mu,Tc

!!c    mu = 0.000171*(d_T/273.0)**0.7
    Tc=d_T-273.15_PS
    ! from Pruppacker's book 
    if(Tc>=0.0_PS) then
       mu=(1.718_PS+0.0049_PS*Tc)*1.0e-4_PS
    else
       mu=(1.718_PS+0.0049_PS*Tc-1.2e-5_PS*Tc*Tc)*1.0e-4_PS
    end if
  end function get_dynamic_viscosity
  function get_sfc_tension( d_T) result (sig_wa)
    ! calculate the surface tension of water: erg/cm^2
    real (PS), intent(in) :: d_T
    real (PS)             :: sig_wa,Tc
    ! (5-12) of Prupacker and Klett 97
    real(PS),parameter,dimension(1:7) :: &
        an=(/75.93,0.115,6.818e-2,6.511e-3,2.933e-4,6.283e-6,5.285e-8/)
!!!    ! Tc range from -20 to 40 C. 
!!!    Tc=d_T-273.15_PS
!!!!    sig_wa = min(76.01-0.155*Tc,79.2_PS)
    ! valid between 40 and -40 Celisus
    Tc=max(-45.0_RP,min(40.0_RP,d_T-273.15_RP))
    sig_wa=an(1)+an(2)*Tc+an(3)*Tc*Tc+an(4)*Tc*Tc*Tc+an(5)*Tc*Tc*Tc*Tc+&
           an(6)*Tc*Tc*Tc*Tc*Tc+an(7)*Tc*Tc*Tc*Tc*Tc*Tc
  end function get_sfc_tension

  function get_GTP( d_T, d_P, d_D_v, d_k_a, d_es, iphase) result (gtp)
    ! calculate the diffusivity function.
    ! ambient temperature (K), pressure (g/s^2/cm) 
    real (PS), intent(in)            :: d_T, d_P
    ! diffusivity of water vapor in air
    real (PS), intent(in)            :: d_D_v
    ! thermal conductivity of dry air
    real (PS), intent(in)            :: d_k_a
    ! saturation vapor pressure ove water or ice.
    real (PS), intent(in)            :: d_es
    ! phase of water, 1: liquid, 2: solid
    integer, intent(in)              :: iphase
    !
    ! diffusivity function
    real (PS)                        :: gtp
    !

    if( iphase == 1 ) then
       gtp = R_v*d_T/(d_es*d_D_v)+L_e*( L_e/(R_v*d_T) - 1.0_PS) /(d_k_a*d_T)
       gtp = 1.0_PS/gtp
    else if( iphase == 2 ) then
       gtp = R_v*d_T/(d_es*d_D_v)+L_s*( L_s/(R_v*d_T) - 1.0_PS) /(d_k_a*d_T)
       gtp = 1.0_PS/gtp
    endif
  end function get_GTP


  function get_rv( d_T, d_e, d_den) result (rv)
    ! calculate the diffusivity function.
    ! ambient temperature (K)
    real (PS), intent(in)            :: d_T
    ! vapor pressure ove water or ice.
    real (PS), intent(in)            :: d_e
    ! density of ambient air
    real (PS), intent(in)            :: d_den
    !
    ! mixing ratio of vapor
    real (PS)                        :: rv
    !

    rv = d_e/(R_v*d_T)/d_den
  end function get_rv

  subroutine renew_rv_var(th_var, rv)
    type (Thermo_Var)  :: th_var
    real(PS), intent(in)    :: rv
    integer   :: i
    !
    th_var%rv = rv
    th_var%e = (th_var%rv*th_var%den*R_v)*th_var%T

    ! for vectorization purpose.
    th_var%s_v(1) = th_var%e/th_var%e_sat(1) - 1.0_PS
    th_var%s_v_n(1)=th_var%s_v(1)
    th_var%s_v(2) = th_var%e/th_var%e_sat(2) - 1.0_PS
    th_var%s_v_n(2)=th_var%s_v(2)

  end subroutine renew_rv_var

  function get_T_fesv_lk(phase, T, des,estbar, esitbar ) result(T_n)
    !_______________________________________________________________________
    !     THIS ROUTINE COMPUTES SATURATION MIXING RATIO OVER ICE
    !     WATER BASED ON A TABLE LOOK UP PROCEEDURE DEVELOPED BY
    !     DERICKSON AND COTTON (1977) FOR THE LIQUID PHASE
    !_______________________________________________________________________
    implicit none
    integer,intent(in)  :: phase
    real(PS),intent(in) :: T
    real(DS),intent(in) :: estbar(150),esitbar(111)
    integer :: I,J,k
    real(PS) :: wt,des,es,T_n
    es=des/10.0_PS
    ! g/s^2/cm
    if( phase==1 ) then
       I=MAX(1,MIN(int(T)-163,149))
       wt=MAX(MIN(T-real(I+163,PS_KIND),1.0_RP),0.0_RP)
!!c       e_sat=(estbar(I)*(1.0-wt)+estbar(I+1)*wt)*10.0_RP

       do k=max(1,I-5),min(149,I+5)
          if(estbar(k)<=es.and.estbar(k+1)>es) then
             wt=(es-estbar(k))/(estbar(k+1)-estbar(k))
             T_n=real(k+163,PS_KIND)+wt
!!c             write(*,*) 'here',wt,k,T_n
             goto 10
          end if
       end do
       do k=1,I-4
          if(estbar(k)<=es.and.estbar(k+1)>es) then
             wt=(es-estbar(k))/(estbar(k+1)-estbar(k))
             T_n=real(k+163,PS_KIND)+wt
             goto 10
          end if
       end do
       do k=I+5,149
          if(estbar(k)<=es.and.estbar(k+1)>es) then
             wt=(es-estbar(k))/(estbar(k+1)-estbar(k))
             T_n=real(k+163,PS_KIND)+wt
             goto 10
          end if
       end do
    else if( phase==2 ) then
       J=MAX(1,MIN(int(T)-163,110))
       wt=MAX(MIN(T-real(J+163,PS_KIND),1.0_RP),0.0_RP)
!!c       e_sat=(esitbar(J)*(1.0-wt)+esitbar(J+1)*wt)*10.0
       do k=max(1,J-5),min(110,J+5)
          if(esitbar(k)<=es.and.esitbar(k+1)>es) then
             wt=(es-esitbar(k))/(esitbar(k+1)-esitbar(k))
             T_n=real(k+163,PS_KIND)+wt
             goto 10
          end if
       end do
       do k=1,J-4
          if(esitbar(k)<=es.and.esitbar(k+1)>es) then
             wt=(es-esitbar(k))/(esitbar(k+1)-esitbar(k))
             T_n=real(k+163,PS_KIND)+wt
             goto 10
          end if
       end do
       do k=J+5,110
          if(esitbar(k)<=es.and.esitbar(k+1)>es) then
             wt=(es-esitbar(k))/(esitbar(k+1)-esitbar(k))
             T_n=real(k+163,PS_KIND)+wt
             goto 10
          end if
       end do
    end if
10 continue
    
  end function get_T_fesv_lk

  subroutine update_DvKa(th_var,radius)
    type (Thermo_Var)  :: th_var
    real(PS), intent(in)   :: radius ! precision test (MP_KIND)
    integer :: i
    real (PS) :: mD_v,mk_a
    ! calculate modified diffusivity and conducitivity
    mk_a = get_mod_thermal_cond(radius,th_var%K_a,th_var%T,th_var%den)
    
    do i=1, 2
       mD_v = get_mod_diffusivity(i,radius,th_var%D_v,th_var%T)
       th_var%MGTP(i) = get_GTP( th_var%T, th_var%P, mD_v, mk_a, &
         th_var%e_sat(i), i)
    end do

  end subroutine update_DvKa

  function get_fkn(th_var,phase,r0) result(fkn)
    ! NOTE:
    ! this should be the same as get_mod_diffusivity.
    !
    type (Thermo_Var)  :: th_var
    real(PS), intent(in)   :: r0
    integer :: phase
    real(PS) :: beta
    ! deposition coefficient, thickness of a bc
    real(PS),parameter :: beta_w=0.036
    real(PS),parameter :: beta_i1=0.5,beta_i2=0.006
    real(PS),parameter :: delta=1.0e-5    
    real(PS) :: fkn

    if(phase==1) then
       beta=beta_w
    else
       beta=beta_i1
    end if
    ! +++ calculate kinetic effect +++
    fkn = r0/(r0+delta)+(th_var%D_v/beta/r0)*&
         sqrt(2.0_PS*PI/R_v/th_var%T)
    fkn = 1.0_PS/fkn
  end function get_fkn
  function get_esi_fesw(esw) result (esi)
    ! **********************************************************************
    ! Calculate the saturation vapor pressure over water and ice
    ! according to Lowe and Ficke (1974)
    ! **********************************************************************
    real (PS),intent(in) :: esw
    real (PS) :: esi
    ! ambient temperature (K)
    real (PS)  :: d_T
    ! coefficient
    real(PS),parameter :: a=6.1070,b=17.15,c=273.16,d=38.25 

    d_T=(d*log(esw/1000.0/a)/b-c)/(log(esw/1000.0/a)/b-1.0_PS)
    esi = 6.1064_PS*exp(21.88_PS*(d_T-273.16_PS)/(d_T-7.65_PS))
    
  end function get_esi_fesw
  function get_T_fesw(esw) result (T)
    ! **********************************************************************
    ! Calculate the saturation vapor pressure over water and ice
    ! according to Lowe and Ficke (1974)
    ! **********************************************************************
    real (PS),intent(in) :: esw
    ! ambient temperature (K)
    real (PS)  :: T
    ! coefficient
    real(PS),parameter :: a=6.1070,b=17.15,c=273.16,d=38.25 

    T=(d*log(esw/1000.0/a)/b-c)/(log(esw/1000.0/a)/b-1.0_PS)
  end function get_T_fesw
END MODULE CLASS_THERMO_VAR
