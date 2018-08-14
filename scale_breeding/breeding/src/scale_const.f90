!-------------------------------------------------------------------------------
!> module CONSTANT
!!
!! @par Description
!!          Physical constants module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro)  [new]
!! @li      2013-08-31 (T.Yamaura)  [add] Stefan-Boltzman constant
!!
!<
module scale_const
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use common
!  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  public 
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(r_size), public            :: CONST_PI      = 3.14159265358979 !< pi
  real(r_size), public            :: CONST_D2R                           !< degree to radian
  real(r_size), public            :: CONST_EPS     = 1.E-16           !< small number
  real(r_size), public            :: CONST_EPS1    = 0.99999999999999 !< small number
  real(r_size), public            :: CONST_HUGE    = 1.E+30           !< huge  number

  integer,  public, parameter :: CONST_UNDEF2  = -32768              !< undefined value (INT2)
  real(r_sngl), public, parameter :: CONST_UNDEF4  = -9.9999E30          !< undefined value (REAL4)
  real(r_sngl), public, parameter :: CONST_UNDEF8  = -9.9999D30          !< undefined value (REAL8)
  real(r_size), public            :: CONST_UNDEF

  ! adopted constants
  real(r_size), public            :: CONST_RADIUS  = 6.37122E+6       !< radius of the planet [m]
  real(r_size), public            :: CONST_OHM     = 7.2920E-5        !< angular velocity of the planet [1/s]
  real(r_size), public            :: CONST_GRAV    = 9.80665          !< standard acceleration of gravity [m/s2]

  ! physical constants
  real(r_size), public, parameter :: CONST_STB     = 5.670373E-8      !< Stefan-Boltzman constant [W/m2/K4]
  real(r_size), public, parameter :: CONST_KARMAN  = 0.4              !< von Karman constant
  real(r_size), public, parameter :: CONST_R       = 8.3144621        !< universal gas constant [J/mol/K]

  ! dry air constants
  real(r_size), public            :: CONST_Mdry    =   28.97          !< mass weight (dry air)                     [g/mol]
  real(r_size), public            :: CONST_Rdry    =  287.04          !< specific gas constant (dry air)           [J/kg/K]
  real(r_size), public            :: CONST_CPdry   = 1004.64          !< specific heat (dry air,constant pressure) [J/kg/K]
  real(r_size), public            :: CONST_CVdry                      !< specific heat (dry air,constant volume)   [J/kg/K]
  real(r_size), public            :: CONST_LAPS    = 6.5E-3           !< lapse rate of ISA                         [K/m]
  real(r_size), public            :: CONST_LAPSdry                    !< dry adiabatic lapse rate                  [K/m]

  ! water constants
  real(r_size), public            :: CONST_Mvap    =  18.02           !< mass weight (water vapor)                      [g/mol]
  real(r_size), public, parameter :: CONST_Rvap    = 461.46           !< specific gas constant (water vapor)            [J/kg/K]
  real(r_size), public, parameter :: CONST_CPvap   = 1845.60          !< specific heat (water vapor, constant pressure) [J/kg/K]
  real(r_size), public            :: CONST_CVvap                      !< specific heat (water vapor, constant volume)   [J/kg/K]
  real(r_size), public, parameter :: CONST_CL      = 4218.0           !< specific heat (liquid water)                   [J/kg/K]
  real(r_size), public, parameter :: CONST_CI      = 2006.0           !< specific heat (ice)                            [J/kg/K]

  real(r_size), public            :: CONST_EPSvap                        !< Rdry / Rvap
  real(r_size), public            :: CONST_EPSTvap                       !< 1 / epsilon - 1

  real(r_size), public, parameter :: CONST_EMELT   = 3.4E5
  real(r_size), public, parameter :: CONST_TMELT   = 273.15

  real(r_size), public            :: CONST_LHV                           !< latent heat of vaporizaion for use
  real(r_size), public            :: CONST_LHS                           !< latent heat of sublimation for use
  real(r_size), public            :: CONST_LHF                           !< latent heat of fusion      for use
  real(r_size), public, parameter :: CONST_LHV0    = 2.5008E6         !< latent heat of vaporizaion at 0C [J/kg]
  real(r_size), public            :: CONST_LHV00                         !< latent heat of vaporizaion at 0K [J/kg]
  real(r_size), public, parameter :: CONST_LHS0    = 2.8342E6         !< latent heat of sublimation at 0C [J/kg]
  real(r_size), public            :: CONST_LHS00                         !< latent heat of sublimation at 0K [J/kg]
  real(r_size), public            :: CONST_LHF0                          !< latent heat of fusion      at 0C [J/kg]
  real(r_size), public            :: CONST_LHF00                         !< latent heat of fusion      at 0K [J/kg]
  real(r_size), public, parameter :: CONST_PSAT0   =  610.7           !< saturate pressure of water vapor at 0C [Pa]
  real(r_size), public, parameter :: CONST_DWATR   = 1000.0           !< density of water [kg/m3]
  real(r_size), public, parameter :: CONST_DICE    =  916.8          !< density of ice   [kg/m3]

  real(r_size), public            :: CONST_SOUND                         !< speed of sound (dry air at 0C) [m/s]

  real(r_size), public            :: CONST_Pstd    = 101325.0         !< standard pressure [Pa]
  real(r_size), public            :: CONST_PRE00   = 100000.0         !< pressure reference [Pa]
  real(r_size), public            :: CONST_Tstd    = 288.15           !< standard temperature (15C) [K]
  real(r_size), public, parameter :: CONST_TEM00   = 273.15           !< temperature reference (0C) [K]
  real(r_size), public, parameter :: CONST_PPM     = 1.E-6            !< parts par million

  integer,  public            :: CONST_I_LW    = 1                   !< long-wave radiation index
  integer,  public            :: CONST_I_SW    = 2                   !< short-wave radiation index

contains 

  SUBROUTINE const_init 


    CONST_UNDEF = real(CONST_UNDEF8,kind=r_size)
    CONST_PI      = 4.D0 * atan( 1.0D0 )
    CONST_D2R     = CONST_PI / 180.0D0
    CONST_EPS     =          epsilon(0.0D0)
    CONST_EPS1    = 1.0D0 - epsilon(0.0D0)
    CONST_HUGE    =             huge(0.0D0)

    !CONST_RADIUS  = CONST_RADIUS / CONST_SmallPlanetFactor
    !CONST_OHM     = CONST_OHM    * CONST_SmallPlanetFactor

    CONST_CVdry   = CONST_CPdry - CONST_Rdry
    CONST_LAPSdry = CONST_GRAV / CONST_CPdry

    CONST_CVvap   = CONST_CPvap - CONST_Rvap
    CONST_EPSvap  = CONST_Rdry / CONST_Rvap
    CONST_EPSTvap = 1.0D0 / CONST_EPSvap - 1.0D0

    CONST_LHF0    = CONST_LHS0 - CONST_LHV0

    CONST_LHV00   = CONST_LHV0 - ( CONST_CPvap - CONST_CL ) * CONST_TEM00
    CONST_LHS00   = CONST_LHS0 - ( CONST_CPvap - CONST_CI ) * CONST_TEM00
    CONST_LHF00   = CONST_LHF0 - ( CONST_CL    - CONST_CI ) * CONST_TEM00

    !if ( CONST_THERMODYN_TYPE == 'EXACT' ) then
       CONST_LHV = CONST_LHV00
       CONST_LHS = CONST_LHS00
       CONST_LHF = CONST_LHF00
    !elseif( CONST_THERMODYN_TYPE == 'SIMPLE' ) then
    !   CONST_LHV = CONST_LHV0
    !   CONST_LHS = CONST_LHS0
    !   CONST_LHF = CONST_LHF0
    !else
    !   write(*,*) 'xxx Not appropriate ATMOS_THERMODYN_ENERGY_TYPE. Check!', trim(CONST_THERMODYN_TYPE)
    !   call PRC_MPIstop
    !endif

    CONST_SOUND = sqrt( CONST_CPdry * CONST_Rdry / ( CONST_CPdry - CONST_Rdry ) * CONST_TEM00 )



  END SUBROUTINE const_init

end module scale_const
