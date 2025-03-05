Module acc_amps
  use scale_precision, only: &
     DP, &
     SP, &
     RP_SCALE => RP
  implicit none
  public
  ! just double precision 15 significant digits, and -10^307 to 10^307
  !integer, parameter  :: DS = SELECTED_REAL_KIND(15, 307)
  integer, parameter  :: DS = DP
  ! just real(4)
  !integer, parameter  :: PS = SELECTED_REAL_KIND(6, 37)
  !integer, parameter  :: PS = SELECTED_REAL_KIND(15, 307)

  integer, parameter  :: PS = RP_SCALE
  integer, parameter  :: MP_KIND = PS
  integer, parameter  :: PS_KIND = PS

  integer, parameter  :: RP = PS

  !---------------------------------------------------------------------------------------------
  ! precalculated coefficients
  real(MP_KIND), parameter :: PI = 3.141592653589793238462643_PS
  real(MP_KIND), parameter :: sq_three = 1.7320508075688772935_PS
  ! 3 sqrt(3), 4pi/3, 2.0*pi
  real(MP_KIND),parameter :: coef3s=5.196152423_MP_KIND, coef4pi3=4.18879020478639_MP_KIND,coef2p=6.28318530717959_MP_KIND
  ! pi/6.0, sqrt(2*pi), 3*sqrt(3)
  real(MP_KIND),parameter :: coefpi6=0.523598776_MP_KIND,coefsq2p=2.506628274631_MP_KIND,coef3sq3=5.19615242270663_MP_KIND
  ! 4.0*pi,(3.0/4.0/pi)**(1.0/3.0)
  real(MP_KIND),parameter :: coef4p=1.25663706144E+01_MP_KIND,coef3i4p1i3=0.62035049089940_MP_KIND
  ! pi/6.0, sqrt(2*pi), 3*sqrt(3)
  real(DS),parameter :: coedpi6=0.523598775598299_DS,coedsq2p=2.506628274631_DS,coed3sq3=5.19615242270663_DS
  ! pi/180.0
  real(MP_KIND),parameter :: coefpi180=0.0174532925199433_MP_KIND
end Module acc_amps
