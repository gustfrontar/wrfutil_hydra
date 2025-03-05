module mod_amps_const
  use acc_amps
  implicit none

  public

  ! reference pressure (g/cm/s^2)
  real(PS), parameter :: p00 = 1.0e6_PS

  ! temperature at triple point
  real(PS), parameter :: T_0 = 273.16_PS

  ! gravity (cm/s^2)
  real(PS), parameter :: gg = 980.0_PS

  ! gas constant for vapor (ergs/deg/g)
  real(PS), parameter :: M_w = 18.016_PS

  real(PS), parameter :: R_u = 8.31436e7_PS
  real(PS), parameter :: MR = M_w / R_u

  ! density of water (g/cm^3)
  real(PS), parameter :: den_w = 1.0_PS
  ! density of ice (bulk density) (g/cm^3)
  real(PS), parameter :: den_i = 0.91668_PS

  ! gas constant for dry air (ergs/deg/g)
  real(PS), parameter :: R_d = 287.04e4_PS
  ! gas constant for vapor (ergs/deg/g)
  real(PS), parameter :: R_v = 461.5e4_PS

  ! heat capacity of dry air
  real(PS), parameter :: C_pa = 1004.64e4_PS

  real(PS), parameter :: Racp = R_d / C_pa
  real(PS), parameter :: Rdvchiarui = R_d / R_v
  real(PS), parameter :: M_a = M_w / Rdvchiarui



  ! heat content of water
  real(PS), parameter :: c_w = 4.187e5_PS

  ! latent heat of condensation in ergs/g
  real(PS), parameter :: L_e = 2.5e10_PS
  ! latent heat of freeze in ergs/g
  real(PS), parameter :: L_f = 0.3337e10_PS
  ! latent heat of sublimation in ergs/g
  real(PS), parameter :: L_s = 2.8337e10_PS

  ! condensation coefficient  for liquid
  real(PS), parameter :: a_cliq = 0.036_PS


  ! thermal conductivity
  real(PS), parameter :: k_w = 0.58e5_PS


  real(PS), parameter :: undef = 999.9e30_PS
  integer,  parameter :: iundef = 999

end module mod_amps_const
