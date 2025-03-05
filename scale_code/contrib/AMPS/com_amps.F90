Module com_amps
  use acc_amps
  implicit none
  public
  integer :: IRLEVEL,ISLEVEL,IALEVEL
  common/MLVL/IRLEVEL,ISLEVEL,IALEVEL
  real(PS) :: drpdrp, bbcdrp, hexdrp, coldrp, &
       gp1drp, gp4drp, gp8drp,rvtm,bu_d!,bu_fd,bu_tmass
  real(PS), pointer :: bu_fd(:,:),bu_tmass(:)
  common/ECTBL/drpdrp(201,201), bbcdrp(64,71), hexdrp(64,71), &
       coldrp(62,71),&
       gp1drp(37,125), gp4drp(27,125), gp8drp(21,125),rvtm(4,3,15,9),&
       bu_fd, bu_tmass
!!!       bu_fd(2,17835),bu_tmass(435)
!!!       bu_fd(2,62400),bu_tmass(780)
!!!       bu_fd(2,956960),bu_tmass(7800)

  integer :: nr_drpdrp,nc_drpdrp,nr_gp1drp,nc_gp1drp,nr_gp4drp,nc_gp4drp,&
       nr_gp8drp,nc_gp8drp,nr_hexdrp,nc_hexdrp,nr_bbcdrp,nc_bbcdrp,nr_coldrp,nc_coldrp
  common/IECTBLH/nr_drpdrp,nc_drpdrp,nr_gp1drp,nc_gp1drp,nr_gp4drp,nc_gp4drp,&
       nr_gp8drp,nc_gp8drp,nr_hexdrp,nc_hexdrp,nr_bbcdrp,nc_bbcdrp,nr_coldrp,nc_coldrp
  real(PS) :: xs_drpdrp,dx_drpdrp,ys_drpdrp,dy_drpdrp,xs_gp1drp,dx_gp1drp,ys_gp1drp,dy_gp1drp,&
       xs_gp4drp,dx_gp4drp,ys_gp4drp,dy_gp4drp,xs_gp8drp,dx_gp8drp,ys_gp8drp,dy_gp8drp,&
       xs_hexdrp,dx_hexdrp,ys_hexdrp,dy_hexdrp,xs_bbcdrp,dx_bbcdrp,ys_bbcdrp,dy_bbcdrp,&
       xs_coldrp,dx_coldrp,ys_coldrp,dy_coldrp
  common/RECTBLH/xs_drpdrp,dx_drpdrp,ys_drpdrp,dy_drpdrp,xs_gp1drp,dx_gp1drp,ys_gp1drp,dy_gp1drp,&
       xs_gp4drp,dx_gp4drp,ys_gp4drp,dy_gp4drp,xs_gp8drp,dx_gp8drp,ys_gp8drp,dy_gp8drp,&
       xs_hexdrp,dx_hexdrp,ys_hexdrp,dy_hexdrp,xs_bbcdrp,dx_bbcdrp,ys_bbcdrp,dy_bbcdrp,&
       xs_coldrp,dx_coldrp,ys_coldrp,dy_coldrp


  integer :: imin_bk,imax_bk,jmin_bk,jmax_bk,npair_bk
  common/BKUPFD/imin_bk,imax_bk,jmin_bk,jmax_bk,npair_bk

  integer :: level_comp, debug_level, coll_level, out_type, T_print_period, &
       token_c, token_r, token_s, token_a,&
       dtype_c, dtype_r, dtype_s, dtype_a, hbreak_c, hbreak_r, hbreak_s, &
       flagp_c, flagp_r, flagp_s, flagp_a, act_type,&
       nrmic, nrtime,nbin_h
  common/BMINOP/level_comp, debug_level, out_type, T_print_period, &
       token_c, token_r, token_s, token_a, &
       dtype_c, dtype_r, dtype_s, dtype_a(4), hbreak_c, hbreak_r, hbreak_s, &
       flagp_c, flagp_r, flagp_s, flagp_a, act_type,&
       nrmic, nrtime(17),nbin_h

  real(PS) :: srat_c, srat_r, srat_s, srat_a, sadd_c, sadd_r, sadd_s, sadd_a, sth_r,&
       fcon_c, minmass_c, minmass_r, minmass_s, minmass_a, binbi,binbr
  
  common/BMREOP/srat_c, srat_r, srat_s, srat_a, sadd_c, sadd_r, sadd_s, sadd_a, sth_r,&
       fcon_c,minmass_c,minmass_r,minmass_s, minmass_a, binbi(101),binbr(101)
!       fcon_c,minmass_c,minmass_r,minmass_s, minmass_a, binbi(101),binbr(201)

  real(PS) :: pol_frq,pla_frq,col_frq,ros_frq,ppo_frq
  common/FRQTBL/pol_frq(51,101),pla_frq(51,101),col_frq(51,101),&
       ros_frq(51,101),ppo_frq(51,101)

  ! for retrieving geometry given temperature and mass
  real(PS) :: mtac_map_col,lmt_mass_col,mtac_map_pla,lmt_mass_pla
  common/RMTACTBL/mtac_map_col(50,101,2),lmt_mass_col(50),mtac_map_pla(50,101,2),lmt_mass_pla(50)
  integer :: i_lmt_mass_col,i_lmt_mass_pla
  common/IMTACTBL/i_lmt_mass_col(50),i_lmt_mass_pla(50)

  real(PS) :: ap_lnsig,ap_mean,apt_st,apt_add,apt_max,frac_apact,&
       M_aps,M_api,den_aps,den_api,den_apt,nu_aps,phi_aps,eps_ap,ap_mass,coef_ap, &
       N_ap_ini,ap_mean_cp,ap_sig_cp,cdf_cp_0,cdf_cp_180m0
  common/RAPTBL/ap_lnsig(4),ap_mean(4),apt_st(4),apt_add(4),apt_max(4),frac_apact(2,8320),&
       M_aps(4),M_api(4),den_aps(4),den_api(4),den_apt(4),nu_aps(4),phi_aps(4),eps_ap(4),ap_mass(4),coef_ap(4), &
       N_ap_ini(4),ap_mean_cp(4),ap_sig_cp(4),cdf_cp_0(4),cdf_cp_180m0(4)
  integer :: napt
  common/IAPTBL/napt(4)
  character (len=16) :: APSNAME
  common/CAPTBL/APSNAME(4)

  character (len=300) :: DRCETB,DRAPTB,DRSTTB
  common/BMCHOP/DRCETB,DRAPTB,DRSTTB

  character (len=16) :: output_format
  common/CIOOPT/output_format
!!c  integer :: N_mc2, N_sat2, is_1Dad
!!c  real :: s_start2, ds2, s_end2, beta2
!!c  common/GMTBL2/N_mc2, N_sat2, s_start2, ds2, s_end2, beta2(5,6,2), &
!!c       is_1Dad(2)
!!c  integer :: N_mc3, N_sat3
!!c  real :: s_start3, ds3, s_end3, beta3
!!c  common/GMTBL3/N_mc3, N_sat3, s_start3, ds3, s_end3, beta3(4,6,3)

  real(PS) :: znorm
  common/STTBL/znorm(4,451)

  real(PS) :: FRHW,FRHI,FTMP

  ! for AP bugdet analysis
  real(PS) :: total_ap0
  common/RAPBGD/total_ap0(3,2,5)
  integer :: iap_z,iap_zbnd,nbglev
  common/IAPBGD/iap_z(100),iap_zbnd(4),nbglev


  ! execution flag for microphysics
  integer :: micexfg
  common/IMICEXFG/micexfg(20)

  ! smallest and largest bin numbers for fragment distribution
  integer :: ibinr_frg0,ibinr_frg1
  common/IFRGDIS/ibinr_frg0,ibinr_frg1
  real(PS) :: r_frg0,r_frg1
  common/RFRGDIS/r_frg0,r_frg1

  ! for random number generator
  integer :: sizeseed,seed,seed_old
!sudachi  common/IRNG/sizeseed,seed(2),seed_old(2)
  common/IRNG/sizeseed,seed(256),seed_old(256)
!uzushio  common/IRNG/sizeseed,seed(8),seed_old(8)

  ! for time steps
  integer :: n_step_cl,n_step_vp
  common/AMPSTSTEP/n_step_cl,n_step_vp

!!!  integer :: seed_priv,isec_seed
!!!  common/IRNGPRIV/seed_priv(256),isec_seed

  ! for inherent growth parameterization
  integer :: nok_igp=23
  real(PS) :: x_igp,a_igp,b_igp
  common/VAPIGP/x_igp(23),a_igp(23,4),b_igp(23,4)

  ! for osmotic coef.
  integer :: n_osm_nh42so4,n_osm_sodchl
  real(PS) :: xs_osm_nh42so4,xs_osm_sodchl
  real(PS) :: dx_osm_nh42so4,dx_osm_sodchl
  real(PS),dimension(100) :: y_osm_nh42so4
  real(PS),dimension(100) :: y_osm_sodchl

  ! growth diagnosis for habit 
  !   default: random,yes
  integer :: ihabit_gm_random=1

  ! for standard normal distribution
  integer :: n_snrml
  real(PS) :: xs_snrml
  real(PS) :: dx_snrml
  real(PS),dimension(501) :: y_snrml
  ! for inverse of standard normal distribution
  integer :: n_isnrml
  real(PS) :: xs_isnrml
  real(PS) :: dx_isnrml
  real(PS),dimension(501) :: y_isnrml

  ! for file id of log
  integer :: fid_alog
  character(len=100) :: fname_ampslog
  ! file id for each process
  character(len=6)   :: string_fid_amps

  !
  ! max activated CCN number
  real(PS) :: CCNMAX=230.0_PS
  real(PS) :: frac_dust=0.01_PS
  real(PS) :: CRIC_RN_IMM=0.25e-4_PS

  logical :: IsMaster

  logical :: debug = .false.

end Module com_amps
