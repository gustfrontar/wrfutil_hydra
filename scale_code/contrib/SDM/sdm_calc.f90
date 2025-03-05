!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Set the number of super-droplets
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2021-03-06 (S.Shima) [new] separated sdm_calc from scale_atmos_phy_mp_sdm
!! @li      2022-05-21 (S.Shima) [mod] support terrain following and vertically stretched coordinate
!! @li      2022-09-17 (S.Shima) [add] subroutine sdm_fillboundaries and inflow boundary condition option
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module m_sdm_calc

  implicit none
  private
  public :: sdm_calc

contains
  subroutine sdm_calc(  &
                      KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                      MOMX,MOMY,MOMZ,DENS,RHOT,QTRC,              & 
                      sdm_calvar,sdm_rdnc,sdm_aslset,             &
                      sdm_inisdnc,sdm_sdnmlvol,                   &
                      sdm_mvexchg,   &
                      prec_crs,zph_crs,   &
                      lsdmup,ni_sdm,nj_sdm,nk_sdm,                &
                      sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,      &
                      sd_u,sd_v,sd_vz,sd_r,sd_asl,sdi,            &
                      sd_fmnum,        sd_fmn,sd_fmliqice,sd_fmx,sd_fmy,sd_fmz,    &
                      sd_fmri,sd_fmrj,sd_fmrk,sd_fmr,sd_fmasl,sd_fmi, &
                      sd_rkl,sd_rku,  &
                      sd_rng,sd_rand,sort_id,sort_key,sort_freq,  &
                      sort_tag,                                   &
                      bufsiz1,                                    & 
                      bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4,      &
                      sdm_itmp1,sdm_itmp2,        &
                      sd_itmp1,sd_itmp2,sd_itmp3,sd_dtmp1,sd_dtmp2,sd_dtmp3,sd_dtmp4,   &
                      crs_val1p,crs_val1c,crs_val2p,crs_val2c,    &
                      rbuf_r8,sbuf_r8,rbuf_i8,sbuf_i8,rbuf_i2,sbuf_i2,rbuf_i4,sbuf_i4, &
                      dt )
   use scale_precision
   use scale_comm_cartesC, only: &
        COMM_vars8, &
        COMM_wait
   use scale_const, only: &
        GRAV   => CONST_GRAV

   use rng_uniform_mt, only: &
        c_rng_uniform_mt

   use m_sdm_common, only: &
        VALID2INVALID, &
        sdm_dtcmph,sdm_cold,sdm_colkrnl,sdm_colbrwn, &
        sdm_aslmw,sdm_aslion,sdm_aslrho,sdm_wbc,sdm_ebc,sdm_sbc,sdm_nbc, &
        nclstp,sd_dtmp5,sd_dtmp6, &
        i2,sdicedef,I_QV_sdm,QA_MP_sdm,tag
   use m_sdm_fluidconv, only: &
        sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw
   use m_sdm_sd2fluid, only: &
        sdm_sd2rhow, sdm_sd2prec, sdm_sd2rhosol
   use m_sdm_boundary, only: &
        sdm_jdginvdv, sdm_boundary
   use m_sdm_motion, only: &
        sdm_getvz_liq, sdm_getvz_ice, sdm_getvel, sdm_move
   use m_sdm_coalescence, only: &
        sdm_coales
   use m_sdm_coalescence_cold, only: &
        sdm_coales_cold
   use m_sdm_condensation_water, only: &
        sdm_condevp, sdm_condevp_updatefluid
   use m_sdm_meltfreeze, only: &
        sdm_meltfreeze, sdm_meltfreeze_updatefluid
   use m_sdm_subldep, only: &
        sdm_subldep, sdm_subldep_updatefluid
   use m_sdm_ptlformation, only: &
        sdm_fillboundaries
   
   integer, intent(in) :: KA, KS, KE
   integer, intent(in) :: IA, IS, IE
   integer, intent(in) :: JA, JS, JE
   real(RP), intent(inout) :: DENS(KA,IA,JA)        !! Density [kg/m3]
   real(RP), intent(inout) :: MOMZ(KA,IA,JA)        !! Momentum [kg/s/m2]
   real(RP), intent(inout) :: MOMX(KA,IA,JA)
   real(RP), intent(inout) :: MOMY(KA,IA,JA)
   real(RP), intent(inout) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
   real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_MP_sdm)    !! Ratio of mass of tracer to total mass[kg/kg]
   ! Input variables
   logical,  intent(in) :: sdm_calvar(5)
   real(RP),intent(in) :: sdm_rdnc     ! Number concentration of real droplets
   real(RP),intent(in) :: sdm_sdnmlvol ! Normal volume for number concentration of super droplets
   integer, intent(in) :: sdm_aslset   ! Option for aerosol species
   real(RP),intent(in) :: sdm_inisdnc  ! Initial number of super droplets per sdm_sdnmlvol
   integer,  intent(in) :: sdm_mvexchg
   real(RP), intent(in) :: zph_crs(KA,IA,JA)  ! z physical coordinates
   integer, intent(in) :: ni_sdm                    ! SDM model dimension in x direction
   integer, intent(in) :: nj_sdm                    ! SDM model dimension in y direction
   integer, intent(in) :: nk_sdm                    ! SDM model dimension in z direction
   integer, intent(in) :: sd_num                    ! number of super-droplets
   integer, intent(in) :: sd_numasl     ! number of kind of chemical material contained as water-soluble aerosol in super droplets
   integer, intent(in) :: sd_fmnum                  ! size of array for super-droplets creation
   real(RP), intent(in) :: sd_rkl(IA,JA) ! index[k/real] at lower boundary in SDM
   real(RP), intent(in) :: sd_rku(IA,JA) ! index[k/real] at upper boundary in SDM
   integer, intent(in) :: bufsiz1       ! buffer size for MPI
   integer, intent(in) :: bufsiz2_r8 ! buffer size for MPI (real8)
   integer, intent(in) :: bufsiz2_i8 ! buffer size for MPI (int8)
   integer, intent(in) :: bufsiz2_i2 ! buffer size for MPI (int2)
   integer, intent(in) :: bufsiz2_i4 ! buffer size for MPI (int4)
   ! Input and output variables
   integer(DP), intent(inout) :: sd_n(1:sd_num)    ! multiplicity of super-droplets
   integer(i2), intent(inout) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
   real(RP), intent(inout) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
   real(RP), intent(inout) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
   real(RP), intent(inout) :: sd_ri(1:sd_num) ! face index-i(real) of super-droplets
   real(RP), intent(inout) :: sd_rj(1:sd_num) ! face index-j(real) of super-droplets
   real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets in vertical
   real(RP), intent(inout) :: sd_u(1:sd_num)  ! x-components velocity of super-droplets
   real(RP), intent(inout) :: sd_v(1:sd_num)  ! y-components velocity of super-droplets
   real(RP), intent(inout) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets /z velocity of super-droplets
   real(RP), intent(inout) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
   real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
   type(sdicedef), intent(inout) :: sdi   ! ice phase super-droplets
   integer(DP), intent(out) :: sd_fmn(1:sd_fmnum) ! multiplicity of super-droplets at aerosol formation
   integer(i2), intent(out) :: sd_fmliqice(1:sd_fmnum)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
   real(RP), intent(out) :: sd_fmx(1:sd_fmnum)  ! x-coordinate of super-droplets at aerosol formation
   real(RP), intent(out) :: sd_fmy(1:sd_fmnum)  ! y-coordinate of super-droplets at aerosol formation
   real(RP), intent(out) :: sd_fmz(1:sd_fmnum)  ! z-coordinate of super-droplets at aerosol formation
   real(RP), intent(out) :: sd_fmri(1:sd_fmnum) ! index[i/real] of super-droplets at aerosol formation
   real(RP), intent(out) :: sd_fmrj(1:sd_fmnum) ! index[j/real] of super-droplets at aerosol formation
   real(RP), intent(out) :: sd_fmrk(1:sd_fmnum) ! index[k/real] of super-droplets at aerosol formation
   real(RP), intent(out) :: sd_fmr(1:sd_fmnum)  ! equivalent radius of super-droplets at aerosol formation
   real(RP), intent(inout) :: sd_fmasl(1:sd_fmnum,1:sd_numasl)  ! aerosol mass of super-droplets at aerosol formation
   type(sdicedef), intent(inout) :: sd_fmi           ! ice phase super-droplets
   type(c_rng_uniform_mt), intent(inout) :: sd_rng ! random number generator
   real(DP), intent(inout) :: sd_rand(1:sd_num) ! random numbers
   integer, intent(inout) :: sort_id(1:sd_num)  ! id that super-droplets sorted by grids
   integer, intent(inout) :: sort_key(1:sd_num) ! sort key
   integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1) ! number of super-droplets in each grid
   integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)  ! accumulated number of super-droplets in each grid
   real(RP), intent(inout) :: prec_crs(IA,JA,1:6)! precipitation rate and accumlation [m]
                                               ! 1:rain rate, 2:rain accumulation
                                               ! 3:snow rate, 4:snow accumulation
                                               ! 5:precipitation rate, 6:precipitation accumulation!
   real(DP), intent(inout) ::                                   &
         &                             rbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Receiving buffer (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+4) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    real(DP), intent(inout) ::                                   &
         &                             sbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Sending buffer (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+4) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    integer(DP), intent(inout) ::                                &
         &                             rbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2)
                       ! reciving buffer for MPI (int8)
                       ! dim02 = 1 (multiplicity of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    integer(DP), intent(inout) ::                                &
         &                             sbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2)
                       ! sending buffer for MPI (int8)
                       ! dim02 = 1 (multiplicity of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    integer(i2), intent(inout) ::                                &
         &                             rbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2)
                       ! reciving buffer for MPI (int2)
                       ! dim02 = 1 (status(liq/ice) of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    integer(i2), intent(inout) ::                                &
         &                             sbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2)
                       ! sending buffer for MPI (int2)
                       ! dim02 = 1 (status(liq/ice) of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    integer, intent(inout) ::                                &
         &                             rbuf_i4(1:,1:,1:)
                       ! rbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2) or rbuf_i4(1:1,1:1,1:2) 
                       ! reciving buffer for MPI (int4)
                       ! dim02 = 1 (sdi_nmono of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    integer, intent(inout) ::                                &
         &                             sbuf_i4(1:,1:,1:)
                       ! sbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2) or sbuf_i4(1:1,1:1,1:2) 
                       ! sending buffer for MPI (int4)
                       ! dim02 = 1 (sdi_nmono of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
   ! Output variables
   logical, intent(out) :: lsdmup  ! flag for updating water hydrometeor by SDM
   integer, intent(out) :: sdm_itmp1(1:ni_sdm*nj_sdm*nk_sdm+2) ! temporary array of SDM dimension
   integer, intent(out) :: sdm_itmp2(1:ni_sdm*nj_sdm*nk_sdm+2) ! temporary array of SDM dimension
   integer, intent(out) :: sd_itmp1(1:sd_num) ! temporary array of the size of the number of super-droplets.
   integer, intent(out) :: sd_itmp2(1:sd_num) ! temporary array of the size of the number of super-droplets.
   integer, intent(out) :: sd_itmp3(1:sd_num) ! temporary array of the size of the number o fsuper-droplets.
   real(RP), intent(out) :: sd_dtmp1(1:sd_num) ! temporary array of the size of the number of super-droplets.
   real(RP), intent(out) :: sd_dtmp2(1:sd_num) ! temporary array of the size of the number of super-droplets.
   real(RP), intent(out) :: sd_dtmp3(1:sd_num) ! temporary array of the size of the number of super-droplets.
   real(RP), intent(out) :: sd_dtmp4(1:sd_num) ! temporary array of the size of the number of super-droplets.
   real(RP), intent(out) :: crs_val1p(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val1c(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val2p(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val2c(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(DP), intent(in)  :: dt

   ! Internal shared variables
   integer :: istep_sdm        ! step number of SDM
   integer :: istep_evl        ! step number of {condensation/evaporation} process
   integer :: istep_col        ! step number of {stochastic coalescence} process
   integer :: istep_adv        ! step number of {motion of super-droplets} process
   integer :: istep_mlt        ! step number of {melt/freeze of super-droplets} process
   integer :: istep_sbl        ! step number of {sublimation/deposition of super-droplets} process
   ! Work variables
   integer :: t         ! index
   integer :: k,i,j,n   ! index
   real(RP) :: u_scale(KA,IA,JA)   ! u components of velocity
   real(RP) :: v_scale(KA,IA,JA)   ! v components of velocity
   real(RP) :: w_scale(KA,IA,JA)   ! w components of velocity
   real(RP) :: pres_scale(KA,IA,JA)  ! Pressure
   real(RP) :: t_scale(KA,IA,JA)    ! Temperature 
   real(RP) :: sdm_dtevl  ! time step of {condensation/evaporation} process
   real(RP) :: sdm_dtcol  ! time step of {stochastic coalescence} process
   real(RP) :: sdm_dtadv  ! time step of {motion of super-droplets} process
   real(RP) :: sdm_dtmlt  ! time step of {melt/freeze of super-droplets} process
   real(RP) :: sdm_dtsbl  ! time step of {sublimation/deposition of super-droplets} process
   real(RP) :: tmp_mink

   integer :: IS_bnd, IE_bnd, JS_bnd, JE_bnd
   !---------------------------------------------------------------------

      IS_bnd = IS
      IE_bnd = IE
      JS_bnd = JS
      JE_bnd = JE
      if ( sdm_wbc .eq. 2 ) IS_bnd = max(IS-1,1)
      if ( sdm_ebc .eq. 2 ) IE_bnd = min(IE+1,IA)
      if ( sdm_sbc .eq. 2 ) JS_bnd = max(JS-1,1)
      if ( sdm_nbc .eq. 2 ) JE_bnd = min(JE+1,JA)


      ! Initialize and rename variables
      sdm_dtevl = real( sdm_dtcmph(1),kind=RP )  !! condensation/evaporation
      sdm_dtcol = real( sdm_dtcmph(2),kind=RP )  !! stochastic coalescence
      sdm_dtadv = real( sdm_dtcmph(3),kind=RP )  !! motion of super-droplets
      sdm_dtmlt = real( sdm_dtcmph(4),kind=RP )  !! melting/freezing
      sdm_dtsbl = real( sdm_dtcmph(5),kind=RP )  !! sublimation/depsition

      istep_sdm = nclstp(0)                !! maximum step
      istep_evl = nclstp(1)                !! condensation/evaporation
      istep_col = nclstp(2)                !! stochastic coalescence
      istep_adv = nclstp(3)                !! motion of super-droplets
      istep_mlt = nclstp(4)                !! motion of super-droplets
      istep_sbl = nclstp(5)                !! motion of super-droplets

      lsdmup = .false.

      ! Calculate super-droplets process.
      !   1 : motion of super-droplets (advection, terminal velocity)
      !   2 : melting / freezing [cold-SDM only]
      !   3 : condensation / evaporation
      !   4 : sublimation / deposition [cold-SDM only]
      !   5 : stochastic coalescence
      do t=1,istep_sdm
         !### Run SDM  ###!

         !=== 1 : motion of super-droplets ===!
#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_start("sdm_motion",0,0)
#endif

         if( sdm_calvar(3) .and.                                 &
             mod(t,istep_sdm/istep_adv)==0 .and. sdm_dtadv>0.d0 ) then

            lsdmup = .true.

            ! fill the boundaries with new super-droplets for inflow BC
            if ( (sdm_wbc==2) .or. (sdm_ebc==2) .or. (sdm_sbc==2) .or. (sdm_nbc==2) ) then
               call sdm_fillboundaries(&
                         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                         DENS, RHOT, QTRC,                   &
                         sd_rng,             &
                         sdm_dtcmph,sdm_calvar,      &
                         sdm_rdnc,sdm_aslset,        &
                         sdm_inisdnc,sdm_sdnmlvol,   &
                         sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_dtmp3,    & ! sd_dtmp3 is for sd_z
                         sd_ri,sd_rj,sd_rk,sd_r,sd_asl,sdi, &
                         sd_fmnum,        sd_fmn,sd_fmliqice,sd_fmx,sd_fmy,sd_fmz,    &
                         sd_fmri,sd_fmrj,sd_fmrk,sd_fmr,sd_fmasl,sd_fmi, &
                         sd_dtmp1, sd_dtmp2, sd_itmp1)
            end if

            ! get the terminal velocity of super-droplets
            !! diagnose necessary fluid variables
            call sdm_rhot_qtrc2p_t( &
                 KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                 RHOT,QTRC,DENS,pres_scale,t_scale)
            !! evaluate the terminal velocity
            call sdm_getvz_liq( KA, IA,IS,IE, JA,JS,JE, &
                           pres_scale,DENS,t_scale,            &
                           sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_r,sd_vz,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'no_interpolation' )
            if( sdm_cold )then
               call sdm_getvz_ice( KA, IA,IS,IE, JA,JS,JE, &
                           DENS,t_scale,            &
                           sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sdi,sd_vz,  &
                           sd_itmp1,'no_interpolation' )
            end if

            ! for predictor-corrector
            !$omp parallel do simd schedule(static,256)
            do n = 1, sd_num
               sd_dtmp3(n) = sd_rk(n)
               if ( sd_rk(n) > VALID2INVALID ) then
                  sd_dtmp1(n) = sd_x(n)
                  sd_dtmp2(n) = sd_y(n)
                  sd_dtmp6(n) = sd_vz(n)
               end if
            end do

            ! get the moving velocity of super-droplets
            !! diagnose necessary fluid variables
            call sdm_rho_mom2uvw( KA, KS, KE, IA, IS_bnd, IE_bnd, JA, JS_bnd, JE_bnd, &
                                  DENS,MOMX,MOMY,MOMZ,u_scale,v_scale,w_scale)
            !! update the HALO region of the fluid variables
            call COMM_vars8( u_scale(:,:,:), 1 )
            call COMM_vars8( v_scale(:,:,:), 2 )
            call COMM_vars8( w_scale(:,:,:), 3 )
            call COMM_wait ( u_scale(:,:,:), 1, FILL_BND=.false. )
            call COMM_wait ( v_scale(:,:,:), 2, FILL_BND=.false. )
            call COMM_wait ( w_scale(:,:,:), 3, FILL_BND=.false. )

            !! evaluate the velocity of each SD
            call sdm_getvel(KA,KS,KE, IA,IS,IE, JA,JS,JE, &
                            u_scale,v_scale,w_scale,                   &
                            sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_u,sd_v,sd_vz)

            ! momentum coupling
            !! rhow coupling
            if( sdm_mvexchg>=1 ) then

               call sdm_sd2rhow( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             zph_crs,crs_val1p,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2p,sd_itmp1)

               do j = JS, JE
                  do i = IS, IE
                     MOMZ(KS-1,i,j) = 0
                     do k = KS, KE-1
                        MOMZ(k,i,j) = MOMZ(k,i,j) - ( crs_val1p(k,i,j) + crs_val1p(k+1,i,j) )*0.5_RP*GRAV*sdm_dtadv
                     enddo
                     MOMZ(KE,i,j) = 0
                  enddo
               enddo

               call COMM_vars8( MOMZ(:,:,:), 1 )
               call COMM_wait ( MOMZ(:,:,:), 1 )

            end if

            !! calculate the current SD momentum field and update fluid
            !!(not supported yet)
            if( sdm_mvexchg>=2 ) then
               ! calculate the current SD momentum field (a)

               ! update fluid using (a) current and (b) previous SD momentum fields

            end if

            ! { motion } in SDM
            !! evaluate the motion eq.
!!$            !!!! explicit Euler scheme
!!$            call sdm_move(sdm_dtadv,                         &
!!$                          sd_num,sd_u,sd_v,sd_vz,sd_x,sd_y,sd_rk)
!!$
            !!!! predictor-corrector scheme
            !!!!!! predictor step
            call sdm_move(IA,IS,IE,JA,JS,JE,sdm_dtadv,                         &
                          sd_num,sd_u,sd_v,sd_vz,sd_dtmp1,sd_dtmp2,sd_ri,sd_rj,sd_dtmp3)
            !!!!!! corrector step
            tmp_mink = real(KS-1,kind=RP)
            !$omp  parallel do simd schedule(static,256)
            do n = 1, sd_num
               if( sd_dtmp3(n)<VALID2INVALID ) cycle
               sd_dtmp3(n) = max(tmp_mink,sd_dtmp3(n))
            end do
            !$omp end parallel do simd
            call sdm_getvel(KA,KS,KE, IA,IS,IE, JA,JS,JE, &
                            u_scale,v_scale,w_scale,                   &
                            sd_num,sd_dtmp1,sd_dtmp2,sd_ri,sd_rj,sd_dtmp3,sd_dtmp4,sd_dtmp5,sd_dtmp6)
            !$omp parallel do simd schedule(static,256)
            do n = 1, sd_num
               if ( sd_dtmp3(n) > VALID2INVALID ) then
                  sd_dtmp4(n) = 0.5_RP*(sd_u(n) +sd_dtmp4(n))
                  sd_dtmp5(n) = 0.5_RP*(sd_v(n) +sd_dtmp5(n))
                  sd_dtmp6(n) = 0.5_RP*(sd_vz(n)+sd_dtmp6(n))
               end if
            end do
            call sdm_move(IA,IS,IE,JA,JS,JE,sdm_dtadv,                         &
                          sd_num,sd_dtmp4,sd_dtmp5,sd_dtmp6,sd_x,sd_y,sd_ri,sd_rj,sd_rk)

            ! lateral boundary routine in SDM
            !! judge super-droplets as invalid or valid in horizontal
            !! do MPI communication to send/receiv SDs
            tag = 0
            call sdm_boundary( IS, IE, JS, JE, &
                             sdm_wbc,sdm_ebc,sdm_sbc,sdm_nbc,                           &
                             sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_rk,     &
                             sd_u,sd_v,sd_vz,sd_r,sd_asl,sdi,               &
                             bufsiz1,                                    &
                             bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4,      &
                             sd_itmp1,                              &
                             rbuf_r8,sbuf_r8,rbuf_i8,sbuf_i8,rbuf_i2,sbuf_i2,rbuf_i4,sbuf_i4) 

            ! judge super-droplets as invalid or valid in vartical
            call sdm_jdginvdv( IA,IS,IE, JA,JS,JE, &
                             sd_rkl,sd_rku,sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk)

            !! save the SDM momentum field at new position (b)
            !!(not supported yet)
            if( sdm_mvexchg>=2 ) then

            end if

         end if

#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_stop("sdm_motion",0,0)
#endif

         !=== 2 : melt / freeze ===!
#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_start("sdm_meltfreeze",0,0)
#endif
         if( sdm_cold .and. sdm_calvar(4) .and.                                 &
             mod(t,istep_sdm/istep_mlt)==0 .and. sdm_dtmlt>0.d0 ) then

            lsdmup = .true.

            ! get density of solid-water before melt/freeze
            !! cres_val1p is the rhosol before meltfreeze
            call sdm_sd2rhosol( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             zph_crs,crs_val1p,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2p,sd_itmp1)

            ! { melting/freezing } in SDM
            !! diagnose necessary fluid variables
            call sdm_rhot_qtrc2p_t( &
                 KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                 RHOT,QTRC,DENS,pres_scale,t_scale)
            !! update the phase of SDs
            call sdm_meltfreeze(                              &
              KA,KS,KE, IA,IS,IE, JA,JS,JE, &
              t_scale,pres_scale,QTRC(:,:,:,I_QV_sdm),         &
              sd_num,sd_liqice,sd_x,sd_y,      &
              sd_r,sd_ri,sd_rj,sd_rk,sdi )

            ! get density of solid-water after melt/freeze
            !! here cres_val1c is the rhosol after meltfreeze
            call sdm_sd2rhosol( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             zph_crs,crs_val1c,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2c,sd_itmp1)

            ! exchange the heat to fluid variables
            call sdm_meltfreeze_updatefluid( &
                             KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             RHOT,QTRC,DENS,crs_val1p,crs_val1c)
            !! update the HALO region of the fluid variables
            call COMM_vars8( RHOT(:,:,:), 1 )
            call COMM_wait ( RHOT(:,:,:), 1 )

         end if

#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_stop("sdm_meltfreeze",0,0)
#endif

         !=== 3 : condensation / evaporation ===!
#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_start("sdm_condevp",0,0)
#endif
         if( sdm_calvar(1) .and.                                 &
             mod(t,istep_sdm/istep_evl)==0 .and. sdm_dtevl>0.d0 ) then

            lsdmup = .true.

            ! get density of liquid-water(qw) before process-1
            !! cres_val1p is the rhow before condevp
            call sdm_sd2rhow( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             zph_crs,crs_val1p,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2p,sd_itmp1)

            ! { condensation/evaporation } in SDM
            !! diagnose necessary fluid variables
            call sdm_rhot_qtrc2p_t( &
                 KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                 RHOT,QTRC,DENS,pres_scale,t_scale)
            !! update the equivalent radius of SDs
            call sdm_condevp(KA,KS,KE, IA,IS,IE, JA,JS,JE, &
                             sdm_aslset,            &
                             sdm_aslmw,sdm_aslion,sdm_dtevl,      &
                             pres_scale,t_scale,QTRC(:,:,:,I_QV_sdm), &
                             sd_num,sd_numasl,sd_liqice,sd_x,sd_y,sd_r,sd_asl,&
                             sd_ri,sd_rj,sd_rk)

            ! get density of liquid-water(qw) after process-1
            !! here cres_val1c is the rhow after condevp
            call sdm_sd2rhow( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             zph_crs,crs_val1c,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2c,sd_itmp1)

            ! exchange the vapor and heat to fluid variables
            call sdm_condevp_updatefluid( &
                             KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             RHOT,QTRC,DENS,crs_val1p,crs_val1c)
            !! update the HALO region of the fluid variables
            call COMM_vars8( RHOT(:,:,:), 1 )
            call COMM_vars8( QTRC(:,:,:,I_QV_sdm), 2 )
            call COMM_vars8( DENS(:,:,:), 3 )
            call COMM_wait ( RHOT(:,:,:), 1 )
            call COMM_wait ( QTRC(:,:,:,I_QV_sdm), 2 )
            call COMM_wait ( DENS(:,:,:), 3 )

         end if

#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_stop("sdm_condevp",0,0)
#endif

         !=== 4 : sublimation / deposition ===!
#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_start("sdm_subldep",0,0)
#endif
         if( sdm_cold .and. sdm_calvar(5) .and.                                 &
             mod(t,istep_sdm/istep_sbl)==0 .and. sdm_dtsbl>0.d0 ) then

            lsdmup = .true.

            ! get density of solid-water before sublimation/deposition
            !! cres_val1p is the rhosol before sublimation/deposition
            call sdm_sd2rhosol( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             zph_crs,crs_val1p,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2p,sd_itmp1)

            ! { sublimation/deposition } in SDM
            !! diagnose necessary fluid variables
            call sdm_rhot_qtrc2p_t(&
                 KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                 RHOT,QTRC,DENS,pres_scale,t_scale)
            !! evaluate the terminal velocity
            call sdm_getvz_ice( KA, IA,IS,IE, JA,JS,JE, &
                           DENS,t_scale,            &
                           sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sdi,sd_vz,  &
                           sd_itmp1,'no_interpolation' )
            !! update the equivalent radius of SDs
            call sdm_subldep(KA,KS,KE, IA,IS,IE, JA,JS,JE,             &
                             sdm_dtsbl,      &
                             pres_scale,t_scale,QTRC(:,:,:,I_QV_sdm),DENS,&
                             sd_num,sd_liqice,sd_x,sd_y,sdi,sd_vz,&
                             sd_ri,sd_rj,sd_rk,sd_itmp1)

            ! get density of solid-water after sublimation/deposition
            !! here cres_val1c is the rhosol after sublimation/deposition
            call sdm_sd2rhosol( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             zph_crs,crs_val1c,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2c,sd_itmp1)

            ! exchange the vapor and heat to fluid variables
            call sdm_subldep_updatefluid( &
                             KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             RHOT,QTRC,DENS,crs_val1p,crs_val1c)
            !! update the HALO region of the fluid variables
            call COMM_vars8( RHOT(:,:,:), 1 )
            call COMM_vars8( QTRC(:,:,:,I_QV_sdm), 2 )
            call COMM_vars8( DENS(:,:,:), 3 )
            call COMM_wait ( RHOT(:,:,:), 1 )
            call COMM_wait ( QTRC(:,:,:,I_QV_sdm), 2 )
            call COMM_wait ( DENS(:,:,:), 3 )

         end if

#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_stop("sdm_subldep",0,0)
#endif

         !=== 5 : stochastic coalescence ===!
#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_start("sdm_coales",0,0)
#endif

         if( sdm_calvar(2) .and.                                 &
             mod(t,istep_sdm/istep_col)==0 .and. sdm_dtcol>0.d0 ) then

            lsdmup = .true.

            ! get the terminal velocity of super-droplets
            !! diagnose necessary fluid variables
            call sdm_rhot_qtrc2p_t(&
                 KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                 RHOT,QTRC,DENS,pres_scale,t_scale)
            !! evaluate the terminal velocity
            call sdm_getvz_liq( KA, IA,IS,IE, JA,JS,JE, &
                           pres_scale,DENS,t_scale,            &
                           sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_r,sd_vz,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'no_interpolation' )
            if( sdm_cold )then
               call sdm_getvz_ice( KA, IA,IS,IE, JA,JS,JE, &
                           DENS,t_scale,            &
                           sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sdi,sd_vz,  &
                           sd_itmp1,'no_interpolation' )
            end if

            ! { coalescence } in SDM
            if( sdm_cold )then
               ! get density of solid-water before riming
               !! cres_val1p is the rhosol before riming
               call sdm_sd2rhosol( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             zph_crs,crs_val1p,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2p,sd_itmp1)

               call sdm_coales_cold(KA,KS,KE, IA,IS,IE, JA,JS,JE, &
                            sdm_colkrnl,sdm_colbrwn,sdm_aslset,         &
                            sdm_aslrho,sdm_dtcol,                       &
                            pres_scale,t_scale,QTRC(:,:,:,I_QV_sdm),DENS,&
                            zph_crs,                                    &
                            ni_sdm,nj_sdm,nk_sdm,sd_num,sd_numasl,      &
                            sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk,     &
                            sdi,                                        & 
                            sort_id,sort_key,sort_freq,sort_tag,        &
                            sd_rng,sd_rand,                             &
                            sdm_itmp1,sdm_itmp2,                        &
                            sd_itmp1(1:sd_num),sd_itmp2(1:sd_num),  &
                            sd_dtmp1)

               ! get density of solid-water after riming
               !! here cres_val1c is the rhosol after riming
               call sdm_sd2rhosol( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             zph_crs,crs_val1c,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2c,sd_itmp1)

               ! exchange the heat to fluid variables
               call sdm_meltfreeze_updatefluid( &
                             KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             RHOT,QTRC,DENS,crs_val1p,crs_val1c)
               !! update the HALO region of the fluid variables
               call COMM_vars8( RHOT(:,:,:), 1 )
               call COMM_wait ( RHOT(:,:,:), 1 )

            else
               call sdm_coales(KA, IA,IS,IE, JA,JS,JE, &
                            sdm_colkrnl,sdm_colbrwn,sdm_aslset,         &
                            sdm_aslrho,sdm_dtcol,                       &
                            pres_scale, t_scale,                        &
                            zph_crs,                                    &
                            ni_sdm,nj_sdm,nk_sdm,sd_num,sd_numasl,      &
                            sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk,     &
                            sort_id,sort_key,sort_freq,sort_tag,        &
                            sd_rng,sd_rand,                             &
                            sdm_itmp1,sdm_itmp2,                        &
                            sd_itmp1(1:sd_num),sd_itmp2(1:sd_num),  &
                            sd_dtmp1)
            end if

         end if

#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_stop("sdm_coales",0,0)
#endif

      end do

    ! Convert super-droplets to precipitation

#ifdef _FAPP_
      ! Section specification for fapp profiler
      call fapp_start("sdm_sd2prec",0,0)
#endif
      call sdm_sd2prec(KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                       real(dt,kind=RP),                  &
                       prec_crs,                          &
                       sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_r,sdi,sd_ri,sd_rj,sd_rk,  &
                       sd_itmp1,sd_itmp2,crs_val1c(1,1:IA,1:JA))

#ifdef _FAPP_
      ! Section specification for fapp profiler
      call fapp_stop("sdm_sd2prec",0,0)
#endif

    return
  end subroutine sdm_calc
!---------------------------------------------------------------------------------------
end module m_sdm_calc
