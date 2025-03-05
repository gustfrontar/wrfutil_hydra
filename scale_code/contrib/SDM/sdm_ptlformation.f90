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
!! @li      2014-06-14 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2016-07-14 (S.Shima) [mod] Evaluation of bufsiz1, bufsiz2_r8, etc. moved from sdm_allocinit to sdm_numset
!! @li      2018-06-30 (S.Shima) [add] rime mass and number of monomers as SD attributes
!! @li      2019-10-07 (S.Shima) [add] aslset=5 for DYCOMSII(RF02) (Ackerman et al. 2009)
!! @li      2021-03-06 (S.Shima) [mod] m_sdm_numset was renamed as m_sdm_ptlformation. sdm_iniset was moved from scale_atmos_phy_mp_sdm to m_sdm_ptlformation.
!! @li      2022-05-21 (S.Shima) [mod] support terrain following and vertically stretched coordinate
!! @li      2022-09-14 (S.Shima) [fix] MPI_REAL -> MPI_REAL8
!! @li      2022-09-14 (S.Shima) [mod] sdnum_s2c -> nint(sdininum_s2c)
!! @li      2022-09-15 (S.Shima) [add] subroutine sdm_sdcreate
!! @li      2022-09-17 (S.Shima) [add] subroutine sdm_fillboundaries and inflow boundary condition option
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module m_sdm_ptlformation
  use scale_io

  implicit none
  private
  public :: sdm_numset,sdm_iniset,sdm_fillboundaries

  interface gen_rand_array
     module procedure gen_rand_array_DP
     module procedure gen_rand_array_SP
  end interface gen_rand_array

contains
  subroutine sdm_numset(  &
      KA, KS, KE, IA, IS, IE, JA, JS, JE, &
      sdm_extbuf,                &
      sdm_aslset, sdm_sdnmlvol,  &
      sdm_inisdnc,               &
      sdm_aslfmsdnc,             &
      sdm_zlower, sdm_zupper,    &
      minzph, sdininum_s2c,      &
      sdfmnum_s2c, sdnum_s2c,    &
      ni_s2c, nj_s2c, nk_s2c,zph_crs )

    use mpi
    use scale_prc, only: &
         PRC_LOCAL_COMM_WORLD
    use scale_prc_cartesC, only: &
         PRC_HAS_W,   &
         PRC_HAS_E,   &
         PRC_HAS_S,   &
         PRC_HAS_N
    use scale_precision
    use scale_topography, only: &
         TOPO_Zsfc => TOPOGRAPHY_Zsfc
    use scale_atmos_grid_cartesC_real, only: &
!         REAL_CZ
         REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
    use scale_atmos_grid_cartesC, only: &
         GRID_FX => ATMOS_GRID_CARTESC_FX,    &
         GRID_FY => ATMOS_GRID_CARTESC_FY
    use m_sdm_common, only: &
         sdm_cold,     &
         sdnumasl_s2c, & 
         bufsiz1,      &      ! Buffer size for MPI
         bufsiz2_r8,   &      ! buffer size for MPI (real8)
         bufsiz2_i8,   &      ! buffer size for MPI (int8)
         bufsiz2_i2,   &      ! buffer size for MPI (int2)
         bufsiz2_i4,   &      ! buffer size for MPI (int4)
         sdm_aslmw,    &
         xmax_sdm, &    ! Width of the domain. xmax_sdm = FX(IE)-FX(IS-1)
         ymax_sdm, &    ! Depth of the domain. ymax_sdm = FY(JE)-FY(JS-1)
         sdm_wbc,sdm_ebc,sdm_nbc,sdm_sbc
    
    ! These argument variables are all defined in m_sdm_common. Perhaps we should make it arguments only or common module only?
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE ! S:start, E: end of active grids. A: num of grids including HALO.
    integer, intent(in) :: sdm_extbuf    ! Rate of increase for extra buffer of super droplets
    integer, intent(in) :: sdm_aslset    ! Option for aerosol species
    real(RP),intent(in) :: sdm_sdnmlvol  ! Normal volume for number concentration
    real(RP),intent(in) :: sdm_inisdnc   ! Number of super droplets at initial per sdm_sdnmlvol
    real(RP),intent(in) :: sdm_aslfmsdnc ! Number of super droplets at aerosol formation per sdm_sdnmlvol
    real(RP),intent(in) :: sdm_zlower    ! Lowest limitaion of initial super droplets distribution
    real(RP),intent(in) :: sdm_zupper    ! Highest limitaion of initial super droplets distribution
    real(RP),intent(out) :: minzph       ! Lowest surface height of the domain
    real(RP),intent(out) :: sdininum_s2c ! Initial number of SDs used 
    integer,intent(out) :: sdfmnum_s2c   ! Number of temp SD arrays used for filling the boundary or aerosol formation
    integer, intent(out) :: sdnum_s2c    ! Total number of SDs allocated
    integer, intent(out) :: ni_s2c, nj_s2c, nk_s2c ! Num of grids used by SDs
    real(RP),intent(out) :: zph_crs(KA,IA,JA) ! This is Real_CZ. Waste of memory. Will be removed in the near future.
    
    ! Work variables
    real(RP) :: dtmp,dtmp_wbc,dtmp_ebc,dtmp_nbc,dtmp_sbc       ! Temporary
    real(RP) :: xmin,xmax,ymin,ymax,zmin,zmax                  ! Temporary
    integer :: i, j, k, n  ! index
    integer :: ierr        ! Error descriptor
    logical :: corner_wn, corner_ws, corner_en, corner_es

    ! Num of grids used by SDs
    ni_s2c = IE-IS+1
    nj_s2c = JE-JS+1
    nk_s2c = KE-KS+1

    minzph = TOPO_Zsfc(IS,JS)
    do j = 1, JA
       do i = 1, IA
          minzph = min( minzph, TOPO_Zsfc(i,j) )
       enddo
    enddo

    call mpi_allreduce(minzph,dtmp,1,mpi_real8,MPI_MIN,              &
         &                     PRC_LOCAL_COMM_WORLD,ierr)
    minzph = dtmp

    ! Waste of memory. zph_crs will be removed in the near future. At least we can use POINTER.
    ! Note also that CReSS uses (i,j,k) but SCALE uses (k,i,j). Be careful..
    do k = 1, KA
       do j = 1, JA
          do i = 1, IA
!             zph_crs(k,i,j) = REAL_CZ(k,i,j)
             zph_crs(k,i,j) = REAL_FZ(k,i,j)
          enddo
       enddo
    enddo
    
    !### setting number of super-droplets at initial ###!
    sdininum_s2c = real(xmax_sdm,kind=RP) * real(ymax_sdm,kind=RP)         &
         * real(sdm_zupper-(sdm_zlower+minzph),kind=RP)       &
         * real(sdm_inisdnc,kind=RP)/real(sdm_sdnmlvol,kind=RP)
    
    !### setting number of super-droplets at aerosol formation ###!
    if( abs(sdm_aslset)>10 ) then
       dtmp = real(xmax_sdm,kind=RP) * real(ymax_sdm,kind=RP)              &
            * real(sdm_zupper-(sdm_zlower+minzph),kind=RP)            &
            * real(sdm_aslfmsdnc,kind=RP)/real(sdm_sdnmlvol,kind=RP)
       sdfmnum_s2c = nint(dtmp)
    else
       sdfmnum_s2c = 1
    end if

    !### setting number of super-droplets to fill the surrounding halo for inflow boundary ###!
    zmin = sdm_zlower + minzph
    zmax = sdm_zupper
    corner_ws = .true. ! if true, west south corner is not taken into account yet
    corner_wn = .true.
    corner_es = .true.
    corner_en = .true.
    dtmp_wbc=0.0_RP
    dtmp_ebc=0.0_RP
    dtmp_nbc=0.0_RP
    dtmp_sbc=0.0_RP

    !! west boundary
    if( (sdm_wbc==2) .and. (.NOT. PRC_HAS_W ) ) then
       xmin = GRID_FX(IS-2)
       xmax = GRID_FX(IS-1)
       if( (.NOT. PRC_HAS_S ) .and. corner_ws ) then ! west south corner
          ymin = GRID_FY(JS-2)
          corner_ws = .false.
       else
          ymin = GRID_FY(JS-1)
       end if
       if( (.NOT. PRC_HAS_N ) .and. corner_wn ) then ! west north corner
          ymax = GRID_FY(JE+1)
          corner_wn = .false.
       else
          ymax = GRID_FY(JE)
       end if

       dtmp_wbc = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)*sdm_inisdnc/sdm_sdnmlvol
    end if

    !! east boundary
    if( (sdm_ebc==2) .and. (.NOT. PRC_HAS_E ) ) then
       xmin = GRID_FX(IE)
       xmax = GRID_FX(IE+1)
       if( (.NOT. PRC_HAS_S ) .and. corner_es ) then ! east south corner
          ymin = GRID_FY(JS-2)
          corner_es = .false.
       else
          ymin = GRID_FY(JS-1)
       end if
       if( (.NOT. PRC_HAS_N ) .and. corner_en ) then ! east north corner
          ymax = GRID_FY(JE+1)
          corner_en = .false.
       else
          ymax = GRID_FY(JE)
       end if

       dtmp_ebc = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)*sdm_inisdnc/sdm_sdnmlvol
    end if

    !! south boundary
    if( (sdm_sbc==2) .and. (.NOT. PRC_HAS_S ) ) then
       ymin = GRID_FY(JS-2)
       ymax = GRID_FY(JS-1)
       if( (.NOT. PRC_HAS_W ) .and. corner_ws ) then ! west south corner
          xmin = GRID_FX(IS-2)
          corner_ws = .false.
       else
          xmin = GRID_FX(IS-1)
       end if
       if( (.NOT. PRC_HAS_E ) .and. corner_es ) then ! east south corner
          xmax = GRID_FX(IE+1)
          corner_es = .false.
       else
          xmax = GRID_FX(IE)
       end if

       dtmp_sbc = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)*sdm_inisdnc/sdm_sdnmlvol
    end if

    !! north boundary
    if( (sdm_nbc==2) .and. (.NOT. PRC_HAS_N ) ) then
       ymin = GRID_FY(JE)
       ymax = GRID_FY(JE+1)
       if( (.NOT. PRC_HAS_W ) .and. corner_wn ) then ! west north corner
          xmin = GRID_FX(IS-2)
          corner_wn = .false.
       else
          xmin = GRID_FX(IS-1)
       end if
       if( (.NOT. PRC_HAS_E ) .and. corner_en ) then ! east north corner
          xmax = GRID_FX(IE+1)
          corner_en = .false.
       else
          xmax = GRID_FX(IE)
       end if

       dtmp_nbc = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)*sdm_inisdnc/sdm_sdnmlvol
    end if

    sdfmnum_s2c = max(sdfmnum_s2c,nint(dtmp_wbc),nint(dtmp_ebc),nint(dtmp_sbc),nint(dtmp_nbc))

    !### setting buffer size of super-droplets ###!
!!$    bufsiz = nint(sdininum_s2c)
!!$    bufsiz = nint( real(bufsiz) * (real(sdm_extbuf)*1.E-2_RP) )
!!$    sdnum_s2c = nint(sdininum_s2c) + bufsiz

    ! Get numbers of kind of chemical material contained as

    ! water-soluble aerosol in super droplets.
    if( abs(mod(sdm_aslset,10))==1 ) then
       !### init+rest : (NH4)2SO4 ###!
       sdnumasl_s2c = 1
    else if( abs(mod(sdm_aslset,10))==2 ) then
       if( abs(sdm_aslset)==2 ) then
          !### init : NaCl ###!
          sdnumasl_s2c = 1
       else if( abs(sdm_aslset)==12 ) then
          !### init : NaCl, rest : (NH4)2SO4 ###!
          sdnumasl_s2c = 2
       end if
    else if( abs(mod(sdm_aslset,10))==3 ) then
       !### init+rest : (NH4)2SO4, NaCl, ... ###!
       sdnumasl_s2c = 2   !! default
       do n=1,20
          if( sdm_aslmw(n)>0.0_RP ) then
             sdnumasl_s2c = n + 2
          end if
       end do
    else if( abs(mod(sdm_aslset,10))==5 ) then
       !### init+rest : (NH4)HSO4 ###! 
       sdnumasl_s2c = 1
    end if

    bufsiz1 = nint( sdininum_s2c*(real(sdm_extbuf)*1.E-2_RP) )
    sdnum_s2c = nint(sdininum_s2c) + bufsiz1

    bufsiz2_r8 = 7 + sdnumasl_s2c    !! x,y,rk,u,v,wc(vz),r,asl
    bufsiz2_i8 = 1                   !! n
    bufsiz2_i2 = 1                   !! liqice
    bufsiz2_i4 = 1                   !! nmono (cold)

    if( sdm_cold ) then
       bufsiz2_r8 = bufsiz2_r8 + 5   !! re,ro,rho,tf,mrime
    end if

    return
  end subroutine sdm_numset
  !-----------------------------------------------------------------------------
  subroutine sdm_iniset(  &
                         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                         DENS, RHOT, QTRC,                   &
                         sdm_dtcmph,         &
                         sdm_rdnc,sdm_sdnmlvol,sdm_aslset,   &
                         sdm_inisdnc,sdm_zlower,             &
                         sdm_zupper,sdm_calvar)
  !***********************************************************************
  ! Input variables
    use scale_precision
    use scale_prc, only: &
         mype => PRC_myrank
    use scale_io, only: &
         IO_L,IO_FID_LOG
    use scale_atmos_grid_cartesC, only: &
         GRID_FX => ATMOS_GRID_CARTESC_FX,    &
         GRID_FY => ATMOS_GRID_CARTESC_FY
    use m_sdm_common, only: &
         rng_s2c,RANDOM_IN_BASENAME,minzph, & 
         sdininum_s2c,sdnumasl_s2c,sdnum_s2c,INVALID, &
         sdn_s2c,sdliqice_s2c,sdri_s2c,sdrj_s2c,sdice_s2c, &
         sdasl_s2c,sdx_s2c,sdy_s2c,sdz_s2c,sdr_s2c,sdrk_s2c, &
         sd_dtmp1,sd_dtmp2,QA_MP_sdm

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: DENS(KA,IA,JA) ! Density     [kg/m3]
    real(RP), intent(in) :: RHOT(KA,IA,JA) ! DENS * POTT [K*kg/m3]
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_MP_sdm) ! ratio of mass of tracer to total mass[kg/kg]
    real(RP),intent(in) :: sdm_dtcmph(5)  ! Time interval of cloud micro physics
    logical, intent(in) :: sdm_calvar(5)! Flag for cond./coll/move calculation
    real(RP),intent(in) :: sdm_rdnc     ! Number concentration of real droplets
    real(RP),intent(in) :: sdm_sdnmlvol ! Normal volume for number concentration of super droplets
    integer, intent(in) :: sdm_aslset   ! Option for aerosol species
    real(RP),intent(in) :: sdm_inisdnc  ! Initial number of super droplets per sdm_sdnmlvol
    real(RP),intent(in) :: sdm_zlower   ! Lower limitaion of initial SDs position
    real(RP),intent(in) :: sdm_zupper   ! Upper limitaion of initial SDs position

    real(RP):: xmin,xmax,ymin,ymax,zmin,zmax
    integer :: n                        ! index
    !---------------------------------------------------------------------
 
    !### Initialize super-droplets ###!
    !! initialize super-droplets(1:nint(sdininum_s2c))
    zmin = sdm_zlower + minzph
    zmax = sdm_zupper
    xmin = GRID_FX(IS-1)
    xmax = GRID_FX(IE)
    ymin = GRID_FY(JS-1)
    ymax = GRID_FY(JE)
    call sdm_sdcreate(   KA,KS,KE, IA,IS,IE, JA,JS,JE, &
                         DENS, RHOT, QTRC,                   &
                         rng_s2c,                &
                         xmin,xmax,ymin,ymax,zmin,zmax,sdm_dtcmph,sdm_calvar,    &
                         sdm_rdnc,sdm_aslset,   &
                         sdm_inisdnc,sdm_sdnmlvol,              &
                         nint(sdininum_s2c),sdnumasl_s2c,sdn_s2c,sdliqice_s2c,sdasl_s2c, sdx_s2c, sdy_s2c,        &
                         sdz_s2c, sdr_s2c,                   &
                         sdri_s2c,sdrj_s2c,sdrk_s2c,sdice_s2c,&
                         sd_dtmp1,sd_dtmp2)

    !! set invalid flags to super-droplets(nint(sdininum_s2c)+1:sdnum_s2c)
    do n=nint(sdininum_s2c)+1,sdnum_s2c
       sdrk_s2c(n) = INVALID
    end do

    !### Output log message ###!
    if( mype==0 ) then
       if( IO_L ) then
          write(IO_FID_LOG,*)
          write(IO_FID_LOG,'(a)')"  ### [SDM] : super-droplets intialized  ###"
          write(IO_FID_LOG,*)
       endif
    end if
      
    return

  end subroutine sdm_iniset
  !-----------------------------------------------------------------------------
  subroutine sdm_fillboundaries(  &
                         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                         DENS, RHOT, QTRC,                   &
                         sd_rng,             & 
                         sdm_dtcmph,sdm_calvar,      &
                         sdm_rdnc,sdm_aslset,        &
                         sdm_inisdnc,sdm_sdnmlvol,   &
                         sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_z,    &
                         sd_ri,sd_rj,sd_rk,sd_r,sd_asl,sdi, &
                         sd_fmnum,        sd_fmn,sd_fmliqice,sd_fmx,sd_fmy,sd_fmz,    &
                         sd_fmri,sd_fmrj,sd_fmrk,sd_fmr,sd_fmasl,sd_fmi, &
                         sd_dtmp1, sd_dtmp2, ilist)
  !***********************************************************************
  ! Input variables
    use scale_precision
    use scale_prc, only: &
         PRC_abort
    use scale_prc_cartesC, only: &
         PRC_HAS_W,   &
         PRC_HAS_E,   &
         PRC_HAS_S,   &
         PRC_HAS_N
    use scale_io, only: &
         IO_L,IO_FID_LOG
    use scale_atmos_grid_cartesC, only: &
         GRID_FX => ATMOS_GRID_CARTESC_FX,    &
         GRID_FY => ATMOS_GRID_CARTESC_FY
    use m_sdm_common, only: &
         minzph,sdm_zlower,sdm_zupper,& 
         sdm_wbc,sdm_ebc,sdm_nbc,sdm_sbc,i2,sdicedef,VALID2INVALID,sdm_cold,QA_MP_sdm
    use rng_uniform_mt, only: &
         c_rng_uniform_mt

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: DENS(KA,IA,JA) ! Density     [kg/m3]
    real(RP), intent(in) :: RHOT(KA,IA,JA) ! DENS * POTT [K*kg/m3]
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_MP_sdm) ! ratio of mass of tracer to total mass[kg/kg]
    type(c_rng_uniform_mt), intent(inout) :: sd_rng     ! random number generator
    real(RP),intent(in) :: sdm_dtcmph(5)  ! Time interval of cloud micro physics
    logical, intent(in) :: sdm_calvar(5)! Flag for cond./coll/move calculation
    real(RP),intent(in) :: sdm_rdnc     ! Number concentration of real droplets
    real(RP),intent(in) :: sdm_sdnmlvol ! Normal volume for number concentration of super droplets
    integer, intent(in) :: sdm_aslset   ! Option for aerosol species
    real(RP),intent(in) :: sdm_inisdnc  ! Initial number of super droplets per sdm_sdnmlvol
    integer, intent(in) :: sd_num   ! number of super-droplets
    integer, intent(in) :: sd_numasl! number of kind of chemical material contained as water-soluble aerosol in super droplets
    integer, intent(in) :: sd_fmnum ! size of array for super-droplets creation
    integer(DP), intent(inout) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    integer(i2), intent(inout) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
    real(RP), intent(inout) :: sd_x(1:sd_num)      ! x-coordinate of super-droplets
    real(RP), intent(inout) :: sd_y(1:sd_num)      ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_z(1:sd_num)      ! z-coordinate of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)     ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)     ! index[j/real] of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num)     ! index[k/real] of super-droplets
    real(RP), intent(inout) :: sd_r(1:sd_num)      ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl)  ! aerosol mass of super-droplets
    type(sdicedef), intent(inout) :: sdi           ! ice phase super-droplets
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

    real(RP),intent(out) :: sd_dtmp1(1:sd_num) ! temp 
    real(RP),intent(out) :: sd_dtmp2(1:sd_num) ! temp
    integer, intent(out) :: ilist(1:sd_num)

    real(RP):: dtmp,dtmp_wbc,dtmp_ebc,dtmp_nbc,dtmp_sbc       ! Temporary
    real(RP):: xmin,xmax,ymin,ymax,zmin,zmax
    integer :: n,k,cnt,available_sdnum,id_invd
    logical :: corner_wn, corner_ws, corner_en, corner_es

    !---------------------------------------------------------------------

    ! check the size of sd_dtmp and sd_fm
    if(sd_num<sd_fmnum)then
       if( IO_L ) then
          LOG_ERROR("sdm_fillboundaries",*) "sd_num<sd_fmnum. Sizes of SD arrays are not consistent."
          call PRC_abort
       endif
    end if

    ! Count the number of invalid SDs and make a list
    cnt = 0
    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) then
          cnt = cnt + 1
          ilist(cnt) = n
       end if
    end do
    available_sdnum = cnt
    
    ! Fill the boundaries with new SDs
    zmin = sdm_zlower + minzph
    zmax = sdm_zupper
    corner_ws = .true. ! if true, west south corner is not taken into account yet
    corner_wn = .true.
    corner_es = .true.
    corner_en = .true.
    dtmp_wbc = 0.0_RP
    dtmp_ebc = 0.0_RP
    dtmp_nbc = 0.0_RP
    dtmp_sbc = 0.0_RP
    cnt = 0

    !! west boundary
    if( (sdm_wbc==2) .and. (.NOT. PRC_HAS_W ) ) then
       xmin = GRID_FX(IS-2)
       xmax = GRID_FX(IS-1)
       if( (.NOT. PRC_HAS_S ) .and. corner_ws ) then ! west south corner
          ymin = GRID_FY(JS-2)
          corner_ws = .false.
       else
          ymin = GRID_FY(JS-1)
       end if
       if( (.NOT. PRC_HAS_N ) .and. corner_wn ) then ! west north corner
          ymax = GRID_FY(JE+1)
          corner_wn = .false.
       else
          ymax = GRID_FY(JE)
       end if

       dtmp_wbc = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)*sdm_inisdnc/sdm_sdnmlvol

       dtmp = dtmp_wbc
       call sdm_sdcreate(KA,KS,KE, IA,IS,IE, JA,JS,JE, &
                         DENS, RHOT, QTRC,                   &
                         sd_rng,                &
                         xmin,xmax,ymin,ymax,zmin,zmax,sdm_dtcmph,sdm_calvar,    &
                         sdm_rdnc,sdm_aslset,   &
                         sdm_inisdnc,sdm_sdnmlvol,              &
                         nint(dtmp),sd_numasl,sd_fmn,sd_fmliqice,sd_fmasl, sd_fmx, sd_fmy,        &
                         sd_fmz, sd_fmr,                   &
                         sd_fmri,sd_fmrj,sd_fmrk,sd_fmi,&
                         sd_dtmp1,sd_dtmp2)
       
       if( (available_sdnum-cnt) < nint(dtmp) ) then
          LOG_ERROR("sdm_fillboundaries",*) 'Number of SDs not sufficient for filling the boundaries'
          LOG_ERROR_CONT(*) 'Num of SDs needed:',nint(dtmp),', Remaining SDs: ', (available_sdnum-cnt)
          call PRC_abort
       end if

       ! Add new super-droplets to the array
       do n=1,nint(dtmp)
          
          id_invd = ilist(n+cnt)

          sd_liqice(id_invd)  = sd_fmliqice(n)

          sd_n(id_invd)  = sd_fmn(n)

          sd_x(id_invd)  = sd_fmx(n)
          sd_y(id_invd)  = sd_fmy(n)
          sd_z(id_invd)  = sd_fmz(n)
          sd_rk(id_invd) = sd_fmrk(n)
          
!!$          sd_u(id_invd)  = 0.0_RP
!!$          sd_v(id_invd)  = 0.0_RP
!!$          sd_vz(id_invd) = sd_fmvz(n)

          sd_r(id_invd)  = sd_fmr(n)

       end do

       do k=1,sd_numasl
          do n=1,nint(dtmp)

            id_invd = ilist(n+cnt)

            sd_asl(id_invd,k) = sd_fmasl(n,k)

         end do
      end do

      if( sdm_cold ) then
         do n=1,nint(dtmp)

            id_invd = ilist(n+cnt)

            sdi%re(id_invd)  = sd_fmi%re(n)
            sdi%rp(id_invd)  = sd_fmi%rp(n)
            sdi%rho(id_invd) = sd_fmi%rho(n)
            sdi%tf(id_invd)  = sd_fmi%tf(n)
            sdi%mrime(id_invd) = sd_fmi%mrime(n)
            sdi%nmono(id_invd) = sd_fmi%nmono(n)
            
         end do
      end if

      cnt = cnt + nint(dtmp)
    end if

    !! east boundary
    if( (sdm_ebc==2) .and. (.NOT. PRC_HAS_E ) ) then
       xmin = GRID_FX(IE)
       xmax = GRID_FX(IE+1)
       if( (.NOT. PRC_HAS_S ) .and. corner_es ) then ! east south corner
          ymin = GRID_FY(JS-2)
          corner_es = .false.
       else
          ymin = GRID_FY(JS-1)
       end if
       if( (.NOT. PRC_HAS_N ) .and. corner_en ) then ! east north corner
          ymax = GRID_FY(JE+1)
          corner_en = .false.
       else
          ymax = GRID_FY(JE)
       end if

       dtmp_ebc = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)*sdm_inisdnc/sdm_sdnmlvol

       dtmp = dtmp_ebc
       call sdm_sdcreate(KA,KS,KE, IA,IS,IE, JA,JS,JE, &
                         DENS, RHOT, QTRC,                   &
                         sd_rng,                &
                         xmin,xmax,ymin,ymax,zmin,zmax,sdm_dtcmph,sdm_calvar,    &
                         sdm_rdnc,sdm_aslset,   &
                         sdm_inisdnc,sdm_sdnmlvol,              &
                         nint(dtmp),sd_numasl,sd_fmn,sd_fmliqice,sd_fmasl, sd_fmx, sd_fmy,        &
                         sd_fmz, sd_fmr,                   &
                         sd_fmri,sd_fmrj,sd_fmrk,sd_fmi,&
                         sd_dtmp1,sd_dtmp2)

       if( (available_sdnum-cnt) < nint(dtmp) ) then
          LOG_ERROR("sdm_fillboundaries",*) 'Number of SDs not sufficinet for filling the boundaries'
          LOG_ERROR_CONT(*) 'Num of SDs needed:',nint(dtmp),', Remaining SDs: ', (available_sdnum-cnt)
          call PRC_abort
       end if

       ! Add new super-droplets to the array
       do n=1,nint(dtmp)
          
          id_invd = ilist(n+cnt)

          sd_liqice(id_invd)  = sd_fmliqice(n)

          sd_n(id_invd)  = sd_fmn(n)

          sd_x(id_invd)  = sd_fmx(n)
          sd_y(id_invd)  = sd_fmy(n)
          sd_z(id_invd)  = sd_fmz(n)
          sd_rk(id_invd) = sd_fmrk(n)
          
!!$          sd_u(id_invd)  = 0.0_RP
!!$          sd_v(id_invd)  = 0.0_RP
!!$          sd_vz(id_invd) = sd_fmvz(n)

          sd_r(id_invd)  = sd_fmr(n)

       end do

       do k=1,sd_numasl
          do n=1,nint(dtmp)

            id_invd = ilist(n+cnt)

            sd_asl(id_invd,k) = sd_fmasl(n,k)

         end do
      end do

      if( sdm_cold ) then
         do n=1,nint(dtmp)

            id_invd = ilist(n+cnt)

            sdi%re(id_invd)  = sd_fmi%re(n)
            sdi%rp(id_invd)  = sd_fmi%rp(n)
            sdi%rho(id_invd) = sd_fmi%rho(n)
            sdi%tf(id_invd)  = sd_fmi%tf(n)
            sdi%mrime(id_invd) = sd_fmi%mrime(n)
            sdi%nmono(id_invd) = sd_fmi%nmono(n)
            
         end do
      end if

      cnt = cnt + nint(dtmp)
    end if

    !! south boundary
    if( (sdm_sbc==2) .and. (.NOT. PRC_HAS_S ) ) then
       ymin = GRID_FY(JS-2)
       ymax = GRID_FY(JS-1)
       if( (.NOT. PRC_HAS_W ) .and. corner_ws ) then ! west south corner
          xmin = GRID_FX(IS-2)
          corner_ws = .false.
       else
          xmin = GRID_FX(IS-1)
       end if
       if( (.NOT. PRC_HAS_E ) .and. corner_es ) then ! east south corner
          xmax = GRID_FX(IE+1)
          corner_es = .false.
       else
          xmax = GRID_FX(IE)
       end if

       dtmp_sbc = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)*sdm_inisdnc/sdm_sdnmlvol

       dtmp = dtmp_sbc
       call sdm_sdcreate(KA,KS,KE, IA,IS,IE, JA,JS,JE, &
                         DENS, RHOT, QTRC,                   &
                         sd_rng,                &
                         xmin,xmax,ymin,ymax,zmin,zmax,sdm_dtcmph,sdm_calvar,    &
                         sdm_rdnc,sdm_aslset,   &
                         sdm_inisdnc,sdm_sdnmlvol,              &
                         nint(dtmp),sd_numasl,sd_fmn,sd_fmliqice,sd_fmasl, sd_fmx, sd_fmy,        &
                         sd_fmz, sd_fmr,                   &
                         sd_fmri,sd_fmrj,sd_fmrk,sd_fmi,&
                         sd_dtmp1,sd_dtmp2)

       if( (available_sdnum-cnt) < nint(dtmp) ) then
          LOG_ERROR("sdm_fillboundaries",*) 'Number of SDs not sufficinet for filling the boundaries'
          LOG_ERROR_CONT(*) 'Num of SDs needed:',nint(dtmp),', Remaining SDs: ', (available_sdnum-cnt)
          call PRC_abort
       end if

       ! Add new super-droplets to the array
       do n=1,nint(dtmp)
          
          id_invd = ilist(n+cnt)

          sd_liqice(id_invd)  = sd_fmliqice(n)

          sd_n(id_invd)  = sd_fmn(n)

          sd_x(id_invd)  = sd_fmx(n)
          sd_y(id_invd)  = sd_fmy(n)
          sd_z(id_invd)  = sd_fmz(n)
          sd_rk(id_invd) = sd_fmrk(n)
          
!!$          sd_u(id_invd)  = 0.0_RP
!!$          sd_v(id_invd)  = 0.0_RP
!!$          sd_vz(id_invd) = sd_fmvz(n)

          sd_r(id_invd)  = sd_fmr(n)

       end do

       do k=1,sd_numasl
          do n=1,nint(dtmp)

            id_invd = ilist(n+cnt)

            sd_asl(id_invd,k) = sd_fmasl(n,k)

         end do
      end do

      if( sdm_cold ) then
         do n=1,nint(dtmp)

            id_invd = ilist(n+cnt)

            sdi%re(id_invd)  = sd_fmi%re(n)
            sdi%rp(id_invd)  = sd_fmi%rp(n)
            sdi%rho(id_invd) = sd_fmi%rho(n)
            sdi%tf(id_invd)  = sd_fmi%tf(n)
            sdi%mrime(id_invd) = sd_fmi%mrime(n)
            sdi%nmono(id_invd) = sd_fmi%nmono(n)
            
         end do
      end if

      cnt = cnt + nint(dtmp)
    end if

    !! north boundary
    if( (sdm_nbc==2) .and. (.NOT. PRC_HAS_N ) ) then
       ymin = GRID_FY(JE)
       ymax = GRID_FY(JE+1)
       if( (.NOT. PRC_HAS_W ) .and. corner_wn ) then ! west north corner
          xmin = GRID_FX(IS-2)
          corner_wn = .false.
       else
          xmin = GRID_FX(IS-1)
       end if
       if( (.NOT. PRC_HAS_E ) .and. corner_en ) then ! east north corner
          xmax = GRID_FX(IE+1)
          corner_en = .false.
       else
          xmax = GRID_FX(IE)
       end if

       dtmp_nbc = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)*sdm_inisdnc/sdm_sdnmlvol

       dtmp = dtmp_nbc
       call sdm_sdcreate(KA,KS,KE, IA,IS,IE, JA,JS,JE, &
                         DENS, RHOT, QTRC,                   &
                         sd_rng,                &
                         xmin,xmax,ymin,ymax,zmin,zmax,sdm_dtcmph,sdm_calvar,    &
                         sdm_rdnc,sdm_aslset,   &
                         sdm_inisdnc,sdm_sdnmlvol,              &
                         nint(dtmp),sd_numasl,sd_fmn,sd_fmliqice,sd_fmasl, sd_fmx, sd_fmy,        &
                         sd_fmz, sd_fmr,                   &
                         sd_fmri,sd_fmrj,sd_fmrk,sd_fmi,&
                         sd_dtmp1,sd_dtmp2)

       if( (available_sdnum-cnt) < nint(dtmp) ) then
          LOG_ERROR("sdm_fillboundaries",*) 'Number of SDs not sufficinet for filling the boundaries'
          LOG_ERROR_CONT(*) 'Num of SDs needed:',nint(dtmp),', Remaining SDs: ', (available_sdnum-cnt)
          call PRC_abort
       end if

       ! Add new super-droplets to the array
       do n=1,nint(dtmp)
          
          id_invd = ilist(n+cnt)

          sd_liqice(id_invd)  = sd_fmliqice(n)

          sd_n(id_invd)  = sd_fmn(n)

          sd_x(id_invd)  = sd_fmx(n)
          sd_y(id_invd)  = sd_fmy(n)
          sd_z(id_invd)  = sd_fmz(n)
          sd_rk(id_invd) = sd_fmrk(n)
          
!!$          sd_u(id_invd)  = 0.0_RP
!!$          sd_v(id_invd)  = 0.0_RP
!!$          sd_vz(id_invd) = sd_fmvz(n)

          sd_r(id_invd)  = sd_fmr(n)

       end do

       do k=1,sd_numasl
          do n=1,nint(dtmp)

            id_invd = ilist(n+cnt)

            sd_asl(id_invd,k) = sd_fmasl(n,k)

         end do
      end do

      if( sdm_cold ) then
         do n=1,nint(dtmp)

            id_invd = ilist(n+cnt)

            sdi%re(id_invd)  = sd_fmi%re(n)
            sdi%rp(id_invd)  = sd_fmi%rp(n)
            sdi%rho(id_invd) = sd_fmi%rho(n)
            sdi%tf(id_invd)  = sd_fmi%tf(n)
            sdi%mrime(id_invd) = sd_fmi%mrime(n)
            sdi%nmono(id_invd) = sd_fmi%nmono(n)

         end do
      end if

      cnt = cnt + nint(dtmp)
    end if
   
    return

  end subroutine sdm_fillboundaries
  !-----------------------------------------------------------------------------
  subroutine sdm_sdcreate(    &
                         KA,KS,KE, IA,IS,IE, JA,JS,JE, &
                         DENS, RHOT, QTRC,                   &
                         sd_rng,                &
                         xmin,xmax,ymin,ymax,zmin,zmax,sdm_dtcmph,sdm_calvar,   &
                         sdm_rdnc,sdm_aslset,                    &
                         sdm_inisdnc,sdm_sdnmlvol,               &
                         sd_num,sd_numasl,sd_n,sd_liqice,sd_asl, &
                         sd_x, sd_y, sd_z, sd_r,                 &
                         sd_ri,sd_rj,sd_rk,sdi,          &
                         sd_dtmp1,sd_dtmp2                    )
  !***********************************************************************
  ! Input variables
    use scale_precision
    use scale_const, only: &
         Rvap => CONST_Rvap
    use scale_prc, only: &
         PRC_abort
    use scale_io, only: &
         IO_L,IO_FID_LOG,H_LONG
    use m_sdm_common, only: &
         n1_amsul,n1_amsul_derksn,n2_amsul,n2_amsul_derksn, &
         rb1_amsul,rb1_amsul_derksn,rb2_amsul,rb2_amsul_derksn, &
         sgm1_amsul,sgm1_amsul_derksn,sgm2_amsul,sgm2_amsul_derksn, &
         rmax_amsul,rmax_amsul_derksn,rmin_amsul,rmin_amsul_derksn, &
         n3_nacl,rb3_nacl,sgm3_nacl,rmax_nacl,rmin_nacl, &
         n1_amsul_zanten,n2_amsul_zanten,rb1_amsul_zanten,rb2_amsul_zanten, &
         sgm1_amsul_zanten,sgm2_amsul_zanten,rmax_amsul_zanten,rmin_amsul_zanten, &
         n3_nacl_derksn,rb3_nacl_derksn,sgm3_nacl_derksn,rmax_nacl_derksn,rmin_nacl_derksn, &
         n1_amsul_dycoms,n2_amsul_dycoms,rb1_amsul_dycoms,rb2_amsul_dycoms, &
         sgm1_amsul_dycoms,sgm2_amsul_dycoms,rmax_amsul_dycoms,rmin_amsul_dycoms, &
         rho_amsul,rho_nacl,sdm_aslmw,sdm_aslion, &
         n_mdust,tfmax,tfmin,mdust_dia,&
         sdnumratio_soluble,INIAsd_ratio,a0_N12,a1_N12,a2_N12, &
         sdm_cold,sdm_fctr2multi, &
         STAT_LIQ,F_THRD,ONE_PI,O_THRD,INVALID, &
         i2,sdicedef,I_QV_sdm,QA_MP_sdm, &
         sdm_wbc,sdm_ebc,sdm_sbc,sdm_nbc
    use m_sdm_coordtrans, only: &
         sdm_z2rk
    use m_sdm_fluidconv, only: &
         sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw
    use m_sdm_sd2fluid, only: &
         sdm_sd2qcqr,sdm_sd2qiqsqg
    use m_sdm_condensation_water, only: &
         sdm_condevp
    use m_sdm_meltfreeze, only: &
         sdm_meltfreeze
    use rng_uniform_mt, only: &
         c_rng_uniform_mt
    integer, intent(in) :: KA,KS,KE
    integer, intent(in) :: IA,IS,IE
    integer, intent(in) :: JA,JS,JE
    real(RP), intent(in) :: DENS(KA,IA,JA) ! Density     [kg/m3]
    real(RP), intent(in) :: RHOT(KA,IA,JA) ! DENS * POTT [K*kg/m3]
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_MP_sdm) ! ratio of mass of tracer to total mass[kg/kg]
    type(c_rng_uniform_mt), intent(inout) :: sd_rng ! random number generator
    real(RP),intent(in) :: sdm_dtcmph(5)  ! Time interval of cloud micro physics
    logical, intent(in) :: sdm_calvar(5)! Flag for cond./coll/move calculation
    real(RP),intent(in) :: sdm_rdnc     ! Number concentration of real droplets
    real(RP),intent(in) :: sdm_sdnmlvol ! Normal volume for number concentration of super droplets
    integer, intent(in) :: sdm_aslset   ! Option for aerosol species
    real(RP),intent(in) :: sdm_inisdnc  ! Initial number of super droplets per sdm_sdnmlvol
    real(RP),intent(in) :: xmin   ! Lower limit of new SD x-positions
    real(RP),intent(in) :: xmax   ! Upper limit of new SD x-positions
    real(RP),intent(in) :: ymin   ! Lower limit of new SD y-positions
    real(RP),intent(in) :: ymax   ! Upper limit of new SD y-positions
    real(RP),intent(in) :: zmin   ! Lower limit of new SD z-positions
    real(RP),intent(in) :: zmax   ! Upper limit of new SD z-positions

    integer, intent(in) :: sd_num  ! number of super-droplets
    integer, intent(in) :: sd_numasl ! Number of kind of chemical material contained as water-soluble aerosol in super droplets
    integer(DP), intent(out) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer(i2), intent(out) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP),intent(out) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets [kg]
    real(RP),intent(out) :: sd_x(1:sd_num)               ! x-coordinate of super-droplets [m]
    real(RP),intent(out) :: sd_y(1:sd_num)               ! y-coordinate of super-droplets [m]
    real(RP),intent(out) :: sd_z(1:sd_num)               ! z-coordinate of super-droplets [m]
    real(RP),intent(out) :: sd_r(1:sd_num)               ! equivalent radius of super-droplets [m]
    real(RP),intent(out) :: sd_ri(1:sd_num)              ! index[i/real] of super-droplets
    real(RP),intent(out) :: sd_rj(1:sd_num)              ! index[j/real] of super-droplets
    real(RP),intent(out) :: sd_rk(1:sd_num)              ! index[k/real] of super-droplets
    type(sdicedef), intent(inout) :: sdi                 ! ice phase super-droplets

    real(RP),intent(out) :: sd_dtmp1(1:sd_num) ! temp 
    real(RP),intent(out) :: sd_dtmp2(1:sd_num) ! temp

    ! Work variables
    real(RP) :: n0                            ! number of real droplets per unit volume and per aerosol radius
    real(RP) :: dry_r                         ! aerosol radius
    real(RP) :: delta1, delta2, sdn_tmp       ! temporary
    logical  :: lexced                        ! temporary
    integer :: i, k, n                        ! index

    real(RP) :: pres_scale(KA,IA,JA)  ! Pressure
    real(RP) :: t_scale(KA,IA,JA)    ! Temperature

    real(RP) :: sdm_dtevl  ! time step of {condensation/evaporation} process
    real(RP) :: sdm_dtcol  ! time step of {stochastic coalescence} process
    real(RP) :: sdm_dtadv  ! time step of {motion of super-droplets} process
    real(RP) :: sdm_dtmlt  ! time step of {melt/freeze of super-droplets} process
    real(RP) :: sdm_dtsbl  ! time step of {sublimation/deposition of super-droplets} process

    !
    real(RP) :: area, INAS_max, prob_INIA, INAS_tf, probdens_tf

    integer :: IS_bnd, IE_bnd, JS_bnd, JE_bnd

    !---------------------------------------------------------------------

    ! Initialize and rename variables
    sdm_dtevl = real( sdm_dtcmph(1),kind=RP )  !! condensation/evaporation
    sdm_dtcol = real( sdm_dtcmph(2),kind=RP )  !! stochastic coalescence
    sdm_dtadv = real( sdm_dtcmph(3),kind=RP )  !! motion of super-droplets
    sdm_dtmlt = real( sdm_dtcmph(4),kind=RP )  !! melting/freezing
    sdm_dtsbl = real( sdm_dtcmph(5),kind=RP )  !! sublimation/depsition

    if( .not. sdm_calvar(1) .and. &
         .not. sdm_calvar(2) .and. &
         .not. sdm_calvar(3) .and. &
         .not. sdm_calvar(4) .and. &
         .not. sdm_calvar(5)        ) return

    ! Initialized super-droplets.
    !### Get parameter for 3mode log-nomiral distribution ###!
    
    if( mod(sdm_aslset,10)==1 .or. mod(sdm_aslset,10)==3 ) then
       
       !! Derksen(2009)
!!$       if( IO_L ) then
!!$          write(IO_FID_LOG,*)  "    number concentration of (NH4)2SO4    : Derksen"
!!$       endif
       n1_amsul   = n1_amsul_derksn
       n2_amsul   = n2_amsul_derksn
       rb1_amsul  = rb1_amsul_derksn
       rb2_amsul  = rb2_amsul_derksn
       sgm1_amsul = sgm1_amsul_derksn
       sgm2_amsul = sgm2_amsul_derksn
       rmax_amsul = rmax_amsul_derksn
       rmin_amsul = rmin_amsul_derksn

       if( mod(sdm_aslset,10)==1 ) then
          n3_nacl    = 1.E-10_RP     !! temporary initialize
          rb3_nacl   = 1.E-10_RP
          sgm3_nacl  = 1.E-10_RP
          rmax_nacl  = 1.E-10_RP
          rmin_nacl  = 1.E-10_RP
       end if
       
    else if( mod(sdm_aslset,10)==-1 .or. mod(sdm_aslset,10)==-3 ) then
       
       !! vanZanten(2010)
!!$       if( IO_L ) then
!!$          write(IO_FID_LOG,*) "    number concentration of (NH4)2SO4    : vanZanten"
!!$       endif
       
       n1_amsul   = n1_amsul_zanten
       n2_amsul   = n2_amsul_zanten
       rb1_amsul  = rb1_amsul_zanten
       rb2_amsul  = rb2_amsul_zanten
       sgm1_amsul = sgm1_amsul_zanten
       sgm2_amsul = sgm2_amsul_zanten
       rmax_amsul = rmax_amsul_zanten
       rmin_amsul = rmin_amsul_zanten

       if( mod(sdm_aslset,10)==-1 ) then
          n3_nacl    = 1.E-10_RP     !! temporary initialize
          rb3_nacl   = 1.E-10_RP
          sgm3_nacl  = 1.E-10_RP
          rmax_nacl  = 1.E-10_RP
          rmin_nacl  = 1.E-10_RP
       end if
       
    end if

    if( mod(sdm_aslset,10)==2 .or. abs(mod(sdm_aslset,10))==3 ) then

       !! Derksen(2009)
       if( IO_L ) then
          write(IO_FID_LOG,*) "number concentration of NaCl aerosol : Derksen"
       endif
       if( mod(sdm_aslset,10)==2 ) then
          n1_amsul   = 1.E-10_RP     !! temporary initialize
          n2_amsul   = 1.E-10_RP
          rb1_amsul  = 1.E-10_RP
          rb2_amsul  = 1.E-10_RP
          sgm1_amsul = 1.E-10_RP
          sgm2_amsul = 1.E-10_RP
          rmax_amsul = 1.E-10_RP
          rmin_amsul = 1.E-10_RP
       end if
       n3_nacl   = n3_nacl_derksn
       rb3_nacl  = rb3_nacl_derksn
       sgm3_nacl = sgm3_nacl_derksn
       rmax_nacl = rmax_nacl_derksn
       rmin_nacl = rmin_nacl_derksn
       
    end if

    if( mod(sdm_aslset,10)==5 ) then
       
       !! Ackerman et al. (2009)
!!$       if( IO_L ) then
!!$          write(IO_FID_LOG,*) "    number concentration of (NH4)2SO4    : Ackerman"
!!$       endif

       n1_amsul   = n1_amsul_dycoms
       n2_amsul   = n2_amsul_dycoms
       rb1_amsul  = rb1_amsul_dycoms
       rb2_amsul  = rb2_amsul_dycoms
       sgm1_amsul = sgm1_amsul_dycoms
       sgm2_amsul = sgm2_amsul_dycoms
       rmax_amsul = rmax_amsul_dycoms
       rmin_amsul = rmin_amsul_dycoms

       n3_nacl    = 1.E-10_RP     !! temporary initialize
       rb3_nacl   = 1.E-10_RP
       sgm3_nacl  = 1.E-10_RP
       rmax_nacl  = 1.E-10_RP
       rmin_nacl  = 1.E-10_RP
       
    end if

    !### Get random number ###!
    call gen_rand_array( sd_rng, sd_x(1:sd_num) )
    call gen_rand_array( sd_rng, sd_y(1:sd_num) )
    call gen_rand_array( sd_rng, sd_z(1:sd_num) )
    call gen_rand_array( sd_rng, sd_r(1:sd_num) )

    lexced = .false.     !! check for memory size of int*8

    do k=1,sd_numasl
       call gen_rand_array( sd_rng, sd_dtmp1(1:sd_num) )
       do n=1,sd_num
          sd_asl(n,k) = sd_dtmp1(n)
       end do
    end do
   
    ! Initialized all super-droplets as water droplet.
    !### status(liquid/ice) of super-droplets ###!
    sd_liqice(1:sd_num) = STAT_LIQ
      
    !### Aerosol mass, muliplicity ###!
    !$omp parallel
    do k=1,sd_numasl
       !$omp do private(i,delta1,delta2,dry_r,n0,sdn_tmp) reduction(.or.:lexced)
       do n=1,sd_num
          i = mod(n-1,sd_numasl) + 1    !! select aerosol index
          if( abs(sdm_aslset)==12 ) then
             i = 2                         !! 1 : (NH4)2SO4, 2: NaCl
          end if
          if( k==i ) then
             !! match aerosol type
             if( abs(mod(sdm_aslset,10))==1 .or.                      &
                  ( k==1 .and. abs(mod(sdm_aslset,10))==3 ) .or.     &
                  ( mod(sdm_aslset,10)==5 ) ) then

                !### (NH4)2SO4 [g] ###!
                delta1 = log(rmax_amsul) - log(rmin_amsul)
                dry_r  = exp( log(rmin_amsul)+delta1*sd_asl(n,k) )
                sd_asl(n,k) = F_THRD * ONE_PI                 &
                     * (dry_r*dry_r*dry_r) * rho_amsul
                !! n0(log(dry_r)) [m-3]
                !! 2-mode log-noraml distribution for ammonium sulfate
                delta1 = log(dry_r) - log(rb1_amsul)
                delta1 = -(delta1*delta1)                             &
                     /(2.0_RP*log(sgm1_amsul)*log(sgm1_amsul))
                delta2 = log(dry_r) - log(rb2_amsul)
                delta2 = -(delta2*delta2)                             &
                     /(2.0_RP*log(sgm2_amsul)*log(sgm2_amsul))
                n0 = (n1_amsul*exp(delta1))                           &
                     /(sqrt(2.0_RP*ONE_PI)*log(sgm1_amsul))  &
                     + (n2_amsul*exp(delta2))                           &
                     /(sqrt(2.0_RP*ONE_PI)*log(sgm2_amsul))
                !! number per unit volume and per aerosol species
                delta1 = real(sdm_inisdnc,kind=RP)                    &
                     /real(sdm_sdnmlvol,kind=RP)
                delta1 = delta1/real(sd_numasl,kind=RP)
                !! continuous uniform distribution
                delta2 = 1.0_RP/(log(rmax_amsul)-log(rmin_amsul))
                !! muliplicity
                sdn_tmp = n0/(delta1*delta2)
             else if( abs(mod(sdm_aslset,10))==2 .or.                 &
                  ( k==2 .and. abs(mod(sdm_aslset,10))==3 ) ) then
                !### NaCl(seasalt,[g]) ###!
                delta1 = log(rmax_nacl) - log(rmin_nacl)
                dry_r  = exp( log(rmin_nacl) + delta1*sd_asl(n,k) )
                sd_asl(n,k) = F_THRD * ONE_PI                 &
                     * (dry_r*dry_r*dry_r) * rho_nacl
                !! n0(log(dry_r)) [m-3]
                !! 1-mode log-noraml distribution for seasalt
                delta1 = log(dry_r) - log(rb3_nacl)
                delta1 = -(delta1*delta1)                             &
                     /(2.0_RP*log(sgm3_nacl)*log(sgm3_nacl))
                n0 = (n3_nacl*exp(delta1))                            &
                       /(sqrt(2.0_RP*ONE_PI)*log(sgm3_nacl))
                !! number per unit volume and per aerosol species
                delta1 = real(sdm_inisdnc,kind=RP)                    &
                     /real(sdm_sdnmlvol,kind=RP)
                delta1 = delta1/real(sd_numasl,kind=RP)
                !! continuous uniform distribution
                delta2 = 1.0_RP/(log(rmax_nacl)-log(rmin_nacl))
                !! muliplicity
                sdn_tmp = n0/(delta1*delta2)
             else
                !### Other ###!
                sd_asl(n,k) = 1.0E-18_RP                           &
                     + 1.0E-14_RP                           &
                     * (log(1.0_RP/(1.0_RP-sd_asl(n,k))))
                sdn_tmp = real(sdm_rdnc,kind=RP)                      &
                     * real(sdm_sdnmlvol,kind=RP)                  &
                     / real(sdm_inisdnc,kind=RP)
             end if

             !! Multiply the initial number density of soluble aerosol particles by sdm_fctr2multi
             !! This only for soluble particles, not for insoluble particles
             sdn_tmp = sdn_tmp * sdm_fctr2multi

             !! check muliplicity
             if( sdn_tmp<(2.0_RP**63.0_RP) ) then
                sd_n(n) = nint( sdn_tmp, kind=DP )
             else
                lexced = .true.
             end if
          else
             sd_asl(n,k) = 0.0_RP
          end if
       end do
    end do
    !$omp end parallel

    ! Initialized SDs for ice phase (freezing temperature and multiplicity)
    if( sdm_cold ) then

       !###### initialization
       do n = 1, sd_num
          sdi%rho(n)   = 0.0_RP
          sdi%re(n)    = 0.0_RP
          sdi%rp(n)    = 0.0_RP
          sdi%mrime(n) = 0.0_RP
          sdi%nmono(n) = 0
       end do

       !###### soluble aerosol ######!
       do n=1,nint(sd_num*sdnumratio_soluble)
          sdi%tf(n) = -38.d0  ! homogeneous freezing limit [degC] 
          sd_n(n) =  nint( real(sd_n(n),kind=RP)/sdnumratio_soluble, kind=DP ) ! adjust multiplicity
       end do

       !###### insluble+soluble aerosol (internally mixed) ######!
       if(nint(sd_num*sdnumratio_soluble)+1 < sd_num) then

          !###### set parameters ######!
          ! set aerosol distribution parameters of soluble component
          n1_amsul   = n_mdust * (n1_amsul_zanten/(n1_amsul_zanten+n2_amsul_zanten))
          n2_amsul   = n_mdust * (n2_amsul_zanten/(n1_amsul_zanten+n2_amsul_zanten))
          rb1_amsul  = rb1_amsul_zanten
          rb2_amsul  = rb2_amsul_zanten
          sgm1_amsul = sgm1_amsul_zanten
          sgm2_amsul = sgm2_amsul_zanten
          rmax_amsul = rmax_amsul_zanten
          rmin_amsul = rmin_amsul_zanten
          
          !###### mass of soluble component ######!
          ! reset all chemical components
          do n=nint(sd_num*sdnumratio_soluble)+1,sd_num
             do k=1,sd_numasl
                sd_asl(n,k) = 0.0_RP
             end do
          end do
          !### (NH4)2SO4 [g] ###!
          k = 1
          ! prepare random numbers
          call gen_rand_array( sd_rng, sd_dtmp1(1:sd_num) )
          do n=nint(sd_num*sdnumratio_soluble)+1,sd_num
             delta1 = log(rmax_amsul) - log(rmin_amsul)
             dry_r  = exp( log(rmin_amsul)+delta1*sd_dtmp1(n) )
             sd_asl(n,k) = F_THRD * ONE_PI                 &
                  * (dry_r*dry_r*dry_r) * rho_amsul
          end do

          !###### mass of insoluble component ######!
          !### Assume a monodisperse distribution of mineral dust with a diameter of mdust_dia ###!
          ! mdust_dia [m] is defined in sdm_common.f90

          !###### freezing temperature ######!
          ! prepare random numbers
          call gen_rand_array( sd_rng, sd_dtmp1(1:sd_num) )
          do n=nint(sd_num*sdnumratio_soluble)+1,sd_num
             delta1 = tfmax-tfmin
             sdi%tf(n) = tfmin+delta1*sd_dtmp1(n)
          end do
          
          !###### multiplicity ######!
          ! prepare random numbers
          call gen_rand_array( sd_rng, sd_dtmp1(1:sd_num) )
          call gen_rand_array( sd_rng, sd_dtmp2(1:sd_num) )
          k=1 ! internally mixed with (NH4)2SO4
          do n=nint(sd_num*sdnumratio_soluble)+1,sd_num
             !### contribution from soluble part ###!
             !! n0(log(dry_r)) [m-3]
             !! 2-mode log-noraml distribution for ammonium sulfate
             dry_r = (sd_asl(n,k)/F_THRD/ONE_PI/rho_amsul)**O_THRD
             delta1 = log(dry_r) - log(rb1_amsul)
             delta1 = -(delta1*delta1)                             &
                  /(2.0_RP*log(sgm1_amsul)*log(sgm1_amsul))
             delta2 = log(dry_r) - log(rb2_amsul)
             delta2 = -(delta2*delta2)                             &
                  /(2.0_RP*log(sgm2_amsul)*log(sgm2_amsul))
             n0 = (n1_amsul*exp(delta1))                           &
                  /(sqrt(2.0_RP*ONE_PI)*log(sgm1_amsul))  &
                  + (n2_amsul*exp(delta2))                           &
                  /(sqrt(2.0_RP*ONE_PI)*log(sgm2_amsul))
             !! number of SDs per unit volume
             delta1 = (1.0_RP-sdnumratio_soluble)*real(sdm_inisdnc,kind=RP)   &
                  /real(sdm_sdnmlvol,kind=RP)
             !! continuous uniform distribution
             delta2 = 1.0_RP/(log(rmax_amsul)-log(rmin_amsul))
             !! muliplicity
             sdn_tmp = n0/(delta1*delta2)
             
             !### contribution from insluble part ###!
             area = ONE_PI*mdust_dia**2
             if(sd_dtmp1(n) < INIAsd_ratio) then
                !! IN inactive mineral dust
                INAS_max = a0_N12*exp(-a1_N12*tfmin+a2_N12) ! INAS of Niemand et al. (2012) 
                prob_INIA = exp(-area*INAS_max)
                sdn_tmp = sdn_tmp * prob_INIA / INIAsd_ratio
                sdi%tf(n) = -38.0_RP  ! homogeneous freezing limit
             else
                !! IN active mineral dust
                INAS_tf = a0_N12*exp(-a1_N12*sdi%tf(n)+a2_N12)  ! INAS of Niemand et al. (2012) 
                probdens_tf = area*a1_N12*INAS_tf*exp(-area*INAS_tf)
                delta1 = tfmax-tfmin
                sdn_tmp = sdn_tmp * probdens_tf * delta1 / (1.0_RP - INIAsd_ratio)
             end if

             !! check multiplicity
             if( sdn_tmp<(2.0_RP**63.0_RP) ) then
                if(sdn_tmp>1.0_RP)then
                   sd_n(n) = nint( sdn_tmp, kind=DP )
                else ! to sample rare event
                   sd_n(n) = 0 
                   if( sd_dtmp2(n) < sdn_tmp )then
                      sd_n(n) = 1
                   end if
                end if
             else
                lexced = .true.
             end if
             
          end do
          
       end if

    end if

    if( lexced ) then
       LOG_ERROR("sdm_sdcreate",*) "overflow of multiplicity"
       call PRC_abort
    end if

    !### position of super-droplets in horizontal ###!
    do n=1,sd_num
       sd_x(n) = (xmax-xmin)*sd_x(n) + xmin
       sd_y(n) = (ymax-ymin)*sd_y(n) + ymin
    end do

    !### position of super-droplets in vertical ###!
    !! valid super-droplets
    do n=1,sd_num
       if( sd_n(n)>0 ) then
          sd_z(n) = (zmax-zmin)*sd_z(n) + zmin
          sd_rk(n) = 0.0 ! initialized tentatively
       else
          sd_z(n) = INVALID     !!! check muliplicity
          sd_rk(n) = INVALID
       end if
    end do
      
    !### index[k/real] of super-droplets               ###!
    !### modify position[z] of invalid super-droplets  ###!
    call sdm_z2rk( &
       KS,KE,IS,IE,JS,JE, &
       sd_num,sd_x,sd_y,sd_z,sd_ri,sd_rj,sd_rk )
    
    !### initial equivalent radius                     ###!
    do n=1,sd_num
       sd_r(n) = 1.0E-15_RP

!ORG     sd_r(n) = 1.0e-5 * ( log(1.0/(1.0-sd_r(n))) )**O_THRD
!ORG     sd_r(n) = 1.0d-8

! temporary for test
!         sd_r(n) = 3.0E-3_RP*sd_r(n)
!         sd_r(n) = exp((log(3.0E-3_RP)-log(1.0E-7_RP))*sd_r(n)+log(1.0E-7_RP))
         
    end do

    !### ( at only condensation/evaporation process ) ###!
    if( sdm_calvar(1) ) then

       IS_bnd = IS
       IE_bnd = IE
       JS_bnd = JS
       JE_bnd = JE
       if ( sdm_wbc .eq. 2 ) IS_bnd = max(IS-1,1)
       if ( sdm_ebc .eq. 2 ) IE_bnd = min(IE+1,IA)
       if ( sdm_sbc .eq. 2 ) JS_bnd = max(JS-1,1)
       if ( sdm_nbc .eq. 2 ) JE_bnd = min(JE+1,JA)

       call sdm_rhot_qtrc2p_t( &
            KA, KS, KE, IA, IS_bnd, IE_bnd, JA, JS_bnd, JE_bnd, &
            RHOT,QTRC,DENS,pres_scale,t_scale)

       call sdm_condevp(    &
            KA, KS, KE, IA, IS_bnd, IE_bnd, JA, JS_bnd, JE_bnd, &
            sdm_aslset,                                   &
            sdm_aslmw,sdm_aslion,sdm_dtevl,               &
            pres_scale,t_scale,QTRC(:,:,:,I_QV_sdm),          &
            sd_num,sd_numasl,sd_liqice,sd_x,sd_y,       &
            sd_r,sd_asl,sd_ri,sd_rj,sd_rk)

    end if

    !### ( melting/freezing process ) ###!
    if( sdm_cold .and. sdm_calvar(4) ) then

       call sdm_rhot_qtrc2p_t( &
            KA, KS, KE, IA, IS_bnd, IE_bnd, JA, JS_bnd, JE_bnd, &
            RHOT,QTRC,DENS,pres_scale,t_scale)

       call sdm_meltfreeze(                              &
            KA,KS,KE, IA,IS_bnd,IE_bnd, JA,JS_bnd,JE_bnd, &
            t_scale,pres_scale,QTRC(:,:,:,I_QV_sdm),         &
            sd_num,sd_liqice,sd_x,sd_y,      &
            sd_r,sd_ri,sd_rj,sd_rk,sdi )

    end if

      !### diagnose terminal velocity (no need evaluate them here) ###!
!!$      call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
!!$      call sdm_getvz_liq(pres_scale,DENS,t_scale,                    &
!!$                     sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_r, &
!!$                     sdvz_s2c,sd_itmp1,sd_itmp2,sd_itmp3,'no_interpolation')
      
    return

  end subroutine sdm_sdcreate

  subroutine gen_rand_array_DP(rng, ary)
    use scale_precision
    use rng_uniform_mt, only: &
         c_rng_uniform_mt, &
         gen_rand_array_mt => gen_rand_array

    type(c_rng_uniform_mt), intent(inout) :: rng
    real(DP), intent(inout) :: ary(:)

    call gen_rand_array_mt(rng, ary)

    return
  end subroutine gen_rand_array_DP

  subroutine gen_rand_array_SP(rng, ary)
    use scale_precision
    use rng_uniform_mt, only: &
         c_rng_uniform_mt, &
         gen_rand_array_mt => gen_rand_array

    type(c_rng_uniform_mt), intent(inout) :: rng
    real(SP), intent(inout) :: ary(:)
    real(DP) :: ary_DP(size(ary))

    call gen_rand_array_mt(rng, ary_DP)
    ary(:) = ary_DP(:)

    return
  end subroutine gen_rand_array_SP

end module m_sdm_ptlformation
