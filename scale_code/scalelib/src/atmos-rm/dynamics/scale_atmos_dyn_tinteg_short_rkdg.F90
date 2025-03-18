!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics Temporal integration
!!
!! @par Description
!!          Temporal scheme selecter for RK short time step in DG dynamica core
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_short_rkdg
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_timeint_rk, only: &
    TimeInt_RK

  use scale_fem_element_base, only: &
    ElementBase3D
  use scale_fem_localmesh_3d, only: &
    LocalMesh3D
  use scale_atmos_dyn_dgm_operator, only: &
    mesh3D, elem3D, &
    dyn_bnd, &
    Dx, Dy, Dz, Lift, &
    MODALFILTER_FLAG, modal_filter_3d, &
    SPONGELAYER_FLAG
  use scale_atmos_dyn_dgm_vars, only: &
    PRGVAR_NUM, &
    DENS_VID => PRGVAR_DDENS_ID, THERM_VID => PRGVAR_THERM_ID,&
    MOMX_VID => PRGVAR_MOMX_ID, MOMY_VID => PRGVAR_MOMY_ID,   &
    MOMZ_VID => PRGVAR_MOMZ_ID, &
    MASSFLX_X_ID, MASSFLX_Y_ID, MASSFLX_Z_ID, &
    TRCVARS3D_DENS0_ID, TRCVARS3D_DENS_ID
     
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_Tinteg_short_rkdg_setup
  public :: ATMOS_DYN_Tinteg_short_rkdg
  public :: ATMOS_DYN_Tinteg_short_rkdg_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
  type(TimeInt_RK), allocatable :: tint(:)
  logical :: hevi_flag


contains
  !-----------------------------------------------------------------------------
  !> Register
!OCL SERIAL
  subroutine ATMOS_DYN_Tinteg_short_rkdg_setup( &
       ATMOS_DYN_Tinteg_short_TYPE, ATMOS_DYN_Tstep_short_TYPE, &
       dt )

    use scale_precision
    use scale_prc, only: &
       PRC_abort           
    implicit none

    character(len=*), intent(in)  :: ATMOS_DYN_Tinteg_short_TYPE
    character(len=*), intent(in)  :: ATMOS_DYN_Tstep_short_TYPE
    real(RP), intent(in) :: dt

    integer :: lcdomID
    !---------------------------------------------------------------------------

    select case( ATMOS_DYN_Tinteg_short_TYPE )
    case( 'OFF', 'NONE' )
       ! do nothing
    case default
      allocate( tint(mesh3D%LOCAL_MESH_NUM) )
      do lcdomID=1, mesh3D%LOCAL_MESH_NUM
        call tint(lcdomID)%Init( ATMOS_DYN_Tinteg_short_TYPE, dt, PRGVAR_NUM, &
          2, (/ elem3D%Np, mesh3D%lcmesh_list(lcdomID)%NeA /) )
      enddo

      if ( tint(1)%imex_flag ) then
        hevi_flag = .true.
      else
        hevi_flag = .false.
      end if
    end select

    return
  end subroutine ATMOS_DYN_Tinteg_short_rkdg_setup

!OCL SERIAL
  subroutine ATMOS_DYN_Tinteg_short_rkdg( &
    dyn_vars,                                    &
    step, nstep, last                            )
    use scale_fem_meshfield_base, only: &
      MeshField3D, MeshField2D
    use scale_atmos_dyn_tstep_short, only: &
      ATMOS_DYN_Tstep_dgm_short_ex, &
      ATMOS_DYN_Tstep_dgm_short_vi
    use scale_atmos_dyn_dgm_vars, only: &
      AtmosDynDGM_vars
    use scale_tracer, only: &
      QA
    use scale_atmos_dyn_tstep_tracer_dgm_heve, only: &
      ATMOS_DYN_Tstep_tracer_dgm_heve_save_massflux, &
      ATMOS_DYN_Tstep_tracer_dgm_heve_calc_fct_coef
    use scale_atmos_dyn_dgm_operator, only: &
      ATMOS_DYN_DGM_operator_apply_modalfilter,   &
      ATMOS_DYN_DGM_operator_addtend_spongelayer, &
      mapprojection
    
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: dyn_vars
    integer, intent(in)  :: step
    integer, intent(in)  :: nstep
    logical,  intent(in) :: last      

    integer :: rkstage
    integer :: tintbuf_ind
    real(RP) :: dt
    real(RP) :: implicit_fac

    real(RP) :: tavg_coef_MFLXZ(tint(1)%nstage)

    class(LocalMesh3D), pointer :: lmesh3D
    integer :: ldomID
    integer :: kelem

    type(MeshField3D), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    type(MeshField3D), pointer :: MFLX_x_tavg,  MFLX_y_tavg,  MFLX_z_tavg
    type(MeshField3D), pointer :: PRES_hyd, PRES_hyd_ref, DENS_hyd
    type(MeshField3D), pointer :: PRES, Rtot, CVtot, CPtot
    type(MeshField3D), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p
    !---------------------------------------------

    call dyn_vars%GetVarsPtr( DDENS, MOMX, MOMY, MOMZ, DRHOT, &
      DENS_hyd, PRES_hyd, PRES_hyd_ref )
    call dyn_vars%GetAuxVarsPtr( Rtot, CVtot, CPtot, PRES=PRES )
    call dyn_vars%GetPhyTendPtr( DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p )
    
    if ( hevi_flag ) then
      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        lmesh3D => mesh3D%lcmesh_list(ldomID)
        call tint(ldomID)%StoreVar0( DDENS%local(ldomID)%val,  DENS_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )
        call tint(ldomID)%StoreVar0( DRHOT%local(ldomID)%val, THERM_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )
        call tint(ldomID)%StoreVar0(  MOMZ%local(ldomID)%val,  MOMZ_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )
        call tint(ldomID)%StoreVar0(  MOMX%local(ldomID)%val,  MOMX_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )
        call tint(ldomID)%StoreVar0(  MOMY%local(ldomID)%val,  MOMY_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )           
      enddo
    end if
    
    do rkstage=1, tint(1)%nstage
    
      !- Calculate the tendency with the implicit part      
      if ( hevi_flag ) then
        do ldomID=1, mesh3D%LOCAL_MESH_NUM
          lmesh3D => mesh3D%lcmesh_list(ldomID)

          call PROF_rapstart( 'ATM_DYN_cal_vi', 2)
          implicit_fac = tint(ldomID)%Get_implicit_diagfac(rkstage)
          tintbuf_ind = tint(ldomID)%tend_buf_indmap(rkstage)
          dt = tint(ldomID)%Get_deltime()

          call ATMOS_DYN_Tstep_dgm_short_vi( &
            tint(ldomID)%tend_buf2D_im(:,:,DENS_VID,tintbuf_ind),                           & ! (out)
            tint(ldomID)%tend_buf2D_im(:,:,MOMX_VID ,tintbuf_ind),                          & ! (out)
            tint(ldomID)%tend_buf2D_im(:,:,MOMY_VID ,tintbuf_ind),                          & ! (out)
            tint(ldomID)%tend_buf2D_im(:,:,MOMZ_VID ,tintbuf_ind),                          & ! (out)
            tint(ldomID)%tend_buf2D_im(:,:,THERM_VID,tintbuf_ind),                          & ! (out)
            DDENS%local(ldomID)%val, MOMX%local(ldomID)%val, MOMY%local(ldomID)%val, MOMZ%local(ldomID)%val, & ! (in)
            DRHOT%local(ldomID)%val,                                                        & ! (in)
            DENS_hyd%local(ldomID)%val, PRES_hyd%local(ldomID)%val,                         & ! (in)
            tint(ldomID)%var0_2D(:,:,DENS_VID), tint(ldomID)%var0_2D(:,:,MOMX_VID),         & ! (in)
            tint(ldomID)%var0_2D(:,:,MOMY_VID ), tint(ldomID)%var0_2D(:,:,MOMZ_VID),        & ! (in)
            tint(ldomID)%var0_2D(:,:,THERM_VID ),                                           & ! (in)
            Rtot%local(ldomID)%val, CVtot%local(ldomID)%val, CPtot%local(ldomID)%val,       & ! (in)
            Dz, Lift,                                                                       & ! (in)
            implicit_fac, dt,                                                               & ! (in)
            lmesh3D, lmesh3D%refElem3D, lmesh3D%lcmesh2D, lmesh3D%lcmesh2D%refElem2D        ) ! (in)          

          call PROF_rapend( 'ATM_DYN_cal_vi', 2)  
          
          call PROF_rapstart( 'ATM_DYN_store_impl', 2)    
          call tint(ldomID)%StoreImplicit( rkstage, DDENS%local(ldomID)%val,  DENS_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )
          call tint(ldomID)%StoreImplicit( rkstage, DRHOT%local(ldomID)%val, THERM_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )
          call tint(ldomID)%StoreImplicit( rkstage,  MOMZ%local(ldomID)%val,  MOMZ_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )
          call tint(ldomID)%StoreImplicit( rkstage,  MOMX%local(ldomID)%val,  MOMX_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )
          call tint(ldomID)%StoreImplicit( rkstage,  MOMY%local(ldomID)%val,  MOMY_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )
          call PROF_rapend( 'ATM_DYN_store_impl', 2) 
        enddo
      end if

      !- Calculate pressure perturbation

      call PROF_rapstart( 'ATM_DYN_cal_pres', 2)
      call dyn_vars%Calcuate_pressure()
      call PROF_rapend( 'ATM_DYN_cal_pres', 2)      

      !- Exchange halo data
      call PROF_rapstart( 'ATM_DYN_exchange_prgv', 2)
      call dyn_vars%Exchange_prgvars()
      call PROF_rapend( 'ATM_DYN_exchange_prgv', 2)

      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        lmesh3D => mesh3D%lcmesh_list(ldomID)
        !- Apply boundary conditions
        call PROF_rapstart( 'ATM_DYN_applyBC_prgv', 2)
        call dyn_bnd%ApplyBC_PROGVARS_lc( ldomID,                                                          & ! (in)
          DDENS%local(ldomID)%val, MOMX%local(ldomID)%val, MOMY%local(ldomID)%val, MOMZ%local(ldomID)%val,            & ! (inout)
          DRHOT%local(ldomID)%val,                                                                                    & ! (inout)
          DENS_hyd%local(ldomID)%val, PRES_hyd%local(ldomID)%val,                                                     & ! (in)
          lmesh3D%Gsqrt(:,:), lmesh3D%GsqrtH(:,:), lmesh3D%GIJ(:,:,1,1), lmesh3D%GIJ(:,:,1,2), lmesh3D%GIJ(:,:,2,2),  & ! (in)
          lmesh3D%GI3(:,:,1), lmesh3D%GI3(:,:,2),                                                                     & ! (in)
          lmesh3D%normal_fn(:,:,1), lmesh3D%normal_fn(:,:,2), lmesh3D%normal_fn(:,:,3),                               & ! (in)
          lmesh3D%vmapM, lmesh3D%vmapP, lmesh3D%vmapB,                                     & ! (in)
          lmesh3D, lmesh3D%refElem3D, lmesh3D%lcmesh2D, lmesh3D%lcmesh2D%refElem2D         ) ! (in)
        call PROF_rapend( 'ATM_DYN_applyBC_prgv', 2)
      enddo

      !- Calculate the tendency with the explicit part
      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        lmesh3D => mesh3D%lcmesh_list(ldomID)

        call PROF_rapstart( 'ATM_DYN_update_caltend_ex', 2)
        tintbuf_ind = tint(ldomID)%tend_buf_indmap(rkstage)
        call ATMOS_DYN_Tstep_dgm_short_ex( &
          tint(ldomID)%tend_buf2D_ex(:,:,DENS_VID,tintbuf_ind),                             & ! (out)
          tint(ldomID)%tend_buf2D_ex(:,:,MOMX_VID ,tintbuf_ind),                            & ! (out)
          tint(ldomID)%tend_buf2D_ex(:,:,MOMY_VID ,tintbuf_ind),                            & ! (out)
          tint(ldomID)%tend_buf2D_ex(:,:,MOMZ_VID ,tintbuf_ind),                            & ! (out)
          tint(ldomID)%tend_buf2D_ex(:,:,THERM_VID,tintbuf_ind),                            & ! (out)
          DDENS%local(ldomID)%val, MOMX%local(ldomID)%val, MOMY%local(ldomID)%val, MOMZ%local(ldomID)%val,               & ! (in)
          DRHOT%local(ldomID)%val, dyn_vars%DPRES%local(ldomID)%val,                                                     & ! (in)
          DENS_hyd%local(ldomID)%val, PRES_hyd%local(ldomID)%val, PRES_hyd_ref%local(ldomID)%val,                        & ! (in)
          dyn_vars%Coriolis%local(ldomID)%val, Rtot%local(ldomID)%val, CVtot%local(ldomID)%val, CPtot%local(ldomID)%val, & ! (in)
          Dx, Dy, Dz, Lift, mapprojection%local(ldomID)%Gam,                                                             & ! (in)
          lmesh3D, lmesh3D%refElem3D, lmesh3D%lcmesh2D, lmesh3D%lcmesh2D%refElem2D                                       ) ! (in)
        call PROF_rapend( 'ATM_DYN_update_caltend_ex', 2)

        !- Sponge layer
        if ( SPONGELAYER_FLAG ) then
          call PROF_rapstart('ATM_DYN_caltend_sponge', 2)
          call ATMOS_DYN_DGM_operator_addtend_spongelayer( &
            tint(ldomID)%tend_buf2D_ex(:,:,MOMX_VID ,tintbuf_ind),   & ! (inout)
            tint(ldomID)%tend_buf2D_ex(:,:,MOMY_VID ,tintbuf_ind),   & ! (inout)
            tint(ldomID)%tend_buf2D_ex(:,:,MOMZ_VID ,tintbuf_ind),   & ! (inout)
            MOMX%local(ldomID)%val, MOMY%local(ldomID)%val, MOMZ%local(ldomID)%val, & ! (in)
            lmesh3D, lmesh3D%refElem3D                                              ) ! (in)
          call PROF_rapend('ATM_DYN_caltend_sponge', 2)
        end if
      enddo

      !- Add physics tendency
      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        lmesh3D => mesh3D%lcmesh_list(ldomID)

        call PROF_rapstart( 'ATM_DYN_update_add_tp', 2)          
        tintbuf_ind = tint(ldomID)%tend_buf_indmap(rkstage)
        call add_phy_tend( &
          tint(ldomID)%tend_buf2D_ex(:,:,:,tintbuf_ind),                                   & ! (inout)
          DENS_tp%local(ldomID)%val, MOMX_tp%local(ldomID)%val, MOMY_tp%local(ldomID)%val, & ! (in)
          MOMZ_tp%local(ldomID)%val, RHOT_tp%local(ldomID)%val,                            & ! (in)
          RHOH_p %local(ldomID)%val, PRES%local(ldomID)%val, Rtot%local(ldomID)%val, CPtot%local(ldomID)%val,  & ! (in)
          ldomID, lmesh3D, lmesh3D%refElem3D                                               ) ! (in)
        call PROF_rapend( 'ATM_DYN_update_add_tp', 2)        
      enddo

      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        lmesh3D => mesh3D%lcmesh_list(ldomID)

        !- Save mass flux
        if ( QA > 0 ) then
          call PROF_rapstart( 'ATM_DYN_tavg_mflx', 2)
          if ( hevi_flag ) then
            tavg_coef_MFLXZ(:) = tint(ldomID)%coef_b_im(:)
          else
            tavg_coef_MFLXZ(:) = tint(ldomID)%coef_b_ex(:)
          end if
          call ATMOS_DYN_Tstep_tracer_dgm_heve_save_massflux( &
            dyn_vars%aux_trcflux(MASSFLX_X_ID)%local(ldomID)%val, & ! (inout)
            dyn_vars%aux_trcflux(MASSFLX_Y_ID)%local(ldomID)%val, & ! (inout)
            dyn_vars%aux_trcflux(MASSFLX_Z_ID)%local(ldomID)%val, & ! (inout)
            dyn_vars%alphaDensM%local(ldomID)%face_val, dyn_vars%alphaDensP%local(ldomID)%face_val,          & ! (inout)
            DDENS%local(ldomID)%val, MOMX%local(ldomID)%val, MOMY%local(ldomID)%val, MOMZ%local(ldomID)%val, & ! (in)
            dyn_vars%DPRES%local(ldomID)%val, DENS_hyd%local(ldomID)%val, PRES_hyd%local(ldomID)%val,        & ! (in)
            Rtot%local(ldomID)%val, CVtot%local(ldomID)%val, CPtot%local(ldomID)%val,                        & ! (in)
            lmesh3D, lmesh3D%refElem3D,                                                                                  & ! (in)
            rkstage==1 .and. step==1, tint(ldomID)%coef_b_ex(rkstage)/real(nstep,kind=RP), tavg_coef_MFLXZ(rkstage)/real(nstep,kind=RP), &
            hevi_flag ) ! (in)
          call PROF_rapend( 'ATM_DYN_tavg_mflx', 2)
        end if        
        
        call PROF_rapstart( 'ATM_DYN_update_advance', 2)      
        call tint(ldomID)%Advance( rkstage, DDENS%local(ldomID)%val,  DENS_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )        
        call tint(ldomID)%Advance( rkstage, DRHOT%local(ldomID)%val, THERM_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )    
        call tint(ldomID)%Advance( rkstage,  MOMZ%local(ldomID)%val,  MOMZ_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )    
        call tint(ldomID)%Advance( rkstage,  MOMX%local(ldomID)%val,  MOMX_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )    
        call tint(ldomID)%Advance( rkstage,  MOMY%local(ldomID)%val,  MOMY_VID, 1, lmesh3D%refElem%Np, lmesh3D%NeS, lmesh3D%NeE )
        call PROF_rapend( 'ATM_DYN_update_advance', 2)
      enddo

    enddo

    call PROF_rapstart( 'ATM_DYN_update_post', 2)
    if ( QA > 0 ) then
      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        lmesh3D => mesh3D%lcmesh_list(ldomID)

        !$omp parallel do
        do kelem=lmesh3D%NeS, lmesh3D%NeE
          dyn_vars%trc_vars(TRCVARS3D_DENS0_ID)%local(ldomID)%val(:,kelem) = tint(ldomID)%var0_2D(:,kelem,DENS_VID)
          dyn_vars%trc_vars(TRCVARS3D_DENS_ID )%local(ldomID)%val(:,kelem) = DDENS%local(ldomID)%val(:,kelem)
        enddo        
     enddo
    end if
    call PROF_rapend( 'ATM_DYN_update_post', 2)   

    if ( MODALFILTER_FLAG ) then
      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        lmesh3D => mesh3D%lcmesh_list(ldomID)

        call PROF_rapstart( 'ATM_DYN_update_modalfilter', 2)
        call ATMOS_DYN_DGM_operator_apply_modalfilter(  & 
          DDENS%local(ldomID)%val, MOMX%local(ldomID)%val, MOMY%local(ldomID)%val, MOMZ%local(ldomID)%val, DRHOT%local(ldomID)%val, & ! (inout)
          lmesh3D, lmesh3D%refElem3D, modal_filter_3d,        & ! (in)
!          do_weight_Gsqrt = .false.                                 ) ! (in)          
          do_weight_Gsqrt = .true.                                   ) ! (in)        
        call PROF_rapend( 'ATM_DYN_update_modalfilter', 2)
      end do
    end if

    return
  end subroutine ATMOS_DYN_Tinteg_short_rkdg
  
!OCL SERIAL
  subroutine ATMOS_DYN_Tinteg_short_rkdg_finalize
    implicit none
    integer :: lcdomid
    !----------------------------------------

    if ( allocated(tint) ) then
      do lcdomid=1, mesh3D%LOCAL_MESH_NUM
        call tint(lcdomid)%Final()
      enddo
      deallocate( tint )
    end if

    return
  end subroutine ATMOS_DYN_Tinteg_short_rkdg_finalize

!-- private ---------------

  !OCL SERIAL
  subroutine add_phy_tend( &
    dyn_tends,                                           & ! (inout)
    DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp,         & ! (in)  
    RHOH_p, PRES, Rtot, CPtot,                           & ! (in)
    domID, lcmesh, elem3D                                ) ! (in)
    use scale_const, only: &
      PRES00 => CONST_PRE00
    use scale_fem_element_base, only: ElementBase3D
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(inout) :: dyn_tends(elem3D%Np,lcmesh%NeA,PRGVAR_NUM)
    real(RP), intent(in) :: DENS_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: RHOT_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: RHOH_p(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem3D%Np,lcmesh%NeA)
    integer, intent(in) :: domID

    integer :: ke
    real(RP) :: EXNER(elem3D%Np)
    real(RP) :: rP0
    !---------------------------------------------------------------------------------

    rP0   = 1.0_RP / PRES00

    !$omp parallel private( EXNER, ke )
    
    !$omp do
    do ke=lcmesh%NeS, lcmesh%NeE
      dyn_tends(:,ke,DENS_VID) = dyn_tends(:,ke,DENS_VID) + DENS_tp(:,ke)
      dyn_tends(:,ke,MOMZ_VID) = dyn_tends(:,ke,MOMZ_VID) + MOMZ_tp(:,ke)
      dyn_tends(:,ke,MOMX_VID) = dyn_tends(:,ke,MOMX_VID) + MOMX_tp(:,ke)
      dyn_tends(:,ke,MOMY_VID) = dyn_tends(:,ke,MOMY_VID) + MOMY_tp(:,ke)
    end do
    !$omp end do

    ! if ( this%ENTOT_CONSERVE_SCHEME_FLAG ) then
    !   !$omp do
    !   do ke=lcmesh%NeS, lcmesh%NeE
    !     EXNER(:) = ( PRES(:,ke) * rP0 )**( Rtot(:,ke) / CPtot(:,ke) )

    !     dyn_tends(:,ke,THERM_VID) = dyn_tends(:,ke,THERM_VID) &
    !                                     + RHOH_p(:,ke)                             &
    !                                     + ( CPtot(:,ke) * EXNER(:) ) * RHOT_tp(:,ke)
    !   end do
    !   !$omp end do
    ! else
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE
        EXNER(:) = ( PRES(:,ke) * rP0 )**( Rtot(:,ke) / CPtot(:,ke) )

        dyn_tends(:,ke,THERM_VID) = dyn_tends(:,ke,THERM_VID) &
                                       + RHOT_tp(:,ke)                             &
                                       + RHOH_p  (:,ke) / ( CPtot(:,ke) * EXNER(:) )
      end do
      !$omp end do
    ! end if
    !$omp end parallel
    return
  end subroutine add_phy_tend
    
end module scale_atmos_dyn_tinteg_short_rkdg
