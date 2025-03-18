!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics Temporal integration
!!
!! @par Description
!!          Temporal integration scheme selecter for dynamical tracer advection based on DGM
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_tracer_rkdg
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer

  use scale_timeint_rk, only: &
    TimeInt_RK

  use scale_fem_localmesh_3d, only: &
    LocalMesh3D
  use scale_atmos_dyn_dgm_operator, only: &
    mesh3D, elem3D, &
    dyn_bnd,                  &
    Dx, Dy, Dz, Lift,         &
    FaceIntMat,               &
    ONLY_TRCADV_FLAG,         &
    TRCADV_MODALFILTER_FLAG,  &
    modal_filter_3d_tracer,   &
    TRCADV_disable_limiter
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_Tinteg_tracer_rkdg_setup
  public :: ATMOS_DYN_Tinteg_tracer_short_rkdg
  public :: ATMOS_DYN_Tinteg_tracer_rkdg_finalize

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

contains
  !-----------------------------------------------------------------------------
  !> Register
  subroutine ATMOS_DYN_Tinteg_tracer_rkdg_setup( &
       ATMOS_DYN_Tinteg_tracer_TYPE, &
       dt )

    use scale_precision
    use scale_atmos_grid_cartesC_index
    use scale_index
    use scale_prc, only: &
       PRC_abort
    implicit none
    character(len=*), intent(in)  :: ATMOS_DYN_Tinteg_tracer_TYPE
    real(RP), intent(in) :: dt

    integer :: lcdomID
    !---------------------------------------------------------------------------

    select case( ATMOS_DYN_Tinteg_tracer_TYPE )
    case( 'OFF', 'NONE' )
       ! do nothing
    case default
      allocate( tint(mesh3D%LOCAL_MESH_NUM) )
      do lcdomID=1, mesh3D%LOCAL_MESH_NUM
        call tint(lcdomID)%Init( ATMOS_DYN_Tinteg_tracer_TYPE, dt, 1, &
          2, (/ elem3D%Np, mesh3D%lcmesh_list(lcdomID)%NeA /) )
      enddo
    end select

    return
  end subroutine ATMOS_DYN_Tinteg_tracer_rkdg_setup

!OCL SERIAL
  subroutine ATMOS_DYN_Tinteg_tracer_short_rkdg( &
      QTRC, QTRC_tmp, FCT_coef, iq,                &
      MFLX_x_tavg, MFLX_y_tavg, MFLX_z_tavg,       &
      RHOQ_tp, alphDensM, alphDensP,               &
      DENS_hyd, DDENS, DDENS_TRC, DDENS0_TRC,      &
      dyn_vars, BND_QA, BND_IQ,                    &
      last                                         )
    use scale_fem_meshfield_base, only: &
      MeshField3D, MeshField2D
    use scale_atmos_dyn_tstep_tracer_dgm_heve, only: &
      ATMOS_DYN_Tstep_tracer_dgm_heve_cal_tend,      &
      ATMOS_DYN_Tstep_tracer_dgm_heve_calc_fct_coef, &
      ATMOS_DYN_Tstep_tracer_dgm_heve_TMAR
    use scale_atmos_dyn_dgm_vars, only: &
      AtmosDynDGM_vars
    implicit none
    type(MeshField3D), intent(inout) :: QTRC
    type(MeshField3D), intent(inout) :: QTRC_tmp
    type(MeshField3D), intent(inout) :: FCT_coef
    integer, intent(in) :: iq
    type(MeshField3D), intent(in) :: MFLX_x_tavg
    type(MeshField3D), intent(in) :: MFLX_y_tavg
    type(MeshField3D), intent(in) :: MFLX_z_tavg
    type(MeshField3D), intent(in) :: RHOQ_tp
    type(MeshField3D), intent(in) :: alphDensM
    type(MeshField3D), intent(in) :: alphDensP
    type(MeshField3D), intent(in) :: DENS_hyd
    type(MeshField3D), intent(in) :: DDENS
    type(MeshField3D), intent(in) :: DDENS_TRC
    type(MeshField3D), intent(in) :: DDENS0_TRC
    class(AtmosDynDGM_vars), intent(inout), target :: dyn_vars
    integer,  intent(in)    :: BND_QA
    integer,  intent(in)    :: BND_IQ(QA)
    logical,  intent(in)    :: last   
    
    integer :: rkstage
    integer :: tintbuf_ind
    real(RP) :: dt
    real(RP) :: dttmp_trc
    
    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: lcdomID
    integer :: kelem
    !----------------------------------------------------------------

    do lcdomID=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(lcdomID)

      !$omp parallel do
      do kelem=lcmesh3D%NeS, lcmesh3D%NeE
        QTRC_tmp%local(lcdomID)%val(:,kelem) = QTRC%local(lcdomID)%val(:,kelem)
      enddo            
    enddo

    do rkstage=1, tint(1)%nstage
      if ( TRACER_ADVC(iq) ) then
        call PROF_rapstart( 'ATM_DYN_exchange_trc', 3)
        call dyn_vars%Exchange_trcvar()
        call PROF_rapend( 'ATM_DYN_exchange_trc', 3)

        do lcdomID=1, mesh3D%LOCAL_MESH_NUM
          lcmesh3D => mesh3D%lcmesh_list(lcdomID)

          if ( BND_IQ(iq) > 0 ) then
            call dyn_bnd%ApplyBC_TRCVAR_lc( lcdomID, iq, & ! (in)
              QTRC_tmp%local(lcdomid)%val,                                                        & ! (inout)
              lcmesh3D%vmapM, lcmesh3D%vmapP, lcmesh3D%vmapB,                                     & ! (in)
              lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D        ) ! (in)
          end if
          dt = tint(lcdomID)%Get_deltime()
          dttmp_trc = dt * tint(lcdomID)%coef_gam_ex(rkstage+1,rkstage) &
                       / tint(lcdomID)%coef_sig_ex(rkstage+1,rkstage)

          call ATMOS_DYN_Tstep_tracer_dgm_heve_calc_fct_coef( &
            FCT_coef%local(lcdomID)%val,                                                                    & ! (out)
            QTRC_tmp%local(lcdomID)%val,                                                                    & ! (in)
            MFLX_x_tavg%local(lcdomID)%val, MFLX_y_tavg%local(lcdomID)%val, MFLX_z_tavg%local(lcdomID)%val, & ! (in) 
            RHOQ_tp%local(lcdomID)%val,                                                                     & ! (in)
            alphDensM%local(lcdomID)%face_val, alphDensP%local(lcdomID)%face_val,                           & ! (in)
            DENS_hyd%local(lcdomID)%val, DDENS_TRC%local(lcdomID)%val, DDENS0_TRC%local(lcdomID)%val,       & ! (in)
            tint(lcdomID)%coef_c_ex(rkstage), dttmp_trc,                                                    & ! (in) 
            Dx, Dy, Dz, Lift, FaceIntMat,                                                                   & ! (in)
            lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D,                   & ! (in)
            TRCADV_disable_limiter                                                                          ) ! (in)

          call PROF_rapstart( 'ATM_DYN_exchange_trc', 3)
          call dyn_vars%Exchange_qtrcauxvars()
          call PROF_rapend( 'ATM_DYN_exchange_trc', 3)            
        end do        
      endif

      do lcdomID=1, mesh3D%LOCAL_MESH_NUM
        lcmesh3D => mesh3D%lcmesh_list(lcdomID)
        tintbuf_ind = tint(lcdomID)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'ATM_DYN_update_caltend_ex_trc', 3)  
        if ( TRACER_ADVC(iq) ) then 
          call ATMOS_DYN_Tstep_tracer_dgm_heve_cal_tend( &        
            tint(lcdomID)%tend_buf2D_ex(:,:,1,tintbuf_ind),                                                 & ! (out)
            QTRC_tmp%local(lcdomID)%val,                                                                    & ! (in)
            MFLX_x_tavg%local(lcdomID)%val, MFLX_y_tavg%local(lcdomID)%val, MFLX_z_tavg%local(lcdomID)%val, & ! (in) 
            alphDensM%local(lcdomID)%face_val, alphDensP%local(lcdomID)%face_val,                           & ! (in)
            FCT_coef%local(lcdomID)%val,                                                                    & ! (out)
            RHOQ_tp%local(lcdomID)%val,                                                                     & ! (in)
            Dx, Dy, Dz, Lift, FaceIntMat,                                                                   & ! (in)
            lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D                    ) ! (in)
        else
          !$omp parallel do
          do kelem=lcmesh3D%NeS, lcmesh3D%NeE
            tint(lcdomID)%tend_buf2D_ex(:,kelem,1,tintbuf_ind) = RHOQ_tp%local(lcdomID)%val(:,kelem)
          end do
        end if
        call PROF_rapend( 'ATM_DYN_update_caltend_ex_trc', 3)

        call PROF_rapstart( 'ATM_DYN_update_advance_trc', 3)                
        call tint(lcdomID)%Advance_trcvar( &
          rkstage, QTRC_tmp%local(lcdomID)%val, 1,                                     &
          1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE,                    &
          DDENS_TRC%local(lcdomID)%val, DDENS0_TRC%local(lcdomID)%val, DENS_hyd%local(lcdomID)%val ) 
        call PROF_rapend( 'ATM_DYN_update_advance_trc', 3)

        if ( rkstage == tint(1)%nstage .and. TRCADV_MODALFILTER_FLAG ) then
          call PROF_rapstart( 'ATM_DYN_update_qtrc_modalfilter', 3)
          call atm_dyn_dgm_tracer_modalfilter_apply( &
            QTRC_tmp%local(lcdomID)%val,                               & ! (inout)
            DENS_hyd%local(lcdomID)%val, DDENS_TRC%local(lcdomID)%val, & ! (in)
            lcmesh3D, lcmesh3D%refElem3D, modal_filter_3d_tracer       ) ! (in)
          call PROF_rapend( 'ATM_DYN_update_qtrc_modalfilter', 3)
        end if

        if ( TRACER_ADVC(iq)                   &
          .and. rkstage == tint(1)%nstage      &
          .and. ONLY_TRCADV_FLAG               &            
          .and. ( .not. TRCADV_disable_limiter ) ) then

          call PROF_rapstart( 'ATM_DYN_update_trc_TMAR', 3)             
          call ATMOS_DYN_Tstep_tracer_dgm_heve_TMAR( &
            QTRC_tmp%local(lcdomID)%val,                                                 & ! (inout)
            DENS_hyd%local(lcdomID)%val, DDENS_TRC%local(lcdomID)%val,                   & ! (in)
            lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D ) ! (in)
          call PROF_rapend( 'ATM_DYN_update_trc_TMAR', 3)  
        end if
      end do      

    enddo ! end for RK loop

    do lcdomID=1, mesh3D%LOCAL_MESH_NUM
      !$omp parallel do
      do kelem=lcmesh3D%NeS, lcmesh3D%NeE
        QTRC%local(lcdomID)%val(:,kelem) = ( DENS_hyd%local(lcdomID)%val(:,kelem) + DDENS_TRC%local(lcdomID)%val(:,kelem) ) &
                                         / ( DENS_hyd%local(lcdomID)%val(:,kelem) + DDENS%local(lcdomID)%val(:,kelem) )     &
                                         * QTRC_tmp%local(lcdomID)%val(:,kelem)
      end do            
    end do

    return
  end subroutine ATMOS_DYN_Tinteg_tracer_short_rkdg

!OCL SERIAL
  subroutine ATMOS_DYN_Tinteg_tracer_rkdg_finalize
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
  end subroutine ATMOS_DYN_Tinteg_tracer_rkdg_finalize

!- Private
!OCL SERIAL
  subroutine atm_dyn_dgm_tracer_modalfilter_apply(  &
    QTRC_,                                   & ! (inout)
    DENS_hyd_, DDENS_,                       & ! (inout)
    lmesh, elem, filter                      ) ! (in)

    use scale_fem_element_base, only: ElementBase
    use scale_fem_localmesh_base, only: LocalMeshBase
    use scale_fem_element_modalfilter, only: ModalFilter  
    implicit none

    class(LocalMeshBase), intent(in) :: lmesh
    class(ElementBase), intent(in) :: elem    
    real(RP), intent(inout)  :: QTRC_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    class(ModalFilter), intent(in) :: filter

    integer :: ke
    real(RP) :: tmp(elem%Np)
    integer :: ii, kk
    real(RP) :: Mik

    real(RP) :: weight(elem%Np)
    !------------------------------------

    !$omp parallel do private( tmp, ii, kk, Mik, weight )
    do ke=lmesh%NeS, lmesh%NeE

      tmp(:) = 0.0_RP
      weight(:) = lmesh%Gsqrt(:,ke) &
                * ( DENS_hyd_(:,ke) + DDENS_(:,ke) )

      do ii=1, elem%Np
      do kk=1, elem%Np
        Mik = filter%FilterMat(ii,kk) * weight(kk)

        tmp(ii) = tmp(ii) + Mik * QTRC_(kk,ke)
      end do
      end do

      QTRC_(:,ke) = tmp(:) / weight(:)
    end do    

    return
  end subroutine atm_dyn_dgm_tracer_modalfilter_apply

end module scale_atmos_dyn_tinteg_tracer_rkdg
