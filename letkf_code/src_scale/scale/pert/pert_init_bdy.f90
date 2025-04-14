PROGRAM pert_init_bdy
!=======================================================================
!
! [PURPOSE:] Main program to compute initial and boundary ensemble perturbations.
!
! [HISTORY:]
!   01/16/2009 Takemasa Miyoshi  created
!   10/16/2023 Juan Ruiz adapted to met em perturbation
!   04/06/2025 Arata Amemiya adapted to SCALE initial and boundary perturbation
!
!=======================================================================

  USE common
  USE common_mpi
  USE common_scale
  USE common_mpi_scale
  USE common_nml
  use common_mtx, only: &
    mtx_setup
  USE pert_tools
  IMPLICIT NONE

  integer :: ifile,  iter , niter, ini_mem , end_mem , mem_iter , nbv_iter
  character(len=7) :: stdoutf = '-000000'
  character(len=6400) :: icmd

!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  call initialize_mpi_scale
  call mpi_timer('', 1)

  if (command_argument_count() >= 2) then
    call get_command_argument(2, icmd)
    if (trim(icmd) /= '') then
      write (stdoutf(2:7), '(I6.6)') myrank
!      write (6,'(3A,I6.6)') 'STDOUT goes to ', trim(icmd)//stdoutf, ' for MYRANK ', myrank
      open (6, file=trim(icmd)//stdoutf)
      write (6,'(A,I6.6,2A)') 'MYRANK=', myrank, ', STDOUTF=', trim(icmd)//stdoutf
    end if
  end if

!-----------------------------------------------------------------------

  if (DET_RUN) then
    call set_mem_node_proc(MEMBER+2)
  else
    call set_mem_node_proc(MEMBER+1)
  end if
  call set_scalelib('PERT')

  if ( LOG_OUT ) then
    WRITE(6,'(A)') '============================================='
    WRITE(6,'(A)') '  MET EM PERTURBATION TOOL                   '
    WRITE(6,'(A)') '============================================='
    WRITE(6,'(A)') '              PARAMETERS                     '
    WRITE(6,'(A)') ' ------------------------------------------- '
    WRITE(6,'(A,I15)')   '  nbv_ori      :',MEMBER
    WRITE(6,'(A,I15)')   '  nbv_tar      :',MEMBER_TAR
    WRITE(6,'(A,I15)')   '  method       :',PERT_METHOD
    WRITE(6,'(A)') ' ------------------------------------------- '
    WRITE(6,'(A)') '              FILES to PROCESS               '
    WRITE(6,'(A)') ' ------------------------------------------- '
    WRITE(6,'(A,I15)')   '  # of init    :',NUM_PERT_RESTART
    WRITE(6,'(A,I15)')   '  # of bdy     :',NUM_PERT_BOUNDARY
    WRITE(6,'(A,I15)')   '  iteration    :',NITER_PERT
    WRITE(6,'(A)') '============================================='
  end if 

  if (myrank_use) then

    call set_common_scale
    call set_common_mpi_scale

    call mtx_setup( MEMBER )

    call mpi_timer('INITIALIZE', 1, barrier=MPI_COMM_a)

    !
    ! LETKF GRID setup
    !
    call set_common_mpi_grid

    allocate (eori3d(nij1,nlev,nens,nv3d))
    allocate (eori2d(nij1,nens,nv2d))
    allocate (eper3d(nij1,nlev,nens,nv3d))
    allocate (eper2d(nij1,nens,nv2d))

    call mpi_timer('SET_GRID', 1, barrier=MPI_COMM_a)

    call get_weights

    call mpi_timer('SET_WEIGHTS', 1, barrier=MPI_COMM_a)

!-----------------------------------------------------------------------
! First guess ensemble
!-----------------------------------------------------------------------


    do ifile = 1, NUM_PERT_RESTART
    !
    ! READ RESTART
    !
    GUES_IN_BASENAME=PERT_RESTART_IN_BASENAME(ifile)
    GUES_MEAN_INOUT_BASENAME=PERT_RESTART_IN_BASENAME(ifile)
    call filename_replace_mem(GUES_MEAN_INOUT_BASENAME, 'mean')
    call read_ens_mpi(eori3d, eori2d)

    if (DET_RUN .and. mmdetin /= mmdet) then
      eori3d(:,:,mmdet,:) = eori3d(:,:,mmdetin,:)
      eori2d(:,mmdet,:) = eori2d(:,mmdetin,:)
    end if

    call mpi_timer('READ_RESTART', 1, barrier=MPI_COMM_a)

!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
    do iter=1, niter_pert
    ini_mem=(iter-1)*MEMBER+1
    end_mem=min((iter-1)*MEMBER+MEMBER,MEMBER_TAR)
    nbv_iter=end_mem-ini_mem+1
    !
    ! LETKF
    !
    call transform_pert( wmatrix(:,ini_mem:end_mem),nbv_iter )
    call mpi_timer('PERT_TRANSFORM', 1, barrier=MPI_COMM_a)

!-----------------------------------------------------------------------
! Analysis ensemble
!-----------------------------------------------------------------------

    !
    ! WRITE ANAL and ENS MEAN
    !
    call write_ens_mpi_ext(eper3d, eper2d, iter-1, PERT_RESTART_OUT_BASENAME(ifile),.false.)

    call mpi_timer('WRITE_PERT_RESTART', 1, barrier=MPI_COMM_a)

    end do !!! iter

    end do !!! ifile

!!-----------------------------------------------------------------------
!! Boundary
!!-----------------------------------------------------------------------

    do ifile = 1, NUM_PERT_BOUNDARY
    !
    ! READ BOUNDARY
    !
    call read_ens_bdy_mpi(eori3d, eori2d, PERT_BOUNDARY_IN_BASENAME(ifile))

    if (DET_RUN .and. mmdetin /= mmdet) then
      eori3d(:,:,mmdet,:) = eori3d(:,:,mmdetin,:)
      eori2d(:,mmdet,:) = eori2d(:,mmdetin,:)
    end if

    call mpi_timer('READ_BOUNDARY', 1, barrier=MPI_COMM_a)

!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
    do iter=1, niter_pert
    ini_mem=(iter-1)*MEMBER+1
    end_mem=min((iter-1)*MEMBER+MEMBER,MEMBER_TAR)
    nbv_iter=end_mem-ini_mem+1
    !
    ! LETKF
    !
    call transform_pert( wmatrix(:,ini_mem:end_mem),nbv_iter )
    call mpi_timer('PERT_TRANSFORM', 1, barrier=MPI_COMM_a)

!-----------------------------------------------------------------------
! Analysis ensemble
!-----------------------------------------------------------------------

    !
    ! WRITE ANAL and ENS MEAN
    !
    call write_ens_mpi_ext(eper3d, eper2d, iter-1, PERT_BOUNDARY_OUT_BASENAME(ifile),.true.)

    call mpi_timer('WRITE_PERT_BOUNDARY', 1, barrier=MPI_COMM_a)

!!-----------------------------------------------------------------------
!! Monitor
!!-----------------------------------------------------------------------

    end do !!! iter

    end do !!! ifile

    deallocate (eper3d, eper2d)
    deallocate (eori3d, eori2d)

    call unset_common_mpi_scale

  end if ! [ myrank_use ]

  call unset_scalelib

  call mpi_timer('FINALIZE', 1, barrier=MPI_COMM_WORLD)

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------

  if ( myrank == 0 ) then
    write(6,'(a)') 'pert_init_bdy finished successfully'
  endif

  call finalize_mpi_scale

  STOP
END PROGRAM pert_init_bdy
