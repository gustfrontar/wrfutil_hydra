PROGRAM pert_met_em
!=======================================================================
!
! [PURPOSE:] Main program to compute me_em ensemble perturbations.
!
! [HISTORY:]
!   01/16/2009 Takemasa Miyoshi  created
!   10/16/2023 Juan Ruiz adapted to met em perturbation
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_met_em
  USE common_mpi_met_em
  USE met_em_tools
  USE common_namelist_met_em

  IMPLICIT NONE
  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr , itime 
  CHARACTER(8) :: stdoutf='NOUT-000'
  CHARACTER(4) :: prefix='____'
  CHARACTER(20):: argument
!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
  CALL CPU_TIME(rtimer00)
  CALL initialize_mpi

  !Get configuration from namelist
  CALL read_namelist()
!
  WRITE(stdoutf(6:8), '(I3.3)') myrank
!  WRITE(6,'(3A,I3.3)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I3.3,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf
!
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '  MET EM PERTURBATION TOOL                   '
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '              PARAMETERS                     '
  WRITE(6,'(A)') ' ------------------------------------------- '
  WRITE(6,'(A,I15)')   '  nbv_ori      :',nbv_ori
  WRITE(6,'(A,I15)')   '  nbv_tar      :',nbv_tar
  WRITE(6,'(A,I15)')   '  nbv_ntimes   :',ntimes  
  WRITE(6,'(A,I15)')   '  method       :',method
  WRITE(6,'(A)') '============================================='
  CALL set_common_met_em('eo0100001')
  CALL set_common_mpi_met_em

  DO itime = 1 , ntimes 
!
    CALL CPU_TIME(rtimer)
    WRITE(6,'(A,2F10.2)') '### TIMER(INITIALIZE):',rtimer,rtimer-rtimer00
    rtimer00=rtimer
!-----------------------------------------------------------------------
! Original ensemble
!-----------------------------------------------------------------------
  !
  ! READ ORIGINAL MET EM FILES
  !

    ALLOCATE(eori3d(nij1,nlev,nbv_ori,nv3d))
    ALLOCATE(eori2d(nij1,nbv_ori,nv2d))
    ALLOCATE(eoris3d(nij1,nlev_soil,nbv_ori,ns3d))

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    WRITE(prefix(1:4),'(A,I2.2)')'eo',itime
    CALL read_ens_mpi_alltoall(prefix,nbv_ori,eori3d,eori2d,eoris3d)

    CALL CPU_TIME(rtimer)
    WRITE(6,'(A,2F10.2)') '### TIMER(READ_GUES):',rtimer,rtimer-rtimer00
    rtimer00=rtimer

!-----------------------------------------------------------------------
! Perturbation computation
!-----------------------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    ALLOCATE(eper3d(nij1,nlev,nbv_tar,nv3d))
    ALLOCATE(eper2d(nij1,nbv_tar,nv2d))
    ALLOCATE(epers3d(nij1,nlev_soil,nbv_tar,ns3d))

    CALL transform_pert()

    DEALLOCATE(eori3d,eori2d,eoris3d)

    CALL CPU_TIME(rtimer)
    WRITE(6,'(A,2F10.2)') '### TIMER(PERTURBATION TRANSFORM):',rtimer,rtimer-rtimer00
    rtimer00=rtimer

!-----------------------------------------------------------------------
! Analysis ensemble
!-----------------------------------------------------------------------
  !
  ! WRITE PERTURBED METEM FILES
  !
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    WRITE(prefix(1:4),'(A,I2.2)')'ep',itime
    CALL write_ens_mpi_alltoall(prefix,nbv_tar,eper3d,eper2d,epers3d)

    CALL CPU_TIME(rtimer)
    WRITE(6,'(A,2F10.2)') '### TIMER(WRITE_ANAL_ens):',rtimer,rtimer-rtimer00
    rtimer00=rtimer

    DEALLOCATE(eper3d,eper2d,epers3d)

  END DO ! End do over times  
!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  

  WRITE(6,'(A)') ' MET EM PERTURBATION COMPLETED SUCCESSFULLY !!!!'
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi


STOP
END PROGRAM pert_met_em
