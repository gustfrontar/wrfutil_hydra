PROGRAM gues_mean
!=======================================================================
!
! [PURPOSE:] Main program of LETKF
!
! [HISTORY:]
!   01/16/2009 Takemasa Miyoshi  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_wrf
  USE common_mpi_wrf
  USE common_namelist
  !USE letkf_tools
  !USE common_letkf

  IMPLICIT NONE
  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr 
  CHARACTER(8) :: stdoutf='NOUT-000'

  CHARACTER(10) :: mem_prefix='wrfarw.mem'
  CHARACTER(14) :: mean_prefix='wrfarw.ensmean'
  CHARACTER(1)  :: mem_sufix='' , mean_sufix=''
  CHARACTER(20) :: nbv_input

  REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:) ! background ensemble 3d variables
  REAL(r_size),ALLOCATABLE :: gues2d(:,:,:)   ! background ensemble 2d variables
!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
  CALL CPU_TIME(rtimer00)
  CALL initialize_mpi

  !Get configuration from namelist
  !CALL read_namelist()
  CALL  getarg(1,nbv_input)

  READ(nbv_input,*)nbv !nbv es un argumento de entrada al programa, no lo toma del namelist.
!
  WRITE(stdoutf(6:8), '(I3.3)') myrank
!  WRITE(6,'(3A,I3.3)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I3.3,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf
!
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A,I15)')   '  nbv          :',nbv
  WRITE(6,'(A)') '============================================='
  CALL set_common_wrf('wrfarw.mem001')
  CALL set_common_mpi_wrf
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(INITIALIZE):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! First guess ensemble
!-----------------------------------------------------------------------
  !
  ! READ GUES
  !
  ALLOCATE(gues3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(gues2d(nij1,nbv,nv2d))

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL read_ens_mpi_gsi(mem_prefix,mem_sufix,nbv,gues3d,gues2d)

  ! WRITE ENS MEAN and SPRD
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL write_ensmspr_mpi_gsi(mean_prefix,mean_sufix,nbv,gues3d,gues2d)

!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_GUES):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

STOP
END PROGRAM gues_mean
