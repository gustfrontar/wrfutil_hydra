PROGRAM letkf
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
  USE common_letkf
  USE letkf_obs
  USE letkf_tools
  USE common_namelist

  IMPLICIT NONE
  REAL(r_size),ALLOCATABLE :: omb(:)
  REAL(r_size),ALLOCATABLE :: oma(:)
  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr
  CHARACTER(10) :: stdoutf='NOUT-00000'
  CHARACTER(4) :: guesf='gs00'
  CHARACTER(7) :: monitobs='obs.dat'
  CHARACTER(20):: argument
!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
  CALL CPU_TIME(rtimer00)
  CALL initialize_mpi

  !Get configuration from namelist
  CALL read_namelist()
!
  WRITE(stdoutf(6:10), '(I5.5)') myrank
!  WRITE(6,'(3A,I3.3)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I3.3,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf
!
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '  LOCAL ENSEMBLE TRANSFORM KALMAN FILTERING  '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '   LL      EEEEEE  TTTTTT  KK  KK  FFFFFF    '
  WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(6,'(A)') '   LL      EEEEE     TT    KKK     FFFFF     '
  WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(6,'(A)') '   LLLLLL  EEEEEE    TT    KK  KK  FF        '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '             WITHOUT LOCAL PATCH             '
  WRITE(6,'(A)') '   WITH RADAR FORWARD OBSERVATION OPERATOR   '
  WRITE(6,'(A)') '   WITH DIRECT NETCDF I/O                    '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '          Coded by Takemasa Miyoshi          '
  WRITE(6,'(A)') '  Based on Ott et al (2004) and Hunt (2005)  '
  WRITE(6,'(A)') '  Tested by Miyoshi and Yamane (2006)        '
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '              LETKF PARAMETERS               '
  WRITE(6,'(A)') ' ------------------------------------------- '
  WRITE(6,'(A,I15)')   '  nbv          :',nbv
  WRITE(6,'(A,I15)')   '  nslots       :',nslots
  WRITE(6,'(A,I15)')   '  nbslot       :',nbslot
  WRITE(6,'(A,F15.2)') '  sigma_obs    :',sigma_obs
  WRITE(6,'(A,F15.2)') '  sigma_obs_radar   :',sigma_obs_radar 
  WRITE(6,'(A,F15.2)') '  sigma_obsv   :',sigma_obsv
  WRITE(6,'(A,F15.2)') '  sigma_obst   :',sigma_obst
  WRITE(6,'(A,L7)')    '  estpar       :',ESTPAR
  WRITE(6,'(A)') '============================================='
  CALL set_common_wrf('gs0100001')
  CALL set_common_mpi_wrf
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(INITIALIZE):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Observations
!-----------------------------------------------------------------------
  !
  ! CONVENTIONAL AND RADAR OBS
  !
  CALL set_letkf_obs

  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_OBS):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! First guess ensemble
!-----------------------------------------------------------------------
  !
  ! READ GUES
  !
  !WRITE(*,*)'CTRL1'
  ALLOCATE(gues3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(gues2d(nij1,nbv,nv2d))

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  WRITE(guesf(3:4),'(I2.2)') nbslot
  !WRITE(*,*)'CTRL2'
  CALL read_ens_mpi_alltoall(guesf,nbv,gues3d,gues2d)
  !CALL read_ens_mpi(guesf,nbv,gues3d,gues2d)
  !WRITE(*,*)'CTRL3'
  IF(ESTPAR)ALLOCATE(guesp2d(nij1,nbv,np2d))
  IF(ESTPAR)CALL read_ensp_mpi(guesf,nbv,guesp2d)
  !WRITE(*,*)'CTRL4'
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_GUES):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
  !
  ! WRITE ENS MEAN and SPRD
  !
  !WRITE(*,*)'CTRL5'
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL write_ensmspr_mpi('gues',nbv,gues3d,gues2d)

!  IF(ESTPAR)CALL write_enspmspr_mpi('gues',nbv,guesp2d)
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(WRITE_GUES_MEAN):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
  !
  ! LETKF
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !WRITE(*,*)'CTRL6'
  ALLOCATE(anal3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(anal2d(nij1,nbv,nv2d))

  IF(ESTPAR)ALLOCATE(analp2d(nij1,nbv,np2d))
  !WRITE(*,*)'CTRL7'
  CALL das_letkf()

  DEALLOCATE(gues3d,gues2d)

  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(DAS_LETKF):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Analysis ensemble
!-----------------------------------------------------------------------
  !
  ! WRITE ANAL
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !WRITE(*,*)'CTRL8'
  CALL write_ens_mpi_alltoall(guesf,nbv,anal3d,anal2d)
  !CALL write_ens_mpi(guesf,nbv,anal3d,anal2d)

  IF(ESTPAR)CALL write_ensp_mpi(guesf,nbv,analp2d)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(WRITE_ANAL_ens):',rtimer,rtimer-rtimer00
  rtimer00=rtimer

  !
  ! WRITE ENS MEAN and SPRD
  !
  !WRITE(*,*)'CTRL9'
  CALL write_ensmspr_mpi('anal',nbv,anal3d,anal2d)

  IF(ESTPAR)CALL write_enspmspr_mpi('anal',nbv,analp2d)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(anal3d,anal2d)
  IF(ESTPAR)DEALLOCATE(analp2d)
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(WRITE_ANAL_mean):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
  !WRITE(*,*)'CTRL11'
  IF(myrank == 0) THEN
    ALLOCATE(omb(nobs),oma(nobs))
!   Compute RMSE and BIAS for the gues mean.
    CALL monit_mean('gues',omb)
!   Compute RMSE and BIAS for the analysis mean. 
    CALL monit_mean('anal',oma)
    CALL monit_obs(monitobs,nobs,obselm,obslon,obslat,obslev,obsdat,&
     & obserr,obstyp,omb,oma,obslot)
  END IF
!
  !WRITE(*,*)'CTRL12'

  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(MONIT_MEAN):',rtimer,rtimer-rtimer00
  rtimer00=rtimer
!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  WRITE(6,'(A)') ' LETKF ends SUCCESFULLY !!!!'
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi



STOP
END PROGRAM letkf
