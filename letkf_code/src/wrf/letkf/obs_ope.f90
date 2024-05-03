PROGRAM obs_ope
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
  USE letkf_obs
  USE common_namelist

  IMPLICIT NONE
  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr
  CHARACTER(10) :: stdoutf='OBSO-00000'
  CHARACTER(4) :: guesf='gs00'
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
  WRITE(6,'(A)') '  OBSERVATION OPERATOR                       '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '          Coded by Takemasa Miyoshi          '
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '              LETKF PARAMETERS               '
  WRITE(6,'(A)') ' ------------------------------------------- '
  WRITE(6,'(A,I15)')   '  nbv          :',nbv
  WRITE(6,'(A,I15)')   '  nslots       :',nslots
  WRITE(6,'(A,I15)')   '  nbslot       :',nbslot
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

  !Call the observation operator
  CALL set_letkf_obs

  !Write the different ensemble members in the observation space
  CALL obsope_io(nbv,nobs,2,obselm,obslon,obslat,obslev,obsdat,obserr,obstyp,obslot  &
                        ,obsi,obsj,obsdep,obshdxf )

  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_OBS):',rtimer,rtimer-rtimer00
  rtimer00=rtimer

  WRITE(6,'(A)') ' OBSOPE ends SUCCESFULLY !!!!'
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi



STOP
END PROGRAM obs_ope
