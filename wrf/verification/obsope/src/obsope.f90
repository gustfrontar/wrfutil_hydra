PROGRAM obsope
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
  USE obsop_tools

  IMPLICIT NONE
  REAL(r_size),ALLOCATABLE :: omb(:)
  REAL(r_size),ALLOCATABLE :: oma(:)
  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr
  CHARACTER(8) :: stdoutf='OBSO-000'
!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
  CALL CPU_TIME(rtimer00)
  CALL initialize_mpi

  !Set default observation list
  !Default configuration for variable list.
  !This can be set throuhg namelist.
  variable_list=0
  variable_list(1:13)=                                                                                               &  !List of variables to be verified
  & (/ id_u_obs , id_v_obs , id_t_obs , id_tv_obs , id_rh_obs , id_q_obs , id_reflectivity_obs , id_radialwind_obs , &
  &    id_ps_obs , id_us_obs , id_vs_obs , id_rhs_obs , id_qs_obs /)



  !Get configuration from namelist
  CALL read_namelist()
!
  WRITE(stdoutf(6:8), '(I3.3)') myrank
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
  WRITE(6,'(A)') '         OBSOP & VERIFICATION TOOL           '
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
  ! CONVENTIONAL OBS
  !
  CALL set_letkf_obs
!  WRITE(*,*)"TERMINE DE PROCESAR LAS OBS"
!
  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_OBS):',rtimer,rtimer-rtimer00
  rtimer00=rtimer

  if( myrank == 0 .and. do_obsgrid )then
    !Get the dimensions of the output grid.
    CALL get_regrid_grid()
    !From observations to grid.
    call obsgrid()
  endif

!
!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

STOP
END PROGRAM obsope
