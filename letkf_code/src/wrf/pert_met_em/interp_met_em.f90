PROGRAM interp_met_em
!=======================================================================
!
! [PURPOSE:] Interpolate in time two met_em files.
!
! [HISTORY:]
!   01/16/2009 Takemasa Miyoshi  created
!   10/16/2023 Juan Ruiz adapted to interp_met_em
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_met_em
  USE common_namelist_met_em

  IMPLICIT NONE
  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr 
  INTEGER :: m
  REAL(r_sngl) :: int_w
  CHARACTER(8) :: stdoutf='NOUT-000'
  CHARACTER(4) :: prefix='____'
  REAL(r_sngl) , ALLOCATABLE :: v3d_ini(:,:,:,:) , v2d_ini(:,:,:) , vs3d_ini(:,:,:,:)
  REAL(r_sngl) , ALLOCATABLE :: v3d_end(:,:,:,:) , v2d_end(:,:,:) , vs3d_end(:,:,:,:)
  CHARACTER(len=1100)        :: os_command
!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------

  !Get configuration from namelist
  CALL read_namelist()
!
  OPEN(6,FILE=stdoutf)
!
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '  MET EM INTERPOLATION TOOL                  '
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '              PARAMETERS                     '
  WRITE(6,'(A)') ' ------------------------------------------- '
  WRITE(6,'(A,F10.2)')   '  time_ini    :',time_ini
  WRITE(6,'(A,F10.2)')   '  time_end    :',time_end
  WRITE(6,'(A,F10.2)')   '  time_tar    :',time_tar
  WRITE(6,'(2A)')      '  date_tar    :',date_tar
  WRITE(6,'(2A)')      '  file_ini    :',TRIM(file_ini)
  WRITE(6,'(2A)')      '  file_end    :',TRIM(file_end)
  WRITE(6,'(2A)')      '  file_tar    :',TRIM(file_tar)
  WRITE(6,'(A)') '============================================='
  CALL set_common_met_em( file_ini )

  IF ( time_ini == time_tar .or. time_end == time_tar )THEN
     WRITE(6,'(A)') ' TARGET TIME IS EQUAL TO INITIAL OR FILAN TIME '
     WRITE(6,'(A)') ' WE DO NOT NEED TO INTERPOLATE ' 
     STOP 0
  ENDIF    

  os_command = 'cp '//file_ini//' '//file_tar
  CALL SYSTEM( os_command )
  CALL write_date_met_em(file_tar,TRIM(date_tar))

  !
  ! READ FILES
  !

  ALLOCATE(v3d_ini(nlon,nlat,nlev,nv3d),v3d_end(nlon,nlat,nlev,nv3d))
  ALLOCATE(v2d_ini(nlon,nlat,nv2d),v2d_end(nlon,nlat,nv2d)) 
  ALLOCATE(vs3d_ini(nlon,nlat,nlev,ns3d),vs3d_end(nlon,nlat,nlev,ns3d))

  CALL read_grd(file_ini,v3d_ini,v2d_ini,vs3d_ini)
  CALL read_grd(file_end,v3d_end,v2d_end,vs3d_end)

  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(READ_FILES):',rtimer,rtimer-rtimer00
  rtimer00=rtimer


  !
  ! INTERPOLATE FILES IN TIME
  !

  int_w = ( time_tar - time_ini ) / ( time_end - time_ini )
  v3d_end = v3d_ini * ( 1.0 - int_w ) + v3d_end * int_w 
  v2d_end = v2d_end * ( 1.0 - int_w ) + v2d_end * int_w 
  vs3d_end= vs3d_end* ( 1.0 - int_w ) + vs3d_end* int_w 


  CALL write_grd(file_tar,v3d_end,v2d_end,vs3d_end)
 
!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
 
  DEALLOCATE( v3d_ini , v2d_ini , vs3d_ini )
  DEALLOCATE( v3d_end , v2d_end , vs3d_end )

  WRITE(6,'(A)') ' MET EM INTERPOLATION COMPLETED SUCCESSFULLY !!!!'

STOP
END PROGRAM interp_met_em
