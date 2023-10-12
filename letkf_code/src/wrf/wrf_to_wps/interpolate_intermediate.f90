PROGRAM read_intermediate
!Read intermediate file in WPS format

  USE met_data_module
  USE common

  IMPLICIT NONE
  INTEGER,PARAMETER :: grdfid1 = 50 , grdfid2 = 51 , grdfido = 52
  REAL(r_size)   :: ini_time , end_time , target_time 
  CHARACTER(500) :: grdin_1 , grdin_2 , grdout 
  CHARACTER(30)  :: int_arg 
  TYPE(met_data) :: slab_data_1 , slab_data_2  
  REAL(r_size)   :: int_w
  INTEGER        :: reading_ierr = 0 , number_of_records = 0

  !Get input file name
  CALL GETARG ( 1, grdin_1  ) !Input  file name 1 (corresponding to initial time)
  CALL GETARG ( 2, grdin_2  ) !Input  file name 2 (corresponding to final time)
  CALL GETARG ( 3, grdout   ) !Output file name   (output file with the interpolated fields)
  CALL GETARG ( 4, int_arg  )
  READ( int_arg , * )ini_time     !Initial time corresponding to input file 1 (in any continous time units, eg, hours, seconds, etc).
  CALL GETARG ( 5, int_arg  )
  READ( int_arg , * )end_time     !Final time corresponding to input file 2 (in any continous time units, eg, hours, seconds, etc).
  CALL GETARG ( 6, int_arg  )
  READ( int_arg , * )target_time  !Target time corresponding to the output file (in any continous time units, eg, hours, seconds, etc).


  !WRITE IN WPS FORMAT
  open(unit=grdfid1, file=trim(grdin_1), status='old', form='unformatted')
  open(unit=grdfid2, file=trim(grdin_2), status='old', form='unformatted')
  open(unit=grdfido, file=trim(grdout), status='unknown', form='unformatted')

  if ( target_time == ini_time .or. target_time == end_time ) then
     !We will not overwrite the existing files. 
     WRITE(*,*)'No interpolation is required for this time'
     STOP 0
  endif 

  do while ( reading_ierr == 0 ) 

     !Read data records from the input files (ini_time and end_time)
     CALL GET_SLAB_RECORD( grdfid1 , slab_data_1 , reading_ierr )   !Read the initial time data
     CALL GET_SLAB_RECORD( grdfid2 , slab_data_2 , reading_ierr )   !Read the final time data

     if ( reading_ierr == 0 ) then 

        number_of_records = number_of_records + 1

        write(*,*)'Found a record ',number_of_records,' ',slab_data_1 % field

        !Do some checks ....  verify that we have the same grid and the same
        !variable at both files. This should be the case if both files where
        !generated with wrf_to_wps from a wrfout file.
        if ( ( slab_data_1 % field /= slab_data_2 % field ) .or. &
          & ( slab_data_1 % units /= slab_data_2 % units ) .or. &  
          & ( slab_data_1 % level /= slab_data_2 % level ) .or. &
          & ( slab_data_1 % nx    /= slab_data_2 % nx    ) .or. &
          & ( slab_data_1 % ny    /= slab_data_2 % ny    ) .or. &
          & ( slab_data_1 % iproj /= slab_data_2 % iproj ) .or. &
          & ( slab_data_1 % startlat /= slab_data_2 % startlat ) .or. &
          & ( slab_data_1 % startlon /= slab_data_2 % startlon ) ) then
          WRITE(*,*)'Error: input fields do not correspond to the same variable'
          STOP 1  
        endif

        !Interpolate the data 
        !Compute the interpolation coefficient.
        int_w = ( target_time - ini_time ) / ( end_time - ini_time )
        slab_data_1 % slab = slab_data_1 % slab * (1.0 - int_w ) + slab_data_2 % slab * int_w 
        !Write the interpolated slab. 
        slab_data_1 % nlev = 1
        CALL PUT_SLAB_RECORD( grdfido , slab_data_1 )

        IF( ALLOCATED( slab_data_1 % slab ) ) DEALLOCATE( slab_data_1 % slab )
        IF( ALLOCATED( slab_data_2 % slab ) ) DEALLOCATE( slab_data_2 % slab )
        IF( ALLOCATED( slab_data_1 % xlvl ) ) DEALLOCATE( slab_data_1 % xlvl )
        IF( ALLOCATED( slab_data_2 % xlvl ) ) DEALLOCATE( slab_data_2 % xlvl )


     endif

  end do 

  close(grdfido)
  close(grdfid2)
  close(grdfid1)

  if ( number_of_records == 0 ) then
     WRITE(*,*)'Error: We reached the end of the file, but no record could be found'
     WRITE(*,*)'No interpolation has been performed'
     STOP 1
  endif

  WRITE(*,*)'We successfully interpolated ',number_of_records,' records'
  STOP 0

END PROGRAM read_intermediate



