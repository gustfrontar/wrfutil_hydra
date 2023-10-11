PROGRAM read_intermediate
!Read intermediate file in WPS format

  USE met_data_module
  USE common

  IMPLICIT NONE
  INTEGER(r_sngl),PARAMETER :: grdfid = 50

  CHARACTER(500) :: grdin

  !Get input file name
  CALL GETARG ( 1, grdin  ) !Input file name (with complete path)



  !WRITE IN WPS FORMAT
  open(unit=grdfid, file=trim(grdin), status='old', form='unformatted')

  CALL READ_SLAB(grdfid)



  CLOSE(grdfid)

  STOP
END PROGRAM read_intermediate



