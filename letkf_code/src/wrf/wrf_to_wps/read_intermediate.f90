PROGRAM read_intermediate
!Read intermediate file in WPS format

  USE met_data_module
  USE common

  IMPLICIT NONE
  INTEGER(r_sngl),PARAMETER :: grdfid = 50

  CHARACTER(10) :: grdout= 'output.grd'


  !WRITE IN WPS FORMAT
  open(unit=grdfid, file=trim(grdout), status='old', form='unformatted')
  !OPEN(grdfid,FILE=grdout,FORM='unformatted',ACCESS='sequential')

  CALL READ_SLAB(grdfid)

  CLOSE(grdfid)

  STOP
END PROGRAM read_intermediate



