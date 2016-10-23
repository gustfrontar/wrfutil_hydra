PROGRAM dec_airsco
  USE common
  USE common_airs
  
  IMPLICIT NONE
  INTEGER        :: i,j,k , ii
  INTEGER        :: iy,im,id,ih,imin
  INTEGER        :: iy1,im1,id1,ih1,imin1
  INTEGER        :: iy2,im2,id2,ih2,imin2
  REAL(r_sngl)   :: sec4
  REAL(r_size)   :: sec
  LOGICAL        :: ex
  CHARACTER(14)  :: outfile='yyyymmddhh.dat'
  CHARACTER(100) :: filename
  REAL(r_size)   :: co_prior_profile(nzco)
  INTEGER        :: data_count(25)
!                           123456789012345
  REAL(r_size) :: tile_lat_mean , tile_lon_mean
!=======================================================================
! READING HDF
!=======================================================================

  CALL GETARG( 1, filename )
  !Read data from hdf file
  CALL read_hdf(filename) 
  !Read first guess data 
  CALL read_first_guess()
  !Read trapezoid limits data
  CALL read_trapezoid_limits()

  data_count=0

  !Write the data
  CALL com_tai2utc(time1,iy1,im1,id1,ih1,imin1,sec)
  IF(imin1 > 30) CALL com_tai2utc(time1+1800.d0,iy1,im1,id1,ih1,imin1,sec)
  WRITE(outfile(1:10),'(I4.4,3I2.2)') iy1,im1,id1,ih1

  INQUIRE(FILE=outfile,EXIST=file_exists)
  IF( file_exists )THEN
    OPEN(nunit+ih1,FILE=outfile,FORM='unformatted',ACCESS='sequential')
    READ(nunit+ih1)data_count(ih1+1)
    CLOSE(nunit+ih1)
    WRITE(*,*)"OPENING THE FILE ",outfile
    OPEN(nunit+ih1,FILE=outfile,FORM='unformatted',ACCESS='sequential',position='APPEND')
  ELSE
    OPEN(nunit+ih1,FILE=outfile,FORM='unformatted',ACCESS='sequential')
    WRITE(nunit+ih1)data_count(ih1+1),id_co_obs,airstype
    WRITE(*,*)"CREATING THE FILE ",outfile
  ENDIF

  CALL com_tai2utc(time2,iy2,im2,id2,ih2,imin2,sec)
  IF(imin2 > 30) CALL com_tai2utc(time2+1800.d0,iy2,im2,id2,ih2,imin2,sec)
  IF(ih2 /= ih1) THEN
    WRITE(outfile(1:10),'(I4.4,3I2.2)') iy2,im2,id2,ih2
    INQUIRE(FILE=outfile,EXIST=file_exists)
    IF( file_exists )THEN
      OPEN(nunit+ih2,FILE=outfile,FORM='unformatted',ACCESS='sequential')
      READ(nunit+ih2)data_count(ih2+1)
      CLOSE(nunit+ih2)
      WRITE(*,*)"OPENING THE FILE ",outfile
      OPEN(nunit+ih2,FILE=outfile,FORM='unformatted',ACCESS='sequential',position='APPEND')
    ELSE
      OPEN(nunit+ih2,FILE=outfile,FORM='unformatted',ACCESS='sequential')
      WRITE(nunit+ih2)data_count(ih2+1),id_co_obs,airstype
      WRITE(*,*)"CREATING THE FILE ",outfile
    ENDIF

  END IF

  DO j=1,ny
    DO i=1,nx
      tile_lat_mean=tile_lat_mean + rlat(i,j)
      tile_lon_mean=tile_lon_mean + rlon(i,j)
      IF(rlon(i,j) < 0) rlon(i,j) = rlon(i,j) + 360.0d0
      IF(rlon(i,j) < lon1) CYCLE
      IF(rlon(i,j) > lon2) CYCLE
      IF(rlat(i,j) < lat1) CYCLE
      IF(rlat(i,j) > lat2) CYCLE

      IF(mod(i,iskip) /= 1) CYCLE
      IF(mod(j,iskip) /= 1) CYCLE

      !Determine in which file we will store the data.
      CALL com_tai2utc(time(i,j),iy,im,id,ih,imin,sec)
      IF(imin > 30) CALL com_tai2utc(time(i,j)+1800.d0,iy,im,id,ih,imin,sec)
      data_count(ih+1)=data_count(ih+1)+1

      !Compute the correct first guess depending on the time of the year
      !and the latitude.
      CALL get_first_guess(rlat(i,j),im,co_prior_profile)

      WRITE(nunit+ih)rlon(i),rlat(i),time(i)
      WRITE(nunit+ih)nzco
 
       CALL invert_order(covmr(:,i,j),nzco)
       CALL invert_order(covmre(:,i,j),nzco)
       CALL invert_order(covmrqc(:,i,j),nzco)

      WRITE(nunit+ih)covmr(:,i,j)
      WRITE(nunit+ih)covmre(:,i,j)
      WRITE(nunit+ih)covmrqc(:,i,j)

      WRITE(nunit+ih)pressure_levels
      WRITE(nunit+ih)surface_pressure

      DO ii=1,n_trapezoids
         CALL invert_order(co_averaging_kernel(ii,:,i,j),n_trapezoids)
         WRITE(nunit+ih)co_averaging_kernel(ii,:,i,j)
      ENDDO
      WRITE(nunit+ih)nzco
      WRITE(nunit+ih)co_prior_profile
      WRITE(nunit+ih)pressure_levels

     ENDDO
    ENDDO
 
  !We open the file corresponding to the previous day
   REWIND(nunit+ih1)
   WRITE(nunit+ih1)data_count(ih1+1)
   IF( ih2 /= ih1)THEN
   REWIND(nunit+ih2)
   WRITE(nunit+ih2)data_count(ih2+1)

  STOP



END PROGRAM dec_airsco


