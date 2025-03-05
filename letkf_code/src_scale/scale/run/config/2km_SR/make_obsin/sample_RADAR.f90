program main
use common_ncio
implicit real(a-h,o-z)

real(4),allocatable::axlon(:),axlat(:),axz(:)

integer,parameter::nelm=3
integer,parameter::elms(nelm)=(/4001,4002,4004/) !! ze, vr, ze(zero)
real(4),parameter::errs(nelm)=(/5.0,3.0,5.0/)     
real(4),parameter::vdist_limit = 160.0e3

real(4),parameter::er=6400.0e3
real(4),parameter::drad=3.141592/180.0

real(4)::wk(8)
character(len=200)::cfile
character(len=200)::ncfile_in='history_merge.pe000000.nc'

integer::ncid, vidlon, vidlat,vidz
 
!!! Cordoba
real(4),parameter::radar_lon = 295.8080749511719
real(4),parameter::radar_lat = -31.44132995605469
real(4),parameter::radar_z   = 476.0 
!!!! Ezeiza
!real(4),parameter::radar_lon = 301.4844360351562
!real(4),parameter::radar_lat = -34.80081939697266
!real(4),parameter::radar_z   = 47.0  


  call ncio_open( trim(ncfile_in), nf90_nowrite, ncid )
  call ncio_check( nf90_inq_dimid(ncid, "lon", vidlon))
  call ncio_check( nf90_inq_dimid(ncid, "lat", vidlat))
  call ncio_check( nf90_inq_dimid(ncid, "z", vidz))
  call ncio_check( nf90_inquire_dimension(ncid, vidlon, len=nlon))
  call ncio_check( nf90_inquire_dimension(ncid, vidlat, len=nlat))
  call ncio_check( nf90_inquire_dimension(ncid, vidz, len=nlev))

  allocate(axlon(nlon))
  allocate(axlat(nlat))
  allocate(axz(nlev))

  call ncio_check( nf90_inq_varid(ncid, "lon", vidlon))
  call ncio_check( nf90_inq_varid(ncid, "lat", vidlat))
  call ncio_check( nf90_inq_varid(ncid, "z", vidz))
  call ncio_check( nf90_get_var(ncid, vidlon, axlon))
  call ncio_check( nf90_get_var(ncid, vidlat, axlat))
  call ncio_check( nf90_get_var(ncid, vidz, axz))
  call ncio_close( ncid ) 

axlon=axlon+360.0

cfile="test_obs_radar.dat"

open (21, file=trim(cfile), form='unformatted', access='sequential', convert='little_endian')

! Radar header
write(21) radar_lon
write(21) radar_lat
write(21) radar_z

do ilon=1,nlon
do ilat=1,nlat
do ilev=1,nlev

  vdist=  (er*drad*(axlat(ilat)-radar_lat))**2  &
        + (er*cos(drad*radar_lat)*drad*(axlon(ilon)-radar_lon))**2 & 
        + (axz(ilev)-radar_z)**2 !!! rough estimate
  vdist=sqrt(vdist)
  if(vdist < vdist_limit)then
    do ie=1,nelm
      wk(1)=elms(ie)
      wk(2)=axlon(ilon)
      wk(3)=axlat(ilat)
      wk(4)=axz(ilev)
      wk(5)=10.0  !!! dat (dummy)
      wk(6)=errs(ie)   !!! err (dummy)
      wk(7)=12.0  !!! typ RADARC
      wk(8)=0.0   !!! dif
      write(21,iostat=ios) wk(1:8)
    end do
  else 
  end if

end do
end do
end do

close(21)

stop
end program main
