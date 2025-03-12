program main
use common_ncio
implicit real(a-h,o-z)

real(4),allocatable::axlon(:,:),axlat(:,:),axz(:),pres(:,:,:)

integer,parameter::nelm=3
integer,parameter::elms(nelm)=(/4001,4002,4004/) !! ze, vr, ze(zero)
real(4),parameter::errs(nelm)=(/5.0,3.0,5.0/)     
real(4),parameter::vdist_limit = 80.0e3

real(4),parameter::er=6400.0e3
real(4),parameter::drad=3.141592/180.0

real(4)::wk(8)
character(len=200)::cfile
character(len=200)::ncfile_in='history_merge.pe000000.nc'

integer::ncid, vidlon, vidlat,vidz

real(4),parameter::radar_lon = 135.232
real(4),parameter::radar_lat = 34.662
real(4),parameter::radar_z   = 0.0  

  call ncio_open( trim(ncfile_in), nf90_nowrite, ncid )
  call ncio_check( nf90_inq_dimid(ncid, "x", vidlon))
  call ncio_check( nf90_inq_dimid(ncid, "y", vidlat))
  call ncio_check( nf90_inq_dimid(ncid, "z", vidz))
  call ncio_check( nf90_inquire_dimension(ncid, vidlon, len=nlon))
  call ncio_check( nf90_inquire_dimension(ncid, vidlat, len=nlat))
  call ncio_check( nf90_inquire_dimension(ncid, vidz, len=nlev))

  allocate(pres(nlon,nlat,nlev))
  allocate(axlon(nlon,nlat))
  allocate(axlat(nlon,nlat))
  allocate(axz(nlev))

  call ncio_check( nf90_inq_varid(ncid, "lon", vidlon))
  call ncio_check( nf90_inq_varid(ncid, "lat", vidlat))
  call ncio_check( nf90_inq_varid(ncid, "z", vidz))
  call ncio_check( nf90_get_var(ncid, vidlon, axlon))
  call ncio_check( nf90_get_var(ncid, vidlat, axlat))
  call ncio_check( nf90_get_var(ncid, vidz, axz))
  call ncio_check( nf90_inq_varid(ncid, "PRES", vidz))
  call ncio_check( nf90_get_var(ncid, vidz, pres))
  call ncio_close( ncid ) 

!cfile="test_obs_3d_xyp.dat"
cfile="test.dat"

open (21, file=trim(cfile), form='unformatted', access='sequential', convert='big_endian')

! Radar header
write(21) radar_lon
write(21) radar_lat
write(21) radar_z

do ilon=1,nlon
do ilat=1,nlat
do ilev=1,nlev

  vdist=  (er*drad*(axlat(ilon,ilat)-radar_lat))**2  &
        + (er*cos(drad*radar_lat)*drad*(axlon(ilon,ilat)-radar_lon))**2 & 
        + (axz(ilev)-radar_z)**2
  vdist=sqrt(vdist)
  if(vdist < vdist_limit)then
    do ie=1,nelm
      wk(1)=elms(ie)!!! elm radar ref
      wk(2)=axlon(ilon,ilat)
      wk(3)=axlat(ilon,ilat)
      wk(4)=axz(ilev)
      wk(5)=10.0  !!! dat (dummy)
      wk(6)=errs(ie)   !!! err
      wk(7)=22.0  !!! typ PHARAD
      wk(8)=0.0   !!! dif
      write(21,iostat=ios) wk(1:7)
    end do
  else 
  end if

end do
end do
end do

close(21)

stop
end program main
