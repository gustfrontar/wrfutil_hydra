program main
use common_ncio
implicit real(a-h,o-z)

real(4),allocatable::axlon(:),axlat(:),topo(:,:)

integer,parameter::nelm=5
integer,parameter::elms(nelm)=(/14593,82819,82820,83073,83331/) !! Ps,U10,V10,T2,RH2 
real(4),parameter::errs(nelm)=(/1.0,3.0,3.0,1.0,10.0/)  !! error std (not used)

real(4)::wk(8)
character(len=200)::cfile
character(len=200)::ncfile_in='history_merge.pe000000.nc'

integer::ncid, vidlon, vidlat,vidz

  call ncio_open( trim(ncfile_in), nf90_nowrite, ncid )
  call ncio_check( nf90_inq_dimid(ncid, "lon", vidlon))
  call ncio_check( nf90_inq_dimid(ncid, "lat", vidlat))
  call ncio_check( nf90_inquire_dimension(ncid, vidlon, len=nlon))
  call ncio_check( nf90_inquire_dimension(ncid, vidlat, len=nlat))

  allocate(topo(nlon,nlat))
  allocate(axlon(nlon))
  allocate(axlat(nlat))

  call ncio_check( nf90_inq_varid(ncid, "lon", vidlon))
  call ncio_check( nf90_inq_varid(ncid, "lat", vidlat))
  call ncio_check( nf90_inq_varid(ncid, "topo", vidz))
  call ncio_check( nf90_get_var(ncid, vidlon, axlon))
  call ncio_check( nf90_get_var(ncid, vidlat, axlat))
  call ncio_check( nf90_get_var(ncid, vidz, topo))
  call ncio_close( ncid ) 

cfile="test_obs_AWS.dat"

open (21, file=trim(cfile), form='unformatted', access='sequential', convert='little_endian')

do ilat=1,nlat
do ilon=1,nlon
do ie=1,nelm
  wk(1)=real(elms(ie))  
  wk(2)=axlon(ilon)
  wk(3)=axlat(ilat)
  wk(4)=topo(ilon,ilat)
  wk(5)=10.0  !!! dat
  wk(6)=errs(ie)   !!! err 
  wk(7)=22.0  !!! typ AWSARG
  wk(8)=0.0   !!! dif
  write(21,iostat=ios) wk(1:8)
end do
end do
end do

close(21)

stop
end program main
