program main
use common_ncio
implicit real(a-h,o-z)

real(4),allocatable::axlon(:,:),axlat(:,:),axz(:),pres(:,:,:)

integer,parameter::nelm=4
integer,parameter::elms(nelm)=(/2819,2820,3073,3331/) !! U,V,T,RH 
real(4),parameter::errs(nelm)=(/3.0,3.0,1.0,5.0/)    !! U,V,T,RH(%) 

real(4)::wk(8)
character(len=200)::cfile
character(len=200)::ncfile_in='history_merge.pe000000.nc'

integer::ncid, vidlon, vidlat,vidz

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

do ilon=1,nlon
do ilat=1,nlat
do ilev=1,nlev
do ie=1,nelm
  wk(1)=real(elms(ie))  
  wk(2)=axlon(ilon,ilat)
  wk(3)=axlat(ilon,ilat)
!  wk(4)=axz(iz)
  wk(4)=pres(ilon,ilat,ilev) * 0.01 !!! hPa
  wk(5)=10.0  !!! dat
  wk(6)=errs(ie)   !!! err 
  wk(7)=1.0  !!! typ ADPUPA
  wk(8)=0.0   !!! dif
  write(21,iostat=ios) wk(1:8)
end do
!  write(*,'(F6.1,5F14.4)') wk(1:6)
end do
end do
end do

close(21)

stop
end program main
