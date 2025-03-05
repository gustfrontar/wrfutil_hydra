program main
use netcdf

implicit double precision(a-h,o-z)

include "netcdf.inc"

real(4),allocatable::vd(:,:,:)
real(4),allocatable::vr(:,:,:)
real(4),allocatable::vg(:,:,:)
real(4),allocatable::vs(:,:,:)
real(4),allocatable::vref(:,:,:)

character(len=200)::ncfile
character(len=100)::ctext

character(len=200)::basename_in,basename_out,vars
integer::nprocs_x_out,nprocs_y_out

character(len=200)::conf

namelist / PARAM_SNO / &
       basename_in,           &
       basename_out,           &
       vars,           &
       nprocs_x_out,          &
       nprocs_y_out

call getarg(1,conf)

open(11,file=trim(conf))
  read(11,nml=PARAM_SNO)
close(11)

ncfile=trim(basename_out)//".pe000000.nc"
istat=NF_OPEN(trim(ncfile),NF_WRITE,idnc)
if (istat.ne.0)then
  write(*,*)'No'
  write(*,*) trim(ncfile)
  stop
end if
  write(*,*) trim(ncfile)

istat=NF_INQ_VARID(idnc,"DENS",idvd)
istat=NF_INQ_VARID(idnc,"QR",idvr)
istat=NF_INQ_VARID(idnc,"QG",idvg)
istat=NF_INQ_VARID(idnc,"QS",idvs)

istat=NF_INQ_DIMID(idnc,"y",iddlat)
if (istat.ne.0)then
  write(*,*)'No'
  write(*,*) istat
  write(*,*) NF_STRERROR(istat)
  stop
end if

istat=NF_INQ_DIMLEN(idnc,iddlat,nlat)
if (istat.ne.0)then
  write(*,*)'No'
  write(*,*) istat
  write(*,*) NF_STRERROR(istat)
  stop
end if

istat=NF_INQ_DIMID(idnc,"x",iddlon)
istat=NF_INQ_DIMLEN(idnc,iddlon,nlon)
istat=NF_INQ_DIMID(idnc,"z",iddz)
istat=NF_INQ_DIMLEN(idnc,iddz,nz)

write(*,*) nlon,nlat,nz

if (.not.allocated(vd))then
allocate(vd(nz,nlon,nlat))
allocate(vr(nz,nlon,nlat))
allocate(vg(nz,nlon,nlat))
allocate(vs(nz,nlon,nlat))
allocate(vref(nz,nlon,nlat))
end if

istat=NF_GET_VAR_REAL(idnc,idvd,vd)
istat=NF_GET_VAR_REAL(idnc,idvr,vr)
istat=NF_GET_VAR_REAL(idnc,idvg,vg)
istat=NF_GET_VAR_REAL(idnc,idvs,vs)

write(*,*) "dens",maxval(vd),minval(vd)
write(*,*) "rain",maxval(vr),minval(vr)
write(*,*) "grap",maxval(vg),minval(vg)
write(*,*) "snow",maxval(vs),minval(vs)

vref=calc_ref(nz,nlon,nlat,vd,vr,vg,vs)
write(*,*) "ref ",maxval(vref),minval(vref)

istat=NF_REDEF(idnc)

istat=NF_DEF_VAR(idnc,"REF",NF_REAL,3,(/iddz,iddlon,iddlat/),idvref)
if (istat /=0) then
  istat=NF_INQ_VARID(idnc,"REF",idvref)
end if
istat=NF_PUT_VAR_REAL(idnc,idvref,vref)
ctext="Radar reflectivity"
istat=NF_PUT_ATT_TEXT(idnc,idvref,"long_name",len(trim(ctext)),trim(ctext))
istat=NF_PUT_ATT_TEXT(idnc,idvref,"units",3,"dbz")
!rmiss=-9.9999e+30
!istat=NF_PUT_ATT_REAL(idnc,idvref,"_FillValue",NF_FLOAT,1,rmiss)
!istat=NF_PUT_ATT_REAL(idnc,idvref,"missing_value",NF_FLOAT,1,rmiss)

istat=NF_ENDDEF(idnc)

istat=NF_CLOSE(idnc)


stop

contains 
function calc_ref(nz,nlon,nlat,vd,vr,vg,vs)
real(4)::vd(nz,nlon,nlat)
real(4)::vr(nz,nlon,nlat)
real(4)::vg(nz,nlon,nlat)
real(4)::vs(nz,nlon,nlat)
real(4)::calc_ref(nz,nlon,nlat)

real(4)::vdbz

real(4),parameter::vmiss=-9.9999d30
real(4),parameter::eps=1.0d-10
real(4),parameter::vref_lbnd=10.0
real(4),parameter::vref_norain=5.0
real(4),parameter::afact_r=2.53d4, bfact_r = 1.84
real(4),parameter::afact_s=3.48d3, bfact_s = 1.66
real(4),parameter::afact_g=5.54d3, bfact_g = 1.70

do iz=1,nz
do ilat=1,nlat
do ilon=1,nlon
 if (vd(iz,ilon,ilat) > 0.0 ) then
   vdbz                  = afact_r * ( 1.0d3* vd(iz,ilon,ilat) * max ( eps, vr(iz,ilon,ilat) ) ) ** bfact_r & 
                         + afact_g * ( 1.0d3* vd(iz,ilon,ilat) * max ( eps, vg(iz,ilon,ilat) ) ) ** bfact_g &
                         + afact_s * ( 1.0d3* vd(iz,ilon,ilat) * max ( eps, vs(iz,ilon,ilat) ) ) ** bfact_s 
   vdbz = 10.0 * log10(vdbz)
   if (vdbz < vref_lbnd ) vdbz = vref_norain
 else
   vdbz = vmiss
 end if
 calc_ref(iz,ilon,ilat) = vdbz
end do
end do
end do
end function

end program
