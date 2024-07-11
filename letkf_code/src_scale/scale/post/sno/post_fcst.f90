program main
use netcdf

implicit double precision(a-h,o-z)

include "netcdf.inc"

real(4),allocatable::vd(:,:,:,:)
real(4),allocatable::vr(:,:,:,:)
real(4),allocatable::vg(:,:,:,:)
real(4),allocatable::vs(:,:,:,:)
real(4),allocatable::vref(:,:,:,:)
real(4),allocatable::vmref(:,:,:)

real(4),allocatable::vp_s(:,:,:)
real(4),allocatable::vt_s(:,:,:)
real(4),allocatable::vqv_s(:,:,:)
real(4),allocatable::vtd_s(:,:,:)
real(4),allocatable::vrh_s(:,:,:)

!real(4)::rmiss(1)

character(len=200)::ncfile
character(len=200)::ncfiles
character(len=100)::ctext

character(len=200)::basename_in,basename_out
character(len=20)::vars(100)
character(len=200)::basename_out_s
integer::nprocs_x_out,nprocs_y_out

character(len=200)::conf

namelist / PARAM_SNO / &
       basename_in,           &
       basename_out,           &
       vars,           &
       nprocs_x_out,          &
       nprocs_y_out

call getarg(1,conf)

if (trim(conf)=="")then
  write(*,*) "input conf"
  stop
end if

open(11,file=trim(conf))
  read(11,nml=PARAM_SNO)
close(11)

indx=index(trim(basename_out),'history')

if (indx==0)then
  write(*,*) "input file is not history"
  stop
else
  basename_out_s=basename_out
  write(basename_out_s(indx:indx+6),'(A7)') "histsfc"
end if

ncfile=trim(basename_out)//".pe000000.nc"
ncfiles=trim(basename_out_s)//".pe000000.nc"
istat=NF_OPEN(trim(ncfile),NF_WRITE,idnc)
istat=NF_OPEN(trim(ncfiles),NF_WRITE,idncs)
if (istat.ne.0)then
  write(*,*)'No'
  write(*,*) trim(ncfile)
  write(*,*) trim(ncfiles)
  stop
end if

istat=NF_INQ_VARID(idnc,"DENS",idvd)
istat=NF_INQ_VARID(idnc,"QR",idvr)
istat=NF_INQ_VARID(idnc,"QG",idvg)
istat=NF_INQ_VARID(idnc,"QS",idvs)

istat=NF_INQ_DIMID(idnc,"lat",iddlat)
istat=NF_INQ_DIMLEN(idnc,iddlat,nlat)
istat=NF_INQ_DIMID(idnc,"lon",iddlon)
istat=NF_INQ_DIMLEN(idnc,iddlon,nlon)
istat=NF_INQ_DIMID(idnc,"z",iddz)
istat=NF_INQ_DIMLEN(idnc,iddz,nz)
istat=NF_INQ_DIMID(idnc,"time",iddt)
istat=NF_INQ_DIMLEN(idnc,iddt,nt)

istat=NF_INQ_VARID(idncs,"SFC_PRES",idvsp)
istat=NF_INQ_VARID(idncs,"T2",idvst)
istat=NF_INQ_VARID(idncs,"Q2",idvsq)

write(*,*) nlon,nlat,nz,nt

allocate(vd(nlon,nlat,nz,nt))
allocate(vr(nlon,nlat,nz,nt))
allocate(vg(nlon,nlat,nz,nt))
allocate(vs(nlon,nlat,nz,nt))
allocate(vref(nlon,nlat,nz,nt))
allocate(vmref(nlon,nlat,nt))
allocate(vp_s(nlon,nlat,nt))
allocate(vt_s(nlon,nlat,nt))
allocate(vtd_s(nlon,nlat,nt))
allocate(vqv_s(nlon,nlat,nt))
allocate(vrh_s(nlon,nlat,nt))

!istat=NF_GET_VAR_DOUBLE(idnc,idvd,vd)
!istat=NF_GET_VAR_DOUBLE(idnc,idvr,vr)
!istat=NF_GET_VAR_DOUBLE(idnc,idvg,vg)
!istat=NF_GET_VAR_DOUBLE(idnc,idvs,vs)

istat=NF_GET_VAR_REAL(idnc,idvd,vd)
istat=NF_GET_VAR_REAL(idnc,idvr,vr)
istat=NF_GET_VAR_REAL(idnc,idvg,vg)
istat=NF_GET_VAR_REAL(idnc,idvs,vs)

istat=NF_GET_VAR_REAL(idncs,idvp,vp_s)
istat=NF_GET_VAR_REAL(idncs,idvt,vt_s)
istat=NF_GET_VAR_REAL(idncs,idvq,vqv_s)



write(*,*) "dens",maxval(vd),minval(vd)
write(*,*) "rain",maxval(vr),minval(vr)
write(*,*) "grap",maxval(vg),minval(vg)
write(*,*) "snow",maxval(vs),minval(vs)

vref=calc_ref(nlon,nlat,nz,nt,vd,vr,vg,vs)
vmref=maxval(vref,3)
write(*,*) "ref ",maxval(vref),minval(vref)
write(*,*) "mref ",maxval(vmref),minval(vmref)

call calc_rh_td_2d(nlon,nlat,nt,vp_s,vt_s,vqv_s,vrh_s,vtd_s)


istat=NF_REDEF(idnc)
istat=NF_REDEF(idncs)

istat=NF_DEF_VAR(idnc,"REF",NF_REAL,4,(/iddlon,iddlat,iddz,iddt/),idvref)
if (istat /=0) then
  istat=NF_INQ_VARID(idnc,"REF",idvref)
end if
istat=NF_DEF_VAR(idncs,"MAXREF",NF_REAL,3,(/iddlon,iddlat,iddt/),idvmref)
if (istat /=0) then
  istat=NF_INQ_VARID(idncs,"MAXREF",idvmref)
end if
istat=NF_PUT_VAR_REAL(idnc,idvref,vref)
istat=NF_PUT_VAR_REAL(idncs,idvmref,vmref)
ctext="Radar reflectivity"
istat=NF_PUT_ATT_TEXT(idnc,idvref,"long_name",len(trim(ctext)),trim(ctext))
istat=NF_PUT_ATT_TEXT(idnc,idvref,"units",3,"dbz")
ctext="Column maximum radar reflectivity"
istat=NF_PUT_ATT_TEXT(idncs,idvmref,"long_name",len(trim(ctext)),trim(ctext))
istat=NF_PUT_ATT_TEXT(idncs,idvmref,"units",3,"dbz")
!!rmiss=-9.9999e+30
!istat=NF_PUT_ATT_REAL(idnc,idvref,"_FillValue",NF_FLOAT,1,rmiss)
!istat=NF_PUT_ATT_REAL(idnc,idvref,"missing_value",NF_FLOAT,1,rmiss)

istat=NF_DEF_VAR(idncs,"Td2",NF_REAL,3,(/iddlon,iddlat,iddt/),idvtd)
if (istat /=0) then
  istat=NF_INQ_VARID(idncs,"Td2",idvtd)
end if
istat=NF_DEF_VAR(idncs,"RH2",NF_REAL,3,(/iddlon,iddlat,iddt/),idvrh)
if (istat /=0) then
  istat=NF_INQ_VARID(idncs,"RH2",idvrh)
end if

istat=NF_PUT_VAR_REAL(idncs,idvrh,vrh_s)
istat=NF_PUT_VAR_REAL(idncs,idvtd,vtd_s)
ctext="2m dew point temperature"
istat=NF_PUT_ATT_TEXT(idncs,idvtd,"long_name",len(trim(ctext)),trim(ctext))
istat=NF_PUT_ATT_TEXT(idncs,idvtd,"units",3,"K")
ctext="2m relative humidity"
istat=NF_PUT_ATT_TEXT(idncs,idvrh,"long_name",len(trim(ctext)),trim(ctext))
istat=NF_PUT_ATT_TEXT(idncs,idvrh,"units",3,"%")
!
ctext="CF-1.6"
ilen=len(trim(ctext))
istat=NF_PUT_ATT_TEXT(idnc,NF_GLOBAL,"Conventions",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idncs,NF_GLOBAL,"Conventions",ilen,trim(ctext))
if(istat.ne.0) write(*,*) "put att",NF_STRERROR(istat)
ctext='degrees_east'
ilen=len(trim(ctext))
istat=NF_INQ_VARID(idnc,"lon",idvlon)
istat=NF_INQ_VARID(idncs,"lon",idvlons)
istat=NF_PUT_ATT_TEXT(idnc,idvlon,"units",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idncs,idvlons,"units",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idnc,idvlon,"axis",1,"X")
istat=NF_PUT_ATT_TEXT(idncs,idvlons,"axis",1,"X")
if(istat.ne.0) write(*,*) "put att",NF_STRERROR(istat)
ctext='degrees_north'
ilen=len(trim(ctext))
istat=NF_INQ_VARID(idnc,"lat",idvlat)
istat=NF_INQ_VARID(idncs,"lat",idvlats)
istat=NF_PUT_ATT_TEXT(idnc,idvlat,"units",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idncs,idvlats,"units",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idnc,idvlat,"axis",1,"Y")
istat=NF_PUT_ATT_TEXT(idncs,idvlats,"axis",1,"Y")
if(istat.ne.0) write(*,*) "put att",NF_STRERROR(istat)
istat=NF_ENDDEF(idnc)
istat=NF_ENDDEF(idncs)

!istat=NF_PUT_ATT_REAL(idnc,idvlon,"actual_range",NF_FLOAT,2,(/-66.7156,-61.6683/))
!istat=NF_PUT_ATT_REAL(idnc,idvlat,"actual_range",NF_FLOAT,2,(/-33.6055,-29.2764/))

istat=NF_CLOSE(idnc)
istat=NF_CLOSE(idncs)

stop

contains 
function calc_ref(nlon,nlat,nz,nt,vd,vr,vg,vs)
real(4)::vd(nlon,nlat,nz,nt)
real(4)::vr(nlon,nlat,nz,nt)
real(4)::vg(nlon,nlat,nz,nt)
real(4)::vs(nlon,nlat,nz,nt)
real(4)::calc_ref(nlon,nlat,nz,nt)

real(4)::vdbz

real(4),parameter::vmiss=-9.9999d30
real(4),parameter::eps=1.0d-10
real(4),parameter::vref_lbnd=10.0
real(4),parameter::vref_norain=5.0
real(4),parameter::afact_r=2.53d4, bfact_r = 1.84
real(4),parameter::afact_s=3.48d3, bfact_s = 1.66
real(4),parameter::afact_g=5.54d3, bfact_g = 1.70

do it=1,nt
do iz=1,nz
do ilat=1,nlat
do ilon=1,nlon
 if (vd(ilon,ilat,iz,it) > 0.0 ) then
   vdbz                  = afact_r * ( 1.0d3* vd(ilon,ilat,iz,it) * max ( eps, vr(ilon,ilat,iz,it) ) ) ** bfact_r & 
                         + afact_g * ( 1.0d3* vd(ilon,ilat,iz,it) * max ( eps, vg(ilon,ilat,iz,it) ) ) ** bfact_g &
                         + afact_s * ( 1.0d3* vd(ilon,ilat,iz,it) * max ( eps, vs(ilon,ilat,iz,it) ) ) ** bfact_s 
   vdbz = 10.0 * log10(vdbz)
   if (vdbz < vref_lbnd ) vdbz = vref_norain
 else
   vdbz = vmiss
 end if
 calc_ref(ilon,ilat,iz,it) = vdbz
end do
end do
end do
end do
end function

subroutine calc_rh_td_2d(nlon,nlat,nt,vp,vt,vq,vrh,vtd)
real(4)::vp(nlon,nlat,nt)
real(4)::vt(nlon,nlat,nt)
real(4)::vq(nlon,nlat,nt)
real(4)::vrh(nlon,nlat,nt)
real(4)::vtd(nlon,nlat,nt)

real(4),parameter::vmiss=-9.9999d30
real(4),parameter::eps=1.0d-10

do it=1,nt
do ilat=1,nlat
do ilon=1,nlon
 if (vt(ilon,ilat,it) > 0.0 .and. vq(ilon,ilat,it) > 0.0 ) then
   call calc_rh_core(vt(ilon,ilat,it),vq(ilon,ilat,it),vp(ilon,ilat,it),vrh(ilon,ilat,it),vtd(ilon,ilat,it))
 else
   vtd(ilon,ilat,it) = vmiss
   vrh(ilon,ilat,it) = vmiss
 end if
end do
end do
end do
end subroutine 

subroutine calc_rh_core(t,q,p,rh,td)
  IMPLICIT NONE
  integer,parameter::r_size=4, RP=4
  REAL(r_size),PARAMETER :: t0=273.15d0
  REAL(r_size),INTENT(IN) :: t,q,p
  REAL(r_size),INTENT(OUT) :: rh,td
  REAL(r_size) :: e,es,tc

  REAL(r_size)::CPvap=1846.0
  REAL(r_size)::CVvap
  REAL(r_size)::Rvap=461.5
  REAL(r_size)::CL  =4218.0
  REAL(r_size)::LHV

  logical :: converged

  real(RP), parameter :: A = 17.625_RP
  real(RP), parameter :: B = 243.04_RP
  real(RP), parameter :: C = 610.94_RP
  integer,  parameter :: itelim = 100
  real(RP), parameter :: criteria = 0.1_RP**(2+RP/2)

  real(RP) :: pvap, psat
  real(RP) :: dpsat_dT
  real(RP) :: Tdi
  real(RP) :: dTdew

  integer :: ite

  e = q * p * 0.01d0 / (0.378d0 * q + 0.622d0)

  es=func_psat(t)

  rh = e/es * 100.0 

  if (rh > 100.0) then
    td=t
    return
  else

  ! first guess is calculated by Alduchov and Eskridge (1996)
  ! See Lawrence (2005) BAMS
    pvap = e  * 100.0 !!! Pa
    Td = B * log( pvap / C ) / ( A - log( pvap / C ) ) + t0
    converged = .false.
    do ite = 1, itelim
       psat=func_psat(td)*100.0

     !  write(*,*) "pvap psat",pvap,psat
  
       lhv=func_lhv(td)

       dpsat_dT = psat * lhv / ( Rvap * td**2 )
       dTdew = ( psat - pvap ) / dpsat_dT
       if ( dTdew < criteria ) then
          converged = .true.
          exit
       end if

       Td = Td - dTdew
    end do
    if( .not. converged ) then
      Tdi = B * log( pvap / C ) / ( A - log( pvap / C ) ) + t0
      write(*,*) "SATURATION not converged ",psat, pvap, T, Td, Tdi, dTdew, dpsat_dT
      stop
    endif

  end if
  RETURN
end subroutine 

  real function func_psat(tk)
  real(4)::tk,tc
  integer,parameter::r_size=4
  REAL(r_size),PARAMETER :: t0=273.15d0
  REAL(r_size),PARAMETER :: e0c=6.11d0
  REAL(r_size),PARAMETER :: al=17.3d0
  REAL(r_size),PARAMETER :: bl=237.3d0
  REAL(r_size),PARAMETER :: e0i=6.1121d0
  REAL(r_size),PARAMETER :: ai=22.587d0
  REAL(r_size),PARAMETER :: bi=273.86d0

  tc = tk-t0
  IF(tc >= 0.0d0) THEN
    func_psat = e0c * exp(al*tc/(bl+tc))
  ELSE IF(tc <= -15.d0) THEN
    func_psat = e0i * exp(ai*tc/(bi+tc))
  ELSE
    func_psat = e0c * exp(al*tc/(bl+tc)) * (15.0d0+tc)/15.0d0 &
       + e0i * exp(ai*tc/(bi+tc)) * (-tc) / 15.0d0
  END IF
    
  end function
  real function func_lhv(tk)
  real(4)::tk,tc
  REAL(4)::LHV0 =2.501E+6
  REAL(4)::CPvap=1846.0
  REAL(4)::CL  =4218.0
  REAL(4),PARAMETER :: t0=273.15d0
  tc = tk-t0
  func_lhv=lhv0 + ( CPvap - CL ) * tc 
  end function

end program
