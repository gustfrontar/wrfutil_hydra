program main
use netcdf

implicit double precision(a-h,o-z)

include "netcdf.inc"

real(4),allocatable::axlon(:)
real(4),allocatable::axlat(:)
real(4),allocatable::axz(:)

real(4),allocatable::vd(:,:,:,:)
real(4),allocatable::vqr(:,:,:,:)
real(4),allocatable::vqg(:,:,:,:)
real(4),allocatable::vqs(:,:,:,:)
real(4),allocatable::vqv(:,:,:,:)
real(4),allocatable::vqi(:,:,:,:)
real(4),allocatable::vqc(:,:,:,:)
real(4),allocatable::vu(:,:,:,:)
real(4),allocatable::vv(:,:,:,:)
real(4),allocatable::vw(:,:,:,:)
real(4),allocatable::vt(:,:,:,:)
real(4),allocatable::vr(:,:,:,:)
real(4),allocatable::vp(:,:,:,:)
real(4),allocatable::vref(:,:,:,:)
real(4),allocatable::vmref(:,:,:)
real(4),allocatable::vt2(:,:,:)
real(4),allocatable::vq2(:,:,:)
real(4),allocatable::vtd2(:,:,:)
real(4),allocatable::vr2(:,:,:)
real(4),allocatable::vu10(:,:,:)
real(4),allocatable::vv10(:,:,:)
real(4),allocatable::vps(:,:,:)
real(4),allocatable::vprec(:,:,:)

real(4)::rmiss(1)

character(len=200)::ncfile
character(len=200)::ncfiles
character(len=200)::ncfilea_3d, ncfileg_3d, ncfilea_2d, ncfileg_2d
character(len=14)::ctime
character(len=14)::ctimef
character(len=100)::ctext
character(len=200)::cdir
character(len=200)::cdiro

character(len=4)::cmem
character(len=3)::cdom

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

cmem=basename_out(indx-5:indx-2)

indx=index(trim(basename_out),'201')
ctime=basename_out(indx:indx+13)

indx_SR=index(trim(basename_out),'2km_SR')
indx_SSD=index(trim(basename_out),'2km_SSD')
if (indx_SR/=0) then
  cdom="SR"
elseif (indx_SSD/=0) then
  cdom="SSD"
else
  write(*,*) "domain not detected. stop."
  stop
end if

cdiro="./output"

read(ctime(7:14),*) itime
!read(ctime,*,iostat=ios) itime
!read(ctime(1:14),*,iostat=ios) itime
ctimef=ctime
write(ctimef(7:14),'(I8.8)') itime+500
if (ctimef(11:12)=="60") then
  ctimef(11:12)="00"
  read(ctimef(9:10),'(I2)',iostat=ios) itimef
  write(ctimef(9:10),'(I2.2)') itimef+1
  if (ctimef(9:10)=="24") then
    ctimef(9:10)="00"
    read(ctimef(7:8),'(I2)',iostat=ios) itimef
    write(ctimef(7:8),'(I2.2)') itimef+1
  end if
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

istat=NF_INQ_VARID(idnc,"DENS",idvdi)
istat=NF_INQ_VARID(idnc,"QR",idvqri)
istat=NF_INQ_VARID(idnc,"QG",idvqgi)
istat=NF_INQ_VARID(idnc,"QS",idvqsi)
istat=NF_INQ_VARID(idnc,"QC",idvqci)
istat=NF_INQ_VARID(idnc,"QI",idvqii)
istat=NF_INQ_VARID(idnc,"QV",idvqvi)
istat=NF_INQ_VARID(idnc,"Umet",idvui)
istat=NF_INQ_VARID(idnc,"Vmet",idvvi)
istat=NF_INQ_VARID(idnc,"W",idvwi)
istat=NF_INQ_VARID(idnc,"T",idvwi)
istat=NF_INQ_VARID(idnc,"PRES",idvpi)
istat=NF_INQ_VARID(idnc,"RH",idvri)

istat=NF_INQ_VARID(idncs,"T2",idvt2i)
istat=NF_INQ_VARID(idncs,"Q2",idvq2i)
istat=NF_INQ_VARID(idncs,"U10",idvu10i)
istat=NF_INQ_VARID(idncs,"V10",idvv10i)
istat=NF_INQ_VARID(idncs,"SFC_PRES",idvpsi)
istat=NF_INQ_VARID(idncs,"PREC",idvpreci)

istat=NF_INQ_DIMID(idnc,"lat",iddlati)
istat=NF_INQ_DIMLEN(idnc,iddlati,nlat)
istat=NF_INQ_DIMID(idnc,"lon",iddloni)
istat=NF_INQ_DIMLEN(idnc,iddloni,nlon)
istat=NF_INQ_DIMID(idnc,"z",iddzi)
istat=NF_INQ_DIMLEN(idnc,iddzi,nz)
istat=NF_INQ_DIMID(idnc,"time",iddti)
istat=NF_INQ_DIMLEN(idnc,iddti,nt)

allocate(axlon(nlon))
allocate(axlat(nlat))
allocate(axz(nz))
istat=NF_INQ_VARID(idnc,"lat",idvlati)
istat=NF_INQ_VARID(idnc,"lon",idvloni)
istat=NF_INQ_VARID(idnc,"z",idvzi)
istat=NF_GET_VAR_REAL(idnc,idvlati,axlat)
istat=NF_GET_VAR_REAL(idnc,idvloni,axlon)
istat=NF_GET_VAR_REAL(idnc,idvzi,axz)

allocate(vd(nlon,nlat,nz,nt))
allocate(vqr(nlon,nlat,nz,nt))
allocate(vqg(nlon,nlat,nz,nt))
allocate(vqs(nlon,nlat,nz,nt))
allocate(vqv(nlon,nlat,nz,nt))
allocate(vqi(nlon,nlat,nz,nt))
allocate(vqc(nlon,nlat,nz,nt))
allocate(vref(nlon,nlat,nz,nt))
allocate(vu(nlon,nlat,nz,nt))
allocate(vv(nlon,nlat,nz,nt))
allocate(vw(nlon,nlat,nz,nt))
allocate(vt(nlon,nlat,nz,nt))
allocate(vp(nlon,nlat,nz,nt))
allocate(vr(nlon,nlat,nz,nt))
allocate(vmref(nlon,nlat,nt))
allocate(vt2(nlon,nlat,nt))
allocate(vq2(nlon,nlat,nt))
allocate(vtd2(nlon,nlat,nt))
allocate(vr2(nlon,nlat,nt))
allocate(vu10(nlon,nlat,nt))
allocate(vv10(nlon,nlat,nt))
allocate(vps(nlon,nlat,nt))
allocate(vprec(nlon,nlat,nt))

ncfilea_2d=trim(cdiro)//"/scale_2d_"//ctime(1:12)//"_"//trim(cdom)//"_anal_M"//cmem//".nc"
ncfileg_2d=trim(cdiro)//"/scale_2d_"//ctimef(1:12)//"_"//trim(cdom)//"_gues_M"//cmem//".nc"
ncfilea_3d=trim(cdiro)//"/scale_3d_"//ctime(1:12)//"_"//trim(cdom)//"_anal_M"//cmem//".nc"
ncfileg_3d=trim(cdiro)//"/scale_3d_"//ctimef(1:12)//"_"//trim(cdom)//"_gues_M"//cmem//".nc"
istat=NF_CREATE(trim(ncfilea_2d),NF_CLOBBER,idnca2)
istat=NF_CREATE(trim(ncfileg_2d),NF_CLOBBER,idncg2)
istat=NF_CREATE(trim(ncfilea_3d),NF_CLOBBER,idnca3)
istat=NF_CREATE(trim(ncfileg_3d),NF_CLOBBER,idncg3)

istat=NF_DEF_DIM(idnca3,"lon",nlon,iddlon)
istat=NF_DEF_DIM(idnca3,"lat",nlat,iddlat)
istat=NF_DEF_DIM(idnca3,"z",nz,iddz)
istat=NF_DEF_VAR(idnca3,"lon",NF_FLOAT,1,(/iddlon/),idvlono3)
istat=NF_DEF_VAR(idnca3,"lat",NF_FLOAT,1,(/iddlat/),idvlato3)
istat=NF_DEF_VAR(idnca3,"z",NF_FLOAT,1,(/iddz/),idvzo3)
istat=NF_DEF_DIM(idncg3,"lon",nlon,iddlon)
istat=NF_DEF_DIM(idncg3,"lat",nlat,iddlat)
istat=NF_DEF_DIM(idncg3,"z",nz,iddz)
istat=NF_DEF_VAR(idncg3,"lon",NF_FLOAT,1,(/iddlon/),idvlono3)
istat=NF_DEF_VAR(idncg3,"lat",NF_FLOAT,1,(/iddlat/),idvlato3)
istat=NF_DEF_VAR(idncg3,"z",NF_FLOAT,1,(/iddz/),idvzo3)

ctext='degrees_east'
ilen=len(trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca3,idvlono3,"units",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idncg3,idvlono3,"units",ilen,trim(ctext))
ctext='longitude'
ilen=len(trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca3,idvlono3,"long_name",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idncg3,idvlono3,"long_name",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca3,idvlono3,"axis",1,"X")
istat=NF_PUT_ATT_TEXT(idncg3,idvlono3,"axis",1,"X")
ctext='degrees_north'
ilen=len(trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca3,idvlato3,"units",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idncg3,idvlato3,"units",ilen,trim(ctext))
ctext='latitude'
ilen=len(trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca3,idvlato3,"long_name",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idncg3,idvlato3,"long_name",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca3,idvlato3,"axis",1,"Y")
istat=NF_PUT_ATT_TEXT(idncg3,idvlato3,"axis",1,"Y")

istat=NF_PUT_ATT_TEXT(idnca3,idvzo3,"long_name",1,"Z")
istat=NF_PUT_ATT_TEXT(idnca3,idvzo3,"units",1,"m")
istat=NF_PUT_ATT_TEXT(idncg3,idvzo3,"long_name",1,"Z")
istat=NF_PUT_ATT_TEXT(idncg3,idvzo3,"units",1,"m")

istat=NF_DEF_VAR(idnca3,"T",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvto3)
istat=NF_DEF_VAR(idnca3,"RH",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvro3)
istat=NF_DEF_VAR(idnca3,"U",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvuo3)
istat=NF_DEF_VAR(idnca3,"V",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvvo3)
istat=NF_DEF_VAR(idnca3,"W",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvwo3)
istat=NF_DEF_VAR(idnca3,"PRES",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvpo3)
istat=NF_DEF_VAR(idnca3,"QV",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvqvo3)
istat=NF_DEF_VAR(idnca3,"QR",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvqro3)
istat=NF_DEF_VAR(idnca3,"QS",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvqso3)
istat=NF_DEF_VAR(idnca3,"QG",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvqgo3)
istat=NF_DEF_VAR(idnca3,"QC",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvqco3)
istat=NF_DEF_VAR(idnca3,"QI",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvqio3)
istat=NF_DEF_VAR(idnca3,"REF",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvrefo3)
istat=NF_DEF_VAR(idncg3,"T",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvto3)
istat=NF_DEF_VAR(idncg3,"RH",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvro3)
istat=NF_DEF_VAR(idncg3,"U",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvuo3)
istat=NF_DEF_VAR(idncg3,"V",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvvo3)
istat=NF_DEF_VAR(idncg3,"W",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvwo3)
istat=NF_DEF_VAR(idncg3,"PRES",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvpo3)
istat=NF_DEF_VAR(idncg3,"QV",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvqvo3)
istat=NF_DEF_VAR(idncg3,"QR",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvqro3)
istat=NF_DEF_VAR(idncg3,"QS",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvqso3)
istat=NF_DEF_VAR(idncg3,"QG",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvqgo3)
istat=NF_DEF_VAR(idncg3,"QC",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvqco3)
istat=NF_DEF_VAR(idncg3,"QI",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvqio3)
istat=NF_DEF_VAR(idncg3,"REF",NF_FLOAT,3,(/iddlon,iddlat,iddz/),idvrefo3)

ctext="Radar reflectivity"
istat=NF_PUT_ATT_TEXT(idnca3,idvrefo3,"long_name",len(trim(ctext)),trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca3,idvrefo3,"units",3,"dbz")
istat=NF_PUT_ATT_TEXT(idncg3,idvrefo3,"long_name",len(trim(ctext)),trim(ctext))
istat=NF_PUT_ATT_TEXT(idncg3,idvrefo3,"units",3,"dbz")

istat=NF_DEF_DIM(idnca2,"lon",nlon,iddlon)
istat=NF_DEF_DIM(idnca2,"lat",nlat,iddlat)
istat=NF_DEF_VAR(idnca2,"lon",NF_FLOAT,1,(/iddlon/),idvlono2)
istat=NF_DEF_VAR(idnca2,"lat",NF_FLOAT,1,(/iddlat/),idvlato2)
istat=NF_DEF_DIM(idncg2,"lon",nlon,iddlon)
istat=NF_DEF_DIM(idncg2,"lat",nlat,iddlat)
istat=NF_DEF_VAR(idncg2,"lon",NF_FLOAT,1,(/iddlon/),idvlono2)
istat=NF_DEF_VAR(idncg2,"lat",NF_FLOAT,1,(/iddlat/),idvlato2)
ctext='degrees_east'
ilen=len(trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca2,idvlono2,"units",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idncg2,idvlono2,"units",ilen,trim(ctext))
ctext='longitude'
ilen=len(trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca2,idvlono2,"long_name",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idncg2,idvlono2,"long_name",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca2,idvlono2,"axis",1,"X")
istat=NF_PUT_ATT_TEXT(idncg2,idvlono2,"axis",1,"X")
ctext='degrees_north'
ilen=len(trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca2,idvlato2,"units",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idncg2,idvlato2,"units",ilen,trim(ctext))
ctext='latitude'
ilen=len(trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca2,idvlato2,"long_name",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idncg2,idvlato2,"long_name",ilen,trim(ctext))
istat=NF_PUT_ATT_TEXT(idnca2,idvlato2,"axis",1,"Y")
istat=NF_PUT_ATT_TEXT(idncg2,idvlato2,"axis",1,"Y")

istat=NF_DEF_VAR(idnca2,"T2",NF_FLOAT,2,(/iddlon,iddlat/),idvt2o2)
istat=NF_DEF_VAR(idnca2,"Td2",NF_FLOAT,2,(/iddlon,iddlat/),idvtd2o2)
istat=NF_DEF_VAR(idnca2,"RH2",NF_FLOAT,2,(/iddlon,iddlat/),idvr2o2)
istat=NF_DEF_VAR(idnca2,"U10",NF_FLOAT,2,(/iddlon,iddlat/),idvu10o2)
istat=NF_DEF_VAR(idnca2,"V10",NF_FLOAT,2,(/iddlon,iddlat/),idvv10o2)
istat=NF_DEF_VAR(idnca2,"PS",NF_FLOAT,2,(/iddlon,iddlat/),idvpso2)
istat=NF_DEF_VAR(idnca2,"PREC",NF_FLOAT,2,(/iddlon,iddlat/),idvpreco2)
istat=NF_DEF_VAR(idnca2,"MREF",NF_FLOAT,2,(/iddlon,iddlat/),idvmrefo2)
istat=NF_DEF_VAR(idncg2,"T2",NF_FLOAT,2,(/iddlon,iddlat/),idvt2o2)
istat=NF_DEF_VAR(idncg2,"Td2",NF_FLOAT,2,(/iddlon,iddlat/),idvtd2o2)
istat=NF_DEF_VAR(idncg2,"RH2",NF_FLOAT,2,(/iddlon,iddlat/),idvr2o2)
istat=NF_DEF_VAR(idncg2,"U10",NF_FLOAT,2,(/iddlon,iddlat/),idvu10o2)
istat=NF_DEF_VAR(idncg2,"V10",NF_FLOAT,2,(/iddlon,iddlat/),idvv10o2)
istat=NF_DEF_VAR(idncg2,"PS",NF_FLOAT,2,(/iddlon,iddlat/),idvpso2)
istat=NF_DEF_VAR(idncg2,"PREC",NF_FLOAT,2,(/iddlon,iddlat/),idvpreco2)
istat=NF_DEF_VAR(idncg2,"MREF",NF_FLOAT,2,(/iddlon,iddlat/),idvmrefo2)

istat=NF_ENDDEF(idnca3)
istat=NF_ENDDEF(idncg3)
istat=NF_ENDDEF(idnca2)
istat=NF_ENDDEF(idncg2)


!!!
istat=NF_GET_VAR_REAL(idnc,idvdi,vd)
istat=NF_GET_VAR_REAL(idnc,idvqri,vqr)
istat=NF_GET_VAR_REAL(idnc,idvqgi,vqg)
istat=NF_GET_VAR_REAL(idnc,idvqsi,vqs)
istat=NF_GET_VAR_REAL(idnc,idvqci,vqc)
istat=NF_GET_VAR_REAL(idnc,idvqii,vqi)
istat=NF_GET_VAR_REAL(idnc,idvqvi,vqv)

istat=NF_GET_VAR_REAL(idnc,idvvi,vv)
istat=NF_GET_VAR_REAL(idnc,idvui,vu)
istat=NF_GET_VAR_REAL(idnc,idvwi,vw)
istat=NF_GET_VAR_REAL(idnc,idvti,vt)
istat=NF_GET_VAR_REAL(idnc,idvpi,vp)
istat=NF_GET_VAR_REAL(idnc,idvri,vr)

write(*,*) "dens",maxval(vd),minval(vd)
write(*,*) "rain",maxval(vqr),minval(vqr)
write(*,*) "grap",maxval(vqg),minval(vqg)
write(*,*) "snow",maxval(vqs),minval(vqs)

vref=calc_ref(nlon,nlat,nz,nt,vd,vqr,vqg,vqs)
vmref=maxval(vref,3)
write(*,*) "ref ",maxval(vref),minval(vref)
write(*,*) "mref ",maxval(vmref),minval(vmref)

istat=NF_GET_VAR_REAL(idncs,idvt2i,vt2)
istat=NF_GET_VAR_REAL(idncs,idvq2i,vq2)
istat=NF_GET_VAR_REAL(idncs,idvu10i,vu10)
istat=NF_GET_VAR_REAL(idncs,idvv10i,vv10)
istat=NF_GET_VAR_REAL(idncs,idvpsi,vps)
istat=NF_GET_VAR_REAL(idncs,idvpreci,vprec)

call calc_rh_td_2d(nlon,nlat,nt,vps,vt2,vq2,vr2,vtd2)
!
!ctext="Column maximum radar reflectivity"
!istat=NF_PUT_ATT_TEXT(idncs,idvmref,"long_name",len(trim(ctext)),trim(ctext))
!istat=NF_PUT_ATT_TEXT(idncs,idvmref,"units",3,"dbz")
!!rmiss=-9.9999e+30
!istat=NF_PUT_ATT_REAL(idnc,idvref,"_FillValue",NF_FLOAT,1,rmiss)
!istat=NF_PUT_ATT_REAL(idnc,idvref,"missing_value",NF_FLOAT,1,rmiss)

istat=NF_PUT_VAR_REAL(idnca3,idvlono3,axlon)
istat=NF_PUT_VAR_REAL(idnca3,idvlato3,axlat)
istat=NF_PUT_VAR_REAL(idnca3,idvzo3,axz)
istat=NF_PUT_VAR_REAL(idncg3,idvlono3,axlon)
istat=NF_PUT_VAR_REAL(idncg3,idvlato3,axlat)
istat=NF_PUT_VAR_REAL(idncg3,idvzo3,axz)

istat=NF_PUT_VAR_REAL(idnca3,idvuo3,vu(:,:,:,1))
istat=NF_PUT_VAR_REAL(idnca3,idvvo3,vv(:,:,:,1))
istat=NF_PUT_VAR_REAL(idnca3,idvwo3,vw(:,:,:,1))
istat=NF_PUT_VAR_REAL(idnca3,idvto3,vt(:,:,:,1))
istat=NF_PUT_VAR_REAL(idnca3,idvro3,vr(:,:,:,1))
istat=NF_PUT_VAR_REAL(idnca3,idvpo3,vp(:,:,:,1))
istat=NF_PUT_VAR_REAL(idnca3,idvqvo3,vqv(:,:,:,1))
istat=NF_PUT_VAR_REAL(idnca3,idvqro3,vqr(:,:,:,1))
istat=NF_PUT_VAR_REAL(idnca3,idvqso3,vqs(:,:,:,1))
istat=NF_PUT_VAR_REAL(idnca3,idvqgo3,vqg(:,:,:,1))
istat=NF_PUT_VAR_REAL(idnca3,idvqco3,vqc(:,:,:,1))
istat=NF_PUT_VAR_REAL(idnca3,idvqio3,vqi(:,:,:,1))
istat=NF_PUT_VAR_REAL(idnca3,idvrefo3,vref(:,:,:,1))

istat=NF_PUT_VAR_REAL(idncg3,idvuo3,vu(:,:,:,2))
istat=NF_PUT_VAR_REAL(idncg3,idvvo3,vv(:,:,:,2))
istat=NF_PUT_VAR_REAL(idncg3,idvwo3,vw(:,:,:,2))
istat=NF_PUT_VAR_REAL(idncg3,idvto3,vt(:,:,:,2))
istat=NF_PUT_VAR_REAL(idncg3,idvro3,vr(:,:,:,2))
istat=NF_PUT_VAR_REAL(idncg3,idvpo3,vp(:,:,:,2))
istat=NF_PUT_VAR_REAL(idncg3,idvqvo3,vqv(:,:,:,2))
istat=NF_PUT_VAR_REAL(idncg3,idvqro3,vqr(:,:,:,2))
istat=NF_PUT_VAR_REAL(idncg3,idvqso3,vqs(:,:,:,2))
istat=NF_PUT_VAR_REAL(idncg3,idvqgo3,vqg(:,:,:,2))
istat=NF_PUT_VAR_REAL(idncg3,idvqco3,vqc(:,:,:,2))
istat=NF_PUT_VAR_REAL(idncg3,idvqio3,vqi(:,:,:,2))
istat=NF_PUT_VAR_REAL(idncg3,idvrefo3,vref(:,:,:,2))

istat=NF_PUT_VAR_REAL(idnca2,idvlono2,axlon)
istat=NF_PUT_VAR_REAL(idnca2,idvlato2,axlat)
istat=NF_PUT_VAR_REAL(idncg2,idvlono2,axlon)
istat=NF_PUT_VAR_REAL(idncg2,idvlato2,axlat)

istat=NF_PUT_VAR_REAL(idnca2,idvt2o2,vt2(:,:,1))
istat=NF_PUT_VAR_REAL(idnca2,idvtd2o2,vtd2(:,:,1))
istat=NF_PUT_VAR_REAL(idnca2,idvr2o2,vr2(:,:,1))
istat=NF_PUT_VAR_REAL(idnca2,idvu10o2,vu10(:,:,1))
istat=NF_PUT_VAR_REAL(idnca2,idvv10o2,vv10(:,:,1))
istat=NF_PUT_VAR_REAL(idnca2,idvpso2,vps(:,:,1))
istat=NF_PUT_VAR_REAL(idnca2,idvpreco2,vprec(:,:,1))
istat=NF_PUT_VAR_REAL(idnca2,idvmrefo2,vmref(:,:,1))

istat=NF_PUT_VAR_REAL(idncg2,idvt2o2,vt2(:,:,2))
istat=NF_PUT_VAR_REAL(idncg2,idvtd2o2,vtd2(:,:,2))
istat=NF_PUT_VAR_REAL(idncg2,idvr2o2,vr2(:,:,2))
istat=NF_PUT_VAR_REAL(idncg2,idvu10o2,vu10(:,:,2))
istat=NF_PUT_VAR_REAL(idncg2,idvv10o2,vv10(:,:,2))
istat=NF_PUT_VAR_REAL(idncg2,idvpso2,vps(:,:,2))
istat=NF_PUT_VAR_REAL(idncg2,idvpreco2,vprec(:,:,2))
istat=NF_PUT_VAR_REAL(idncg2,idvmrefo2,vmref(:,:,2))

istat=NF_CLOSE(idnca3)
istat=NF_CLOSE(idncg3)
istat=NF_CLOSE(idnca2)
istat=NF_CLOSE(idncg2)

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
