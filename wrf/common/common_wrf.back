MODULE common_wrf
!=======================================================================
!
! [PURPOSE:] Common Information for WRF
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   07/14/2010 Takemasa Miyoshi  modified for WRF
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE map_utils
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------

  !MODULE VARIABLES (NOT READ FROM NAMELIST) 
  INTEGER, PARAMETER :: nv3d=15
  INTEGER, PARAMETER :: nv2d=5
  INTEGER, PARAMETER :: np2d=3
  INTEGER, PARAMETER :: nid_obs=10

  REAL(r_sngl) :: tmpatt
  INTEGER,SAVE :: dx,dy ! grid spacing [m]
  REAL(r_size),SAVE :: pdx,pdy
  INTEGER,SAVE :: nlon
  INTEGER,SAVE :: nlat
  INTEGER,SAVE :: nlev
  !3D VARIABLES
  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_w=3
  INTEGER,PARAMETER :: iv3d_t=4
  INTEGER,PARAMETER :: iv3d_p=5
  INTEGER,PARAMETER :: iv3d_ph=6
  INTEGER,PARAMETER :: iv3d_qv=7
  INTEGER,PARAMETER :: iv3d_qc=8
  INTEGER,PARAMETER :: iv3d_qr=9
  INTEGER,PARAMETER :: iv3d_qci=10
  INTEGER,PARAMETER :: iv3d_qs=11
  INTEGER,PARAMETER :: iv3d_qg=12
  !3D CHEM VARS
  INTEGER,PARAMETER :: iv3d_coant=13
  INTEGER,PARAMETER :: iv3d_cobck=14
  INTEGER,PARAMETER :: iv3d_cobbu=15
  !2D VARIABLES
  INTEGER,PARAMETER :: iv2d_ps=1
  INTEGER,PARAMETER :: iv2d_t2=2
  INTEGER,PARAMETER :: iv2d_q2=3
  INTEGER,PARAMETER :: iv2d_u10=4
  INTEGER,PARAMETER :: iv2d_v10=5
  !SPECIAL VARIABLES
  INTEGER,PARAMETER :: iv3d_tv=100
  INTEGER,PARAMETER :: iv2d_mu=200
  !2D PARAMETERS
  INTEGER,PARAMETER :: ip2d_hfx=1
  INTEGER,PARAMETER :: ip2d_qfx=2
  INTEGER,PARAMETER :: ip2d_ust=3

  REAL(4),SAVE      :: p0 != 1.0e+5
  REAL(4),PARAMETER :: t0 = 300.0
  REAL(r_size)      :: p_top
  INTEGER(4),PARAMETER :: imt_grd = 50
  INTEGER(4),PARAMETER :: spec_bdy_width = 5
  INTEGER,SAVE      :: nij0
  INTEGER,SAVE      :: nlevall
  INTEGER,SAVE      :: ngpv
  REAL(r_size),ALLOCATABLE,SAVE :: lon(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lat(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lonu(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: latu(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lonv(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: latv(:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: fcori(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: phi0(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: sinalpha(:,:),cosalpha(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: landmask(:,:)
  CHARACTER(20),SAVE :: element(nv3d+nv2d+np2d)
  REAL(r_sngl),ALLOCATABLE,SAVE :: dnw(:)
  REAL(r_sngl),ALLOCATABLE,SAVE :: znu(:)
  REAL(r_size) :: plat1, plon1, pstdlon , ptruelat1 , ptruelat2
  INTEGER      :: pcode

  CHARACTER(19) :: hdate !Date in input files.

  LOGICAL, PARAMETER :: OPTIONAL_VARIABLE(nv3d+nv2d+np2d)= &  !Some input variables are optional (eg. condensates)
  &(/ .false. , .false. , .true. , .false. , .false. , .false. , .true. , .true. , &
  !   U          V         W         T         P         PH      QV        QC        
  & .true. , .true. , .true. , .true. , .true. , .true. , .true. , .false. , .true. , .true. , .true. , .true. , &
  !   QR       QCI      QS       QG     CO_ANT  CO_BCK    CO_BBU     PS        T2M     Q2M       U10M    V10M
  & .true. , .true. , .true.  /)
  !  HFX_F    QFX_F   UST_F 

  LOGICAL :: PRESENT_VARIABLE(nv3d+nv2d+np2d) = .false.  !For checking input variables.

  TYPE(proj_info)  :: projection

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_wrf(inputfile)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  REAL(r_sngl),ALLOCATABLE :: buf4(:,:)
  REAL(r_sngl) :: buf41d
  INTEGER :: i
  INTEGER(4) :: ncid,varid,dimids(3),shape(3)
  INTEGER(4),ALLOCATABLE :: start(:),count(:)
  CHARACTER(*)           :: inputfile
  INTEGER :: ierr
  


  WRITE(6,'(A)') 'Hello from set_common_wrf'

  !
  ! Elements
  !
  element(iv3d_u)  = 'U   '
  element(iv3d_v)  = 'V   '
  element(iv3d_w)  = 'W   '
  element(iv3d_t)  = 'T   '
  element(iv3d_p)  = 'P   '
  element(iv3d_ph) = 'PH   '
  element(iv3d_qv) = 'QVAPOR'
  element(iv3d_qc) = 'QCLOUD'
  element(iv3d_qr) = 'QRAIN'
  element(iv3d_qci) ='QICE'
  element(iv3d_qg) = 'QGRAUP'
  element(iv3d_qs) = 'QSNOW'
  element(iv3d_coant) = 'CO_ANT'
  element(iv3d_cobck) = 'CO_BCK'
  element(iv3d_cobbu) = 'CO_BBU'
  element(nv3d+iv2d_ps)   = 'PSFC'
  element(nv3d+iv2d_t2)   = 'T2  '
  element(nv3d+iv2d_q2)   = 'Q2  '
  element(nv3d+iv2d_u10)   = 'U10  '
  element(nv3d+iv2d_v10)   = 'V10  '
  element(nv3d+nv2d+ip2d_hfx)   = 'HFX_FACTOR'
  element(nv3d+nv2d+ip2d_qfx)   = 'QFX_FACTOR'
  element(nv3d+nv2d+ip2d_ust)   = 'UST_FACTOR'
  

  !
  ! Lon, Lat, F, phi0
  !
  WRITE(6,'(A)') 'READING: ',inputfile
  CALL check_io(NF_OPEN(inputfile,NF_NOWRITE,ncid))
  CALL check_io(NF_INQ_VARID(ncid,'U',varid))
  CALL check_io(NF_INQ_VARDIMID(ncid,varid,dimids))
  DO i=1,3
    CALL check_io(NF_INQ_DIMLEN(ncid,dimids(i),shape(i)))
  END DO
  nlon = shape(1)
  nlat = shape(2) + 1
  nlev = shape(3) + 1
  WRITE(6,'(A)') '*** grid information ***'
  WRITE(6,'(3(2X,A,I5))') 'nlon =',nlon,'nlat =',nlat,'nlev =',nlev
  nij0=nlon*nlat
  nlevall=nlev*nv3d+nv2d
  ngpv=nij0*nlevall

  CALL check_io(NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DX',tmpatt))
   dx=INT(tmpatt)
   pdx=REAL(tmpatt,r_size)
  CALL check_io(NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DY',tmpatt))
   dy=INT(tmpatt)
   pdy=REAL(tmpatt,r_size)

  WRITE(6,*) '*** MODEL RESOLUTION IS = ',dx

  ALLOCATE(buf4(nlon,nlat))
  ALLOCATE(lon(nlon,nlat))
  ALLOCATE(lat(nlon,nlat))
  ALLOCATE(lonu(nlon,nlat))
  ALLOCATE(latu(nlon,nlat))
  ALLOCATE(lonv(nlon,nlat))
  ALLOCATE(latv(nlon,nlat))
!  ALLOCATE(fcori(nlon,nlat))
  ALLOCATE(phi0(nlon,nlat))
  ALLOCATE(landmask(nlon,nlat))
  ALLOCATE(cosalpha(nlon,nlat))
  ALLOCATE(sinalpha(nlon,nlat))

  !!! Time
  CALL check_io(NF_INQ_VARID(ncid,'Times',varid))
  ALLOCATE(start(3),count(3))
  start = (/ 1,1,1 /)
  count = (/ 19,1,1 /)
  CALL check_io(NF_GET_VARA_TEXT(ncid,varid,start(1:2),count(1:2),hdate))
  !!! P_top
  CALL check_io(NF_INQ_VARID(ncid,'P_TOP',varid))
  start = (/ 1,1,1 /)
  count = (/ 1,1,1 /)
  CALL  check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1,1)))
  p_top=REAL(buf4(1,1),r_size)

  start = (/ 1,1,1 /)
  !!! XLONG
  count = (/ nlon-1,nlat-1,1 /)
  lon=0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'XLONG',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon-1,1:nlat-1)))
  lon(1:nlon-1,1:nlat-1)=REAL(buf4(1:nlon-1,1:nlat-1),r_size)  
  !!! XLAT
  count = (/ nlon-1,nlat-1,1 /)
  lat = 0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'XLAT',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon-1,1:nlat-1)))
  lat(1:nlon-1,1:nlat-1)=REAL(buf4(1:nlon-1,1:nlat-1),r_size)
  !!! XLONG_U
  count = (/ nlon,nlat-1,1 /)
  lonu = 0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'XLONG_U',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon,1:nlat-1)))
  lonu(1:nlon,1:nlat-1)=REAL(buf4(1:nlon,1:nlat-1),r_size)
  !!! XLAT_U
  count = (/ nlon,nlat-1,1 /)
  latu= 0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'XLAT_U',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon,1:nlat-1)))
  latu(1:nlon,1:nlat-1)=REAL(buf4(1:nlon,1:nlat-1),r_size)
  !!! XLONG_V
  count = (/ nlon-1,nlat,1 /)
  lonv = 0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'XLONG_V',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon-1,1:nlat)))
  lonv(1:nlon-1,1:nlat)=REAL(buf4(1:nlon-1,1:nlat),r_size)
  !!! XLAT_V
  count = (/ nlon-1,nlat,1 /)
  latv = 0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'XLAT_V',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon-1,1:nlat)))
  latv(1:nlon-1,1:nlat)=REAL(buf4(1:nlon-1,1:nlat),r_size)
  !!! HGT
  count = (/ nlon-1,nlat-1,1 /)
  phi0 = 0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'HGT',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon-1,1:nlat-1)))
  phi0(1:nlon-1,1:nlat-1)=REAL(buf4(1:nlon-1,1:nlat-1),r_size)
  !!! LANDMASK
  count = (/ nlon-1,nlat-1,1 /)
  landmask = 0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'LANDMASK',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon-1,1:nlat-1)))
  landmask(1:nlon-1,1:nlat-1)=REAL(buf4(1:nlon-1,1:nlat-1),r_size)
  count = (/ nlon-1,nlat-1,1 /)
  !!! COSINE OF WIND ROATATION
  cosalpha = 0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'COSALPHA',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon-1,1:nlat-1)))
  cosalpha(1:nlon-1,1:nlat-1)=REAL(buf4(1:nlon-1,1:nlat-1),r_size)
  count = (/ nlon-1,nlat-1,1 /)
  !!! SINE OF WIND ROTATION
  sinalpha = 0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'SINALPHA',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon-1,1:nlat-1)))
  sinalpha(1:nlon-1,1:nlat-1)=REAL(buf4(1:nlon-1,1:nlat-1),r_size)
  deallocate(start,count)

  !!! P00
  CALL check_io(NF_INQ_VARID(ncid,'P00',varid))
  CALL check_io(NF_GET_VAR_REAL(ncid,varid,p0))


  IF( p0 .eq. 0.0d0 )THEN
    WRITE(6,*)'WARNING : p0 is 0, probably leading with an ideal case a value of 1e5 will be assumed'
    p0=1e5
  ENDIF
  !!! DNW
  ALLOCATE(start(2),count(2),dnw(nlev),znu(nlev))
  start = (/ 1,1 /)
  count = (/ nlev-1,1 /)
  dnw(1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'DNW',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,dnw(1:nlev-1)))
  znu(1:nlev) = 0.0
  CALL check_io(NF_INQ_VARID(ncid,'ZNU',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,znu(1:nlev-1)))
       

  DEALLOCATE(start,count)


  WRITE(6,'(a)') 'DOMAIN CORNERS'
  WRITE(6,'(A4,F7.2,A4,F7.2)')"LON",lon(1,1),"LAT",lat(1,1)
  WRITE(6,'(A4,F7.2,A4,F7.2)')"LON",lon(nlon-1,1),"LAT",lat(nlon-1,1)
  WRITE(6,'(A4,F7.2,A4,F7.2)')"LON",lon(1,nlat-1),"LAT",lat(1,nlat-1)
  WRITE(6,'(A4,F7.2,A4,F7.2)')"LON",lon(nlon-1,nlat-1),"LAT",lat(nlon-1,nlat-1)
  WRITE(6,'(A20,A19)')"THE DATE IS ",hdate

  !Lets load the projection information for fast ll_to_ij transformations using
  !WPS original routines.
  CALL check_io(NF_GET_ATT_INT(ncid,NF_GLOBAL,'MAP_PROJ',pcode))
  WRITE(6,*)'MAP_PROJ=',pcode
  CALL check_io(NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT1',tmpatt))
  ptruelat1=REAL(tmpatt,r_size)
  WRITE(6,*)'TRUELAT1=',ptruelat1

  CALL check_io(NF_GET_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT2',tmpatt))
  ptruelat2=REAL(tmpatt,r_size)
  WRITE(6,*)'TRUELAT2=',ptruelat2

  CALL check_io(NF_GET_ATT_REAL(ncid,NF_GLOBAL,'STAND_LON',tmpatt))
  pstdlon=REAL(tmpatt,r_size)
  WRITE(6,*)'STDLON=',pstdlon

  plat1=lat(1,1)
  plon1=lon(1,1)
  WRITE(6,*)'SW LON=',plon1,' SW LAT=',plat1
   CALL map_set(pcode, projection, lat1=plat1, lon1=plon1, knowni=1.0d0,   &
   &             knownj=1.0d0, dx=pdx, dy=pdy, stdlon=pstdlon, truelat1=ptruelat1,   &
   &             truelat2=ptruelat2, r_earth=re)

  WRITE(6,'(a)') '***  END  : READ NETCDF (CNST) ***'

  !Check input variables
  DO i = 1,nv3d+nv2d+np2d
      !Check if the variable is present.    
      ierr=NF_INQ_VARID(ncid,TRIM(element(i)),varid)
      IF( ierr == nf_noerr )THEN
       PRESENT_VARIABLE(i) = .true.
       WRITE(6,*)"Variable ",element(i)," is present in input file."
      ELSE
       if( OPTIONAL_VARIABLE(i) )THEN
        WRITE(6,*)"Warning: Variable ",element(i)," is not present in input file."
       else
        WRITE(6,*)"Error: Non optioanl variable ",element(i)," is not present in input file -> Aborting execution"
        STOP
       endif
      ENDIF
  ENDDO


  CALL check_io(NF_CLOSE(ncid))



  RETURN
END SUBROUTINE set_common_wrf



!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------
SUBROUTINE read_date_wrf(inputfile,iyyyy,imm,idd,ihh,imn,ss)
IMPLICIT NONE
character(*) , intent(in)  :: inputfile
integer , intent(out) :: iyyyy , imm , idd , ihh , imn
real(r_size) , intent(out) :: ss
integer :: start(2) , count(2) , ncid
character(19) :: tmpdate
integer :: varid , tmpsec
INCLUDE 'netcdf.inc'

  CALL open_wrf_file(inputfile,'ro',ncid)
  CALL check_io(NF_INQ_VARID(ncid,'Times',varid))

  start = (/ 1 ,1 /)
  count = (/ 19,1 /)
  CALL check_io(NF_GET_VARA_TEXT(ncid,varid,start,count,tmpdate))
  READ(tmpdate(1:4  ),'(I4)' )iyyyy
  READ(tmpdate(6:7  ),'(I2)' )imm
  READ(tmpdate(9:10 ),'(I2)' )idd
  READ(tmpdate(12:13),'(I2)' )ihh
  READ(tmpdate(15:16),'(I2)' )imn
  READ(tmpdate(18:19),'(I2)')tmpsec
  ss=REAL(tmpsec,r_size)

  CALL close_wrf_file(ncid)

END SUBROUTINE read_date_wrf

SUBROUTINE write_date_wrf(inputfile,iyyyy,imm,idd,ihh,imn,ss)
IMPLICIT NONE
character(*) , intent(in)  :: inputfile
integer , intent(in) :: iyyyy , imm , idd , ihh , imn
real(r_size) , intent(in) :: ss
integer :: start(2) , count(2) , ncid
character(19) :: tmpdate
integer :: varid
INCLUDE 'netcdf.inc'

  CALL open_wrf_file(inputfile,'rw',ncid)
  CALL check_io(NF_INQ_VARID(ncid,'Times',varid))
  start = (/ 1 ,1 /)
  count = (/ 19,1 /)
  !Read the current date in file to get the date format.
  CALL check_io(NF_GET_VARA_TEXT(ncid,varid,start,count,tmpdate))

  WRITE(tmpdate(1:4  ),'(I4.4)' )iyyyy
  WRITE(tmpdate(6:7  ),'(I2.2)' )imm
  WRITE(tmpdate(9:10 ),'(I2.2)' )idd
  WRITE(tmpdate(12:13),'(I2.2)' )ihh
  WRITE(tmpdate(15:16),'(I2.2)' )imn
  WRITE(tmpdate(18:19),'(I2.2)')INT(ss)

  !Write the new date in the file.
  CALL check_io(NF_PUT_VARA_TEXT(ncid,varid,start,count,tmpdate))

  CALL check_io(NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'START_DATE',19,tmpdate))
  CALL check_io(NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'SIMULATION_START_DATE',19,tmpdate))

  CALL close_wrf_file(ncid)

END SUBROUTINE write_date_wrf

SUBROUTINE read_var_wrf(ncid,iv,slev,elev,fieldg,flag)
IMPLICIT NONE
INCLUDE 'netcdf.inc'
INTEGER(4), INTENT(IN) :: ncid,iv,slev,elev
INTEGER                :: elev2,varid
REAL(r_sngl),INTENT(OUT) :: fieldg(nlon,nlat,elev-slev+1)
REAL(r_sngl)             :: auxfieldg(nlon,nlat,elev-slev+1)
INTEGER, ALLOCATABLE         :: start(:),count(:)
INTEGER                      :: nxvar,nyvar,nzvar,ierr
CHARACTER(2), INTENT(IN)     :: flag
CHARACTER(20)                :: varname,auxvarname
LOGICAL                      :: PAREXIST , READVAR
varname='                    '
auxvarname='                    '
fieldg=0.0e0
auxfieldg=0.0e0


IF ( flag .EQ. '3d' )THEN

   IF(iv <= nv3d)varname=element(iv)
   SELECT CASE (iv)
      CASE(iv3d_u)
       nxvar=nlon
       nyvar=nlat-1
       nzvar=nlev-1
       readvar=present_variable(iv)
      CASE(iv3d_v)
       nxvar=nlon-1
       nyvar=nlat
       nzvar=nlev-1  
       readvar=present_variable(iv)
      CASE(iv3d_w)
       nxvar=nlon-1
       nyvar=nlat-1
       nzvar=nlev
       readvar=present_variable(iv)
      CASE(iv3d_p)
       nxvar=nlon-1
       nyvar=nlat-1
       nzvar=nlev-1
       auxvarname='PB'
       readvar=present_variable(iv)
      CASE(iv3d_ph)
       nxvar=nlon-1
       nyvar=nlat-1
       nzvar=nlev-1
       auxvarname='PHB'
       readvar=present_variable(iv)
      CASE(iv3d_t,iv3d_qv,iv3d_qc,iv3d_qr,iv3d_qci,iv3d_qg,iv3d_qs,iv3d_coant,iv3d_cobck,iv3d_cobbu)
       nxvar=nlon-1
       nyvar=nlat-1
       nzvar=nlev-1
       readvar=present_variable(iv)
      CASE(iv3d_tv)
       varname='T'
       nxvar=nlon-1
       nyvar=nlat-1
       nzvar=nlev-1
       readvar=present_variable(iv)
      CASE DEFAULT
       WRITE(6,*)"3D Variable non-recognized, retrieve none: ",iv
       RETURN
  END SELECT

 IF( readvar )THEN

   IF( slev > nzvar )RETURN
   IF( elev > nzvar )THEN
     elev2=nzvar
   ELSE
     elev2=elev
   END IF
   ALLOCATE(start(4),count(4))
   start = (/ 1,1,slev,1 /)
   count = (/ nxvar,nyvar,elev2-slev+1,1 /) 

   !READ ADDITIONAL DATA
   IF( iv == iv3d_t .OR. iv == iv3d_tv )THEN
    CALL check_io(NF_INQ_VARID(ncid,'P',varid))
    CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1:elev2-slev+1)))
    CALL check_io(NF_INQ_VARID(ncid,'PB',varid))
    CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,auxfieldg(1:nxvar,1:nyvar,1:elev2-slev+1)))
    auxfieldg=auxfieldg+fieldg
    auxfieldg=( (auxfieldg/p0 ) ** (rd/cp) )
   ENDIF 
   IF( iv == iv3d_tv )THEN 
    CALL check_io(NF_INQ_VARID(ncid,'QVAPOR',varid))
    CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1:elev2-slev+1)))
    fieldg=fieldg/(1.0d0 + fieldg)
    auxfieldg=auxfieldg*(1.0d0 + 0.608d0*fieldg)
   ENDIF

   IF( iv == iv3d_p .OR. iv == iv3d_ph )THEN
      CALL check_io(NF_INQ_VARID(ncid,TRIM(auxvarname),varid))
      CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,auxfieldg(1:nxvar,1:nyvar,1:elev2-slev+1)))
      fieldg=fieldg+auxfieldg
   ENDIF
  
   !READ VARIABLE
   CALL check_io(NF_INQ_VARID(ncid,TRIM(varname),varid))
   CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1:elev2-slev+1)))

   !SOME EXTRA OPERATIONS
   IF( iv == iv3d_p .OR. iv == iv3d_ph )THEN
     fieldg = fieldg + auxfieldg
   ENDIF
   IF( iv == iv3d_t .OR. iv == iv3d_tv )THEN
     fieldg = ( fieldg + t0 ) * auxfieldg
   ENDIF


   DEALLOCATE(start,count)
 ELSE
   fieldg=undefs
 END IF  ! [PRESENT_VARIABLE]
END IF ! [3DVARIABLE]

IF ( flag .EQ. '2d' )THEN
   IF(iv <= nv2d)varname=element(nv3d+iv)
   ALLOCATE(start(3),count(3))
      SELECT CASE (iv)
      CASE(iv2d_ps,iv2d_t2,iv2d_q2)
       nxvar=nlon-1
       nyvar=nlat-1 
       readvar=present_variable(nv3d+iv)
      CASE(iv2d_mu)
       nxvar=nlon-1
       nyvar=nlat-1
       varname='MU'
       auxvarname='MUB'
       readvar=.true.
      CASE(iv2d_u10)
       nxvar=nlon-1
       nyvar=nlat-1
       readvar=present_variable(nv3d+iv)
      CASE(iv2d_v10)
       nxvar=nlon-1
       nyvar=nlat-1
       readvar=present_variable(nv3d+iv)
      CASE DEFAULT
       WRITE(6,*)"2D variable not recognized, retrieve none : ",iv
       RETURN
      END SELECT

  IF( readvar )THEN
  
      start = (/ 1,1,1 /)
      count = (/ nxvar,nyvar,1,1 /) !READ ONE SINGLE LEVEL
      CALL check_io(NF_INQ_VARID(ncid,TRIM(varname),varid))
      CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1)))
      IF( iv == iv2d_mu )THEN
          CALL check_io(NF_INQ_VARID(ncid,TRIM(auxvarname),varid))
          CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,auxfieldg(1:nxvar,1:nyvar,1)))
          fieldg=fieldg+auxfieldg
      ENDIF

   DEALLOCATE(start,count)
 ELSE
  fieldg=undefs
 END IF ![PRESENT_VARIABLE]

END IF !End 2d variables

!Parameter are more complicate.
!They can be in the file or not.
!If they are not, then we will assume default values for them.
!If estimated parameters are not present they will be perturbed and created in
!the netcdf file before outupt.
IF ( flag .EQ. 'pa' )THEN
   IF(iv <= np2d)varname=element(nv3d+nv2d+iv)
   ALLOCATE(start(3),count(3))
      SELECT CASE (iv)
      CASE(ip2d_hfx,ip2d_qfx,ip2d_ust)
       nxvar=nlon-1
       nyvar=nlat-1
       readvar=present_variable(nv3d+nv2d+iv)
      CASE DEFAULT
       WRITE(6,*)"2D parameter not recognized, will retrieve none : ",iv
       RETURN
      END SELECT

  IF( readvar  )THEN
      start = (/ 1,1,1 /)
      count = (/ nxvar,nyvar,1,1 /) !READ ONE SINGLE LEVEL

      ierr=NF_INQ_VARID(ncid,varname,varid)
      CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1)))

      fieldg(nlon,:,1)=REAL(UNDEF,r_sngl)
      fieldg(:,nlat,1)=REAL(UNDEF,r_sngl)

     DEALLOCATE(start,count)
  ELSE
     fieldg=undefs
  ENDIF ![PRESENT_VARIABLE]

END IF
RETURN

END SUBROUTINE read_var_wrf


SUBROUTINE write_var_wrf(ncid,iv,slev,elev,fieldg,flag)
IMPLICIT NONE
INCLUDE 'netcdf.inc'
INTEGER, INTENT(IN) :: ncid,iv,slev,elev
INTEGER(4)          :: elev2,varid,varidps,ndims,ierr
INTEGER(4), ALLOCATABLE :: dimsid(:)
REAL(r_sngl),    INTENT(INOUT)  :: fieldg(nlon,nlat,elev-slev+1)
REAL(r_sngl)                 :: auxfieldg(nlon,nlat,elev-slev+1)
INTEGER, ALLOCATABLE         :: start(:),count(:)
INTEGER                      :: nxvar,nyvar,nzvar
CHARACTER(2), INTENT(IN)     :: flag
CHARACTER(20)                :: varname,auxvarname
LOGICAL                      :: PAREXIST , WRITEVAR
varname='                    '
auxvarname='                    '
auxfieldg=0.0e0

IF ( flag .EQ. '3d' )THEN
   IF(iv <= nv3d)varname=element(iv)
   SELECT CASE (iv)
      CASE(iv3d_u)
       nxvar=nlon
       nyvar=nlat-1
       nzvar=nlev-1
       writevar=present_variable(iv)
      CASE(iv3d_v)
       nxvar=nlon-1
       nyvar=nlat
       nzvar=nlev-1
       writevar=present_variable(iv)
      CASE(iv3d_w)
       nxvar=nlon-1
       nyvar=nlat-1
       nzvar=nlev
       writevar=present_variable(iv)
      CASE(iv3d_p)
       nxvar=nlon-1
       nyvar=nlat-1
       nzvar=nlev-1
       auxvarname='PB'
       writevar=present_variable(iv)
      CASE(iv3d_ph)
       nxvar=nlon-1
       nyvar=nlat-1
       nzvar=nlev-1
       auxvarname='PHB'
       writevar=present_variable(iv)
      CASE(iv3d_t,iv3d_qv,iv3d_qc,iv3d_qr,iv3d_qci,iv3d_qg,iv3d_qs,iv3d_coant,iv3d_cobck,iv3d_cobbu)
       nxvar=nlon-1
       nyvar=nlat-1
       nzvar=nlev-1
       writevar=present_variable(iv)
      CASE DEFAULT
       WRITE(6,*)"3D Variable not recognized. Won't be written: ",iv
       RETURN
      END SELECT

 IF( writevar )THEN

  ALLOCATE(start(4),count(4))

   IF( slev > nzvar )RETURN
   IF( elev > nzvar )THEN
     elev2=nzvar
   ELSE
     elev2=elev
   END IF
   start = (/ 1,1,slev,1 /)
   count = (/ nxvar,nyvar,elev2-slev+1,1 /)
   IF( iv == iv3d_p .OR. iv == iv3d_ph )THEN
      CALL check_io(NF_INQ_VARID(ncid,TRIM(auxvarname),varid))
      CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,auxfieldg(1:nxvar,1:nyvar,1:elev2-slev+1)))
      fieldg=fieldg-auxfieldg
   ENDIF
   CALL check_io(NF_INQ_VARID(ncid,TRIM(varname),varid))
   CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1:elev2-slev+1)))

   DEALLOCATE(start,count)
  END IF ![PRESENT_VARIABLE]

END IF

IF ( flag .EQ. '2d' )THEN
   IF(iv <= nv2d)varname=element(nv3d+iv)
   ALLOCATE(start(3),count(3))
      SELECT CASE (iv)
      CASE(iv2d_ps,iv2d_t2,iv2d_q2)
       nxvar=nlon-1
       nyvar=nlat-1
       writevar=present_variable(nv3d+iv)
      CASE(iv2d_mu)
       nxvar=nlon-1
       nyvar=nlat-1
       varname='MU'
       auxvarname='MUB'
       writevar=.true.
      CASE(iv2d_u10)
       nxvar=nlon-1
       nyvar=nlat-1
       writevar=present_variable(nv3d+iv)
      CASE(iv2d_v10)
       nxvar=nlon-1
       nyvar=nlat-1
       writevar=present_variable(nv3d+iv)
      CASE DEFAULT
       WRITE(6,*)"Variable not recognized, won't be written",iv
       RETURN
      END SELECT

 IF( writevar )THEN
      start = (/ 1,1,1 /)
      count = (/ nxvar,nyvar,1,1 /) !READ ONE SINGLE LEVEL
      IF( iv == iv2d_mu )THEN
          CALL check_io(NF_INQ_VARID(ncid,TRIM(auxvarname),varid))
          CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,auxfieldg(1:nxvar,1:nyvar,1)))
          fieldg=fieldg-auxfieldg
      ENDIF
      CALL check_io(NF_INQ_VARID(ncid,TRIM(varname),varid))
      CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1)))


   DEALLOCATE(start,count)
 END IF ![PRESENT_VARIABLE]
END IF

IF ( flag .EQ. 'pa' )THEN
   IF(iv <= np2d)varname=element(nv3d+np2d+iv)
   ALLOCATE(start(3),count(3))
      
      SELECT CASE (iv)
      CASE(ip2d_hfx,ip2d_qfx,ip2d_ust)
       writevar=present_variable(nv3d+nv2d+iv)
       nxvar=nlon-1
       nyvar=nlat-1
      CASE DEFAULT
       WRITE(6,*)"Parameter not recognized, won't be written",iv
       RETURN
      END SELECT

 IF( writevar )THEN
      start = (/ 1,1,1 /)
      count = (/ nxvar,nyvar,1,1 /) !WRITE ONE SINGLE LEVEL

      PAREXIST=.false.
      ierr=NF_INQ_VARID(ncid,varname,varid)
      IF(ierr == nf_noerr)PAREXIST=.true.

      IF(.NOT. PAREXIST)THEN !Create the variable using PS as a template
       WRITE(6,*)'WARNING!!!! : PARAMETERS NOT PRESENT IN INPUT FILE'
       WRITE(6,*)'WILL GENERATE THE VARIABLES IN THE INPUT FILE'
       CALL check_io(NF_INQ_VARID(ncid,'PSFC',varidps))
       CALL check_io(NF_INQ_VARNDIMS(ncid,varidps,ndims))
       ALLOCATE(dimsid(ndims))
       CALL check_io(NF_INQ_VARDIMID(ncid,varidps,dimsid))
       CALL check_io(NF_REDEF(ncid))
       CALL check_io(NF_DEF_VAR(ncid,varname,NF_FLOAT,ndims,dimsid,varid))
       !Copy PSFC attribtues to the parameters.
       CALL check_io(NF_COPY_ATT(ncid,varidps,"units",ncid,varid))
       CALL check_io(NF_COPY_ATT(ncid,varidps,"FieldType",ncid,varid))
       CALL check_io(NF_COPY_ATT(ncid,varidps,"MemoryOrder",ncid,varid))
       CALL check_io(NF_COPY_ATT(ncid,varidps,"description",ncid,varid))
       CALL check_io(NF_COPY_ATT(ncid,varidps,"stagger",ncid,varid))
       CALL check_io(NF_COPY_ATT(ncid,varidps,"coordinates",ncid,varid))
       CALL check_io(NF_ENDDEF(ncid))
      END IF
      CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1)))

   DEALLOCATE(start,count)

  END IF ![PRESENT VARIABLE]

END IF


RETURN

END SUBROUTINE write_var_wrf


!-- Read a grid file ---------------------------------------------------
!SUBROUTINE read_grd(filename,v3d,v2d)
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: filename
!  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
!  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
!  REAL(r_sngl) :: buf4(nlon,nlat)
!  INTEGER :: iunit,iolen
!  INTEGER :: k,n,irec
!
!  iunit=11
!  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='sequential')
!
!  DO n=1,nv3d
!    DO k=1,nlev
!      READ(iunit) buf4
!      v3d(:,:,k,n) = REAL(buf4,r_size)
!    END DO
!  END DO
!
!  DO n=1,nv2d
!    READ(iunit) buf4
!    v2d(:,:,n) = REAL(buf4,r_size)
!  END DO
!
!  CLOSE(iunit)
!
!  RETURN
!END SUBROUTINE read_grd

!SUBROUTINE read_grd4(filename,v3d,v2d)
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: filename
!  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
!  INTEGER :: iunit,iolen
!  INTEGER :: i,j,k,n,irec
!
!  iunit=11
!  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='sequential')
!
!  DO n=1,nv3d
!    DO k=1,nlev
!      READ(iunit) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
!    END DO
!  END DO
!
!  DO n=1,nv2d
!    READ(iunit) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
!  END DO
!
!  CLOSE(iunit)
!
!  RETURN
!END SUBROUTINE read_grd4

!-- Write a grid file -------------------------------------------------
!SUBROUTINE write_grd(filename,v3d,v2d)
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: filename
!  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
!  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
!  REAL(r_sngl) :: buf4(nlon,nlat)
!  INTEGER :: iunit,iolen
!  INTEGER :: k,n,irec

!  iunit=55
!  INQUIRE(IOLENGTH=iolen) iolen
!  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='sequential')

!  DO n=1,nv3d
!    DO k=1,nlev
!      buf4 = REAL(v3d(:,:,k,n),r_sngl)
!      WRITE(iunit) buf4
!    END DO
!  END DO
!
!  DO n=1,nv2d
!    buf4 = REAL(v2d(:,:,n),r_sngl)
!    WRITE(iunit) buf4
!  END DO
!
!  CLOSE(iunit)
!
!  RETURN
!END SUBROUTINE write_grd

!SUBROUTINE write_grd4(filename,v3d,v2d)
!  IMPLICIT NONE
!  CHARACTER(*),INTENT(IN) :: filename
!  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
!  INTEGER :: iunit,iolen
!  INTEGER :: i,j,k,n
!
!  iunit=55
!  INQUIRE(IOLENGTH=iolen) iolen
!  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='sequential')
!
!  DO n=1,nv3d
!    DO k=1,nlev
!      WRITE(iunit) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
!    END DO
!  END DO
!
!  DO n=1,nv2d
!    WRITE(iunit) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
!  END DO
!
!  CLOSE(iunit)
!
!  RETURN
!END SUBROUTINE write_grd4

!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
SUBROUTINE monit_grd(v3d,v2d)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: k,n

  DO k=1,nlev
    WRITE(6,'(I2,A)') k,'th level'
    DO n=1,nv3d
      WRITE(6,'(A,2ES10.2)') element(n),MAXVAL(v3d(:,:,k,n)),MINVAL(v3d(:,:,k,n))
    END DO
  END DO

  DO n=1,nv2d
    WRITE(6,'(A,2ES10.2)') element(nv3d+n),MAXVAL(v2d(:,:,n)),MINVAL(v2d(:,:,n))
  END DO

  RETURN
END SUBROUTINE monit_grd
!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij,member,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: i,k,m,n


  DO n=1,nv3d
!!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO k=1,nlev
      DO i=1,nij
        v3dm(i,k,n) = v3d(i,k,1,n)
        DO m=2,member
          v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
        END DO
       v3dm(i,k,n) = v3dm(i,k,n) / REAL(member,r_size)
      END DO
    END DO
!!$OMP END PARALLEL DO
  END DO

  DO n=1,nv2d
    DO i=1,nij
      v2dm(i,n) = v2d(i,1,n)
      DO m=2,member
        v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
      END DO
      v2dm(i,n) = v2dm(i,n) / REAL(member,r_size)
    END DO
  END DO

  RETURN
END SUBROUTINE ensmean_grd

SUBROUTINE ensmean_ngrd(member,nij,v,vm,ndim)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member,ndim
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v(nij,member,ndim)
  REAL(r_size),INTENT(OUT) :: vm(nij,ndim)
  INTEGER :: i,k,m,n

  DO n=1,ndim
    DO i=1,nij
      vm(i,n) = v(i,1,n)
      DO m=2,member
        vm(i,n) = vm(i,n) + v(i,m,n)
      END DO
      vm(i,n) = vm(i,n) / REAL(member,r_size)
    END DO
  END DO

  RETURN
END SUBROUTINE ensmean_ngrd


SUBROUTINE check_io(status)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER(4),INTENT(IN) :: status

  IF(status /= nf_noerr) THEN
    WRITE(6,*) TRIM(nf_strerror(status))
    STOP 10
  ENDIF

  RETURN
END SUBROUTINE check_io

SUBROUTINE open_wrf_file(filename,mode,ncid)
  IMPLICIT NONE
  INTEGER(4), INTENT(OUT)   :: ncid
  CHARACTER(*) , INTENT(IN) :: filename,mode
  INCLUDE 'netcdf.inc'

  IF(mode .eq. 'ro' )THEN
  CALL check_io(NF_OPEN(TRIM(filename),NF_NOWRITE,ncid))
  ELSEIF(mode .eq. 'rw' )THEN
  CALL check_io(NF_OPEN(TRIM(filename),NF_WRITE,ncid))
  END IF

  RETURN
END SUBROUTINE open_wrf_file

SUBROUTINE close_wrf_file(ncid)
  IMPLICIT NONE
  INTEGER(4), INTENT(IN)   :: ncid
  INCLUDE 'netcdf.inc'

  CALL  check_io(NF_CLOSE(ncid))

  RETURN
END SUBROUTINE close_wrf_file


!-----------------------------------------------------------------------
! Adjust lateral boundary
!-----------------------------------------------------------------------
SUBROUTINE adjust_bnd(nlon,nlat,nlev,var)
  IMPLICIT NONE
  INTEGER(4),INTENT(in) :: nlon,nlat,nlev
  REAL(4),INTENT(inout) :: var(nlon,nlat,nlev)
  INTEGER(4) :: i,j,k
  
  DO k = 1, nlev
    DO i = 1, nlon-1
      var(i,1,k) = var(i,2,k)
    ENDDO
    DO i = 1, nlon-1
      var(i,nlat-1,k) = var(i,nlat-2,k)
    ENDDO  
    DO j = 1, nlat-1
      var(1,j,k) = var(2,j,k)
    ENDDO
    DO j = 1, nlat-1
      var(nlon-1,j,k) = var(nlon-2,j,k)
    ENDDO 
  ENDDO
  
  RETURN
END SUBROUTINE 
!-----------------------------------------------------------------------
! Damp lateral boundary
!-----------------------------------------------------------------------
SUBROUTINE damp_latbnd(spec_bdy_width,nlon,nlat,nlev,var)
  IMPLICIT NONE
  INTEGER(4),INTENT(in) :: spec_bdy_width
  INTEGER(4),INTENT(in) :: nlon,nlat,nlev
  REAL(4),INTENT(inout) :: var(nlon,nlat,nlev)
  INTEGER(4) :: i,j,k
  INTEGER(4) :: ist,ied,jst,jed

  ist = 1
  ied = nlon - 1
  jst = 1
  jed = nlat - 1
  DO k = 1, nlev
    DO i = ist, spec_bdy_width
      var(i,1:nlat,k) = var(i,1:nlat,k) * real(i-1) / real(spec_bdy_width) 
    END DO
    DO i = ied-spec_bdy_width+1, ied
      var(i,1:nlat,k) = var(i,1:nlat,k) * real(ied-i) / real(spec_bdy_width)       
    END DO
    DO j = jst, spec_bdy_width
      var(1:nlon,j,k) = var(1:nlon,j,k) * real(j-1) / real(spec_bdy_width)
    END DO
    DO j = jed-spec_bdy_width+1, jed
      var(1:nlon,j,k) = var(1:nlon,j,k) * real(jed-j) / real(spec_bdy_width)
    END DO      
  END DO
  
  RETURN
END SUBROUTINE

SUBROUTINE read_grd(filename,variables3d,variables2d)
IMPLICIT NONE
CHARACTER(*) :: filename
REAL(r_size) :: variables3d(nlon,nlat,nlev,nv3d) , variables2d(nlon,nlat,nv2d)
REAL(r_sngl) :: tmp(nlon,nlat,nlev)
INTEGER(4)   :: ncid 
INTEGER      :: slev , elev , iv 

     slev=1
     elev=nlev
     !Open the file
     CALL open_wrf_file(filename,'ro',ncid)
     !Get data 
     DO iv=1,nv3d
       !WRITE(6,*)iv
       tmp=0.0e0
       CALL read_var_wrf(ncid,iv,slev,elev,tmp,'3d')
       variables3d(:,:,:,iv)=REAL(tmp,r_size)
     ENDDO
     DO iv=1,nv2d
       !WRITE(6,*)iv
       tmp=0.0e0
       CALL read_var_wrf(ncid,iv,slev,elev,tmp(:,:,1),'2d')
       variables2d(:,:,iv)=REAL(tmp(:,:,1),r_size)
     ENDDO

     CALL close_wrf_file(ncid)

END SUBROUTINE read_grd

SUBROUTINE write_grd(filename,variables3d,variables2d)
IMPLICIT NONE
CHARACTER(*) :: filename
REAL(r_size) :: variables3d(nlon,nlat,nlev,nv3d) , variables2d(nlon,nlat,nv2d)
REAL(r_sngl) :: tmp(nlon,nlat,nlev)
INTEGER(4)   :: ncid
INTEGER      :: slev , elev , iv

     slev=1
     elev=nlev
     !Open the file
     CALL open_wrf_file(filename,'rw',ncid)
     !Get data
     DO iv=1,nv3d
       tmp=REAL(variables3d(nlon,nlat,nlev,iv),r_sngl)
       CALL write_var_wrf(ncid,iv,slev,elev,tmp,'3d')
     ENDDO
     DO iv=1,nv2d
       tmp=REAL(variables2d(nlon,nlat,iv),r_sngl)
       CALL write_var_wrf(ncid,iv,slev,elev,tmp(:,:,1),'2d')
     ENDDO

     CALL close_wrf_file(ncid)

END SUBROUTINE write_grd




END MODULE common_wrf
