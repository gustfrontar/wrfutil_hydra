MODULE common_wrf
!=======================================================================
!
! [PURPOSE:] Common Information for WRF
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   07/14/2010 Takemasa Miyoshi  modified for WRF
!   TODO: This routine is for NCEP data, other input sources may produce
!   different variables inside met_em files. In that case the routine 
!   has to be able to process the data including the new variables and
!   automatically detecting the abssence of some other variables.
!=======================================================================
!$USE OMP_LIB
  USE common
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  REAL(r_sngl) :: tmpatt
  INTEGER,SAVE :: dx,dy ! grid spacing [m]
  INTEGER,SAVE :: nlon
  INTEGER,SAVE :: nlat
  INTEGER,SAVE :: nlev
  INTEGER,PARAMETER :: nv3d=7! u,v,w,t,ght,rh
  INTEGER,PARAMETER :: nv2d=11 ! ps , slp , skintemp , soilvariables
  INTEGER,PARAMETER :: np2d=3 ! number of 2d estimated parameters
  LOGICAL,PARAMETER :: ESTPAR=.FALSE.  !If parameters will be estimated or not.
  REAL(r_size) , PARAMETER :: param_default_value(np2d)= (/1.0,1.0,1.0/)!If parameter is not estimated then this value will be assigned to parameter analysis.
  REAL(r_size) , PARAMETER :: param_max_value(np2d)= (/2.0, 2.0, 2.0/)  !This is the maximum allowed value for the parameters.
  REAL(r_size) , PARAMETER :: param_min_value(np2d)= (/0.5, 0.5, 0.5/)  !This is the minimum allowed value for the parameters.
  INTEGER,PARAMETER :: iv3d_u   =1
  INTEGER,PARAMETER :: iv3d_v   =2
  INTEGER,PARAMETER :: iv3d_w   =3
  INTEGER,PARAMETER :: iv3d_t   =4
  INTEGER,PARAMETER :: iv3d_ght =5
  INTEGER,PARAMETER :: iv3d_rh  =6
  INTEGER,PARAMETER :: iv3d_pres=7
  INTEGER,PARAMETER :: iv2d_psfc    =1
  INTEGER,PARAMETER :: iv2d_pmsl    =2
  INTEGER,PARAMETER :: iv2d_skintemp=3
  INTEGER,PARAMETER :: iv2d_st100200=4
  INTEGER,PARAMETER :: iv2d_st040100=5
  INTEGER,PARAMETER :: iv2d_st010040=6
  INTEGER,PARAMETER :: iv2d_st000010=7
  INTEGER,PARAMETER :: iv2d_sm100200=8
  INTEGER,PARAMETER :: iv2d_sm040100=9
  INTEGER,PARAMETER :: iv2d_sm010040=10
  INTEGER,PARAMETER :: iv2d_sm000010=11


  INTEGER,PARAMETER :: ip2d_hfxfactor=1
  INTEGER,PARAMETER :: ip2d_qfxfactor=2
  INTEGER,PARAMETER :: ip2d_ustfactor=3

  REAL(4),SAVE      :: p0 != 1.0e+5
  REAL(4),PARAMETER :: t0 = 300.0
  INTEGER(4),PARAMETER :: imt_grd = 50
  INTEGER(4),PARAMETER :: spec_bdy_width = 10
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
  REAL(r_size),ALLOCATABLE,SAVE :: landmask(:,:)
  CHARACTER(20),SAVE :: element(nv3d+nv2d)
  REAL(r_sngl),ALLOCATABLE,SAVE :: dnw(:)
  REAL(r_sngl),ALLOCATABLE,SAVE :: znu(:)
!  CHARACTER(19),SAVE  :: DATE

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_wrf(inputfile)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  REAL(r_sngl),ALLOCATABLE :: buf4(:,:)
  INTEGER :: i
  INTEGER(4) :: ncid,varid,dimids(3),shape(3)
  INTEGER(4),ALLOCATABLE :: start(:),count(:)
  CHARACTER(*)           :: inputfile

  WRITE(6,'(A)') 'Hello from set_common_wrf'

  !
  ! Elements
  !
  !3D variables
  element(iv3d_u)    = 'UU'
  element(iv3d_v)    = 'VV'
  element(iv3d_w)    = 'WW'
  element(iv3d_t)    = 'TT'
  element(iv3d_ght)  = 'GHT'
  element(iv3d_rh)   = 'RH'
  element(iv3d_pres) = 'PRES'
  !2D variables
  element(nv3d+iv2d_psfc)     = 'PSFC'
  element(nv3d+iv2d_pmsl)     = 'PMSL'
  element(nv3d+iv2d_skintemp) = 'SKINTEMP'
  element(nv3d+iv2d_st100200) = 'ST100200'
  element(nv3d+iv2d_st040100) = 'ST040100'
  element(nv3d+iv2d_st010040) = 'ST010040'
  element(nv3d+iv2d_st000010) = 'ST100200'
  element(nv3d+iv2d_sm100200) = 'SM100200'
  element(nv3d+iv2d_sm040100) = 'SM040100'
  element(nv3d+iv2d_sm010040) = 'SM010040'
  element(nv3d+iv2d_sm000010) = 'SM100200'

  !2D parameters
  element(nv3d+nv2d+ip2d_hfxfactor)='HFXFACTOR'
  element(nv3d+nv2d+ip2d_hfxfactor)='QFXFACTOR'
  element(nv3d+nv2d+ip2d_hfxfactor)='USTFACTOR'
  !
  ! Lon, Lat, F, phi0
  !
  WRITE(6,'(A)') 'READING: ',inputfile
  CALL check_io(NF_OPEN(inputfile,NF_NOWRITE,ncid))
  CALL check_io(NF_INQ_VARID(ncid,'PRES',varid))
  CALL check_io(NF_INQ_VARDIMID(ncid,varid,dimids))
  DO i=1,3
    CALL check_io(NF_INQ_DIMLEN(ncid,dimids(i),shape(i)))
  END DO
  nlon = shape(1)+1
  nlat = shape(2)+1 
  nlev = shape(3) 
  WRITE(6,'(A)') '*** grid information ***'
  WRITE(6,'(3(2X,A,I5))') 'nlon =',nlon,'nlat =',nlat,'nlev =',nlev
  nlevall=nlev*nv3d+nv2d

  CALL check_io(NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DX',tmpatt))
   dx=INT(tmpatt)
  CALL check_io(NF_GET_ATT_REAL(ncid,NF_GLOBAL,'DY',tmpatt))
   dy=INT(tmpatt)
  

  WRITE(6,*) '*** MODEL RESOLUTION IS = ',dx
  WRITE(6,'(a)') '***  END  : READ NETCDF (CNST) ***'
  CALL check_io(NF_CLOSE(ncid))


  RETURN
END SUBROUTINE set_common_wrf
!-----------------------------------------------------------------------
! Get time nc
!-----------------------------------------------------------------------
SUBROUTINE get_date(inputfile,date)
IMPLICIT NONE
CHARACTER(19),INTENT(OUT) :: date
CHARACTER(*) ,INTENT(IN)  :: inputfile
INTEGER :: ncid , start(2) , count(2) , varid
INCLUDE 'netcdf.inc'

  CALL check_io(NF_OPEN(inputfile,NF_NOWRITE,ncid))

  CALL check_io(NF_INQ_VARID(ncid,'Times',varid))
  start = (/  1,1 /)
  count = (/ 19,1 /)
  CALL check_io(NF_GET_VARA_TEXT(ncid,varid,start,count,date))

  CALL check_io(NF_CLOSE(ncid))

END SUBROUTINE get_date
!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------
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
LOGICAL                      :: PAREXIST
varname='                    '
auxvarname='                    '
fieldg=0.0e0
auxfieldg=0.0e0


IF ( flag .EQ. '3d' )THEN

   varname=element(iv)
   SELECT CASE (iv)
      CASE(iv3d_u)
       nxvar=nlon
       nyvar=nlat-1
       nzvar=nlev
      CASE(iv3d_v)
       nxvar=nlon-1
       nyvar=nlat
       nzvar=nlev
      CASE(iv3d_w,iv3d_t,iv3d_rh,iv3d_ght,iv3d_pres)
       nxvar=nlon-1
       nyvar=nlat-1
       nzvar=nlev

      CASE DEFAULT
       WRITE(6,*)"3D Variable non-recognized, retrieve none: ",iv
       RETURN
  END SELECT

   IF( slev > nzvar )RETURN
   IF( elev > nzvar )THEN
     elev2=nzvar
   ELSE
     elev2=elev
   END IF
   ALLOCATE(start(4),count(4))
   start = (/ 1,1,slev,1 /)
   count = (/ nxvar,nyvar,elev2-slev+1,1 /) 

   !READ VARIABLE
   IF( iv == iv3d_w )THEN
     ierr=NF_INQ_VARID(ncid,varname,varid)
     IF(ierr == nf_noerr)THEN
       CALL check_io(NF_INQ_VARID(ncid,TRIM(varname),varid))
       CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1:elev2-slev+1)))
     ELSE
       WRITE(*,*)"WARNING!: WW not present in input file. Initializing WW as 0.0d0"
       fieldg=0.0d0
     ENDIF
   ELSE 
     CALL check_io(NF_INQ_VARID(ncid,TRIM(varname),varid))
     CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1:elev2-slev+1)))
   ENDIF

   DEALLOCATE(start,count)
END IF

IF ( flag .EQ. '2d' )THEN
   ALLOCATE(start(3),count(3))
      varname=element(nv3d+iv)
      nxvar=nlon-1
      nyvar=nlat-1

      start = (/ 1,1,1 /)
      count = (/ nxvar,nyvar,1,1 /) !READ ONE SINGLE LEVEL
      CALL check_io(NF_INQ_VARID(ncid,TRIM(varname),varid))
      CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1)))

   DEALLOCATE(start,count)

END IF !End 2d variables



RETURN

END SUBROUTINE read_var_wrf


SUBROUTINE write_var_wrf(ncid,iv,slev,elev,fieldg,flag)
IMPLICIT NONE
INCLUDE 'netcdf.inc'
INTEGER, INTENT(IN) :: ncid,iv,slev,elev
INTEGER(4)          :: elev2,varid,varidt,varidu,ndims,ierr
INTEGER(4), ALLOCATABLE :: dimsid(:)
REAL(r_sngl),    INTENT(INOUT)  :: fieldg(nlon,nlat,elev-slev+1)
REAL(r_sngl)                 :: auxfieldg(nlon,nlat,elev-slev+1)
INTEGER, ALLOCATABLE         :: start(:),count(:)
INTEGER                      :: nxvar,nyvar,nzvar
CHARACTER(2), INTENT(IN)     :: flag
CHARACTER(20)                :: varname,auxvarname
LOGICAL                      :: PAREXIST
varname='                    '
auxvarname='                    '
auxfieldg=0.0e0

IF ( flag .EQ. '3d' )THEN

   varname=element(iv)
   SELECT CASE (iv)
      CASE(iv3d_u)
       nxvar=nlon
       nyvar=nlat-1
       nzvar=nlev
      CASE(iv3d_v)
       nxvar=nlon-1
       nyvar=nlat
       nzvar=nlev
      CASE(iv3d_w,iv3d_t,iv3d_ght,iv3d_rh,iv3d_pres)
       nxvar=nlon-1
       nyvar=nlat-1
       nzvar=nlev

      CASE DEFAULT
       WRITE(6,*)"3D Variable not recognized. Won't be written: ",iv
       RETURN
      END SELECT

  ALLOCATE(start(4),count(4))

   IF( slev > nzvar )RETURN
   IF( elev > nzvar )THEN
     elev2=nzvar
   ELSE
     elev2=elev
   END IF

   start = (/ 1,1,slev,1 /)
   count = (/ nxvar,nyvar,elev2-slev+1,1 /)

   IF( iv == iv3d_w )THEN
     ierr=NF_INQ_VARID(ncid,varname,varid)
     IF(ierr /= nf_noerr)THEN
       WRITE(6,*)'WARNING!!!! : WW not present in output file. We will create the variable in the output file.'
       CALL check_io(NF_INQ_VARID(ncid,'TT',varidt))
       CALL check_io(NF_INQ_VARID(ncid,'UU',varidu))
       !Use u-wind as a template for WW and create the variable.
       CALL check_io(NF_INQ_VARNDIMS(ncid,varidt,ndims))
       ALLOCATE(dimsid(ndims))
       CALL check_io(NF_INQ_VARDIMID(ncid,varidt,dimsid))
       CALL check_io(NF_REDEF(ncid))
       CALL check_io(NF_DEF_VAR(ncid,varname,NF_FLOAT,ndims,dimsid,varid))
       CALL check_io(NF_COPY_ATT(ncid,varidu,"units",ncid,varid))
       CALL check_io(NF_COPY_ATT(ncid,varidt,"FieldType",ncid,varid))
       CALL check_io(NF_COPY_ATT(ncid,varidt,"MemoryOrder",ncid,varid))
       CALL check_io(NF_COPY_ATT(ncid,varidt,"description",ncid,varid))
       CALL check_io(NF_COPY_ATT(ncid,varidt,"stagger",ncid,varid))
       CALL check_io(NF_COPY_ATT(ncid,varidt,"sr_x",ncid,varid))
       CALL check_io(NF_COPY_ATT(ncid,varidt,"sr_y",ncid,varid))
       CALL check_io(NF_ENDDEF(ncid))
     ENDIF
   ENDIF
       CALL check_io(NF_INQ_VARID(ncid,TRIM(varname),varid))
       CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1:elev2-slev+1)))


   DEALLOCATE(start,count)
END IF

IF ( flag .EQ. '2d' )THEN

   varname=element(nv3d+iv)
   ALLOCATE(start(3),count(3))
       nxvar=nlon-1
       nyvar=nlat-1

      start = (/ 1,1,1 /)
      count = (/ nxvar,nyvar,1,1 /) !READ ONE SINGLE LEVEL
      CALL check_io(NF_INQ_VARID(ncid,TRIM(varname),varid))
      CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,fieldg(1:nxvar,1:nyvar,1)))


   DEALLOCATE(start,count)
END IF


RETURN

END SUBROUTINE write_var_wrf


!-- Read a grid file ---------------------------------------------------
SUBROUTINE read_bin(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,irec
!
  iunit=11
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='sequential')

  DO n=1,nv3d
    DO k=1,nlev
      WRITE(*,*)n,k
      READ(iunit) buf4
      v3d(:,:,k,n) = REAL(buf4,r_size)
    END DO
  END DO

  DO n=1,nv2d
    READ(iunit) buf4
    v2d(:,:,n) = REAL(buf4,r_size)
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE read_bin


!-- Write a grid file -------------------------------------------------
SUBROUTINE write_bin(filename,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: buf4(nlon,nlat)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,irec

  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='sequential')

  DO n=1,nv3d
    DO k=1,nlev
      buf4 = REAL(v3d(:,:,k,n),r_sngl)
      WRITE(iunit) buf4
    END DO
  END DO

  DO n=1,nv2d
    buf4 = REAL(v2d(:,:,n),r_sngl)
    WRITE(iunit) buf4
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_bin


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
!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO k=1,nlev
      DO i=1,nij
        v3dm(i,k,n) = v3d(i,k,1,n)
        DO m=2,member
          v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
        END DO
       v3dm(i,k,n) = v3dm(i,k,n) / REAL(member,r_size)
      END DO
    END DO
!$OMP END PARALLEL DO
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

!-----------------------------------------------------------------------
! Damp lateral boundary
!-----------------------------------------------------------------------
SUBROUTINE damp_latbnd_double(spec_bdy_width,nlon,nlat,nlev,var)
  IMPLICIT NONE
  INTEGER,INTENT(in) :: spec_bdy_width
  INTEGER,INTENT(in) :: nlon,nlat,nlev
  REAL(r_size),INTENT(inout) :: var(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ist,ied,jst,jed

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
INCLUDE 'netcdf.inc'


     slev=1
     elev=nlev
     !Open the file
     CALL open_wrf_file(filename,'ro',ncid)

     DO iv=1,nv3d
       tmp=0.0e0
       WRITE(*,*)'Reading variable ',element(iv),' iv = ',iv
       CALL read_var_wrf(ncid,iv,slev,elev,tmp,'3d')
       variables3d(:,:,:,iv)=REAL(tmp,r_size)
     ENDDO
     DO iv=1,nv2d
       tmp=0.0e0
       WRITE(*,*)'Reading variable ',element(nv3d+iv),' iv = ',iv
       CALL read_var_wrf(ncid,iv,slev,elev,tmp(:,:,1),'2d')
       variables2d(:,:,iv)=REAL(tmp(:,:,1),r_size)
     ENDDO
     WRITE(*,*)'Finish reading file ',filename

END SUBROUTINE read_grd


SUBROUTINE write_grd(filename,v3d,v2d)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(4) :: sdmd, s1md
  REAL(r_sngl) :: fieldg(nlon,nlat,nlev)
  INTEGER :: l,n,ll,iv,ilev,ncid,ierr,i,j,k,varid
  CHARACTER(*),INTENT(IN) :: filename
  INTEGER :: rstart(2) , rend(2)
  REAL(r_size) ::  mu(nlon,nlat),ps(nlon,nlat),qv(nlon,nlat,nlev),t(nlon,nlat,nlev)
  INCLUDE 'netcdf.inc'


  !OPEN NC FILE
  CALL open_wrf_file(filename,'rw',ncid)

  rstart = (/ 1,1 /)
  rend   = (/ 19,1/)

  !Modify input date (if requested)
  !CALL check_io(NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'START_DATE',19,DATE))
  !CALL check_io(NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'SIMULATION_START_DATE',19,DATE))
  !CALL check_io(NF_INQ_VARID(ncid,'Times',varid))
  !CALL check_io(NF_PUT_VARA_TEXT(ncid,varid,rstart,rend,DATE))

  !WRITE 3D VARIABLES
  !GATHER AND WRITE ONE GRID AT A TIME
  DO iv=1,nv3d
   fieldg=REAL(v3d(:,:,:,iv),r_sngl)

   !WRITE THE DATA    
   CALL write_var_wrf(ncid,iv,1,nlev,fieldg,'3d')

 END DO  !En do over variables.

 !WRITE 2D VARIABLES
 !1 GRID AT A TIME
  DO iv=1,nv2d
    fieldg(:,:,1) = v2d(:,:,iv)
    
    CALL write_var_wrf(ncid,iv,1,1,fieldg(:,:,1),'2d')
  ENDDO

  !CLOSE FILES
  CALL close_wrf_file(ncid)

RETURN

END SUBROUTINE write_grd

END MODULE common_wrf
