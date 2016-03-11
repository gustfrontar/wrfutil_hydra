module verify_tools
!=======================================================================
!
! [PURPOSE:] Common procedure for verification against grided data
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
  INTEGER,SAVE :: dx,dy ! grid spacing [m]
  INTEGER,SAVE :: nlon
  INTEGER,SAVE :: nlat
  INTEGER,SAVE :: nlev
  
  INTEGER,SAVE :: nlonreg,nlatreg,nlevreg !Regrid domain size

  !NAMELIST INPUT

  INTEGER :: nbv=1  !Number of ensemble members.
  real(r_size) :: regrid_res = 1.0d0       ! Regrid output resolution (degree)
  real(r_size) :: regrid_vert_res = 10000.0d0! Regrid output vertical resolution (Pa)
  logical      :: regrid_output = .true.


  integer, parameter :: max_narea = 20
  integer :: narea 
  real(r_size) :: vlon1(max_narea) 
  real(r_size) :: vlon2(max_narea) 
  real(r_size) :: vlat1(max_narea) 
  real(r_size) :: vlat2(max_narea)

  CHARACTER(LEN=20) :: NAMELIST_FILE='./verify.namelist'

  !GRID INFO

  REAL(r_size),ALLOCATABLE,SAVE :: lon(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lat(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lev(:) , levall(:) 
  REAL(r_size),ALLOCATABLE,SAVE :: varall(:)

  REAL(r_size),ALLOCATABLE,SAVE :: lonreg(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: latreg(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: levreg(:) , levallreg(:)
  REAL(r_size),ALLOCATABLE,SAVE :: varallreg(:)
 
 
  REAL(r_size) :: minlon , minlat , maxlon , maxlat
  REAL(r_size) :: minlev , maxlev
  REAL(r_size) :: minlonreg , minlatreg , maxlonreg , maxlatreg
  REAL(r_size) :: minlevreg , maxlevreg
  INTEGER,ALLOCATABLE  :: levregmap(:)  !maps the original levels in the regrid levels.

  REAL(r_sngl) :: undefbin
  
  character(20) , ALLOCATABLE :: var_name(:)
  INTEGER , ALLOCATABLE :: var_lev(:)
  INTEGER :: nvar_ctl

  CHARACTER(19) :: hdate !Date in input files.

  TYPE(proj_info)  :: projection

  integer :: nv !Total number of variables (vars times levs)
  integer :: nvreg !Total number of variables in regrid output ( vars times levels)

  logical :: file_exist
  integer :: ierr

CONTAINS

!-----------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------

SUBROUTINE read_namelist
  !Default values for verifcation regions.
                    ! GLOBE    NH     TROPICS   SH
narea=4
vlon1(1:narea) = (/  0.d0,   0.d0,   0.d0,   0.d0/)
vlon2(1:narea) = (/360.d0, 360.d0, 360.d0, 360.d0/)
vlat1(1:narea) = (/-90.d0,  20.d0, -20.d0, -90.d0/)
vlat2(1:narea) = (/ 90.d0,  90.d0,  20.d0, -20.d0/)

NAMELIST / GENERAL / nbv , narea , vlon1 , vlon2 , vlat1 , vlat2 , regrid_output &
                  ,  regrid_res , regrid_vert_res 

INQUIRE(FILE=NAMELIST_FILE, EXIST=file_exist)

IF( .NOT. file_exist )THEN
  WRITE(*,*)"ERROR: Could not find namelist"
ENDIF

OPEN(54,FILE=NAMELIST_FILE)
READ(54,NML=GENERAL,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at GENERAL section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)

END SUBROUTINE read_namelist

!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------
!Read dimensions from ctl file
SUBROUTINE parse_ctl(ctl_file)
IMPLICIT NONE
character(*) , INTENT(IN) :: ctl_file
character(100)            :: buffer , dummyc
integer                   :: iunit , i , j , counter
real(r_size),allocatable  :: tmplat(:) , tmplon(:) , tmplev(:)

iunit=98

 INQUIRE(FILE=ctl_file,EXIST=file_exist)
  IF(file_exist) THEN
    OPEN(iunit,FILE=ctl_file,FORM='formatted')
  ELSE
    WRITE(*,*)'[Error]: Missing ctl-file '
    STOP
  ENDIF

  nlon=0
  nlat=0
  nlev=0
  undefbin=0.0e0
  nv=0
  

 !This is valid for mercator 
 !TODO: Implementar otras proyecciones como LAMBERT
 DO
   READ(iunit,'(A100)',IOSTAT=IERR)buffer
      IF(IERR /= 0)THEN
       WRITE(*,*)'End of file'
       EXIT
      ENDIF
      IF( INDEX(buffer,'undef') > 0 .or. INDEX(buffer,'UNDEF') > 0 )THEN
       READ(buffer,*)dummyc , undefbin
       WRITE(*,*)'Undef value for this file is ',undefbin
      ENDIF
      
      IF( INDEX(buffer,'xdef') > 0 .or. INDEX(buffer,'XDEF') > 0 )THEN
        READ(buffer,*)dummyc ,  nlon
        allocate(tmplon(nlon))
        DO i=1,nlon 
          READ(iunit,*)tmplon(i)
        ENDDO
      ENDIF
      IF( INDEX(buffer,'ydef') > 0 .or. INDEX(buffer,'YDEF') > 0 )THEN
        READ(buffer,*)dummyc ,  nlat
        allocate(tmplat(nlat))
        DO i=1,nlat
          READ(iunit,*)tmplat(i)
        ENDDO
      ENDIF
      IF( INDEX(buffer,'zdef') > 0 .or. INDEX(buffer,'ZDEF') > 0 )THEN
        READ(buffer,*)dummyc ,  nlev
        allocate(tmplev(nlev))
        DO i=1,nlev
          READ(iunit,*)tmplev(i)
        ENDDO
      ENDIF
      IF( INDEX(buffer,'vars') > 0 .or. INDEX(buffer,'VARS') > 0 .AND. .NOT. &
        ( INDEX(buffer,'endvars') > 0 .or. INDEX(buffer,'ENDVARS') > 0 )  )THEN
        READ(buffer,*)dummyc,nvar_ctl
        allocate(var_name(nvar_ctl) )
        allocate(var_lev(nvar_ctl)  )
        DO i=1,nvar_ctl
          READ(iunit,*)var_name(i),var_lev(i)
        ENDDO
        !Get the total number of fields.
        nv=SUM(var_lev)
      ENDIF
 ENDDO

 minlat=undef
 maxlat=-undef
 minlon=undef
 maxlon=-undef
 
 allocate( lat(nlon,nlat) , lon(nlon,nlat) , lev(nlev) )
 DO i = 1,nlon
  DO j = 1,nlat
     lat(i,j)=tmplat(j)
     lon(i,j)=tmplon(i)

     if ( lat(i,j) < minlat )minlat=lat(i,j)
     if ( lat(i,j) > maxlat )maxlat=lat(i,j)
     if ( lon(i,j) < minlon )minlon=lon(i,j)
     if ( lon(i,j) > maxlon )maxlon=lon(i,j) 
  ENDDO
 ENDDO

 lev=tmplev

 minlev=MINVAL(lev)
 maxlev=MAXVAL(lev)

 !Set levall (the level corresponding to each variable)
 ALLOCATE( levall(nv) , varall(nv) )
 counter=1
 DO i=1,nvar_ctl
   DO j=1,var_lev(i)
     levall(counter)=lev(j)
     varall(counter)=REAL(i,r_size)
     counter=counter+1
   ENDDO
 ENDDO

 IF( nlon == 0 )THEN
  WRITE(*,*)'[Error] Failed to get longitudes.'
  STOP
 ELSEIF( nlat == 0 )THEN
  WRITE(*,*)'[Error] Failed to get latitudes.'
  STOP
 ELSEIF( nlev == 0 )THEN
  WRITE(*,*)'[Error] Failed to get levels.'
  STOP
 ELSEIF( undefbin == 0.0e0 )THEN
  WRITE(*,*)'[Error] Failed to get undef.'
  STOP
 ENDIF

 WRITE(*,*)'===================================================='
 WRITE(*,*)' The total number of individual variables is = ',nvar_ctl
 WRITE(*,*)' Total number of matrices in each file is = ',nv
 WRITE(*,*)'===================================================='
 WRITE(*,*)'MINLON = ',minlon
 WRITE(*,*)'MINLAT = ',minlat
 WRITE(*,*)'MAXLON = ',maxlon
 WRITE(*,*)'MAXLAT = ',maxlat

END SUBROUTINE parse_ctl
!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------

!-- Read a grid file ---------------------------------------------------
SUBROUTINE read_grd(filename,nx,ny,nz,var,undefmask)
  IMPLICIT NONE
  INTEGER,INTENT(IN)      :: nx,ny,nz
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: var(nx,ny,nz) 
  LOGICAL     ,INTENT(OUT) :: undefmask(nx,ny,nz)
  REAL(r_sngl) :: buf4(nx,ny)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,reclength,ios

  undefmask=.true.

  INQUIRE(IOLENGTH=reclength) reclength
  reclength = reclength * nx * ny


 iunit=33
 INQUIRE(FILE=filename,EXIST=file_exist)
  IF(file_exist) THEN
   OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=reclength)
  ELSE
   undefmask=.false.
   var=0.0d0
   RETURN
  ENDIF

  DO n=1,nz
    READ(iunit,rec=n,IOSTAT=ios)buf4
     IF( ios == 0 )THEN
      WHERE( buf4 == undefbin )
       undefmask(:,:,n)=.false.
      ENDWHERE
       var(:,:,n) = REAL(buf4,r_size)
     ELSE !We reached the end of the file.
       IF( n==1)THEN
         WRITE(*,*)"[Warning]: Error reading, this might be a new file "
         WRITE(*,*)filename
         var=0.0d0     
         undefmask=.false.
       ELSE
         WRITE(*,*)"[Error]: Unexpectively reached end of file"
         WRITE(*,*)filename
         STOP
       ENDIF
     ENDIF
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE read_grd

!-- Write a grid file -------------------------------------------------
SUBROUTINE write_grd(filename,nx,ny,nz,var,undefmask)
  IMPLICIT NONE
  INTEGER, INTENT(IN)     :: nx,ny,nz
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: var(nx,ny,nz)
  LOGICAL     ,INTENT(IN) :: undefmask(nx,ny,nz)
  REAL(r_sngl) :: buf4(nx,ny)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,reclength

  INQUIRE(IOLENGTH=reclength) reclength
  reclength = reclength * nx * ny


  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=reclength)

  DO n=1,nz
   buf4 = REAL(var(:,:,n),r_sngl)
   WHERE( .NOT. undefmask(:,:,n) )
    buf4=undefbin
   ENDWHERE     
   
   WRITE(iunit,rec=n) buf4
  ENDDO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_grd

!Compute regrid grid. ------------------------------------------------------
SUBROUTINE get_regrid_grid()
IMPLICIT NONE
INTEGER      :: i,j,ivar,ilev,ilevreg,ilevprev,iregprev

  regrid_vert_res=regrid_vert_res / 100.0d0 !From Pa to hPa.

  nlonreg=CEILING( (maxlon-minlon)/regrid_res ) + 1
  nlatreg=CEILING( (maxlat-minlat)/regrid_res ) + 1

  ALLOCATE( lonreg(nlonreg,nlatreg) , latreg(nlonreg,nlatreg) )

  DO i = 1,nlonreg
   DO j = 1,nlatreg

    lonreg(i,j)=minlon + REAL((i-1),r_size)*regrid_res
    latreg(i,j)=minlat + REAL((j-1),r_size)*regrid_res    

   ENDDO
  ENDDO
 
  minlonreg=minlon
  minlatreg=minlat
  maxlonreg=minlon+regrid_res*REAL(nlonreg-1,r_size)
  maxlatreg=minlat+regrid_res*REAL(nlatreg-1,r_size)
  
  maxlevreg=maxlev

  nlevreg= CEILING( (maxlev-minlev)/ regrid_vert_res ) + 1

  minlevreg= maxlevreg - REAL(nlevreg,r_size)*regrid_vert_res 

   WRITE(*,*)'=================================================='
   WRITE(*,*)'    Regrid output enambled'
   WRITE(*,*)'=================================================='
   WRITE(*,*)'Regrid output grid                                '
   WRITE(*,*)'INILON= ',minlonreg
   WRITE(*,*)'INILAT= ',minlatreg
   WRITE(*,*)'NLON  = ',nlonreg
   WRITE(*,*)'NLAT  = ',nlatreg
   WRITE(*,*)'HRES  = ',regrid_res*100.0d0,' (degree)'
   WRITE(*,*)'MAXLEV= ',maxlevreg
   WRITE(*,*)'NLEV  = ',nlevreg
   WRITE(*,*)'LEVRES= ',regrid_vert_res, ' (Pa)'
   WRITE(*,*)'=================================================='
 
  ALLOCATE( levreg(nlevreg) ) 
  DO i = 1,nlevreg
    levreg(i)=maxlevreg - REAL((i-1),r_size)*regrid_vert_res 
  ENDDO
 
  !How many total levels do we have in the regird domain?
  nvreg=1
  ilevprev=NINT( ( maxlevreg - levall(1) )  / regrid_vert_res ) + 1
  DO i = 2,nv
   ilev=NINT( ( maxlevreg - levall(i) )  / regrid_vert_res ) + 1 
   if ( ilev /= ilevprev .or. varall(i) /= varall(i-1) )then
      nvreg = nvreg+1
   endif 
   ilevprev=ilev
  ENDDO

  ilevreg=1

  ALLOCATE( levallreg(nvreg) )
  ALLOCATE( varallreg(nvreg) )
  ALLOCATE( levregmap(nv)    )
  ilevprev=NINT( ( maxlevreg - levall(1) )  / regrid_vert_res ) + 1
  levallreg(1)=levall(ilevprev) 
  levregmap(1)=ilevprev
  varallreg(1)=varall(1)
  DO i = 2,nv
   ilev=NINT( ( maxlevreg - levall(i) )  / regrid_vert_res ) + 1
   if ( ilev /= ilevprev .or. varall(i) /= varall(i-1)  )then
     ilevreg=ilevreg+1
     levallreg(ilevreg)=levreg(ilev)
     varallreg(ilevreg)=varall(i)
   endif
   levregmap(i)=ilevreg
   if( varall(i-1) /= varall(i) )then !Variable change
    WRITE(*,*)"NLEVS for var ",var_name(varall(i-1))," ",ilevprev
   endif
   ilevprev=ilev
  ENDDO


  WRITE(*,*)'=================================================='


END SUBROUTINE get_regrid_grid

!Regrid var 
!This subroutine regrids -----------------------------------------------------------
SUBROUTINE regrid_var(var,undefmask,varreg,undefmaskreg)
IMPLICIT NONE
REAL(r_size), INTENT(IN) :: var(nlon,nlat,nv) 
REAL(r_size), INTENT(OUT):: varreg(nlonreg,nlatreg,nvreg)
INTEGER                  :: num(nlonreg,nlatreg,nvreg)
LOGICAL                  :: undefmask(nlon,nlev,nv),undefmaskreg(nlonreg,nlatreg,nvreg)
INTEGER                  :: i , j , k , ireg , jreg , kreg

   varreg=0.0d0
   num=0.0d0
   undefmaskreg=.true.

   DO i = 1,nlon
    DO j = 1,nlat
     DO k = 1,nv
      if( undefmask(i,j,k) )then
        ireg= NINT( ( lon(i,j)-minlonreg ) / regrid_res  ) + 1
        jreg= NINT( ( lat(i,j)-minlatreg ) / regrid_res  ) + 1
        kreg= levregmap(k)  !Since k dimension includes variables an levels we
                            !have to define a map to know whicho original level
                            !correspond to the regrid levels.
        varreg(ireg,jreg,kreg)=varreg(ireg,jreg,kreg)+var(i,j,k) 
        num(ireg,jreg,kreg)   =num(ireg,jreg,kreg)   +1
       endif
     ENDDO
    ENDDO
   ENDDO

   WHERE( num > 0 )
     varreg=varreg/REAL(num,r_size)
   ELSEWHERE
     undefmaskreg=.false.
   ENDWHERE


END SUBROUTINE regrid_var

END MODULE verify_tools
