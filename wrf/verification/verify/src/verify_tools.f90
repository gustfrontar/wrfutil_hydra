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
  !INTEGER,SAVE :: nlon
  !INTEGER,SAVE :: nlat
  !INTEGER,SAVE :: nlev
  
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

!  REAL(r_size),ALLOCATABLE,SAVE :: lon(:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: lat(:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: lev(:) , levall(:) 
!  REAL(r_size),ALLOCATABLE,SAVE :: varall(:)

  REAL(r_size),ALLOCATABLE,SAVE :: lonreg(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: latreg(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: levreg(:) , levallreg(:)
  REAL(r_size),ALLOCATABLE,SAVE :: varallreg(:)
 
 
!  REAL(r_size) :: minlon , minlat , maxlon , maxlat
!  REAL(r_size) :: minlev , maxlev
  REAL(r_size) :: minlonreg , minlatreg , maxlonreg , maxlatreg
  REAL(r_size) :: minlevreg , maxlevreg
  INTEGER,ALLOCATABLE  :: levregmap(:)  !maps the original levels in the regrid levels.

  REAL(r_sngl) :: undefbin
  
!  character(20) , ALLOCATABLE :: var_name(:)
!  INTEGER , ALLOCATABLE :: varlev(:)
!  INTEGER :: nvar_ctl

  CHARACTER(19) :: hdate !Date in input files.

  TYPE(proj_info)  :: projection

!  integer :: nv !Total number of variables (vars times levs)
  integer :: nvreg !Total number of variables in regrid output ( vars times levels)

  logical :: file_exist
  integer :: ierr

  TYPE ctl_info
  REAL(r_size),ALLOCATABLE  :: lev(:)
  REAL(r_size),ALLOCATABLE  :: levall(:)
  REAL(r_size)              :: undefbin
  
  REAL(r_size)              :: minlon,minlat,maxlon,maxlat,minlev,maxlev
  CHARACTER(20),ALLOCATABLE :: varname(:),varnameall(:)
  CHARACTER(200)            :: datafile
  REAL(r_size), ALLOCATABLE :: varlev(:)

  LOGICAL     , ALLOCATABLE :: readvar(:)
  INTEGER     , ALLOCATABLE :: varindex(:)
  INTEGER                   :: nlon,nlat,nlev,nvar,nfields
  REAL(r_size), ALLOCATABLE :: lon1d(:),lat1d(:),lev(:)
  REAL(r_size), ALLOCATABLE :: lon(:,:),lat(:,:)
  REAL(r_size)              :: minlat,minlon,minlev
  REAL(r_size)              :: maxlat,maxlon,maxlev

  REAL(r_size), ALLOCATABLE  :: commonlevall(:)
  CHARACTER(20), ALLOCATABLE :: commonvarnameall(:)
  INTEGER                    :: commonnvar
  
  ENDTYPE

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
SUBROUTINE parse_ctl(ctl_file,myctl)
IMPLICIT NONE
character(*) , INTENT(IN) :: ctl_file
character(100)            :: buffer , dummyc
integer                   :: iunit , i , j , counter
TYPE ctl_info  INTENT(OUT) :: myctl

iunit=98

 INQUIRE(FILE=ctl_file,EXIST=file_exist)
  IF(file_exist) THEN
    OPEN(iunit,FILE=ctl_file,FORM='formatted')
  ELSE
    WRITE(*,*)'[Error]: Missing ctl-file '
    STOP
  ENDIF

 !This is valid for mercator 
 DO
   READ(iunit,'(A100)',IOSTAT=IERR)buffer
      IF(IERR /= 0)THEN
       WRITE(*,*)'End of file'
       EXIT
      ENDIF
      IF( INDEX(buffer,'undef') > 0 .or. INDEX(buffer,'UNDEF') > 0 )THEN
       READ(buffer,*)dummyc , myctl%undefbin
       WRITE(*,*)'Undef value for this file is ',myctl%undefbin
      ENDIF
      
      IF( INDEX(buffer,'xdef') > 0 .or. INDEX(buffer,'XDEF') > 0 )THEN
        READ(buffer,*)dummyc ,  myctl%nlon
        allocate(myctl%lon1d(nlon))
        DO i=1,nlon 
          READ(iunit,*)myctl%lon1d(i)
        ENDDO
      ENDIF
      IF( INDEX(buffer,'ydef') > 0 .or. INDEX(buffer,'YDEF') > 0 )THEN
        READ(buffer,*)dummyc ,  myctl%nlat
        allocate(myctl%lat1d(nlat))
        DO i=1,nlat
          READ(iunit,*)myctl%lat1d(i)
        ENDDO
      ENDIF
      IF( INDEX(buffer,'zdef') > 0 .or. INDEX(buffer,'ZDEF') > 0 )THEN
        READ(buffer,*)dummyc ,  myctl%nlev
        allocate(myctl%lev(nlev))
        DO i=1,myctl%nlev
          READ(iunit,*)myctl%lev(i)
        ENDDO
      ENDIF
      IF( INDEX(buffer,'vars') > 0 .or. INDEX(buffer,'VARS') > 0 .AND. .NOT. &
        ( INDEX(buffer,'endvars') > 0 .or. INDEX(buffer,'ENDVARS') > 0 )  )THEN
        READ(buffer,*)dummyc,myctl%nvar
        allocate(myctl%varname(myctl%nvar) )
        allocate(myctl%varlev(myctl%nvar)  )
        DO i=1,myctl%nvarctl
          READ(iunit,*)myctl%varname(i),myctl%varlev(i)
        ENDDO
        !Get the total number of fields.
        myctl%nfields=SUM(myctl%varlev)
      ENDIF
 ENDDO

 myctl%minlat=myctl%undef
 myctl%maxlat=-myctl%undef
 myctl%minlon=myctl%undef
 myctl%maxlon=-myctl%undef
 
 allocate( myctl%lat(nlon,nlat) , myctl%lon(nlon,nlat) , myctl%lev(nlev) )
 DO i = 1,myctl%nlon
  DO j = 1,myctl%nlat
     myctl%lat(i,j)=myctl%lat1d(j)
     myctl%lon(i,j)=myctl%lon1d(i)

     if ( myctl%lat(i,j) < myctl%minlat )myctl%minlat=myctl%lat(i,j)
     if ( myctl%lat(i,j) > myctl%maxlat )myctl%maxlat=myctl%lat(i,j)
     if ( myctl%lon(i,j) < myctl%minlon )myctl%minlon=myctl%lon(i,j)
     if ( myctl%lon(i,j) > myctl%maxlon )myctl%maxlon=myctl%lon(i,j) 
  ENDDO
 ENDDO


 myctl%minlev=MINVAL(myctl%lev)
 myctl%maxlev=MAXVAL(myctl%lev)

 !Set levall (the level corresponding to each variable)
 ALLOCATE( myctl%levall(myctl%nfields) , myctl%varnameall(myctl%nfields) )
 counter=1
 DO i=1,myctl%nvar
   DO j=1,myctl%varlev(i)
     myctl%levall(counter)=myctl%lev(j)
     myctl%varnameall(counter)=myctl%varname(j)
     counter=counter+1
   ENDDO
 ENDDO

 IF( myctl%nlon == 0 )THEN
  WRITE(*,*)'[Error] Failed to get longitudes.'
  STOP
 ELSEIF( myctl%nlat == 0 )THEN
  WRITE(*,*)'[Error] Failed to get latitudes.'
  STOP
 ELSEIF( myctl%nlev == 0 )THEN
  WRITE(*,*)'[Error] Failed to get levels.'
  STOP
 ELSEIF( myctl%undefbin == 0.0e0 )THEN
  WRITE(*,*)'[Error] Failed to get undef.'
  STOP
 ENDIF

 WRITE(*,*)'===================================================='
 WRITE(*,*)' The total number of individual variables is = ',myctl%nvar
 WRITE(*,*)' Total number of matrices in each file is = ',myctl%nfields
 WRITE(*,*)'===================================================='
 WRITE(*,*)'MINLON = ',myctl%minlon
 WRITE(*,*)'MINLAT = ',myctl%minlat
 WRITE(*,*)'MAXLON = ',myctl%maxlon
 WRITE(*,*)'MAXLAT = ',myctl%maxlat

END SUBROUTINE parse_ctl

!-----------------------------------------------------------------------
! Match ctl variables
!-----------------------------------------------------------------------
! Get the info from two ctls and match variables from ctl1 with those in 
! ctl2.

SUBROUTINE match_ctl_variables(ctl1,ctl2)
IMPLICIT NONE
TYPE ctl_info , INTENT(IN) :: ctl1 , ctl2
INTEGER :: ii , jj , varindex

ALLOCATE( ctl1%readvar(ctl1%nfields) )
ALLOCATE( ctl2%readvar(ctl2%nfields) )
ctl1%readvar=.false.
ctl2%readvar=.false.

ALLOCATE( ctl1%varindex(ctl1%nfields) )
ALLOCATE( ctl2%varindex(ctl2%nfields) )

ctl1%varindex=0
ctl2%varindex=0

 !First identify which variables are common among both ctl data.
 varindex=0
 DO ii=1,ctl1%nfields
    DO jj=1,ctl2%nfields

       IF( ctl1%varnameall(ii) == ctl2%varnameall(jj) .AND. ctl1%levall(ii) == ctl2%levall(jj) )THEN
         varindex=varindex+1
         !We have a match between on variable from ctl1 and ctl2.
 
         !We will read this variable from the binary files.
         ctl1%readvar(ii)=.true.
         ctl2%readvar(jj)=.true.

         !Varindex defines the place where the variable will be stored in the common array.
         ctl1%varindex(ii)=varindex
         ctl2%varindex(jj)=varindex
        
       ENDIF

    ENDDO
 ENDDO 

 !Store the total number of fields in the common grid.
 ctl1%commonnfields=varindex
 ctl2%commonnfields=varindex 

 !Get the commonlevall y commonvarall
 ALLOCATE( ctl1%commonlevall(ctl1%nfields) )
 ALLOCATE( ctl1%commonvarnameall(ctl1%nfields) )
 ALLOCATE( ctl2%commonlevall(ctl2%nfields) )
 ALLOCATE( ctl2%commonvarnameall(ctl2%nfields) )
 

 varindex=0
 DO ii=1,ctl1%nfields
   IF( ctl1%readvar(ii) )THEN
     varindex=varindex+1
     ctl1%commonlevall(varindex)=ctl1%levall(ii)
     ctl1%commonvarnameall(varindex)=ctl%varnameall(ii) 
     ctl2%commonlevall(varindex)=ctl1%levall(ii)
     ctl2%commonvarnameall(varindex)=ctl%varnameall(ii)
   ENDIF
 ENDDO

 WRITE(*,*)"----------------------------------------------------"
 WRITE(*,*)"  MATCH CTL VARIABLES SUMMARY "
 WRITE(*,*)" ii, levall, varnameall, readvar,varindex for CTL1 "
 DO ii=1,ctl1%nfields
   WRITE(*,*)ii,ctl1%levall(ii),ctl1%varnameall(ii),ctl1%readvar(ii),ctl1%varindex(ii)
 ENDDO
 WRITE(*,*)" ii, levall, varnameall, readvar,varindex for CTL2 "
 DO ii=1,ctl2%nfields
   WRITE(*,*)ii,ctl2%levall(ii),ctl2%varnameall(ii),ctl2%readvar(ii),ctl2%varindex(ii)
 ENDDO
 WRITE(*,*)" ii, commonlevall , commonvarnameall "
 DO ii=1,ctl1%commonnfields
   WRITE(*,*)ii,ctl1%commonlevall(ii),ctl1%commonvarnameall
 ENDDO


END SUBROUTINE match_ctl_variables 

!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------

!-- Read a grid file ---------------------------------------------------
!This subroutine reads in only the fields where readvar is true. 
!The array varindex in myctl indicates where to store the output.
SUBROUTINE read_grd(filename,myctl,var,undefmask)
  IMPLICIT NONE
  TYPE ctl_info , INTENT(IN) :: myctl
!  INTEGER,INTENT(IN)      :: ,ny,nz
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: var(myctl%nlon,myctl%nlat,myctl%commonnfields) 
  LOGICAL     ,INTENT(OUT) :: undefmask(myctl%nlon,myctl%nlat,myctl%commonnfields)
  REAL(r_sngl) :: buf4(myctl%nlon,myctl%nlat)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,reclength,ios
  
  undefmask=.true.

  INQUIRE(IOLENGTH=reclength) reclength
  reclength = reclength * myctl%nlon * myctl%nlat

 iunit=33
 INQUIRE(FILE=filename,EXIST=file_exist)
  IF(file_exist) THEN
   OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=reclength)
  ELSE
   undefmask=.false.
   var=0.0d0
   RETURN
  ENDIF

  DO n=1,myctl%nfields
    READ(iunit,rec=n,IOSTAT=ios)buf4

    IF( ios /= 0 .AND. n==1 )THEN
      WRITE(*,*)"[Warning]: Error reading, this might be a new file "
      WRITE(*,*)filename
      var=0.0d0
      undefmask=.false.
      exit
    ELSEIF( ios /=0 .AND. n > 1 )THEN
      WRITE(*,*)"[Error]: Unexpected end of file for file ",filename 
      STOP
    ENDIF

    IF( myctl%readvar(n) )THEN
      WHERE( buf4 == myctl%undefbin )
       undefmask(:,:,myctl%varindex)=.false.
      ENDWHERE
      var(:,:,myctl%varindex) = REAL(buf4,r_size)
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
SUBROUTINE get_regrid_grid(myctl)
IMPLICIT NONE
TYPE ctl_info, INTENT(IN) :: myctl 
INTEGER      :: i,j,ivar,ilev,ilevreg,ilevprev,iregprev

  regrid_vert_res=regrid_vert_res / 100.0d0 !From Pa to hPa.

  nlonreg=CEILING( (myctl%maxlon-myctl%minlon)/regrid_res ) + 1
  nlatreg=CEILING( (myctl%maxlat-myctl%minlat)/regrid_res ) + 1

  ALLOCATE( lonreg(nlonreg,nlatreg) , latreg(nlonreg,nlatreg) )

  DO i = 1,nlonreg
   DO j = 1,nlatreg

    lonreg(i,j)=myctl%minlon + REAL((i-1),r_size)*regrid_res
    latreg(i,j)=myctl%minlat + REAL((j-1),r_size)*regrid_res    

   ENDDO
  ENDDO
 
  minlonreg=myctl%minlon
  minlatreg=myctl%minlat
  maxlonreg=myctl%minlon+regrid_res*REAL(nlonreg-1,r_size)
  maxlatreg=myctl%minlat+regrid_res*REAL(nlatreg-1,r_size)
  
  maxlevreg=myctl%maxlev

  nlevreg= CEILING( (myctl%maxlev-myctl%minlev)/ regrid_vert_res ) + 1

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
  ilevprev=NINT( ( maxlevreg - myctl%commonlevall(1) )  / regrid_vert_res ) + 1
  DO i = 2,myctl%commonnfields
   ilev=NINT( ( maxlevreg - myctl%commonlevall(i) )  / regrid_vert_res ) + 1 
   if ( ilev /= ilevprev .or. myctl%commonvarnameall(i) /= myctl%commonvarname(i-1) )then
      nvreg = nvreg+1
   endif 
   ilevprev=ilev
  ENDDO

  ilevreg=1

  ALLOCATE( levallreg(nvreg) )
  ALLOCATE( varallreg(nvreg) )
  ALLOCATE( levregmap(commonnfields)    )
  ilevprev=NINT( ( maxlevreg - myctl%commonlevall(1) )  / regrid_vert_res ) + 1
  levallreg(1)=myctl%commonlevall(ilevprev) 
  levregmap(1)=ilevprev
  varname(1)=myctl%commonvarnameall(1)
  DO i = 2,commonnfields
   ilev=NINT( ( maxlevreg - myctl%commonlevall(i) )  / regrid_vert_res ) + 1
   if ( ilev /= ilevprev .or. myctl%commonvarnameall(i) /= myctl%commonvarnameall(i-1)  )then
     ilevreg=ilevreg+1
     levallreg(ilevreg)=levreg(ilev)
     varnameallreg(ilevreg)=myctl%commonvarnameall(i)
   endif
   levregmap(i)=ilevreg
   if( myctl%commonvarnameall(i-1) /= myctl%commonvarnameall(i) )then !Variable change
    WRITE(*,*)"NLEVS for var ",myctl%commonvarnameall(varall(i-1))," ",ilevprev
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
