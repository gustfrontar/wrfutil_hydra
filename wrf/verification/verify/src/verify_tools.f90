module verify_tools
!=======================================================================
!
! [PURPOSE:] Common procedure for verification against grided data
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE map_utils
  USE common_verification
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
 
  INTEGER :: nlon   ,nlat   ,nfields
  INTEGER :: nlonreg,nlatreg,nfieldsreg,nlevreg !Regrid domain size

  !NAMELIST INPUT

  INTEGER :: nbv=1  !Number of ensemble members.
  real(r_size) :: regrid_res = 1.0d0          ! Regrid output resolution (degree)
  real(r_size) :: regrid_vert_res = 10000.0d0 ! Regrid output vertical resolution (Pa)
  logical      :: regrid_output = .true.
 
  character(50) :: inputendian='big_endian'
  character(50) :: outputendian='big_endian'



  integer, parameter :: max_narea = 20
  integer :: narea 
  real(r_size) :: vlon1(max_narea) 
  real(r_size) :: vlon2(max_narea) 
  real(r_size) :: vlat1(max_narea) 
  real(r_size) :: vlat2(max_narea)

  CHARACTER(LEN=20) :: NAMELIST_FILE='./verify.namelist'

  !GRID INFO
  REAL(r_size),ALLOCATABLE      :: lon(:,:),lat(:,:)
  REAL(r_size),ALLOCATABLE      :: lonreg(:,:)
  REAL(r_size),ALLOCATABLE      :: latreg(:,:)
  REAL(r_size),ALLOCATABLE      :: levreg(:) , levallreg(:)
  CHARACTER(20),ALLOCATABLE     :: varnameallreg(:)
 
  REAL(r_size) :: minlonreg , minlatreg , maxlonreg , maxlatreg
  REAL(r_size) :: minlevreg , maxlevreg
  INTEGER,ALLOCATABLE  :: levregmap(:)  !maps the original levels in the regrid levels.


  REAL(r_size), ALLOCATABLE  :: commonlevall(:)
  CHARACTER(20), ALLOCATABLE :: commonvarnameall(:)
  INTEGER                    :: commonnfields
  INTEGER,ALLOCATABLE        :: commonmap1(:),commonmap2(:) !Maps original var distribution to the common var distribution.

!  integer :: nv !Total number of variables (vars times levs)
  integer :: nvreg !Total number of variables in regrid output ( vars times levels)

  logical :: file_exist
  integer :: ierr

CONTAINS

!-----------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------

SUBROUTINE read_namelist
NAMELIST / GENERAL / nbv , narea , vlon1 , vlon2 , vlat1 , vlat2 , regrid_output &
                  ,  regrid_res , regrid_vert_res , inputendian , outputendian



  !Default values for verifcation regions.
                    ! GLOBE    NH     TROPICS   SH
narea=4
vlon1(1:narea) = (/  0.d0,   0.d0,   0.d0,   0.d0/)
vlon2(1:narea) = (/360.d0, 360.d0, 360.d0, 360.d0/)
vlat1(1:narea) = (/-90.d0,  20.d0, -20.d0, -90.d0/)
vlat2(1:narea) = (/ 90.d0,  90.d0,  20.d0, -20.d0/)


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
! Match ctl variables
!-----------------------------------------------------------------------
! Get the info from two ctls and match variables from ctl1 with those in 
! ctl2.

SUBROUTINE match_ctl_variables(ctl1,ctl2,commonmap1,commonmap2)
IMPLICIT NONE
TYPE(ctl_info) , INTENT(INOUT) :: ctl1 , ctl2
INTEGER , INTENT(OUT)  :: commonmap1(ctl1%nfields) , commonmap2(ctl2%nfields) 
INTEGER :: ii , jj , varindex


ALLOCATE( ctl1%readvar(ctl1%nfields) )
ALLOCATE( ctl2%readvar(ctl2%nfields) )
ctl1%readvar=.false.
ctl2%readvar=.false.

ALLOCATE( ctl1%varindex(ctl1%nfields) )
ALLOCATE( ctl2%varindex(ctl2%nfields) )

commonmap1=0
commonmap2=0

ctl1%varindex=0
ctl2%varindex=0

 !First identify which variables are common among both ctl data.
 commonnfields=0
 DO ii=1,ctl2%nfields
   CALL get_var_index(ctl1,ctl2%varnameall(ii),ctl2%levall(ii),varindex)
   !We will read this variable from the binary files.
     if( varindex > 0 )then
       commonnfields=commonnfields+1
       ctl2%readvar(ii)=.true.
       ctl1%readvar(varindex)=.true.
       commonmap2(ii)=commonnfields
       commonmap1(varindex)=commonnfields
     endif
 ENDDO 

 !Get the commonlevall y commonvarall
 ALLOCATE( commonlevall(commonnfields) )
 ALLOCATE( commonvarnameall(commonnfields) )
 

 varindex=0
 DO ii=1,ctl1%nfields
   IF( ctl1%readvar(ii) )THEN
     varindex=varindex+1
     commonlevall(varindex)=ctl1%levall(ii)
     commonvarnameall(varindex)=ctl1%varnameall(ii) 
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
 DO ii=1,commonnfields
   WRITE(*,*)ii,commonlevall(ii),commonvarnameall
 ENDDO


END SUBROUTINE match_ctl_variables 

SUBROUTINE read_grd(filename,nx,ny,nz,var,undefmask,undefbin)
  IMPLICIT NONE
  INTEGER,INTENT(IN)      :: nx,ny,nz
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: var(nx,ny,nz)
  LOGICAL     ,INTENT(OUT) :: undefmask(nx,ny,nz)
  REAL(r_sngl) :: buf4(nx,ny)
  REAL(r_sngl),INTENT(IN) :: undefbin
  INTEGER :: iunit,iolen
  INTEGER :: k,n,reclength,ios
  LOGICAL :: file_exist


  undefmask=.true.

  INQUIRE(IOLENGTH=reclength) reclength
  reclength = reclength * nx * ny


 iunit=33
 INQUIRE(FILE=filename,EXIST=file_exist)
  IF(file_exist) THEN
   OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=reclength,CONVERT=inputendian)
  ELSE
   WRITE(*,*)"[Warning]: Could not finde file! ",filename
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


SUBROUTINE write_grd(filename,nx,ny,nz,var,undefmask,undefbin)
  IMPLICIT NONE
  INTEGER, INTENT(IN)     :: nx,ny,nz
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: var(nx,ny,nz)
  LOGICAL     ,INTENT(IN) :: undefmask(nx,ny,nz)
  REAL(r_sngl),INTENT(IN) :: undefbin
  REAL(r_sngl) :: buf4(nx,ny)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,reclength

  INQUIRE(IOLENGTH=reclength) reclength
  reclength = reclength * nx * ny


  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=reclength,CONVERT=outputendian)

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
SUBROUTINE get_regrid_grid(minlon,maxlon,minlat,maxlat,minlev,maxlev,varnameall,varlevall,nfields)
IMPLICIT NONE
INTEGER      , INTENT(IN) :: nfields
REAL(r_size) , INTENT(IN) :: minlon,maxlon,minlat,maxlat,minlev,maxlev
CHARACTER(*) , INTENT(IN) :: varnameall(nfields)
REAL(r_size) , INTENT(IN) :: varlevall(nfields)

INTEGER      :: i,j,ivar,ilev,ilevreg,ilevprev,iregprev

  regrid_vert_res=regrid_vert_res / 100.0d0 !From Pa to hPa.

  nlonreg=CEILING( (maxlon-minlon)/regrid_res ) + 1
  nlatreg=CEILING( (maxlat-minlat)/regrid_res ) + 1

  ALLOCATE( lonreg(nlonreg,nlatreg) , latreg(nlonreg,nlatreg) )

  !Set horizontal grid points
  DO i = 1,nlonreg
   DO j = 1,nlatreg
    lonreg(i,j)=minlon + REAL((i-1),r_size)*regrid_res
    latreg(i,j)=minlat + REAL((j-1),r_size)*regrid_res    
   ENDDO
  ENDDO
 
  minlonreg=MINVAL(lonreg(:,1))
  minlatreg=MINVAL(latreg(1,:))
  maxlonreg=MAXVAL(lonreg(:,1))
  maxlatreg=MINVAL(latreg(1,:))
  
  !Set vertical levels
  nlevreg= CEILING( (maxlev-minlev)/ regrid_vert_res ) + 1
  ALLOCATE( levreg(nlevreg) ) 
  DO i = 1,nlevreg
    levreg(i)=maxlevreg - REAL((i-1),r_size)*regrid_vert_res 
  ENDDO

  minlevreg=MINVAL(levreg)
  maxlevreg=MAXVAL(levreg)

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
 
  !How many total levels do we have in the regird domain?
  !Get the number of total fields after regrid in z.
  nvreg=1
  ilevprev=NINT( ( maxlevreg - varlevall(1) )  / regrid_vert_res ) + 1
  DO i = 2,nfields
   ilev=NINT( ( maxlevreg - varlevall(i) )  / regrid_vert_res ) + 1 
   if ( ilev /= ilevprev .or. varnameall(i) /= varnameall(i-1) )then
      nvreg = nvreg+1
   endif 
   ilevprev=ilev
  ENDDO

  ALLOCATE( levallreg(nvreg) )
  ALLOCATE( varnameallreg(nvreg) )
  ALLOCATE( levregmap(nfields)  )

  !Define the map between the original grid and the regrided grid in the variable dimension.
  !Find the name and lev for each field after regrid.
  ilevreg=1
  ilevprev=NINT( ( maxlevreg - varlevall(1) )  / regrid_vert_res ) + 1
  levallreg(1)=levreg(ilevprev) 
  levregmap(1)=1
  varnameallreg(1)=varnameall(1)
  DO i = 2,nfields
   ilev=NINT( ( maxlevreg - varlevall(i) )  / regrid_vert_res ) + 1
   if ( ilev /= ilevprev .or. varnameall(i) /= varnameall(i-1)  )then
     ilevreg=ilevreg+1
     levallreg(ilevreg)=levreg(ilev)
     varnameallreg(ilevreg)=varnameall(i)
   endif
   levregmap(i)=ilevreg

   ilevprev=ilev
  ENDDO

END SUBROUTINE get_regrid_grid

!Regrid var 
!This subroutine regrids -----------------------------------------------------------
SUBROUTINE regrid_var(nx,ny,nz,nxreg,nyreg,nzreg,zmap,x,y,xreg,yreg,var,undefmask,varreg,undefmaskreg,num)
IMPLICIT NONE
INTEGER     , INTENT(IN) :: nx,ny,nz,nxreg,nyreg,nzreg
REAL(r_size), INTENT(IN) :: x(nx,ny),y(nx,ny),xreg(nx,ny),yreg(nx,ny)
REAL(r_size), INTENT(IN) :: var(nx,ny,nz) 
LOGICAL     , INTENT(IN) :: undefmask(nx,ny,nz)
LOGICAL     , INTENT(OUT):: undefmaskreg(nxreg,nyreg,nzreg)
INTEGER     , INTENT(IN) :: zmap(nz)
REAL(r_size), INTENT(OUT):: varreg(nxreg,nyreg,nzreg)
INTEGER     , INTENT(OUT):: num(nxreg,nyreg,nzreg)
INTEGER                  :: i , j , k , ireg , jreg , kreg
REAL(r_size)             :: dx , dy , minx , miny 

   varreg=0.0d0
   num=0
   undefmaskreg=.true.

   dx = xreg(2,1) - xreg(1,1) !Assume regular grid for the regrid data.
   dy = yreg(1,2) - yreg(1,1) !

   minx=MINVAL(xreg(:,1))
   miny=MINVAL(yreg(1,:))

   DO i = 1,nx
    DO j = 1,ny
     DO k = 1,nz
      IF( undefmask(i,j,k) )THEN
        ireg= NINT( ( x(i,j)-minx ) / dx  ) + 1
        jreg= NINT( ( y(i,j)-miny ) / dy  ) + 1
        kreg= zmap(k)       !Since k dimension includes variables an levels we
                            !have to define a map to know which original level
                            !correspond to the regrid levels.
        varreg(ireg,jreg,kreg)=varreg(ireg,jreg,kreg)+var(i,j,k) 
        num(ireg,jreg,kreg)   =num(ireg,jreg,kreg)   +1
       ENDIF
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
