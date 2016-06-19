module common_verification
!=======================================================================
! [PURPOSE:] Common procedure for verification codes
! Juan Ruiz 2016 created.
!=======================================================================
!$USE OMP_LIB
  USE common
  USE map_utils
  IMPLICIT NONE
  PUBLIC

  INTEGER , PARAMETER :: clen=20  !Allocatable character lengths.
  
  !Define ctl_info type
  TYPE ctl_info
  REAL(r_sngl)              :: undefbin
  
  CHARACTER(len=:),ALLOCATABLE :: varname(:),varnameall(:)
  REAL(r_size),ALLOCATABLE  :: lev(:),levall(:)

  CHARACTER(200)            :: datafile
  REAL(r_size), ALLOCATABLE :: varlev(:)

  LOGICAL     , ALLOCATABLE :: readvar(:)
  INTEGER     , ALLOCATABLE :: varindex(:)
  INTEGER                   :: nlon,nlat,nlev,nvar,nfields
  REAL(r_size), ALLOCATABLE :: lon1d(:),lat1d(:)
  REAL(r_size), ALLOCATABLE :: lon(:,:),lat(:,:)
  REAL(r_size)              :: minlat,minlon,minlev
  REAL(r_size)              :: maxlat,maxlon,maxlev
  REAL(r_size)              :: lonres,latres,levres
  TYPE(proj_info)           :: proj_data !Projection info (for Lambert only)
  INTEGER                   :: proj_code !Projection type
  ENDTYPE


CONTAINS


!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------
!Read dimensions from ctl file
!Currently working for mercator, latlon and lambert.
SUBROUTINE parse_ctl(ctl_file,myctl)
IMPLICIT NONE
character(*) , INTENT(IN) :: ctl_file
character(300)            :: buffer , dummyc , proj_name
integer                   :: iunit , i , j , counter , ierr
TYPE (ctl_info) ,  INTENT(OUT) :: myctl
REAL(r_size) :: plat1,plon1,pdx,pdy,pstdlon,ptruelat1,ptruelat2,pknowni,pknownj !Projection parameters.

LOGICAL                    :: ignore_nlon_nlat=.false. , file_exist

iunit=98
myctl%nlon=0
myctl%nlat=0
myctl%nlev=0
myctl%undefbin=0

myctl%proj_code = PROJ_MERC  !We start assuming mercartor.

 INQUIRE(FILE=ctl_file,EXIST=file_exist)
  IF(file_exist) THEN
    OPEN(iunit,FILE=ctl_file,FORM='formatted')
  ELSE
    WRITE(*,*)'[Error]: Missing ctl-file '
    STOP
  ENDIF

 DO
   READ(iunit,'(A300)',IOSTAT=IERR)buffer
      IF(IERR /= 0)THEN
       WRITE(*,*)'End of file'
       EXIT
      ENDIF
      IF( INDEX(buffer,'pdef ') > 0 .or. INDEX(buffer,'PDEF ') > 0 )THEN
     
       ignore_nlon_nlat = .true. !nlon and nlat will be provided by pdef line.
       READ(buffer,*)dummyc,myctl%nlon,myctl%nlat,proj_name
                
        IF( INDEX(proj_name,'LCC') > 0 .or. INDEX(proj_name,'lcc') > 0)THEN
          WRITE(*,*)'This is a Lambert conformal grid'
          !Get data for proj info structure.
          myctl%proj_code=PROJ_LC
          READ(buffer,*)dummyc,myctl%nlon,myctl%nlat,proj_name,plat1,plon1,pknowni,pknownj,   &
   &           ptruelat1,ptruelat2,pstdlon,pdx,pdy
          !Populate projection structure with the data provided by the ctl.
          CALL map_init(myctl%proj_data)
          CALL map_set(myctl%proj_code , myctl%proj_data , lat1=plat1, lon1=plon1, knowni=pknowni,   &
   &           knownj=pknownj, dx=pdx, dy=pdy, stdlon=pstdlon, truelat1=ptruelat1,   &
   &           truelat2=ptruelat2, r_earth=re)  
        ENDIF 
        !Add more projections here.
      ENDIF


      IF( INDEX(buffer,'undef') > 0 .or. INDEX(buffer,'UNDEF') > 0 )THEN
       READ(buffer,*)dummyc , myctl%undefbin
       WRITE(*,*)'Undef value for this file is ',myctl%undefbin
      ENDIF
      
      IF( INDEX(buffer,'xdef') > 0 .or. INDEX(buffer,'XDEF') > 0  )THEN
       IF( .not. ignore_nlon_nlat )THEN
        IF( INDEX(buffer,'levels') > 0 .or. INDEX(buffer,'LEVELS') > 0 )THEN
         READ(buffer,*)dummyc ,  myctl%nlon
         allocate(myctl%lon1d(myctl%nlon))
         DO i=1,myctl%nlon
           READ(iunit,*)myctl%lon1d(i)
         ENDDO
        ELSEIF( INDEX(buffer,'linear') > 0 .or. INDEX(buffer,'linear') > 0)THEN
         READ(buffer,*)dummyc, myctl%nlon,dummyc,myctl%minlon,myctl%lonres
         allocate(myctl%lon1d(myctl%nlon))
         DO i=1,myctl%nlon
           myctl%lon1d(i)=myctl%minlon + (i-1)*myctl%lonres
         ENDDO
        ENDIF
       ENDIF
      ENDIF

      IF( INDEX(buffer,'ydef') > 0 .or. INDEX(buffer,'YDEF') > 0 )THEN
       IF( .not. ignore_nlon_nlat )THEN
        IF( INDEX(buffer,'levels') > 0 .or. INDEX(buffer,'LEVELS') > 0 )THEN
         myctl%proj_code=PROJ_MERC   !Irregular grid in y (assume mercator)
         READ(buffer,*)dummyc ,  myctl%nlat
         allocate(myctl%lat1d(myctl%nlat))
         DO i=1,myctl%nlat
           READ(iunit,*)myctl%lat1d(i)
         ENDDO
        ELSEIF( INDEX(buffer,'linear') > 0 .or. INDEX(buffer,'linear') > 0)THEN
         myctl%proj_code=PROJ_LATLON  !Regular grid in y.
         READ(buffer,*)dummyc, myctl%nlat,dummyc,myctl%minlat,myctl%latres
         allocate(myctl%lat1d(myctl%nlat))
         DO i=1,myctl%nlat
           myctl%lat1d(i)=myctl%minlat + (i-1)*myctl%latres
         ENDDO
        ENDIF
       ENDIF
      ENDIF

      IF( INDEX(buffer,'zdef') > 0 .or. INDEX(buffer,'ZDEF') > 0 )THEN
       IF( INDEX(buffer,'levels') > 0 .or. INDEX(buffer,'LEVELS') > 0 )THEN
        READ(buffer,*)dummyc ,  myctl%nlev
        allocate(myctl%lev(myctl%nlev))
        DO i=1,myctl%nlev
          READ(iunit,*)myctl%lev(i)
        ENDDO
       ELSEIF( INDEX(buffer,'linear') > 0 .or. INDEX(buffer,'linear') > 0)THEN
        READ(buffer,*)dummyc, myctl%nlev,dummyc,myctl%minlev,myctl%levres
        allocate(myctl%lev(myctl%nlev))
        DO i=1,myctl%nlev
          myctl%lev(i)=myctl%minlev + (i-1)*myctl%levres
        ENDDO
       ENDIF
      ENDIF

      IF( ( INDEX(buffer,'vars') > 0 .or. INDEX(buffer,'VARS') > 0 ) .AND. .NOT. &
        ( INDEX(buffer,'endvars') > 0 .or. INDEX(buffer,'ENDVARS') > 0 )  )THEN
        READ(buffer,*)dummyc,myctl%nvar
        allocate(character(clen)::myctl%varname(myctl%nvar) )
        allocate(myctl%varlev(myctl%nvar)  )
        DO i=1,myctl%nvar
          READ(iunit,*)myctl%varname(i),myctl%varlev(i)
          IF( myctl%varlev(i) == 0 )myctl%varlev(i)=1
        ENDDO
        !Get the total number of fields.
        myctl%nfields=SUM(myctl%varlev)
      ENDIF
 ENDDO

 myctl%minlat=undef
 myctl%maxlat=-undef
 myctl%minlon=undef
 myctl%maxlon=-undef
 
 allocate( myctl%lat(myctl%nlon,myctl%nlat) , myctl%lon(myctl%nlon,myctl%nlat) )

 IF( myctl%proj_code == PROJ_MERC .OR. myctl%proj_code == PROJ_LATLON )THEN
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
 ENDIF
 IF( myctl%proj_code == PROJ_LC )THEN
  DO i = 1,myctl%nlon
   DO j = 1,myctl%nlat
     CALL ij_to_latlon(myctl%proj_data, REAL(i,r_size), REAL(j,r_size), myctl%lat(i,j), myctl%lon(i,j) )
     if ( myctl%lat(i,j) < myctl%minlat )myctl%minlat=myctl%lat(i,j)
     if ( myctl%lat(i,j) > myctl%maxlat )myctl%maxlat=myctl%lat(i,j)
     if ( myctl%lon(i,j) < myctl%minlon )myctl%minlon=myctl%lon(i,j)
     if ( myctl%lon(i,j) > myctl%maxlon )myctl%maxlon=myctl%lon(i,j)

   ENDDO
  ENDDO

 ENDIF

 myctl%minlev=MINVAL(myctl%lev)
 myctl%maxlev=MAXVAL(myctl%lev)

 !Set levall (the level corresponding to each variable)
 allocate(character(clen)::myctl%varnameall(myctl%nfields) )
 ALLOCATE( myctl%levall(myctl%nfields) )
 counter=1
 DO i=1,myctl%nvar
   DO j=1,myctl%varlev(i)
     myctl%levall(counter)=myctl%lev(j)
     myctl%varnameall(counter)=myctl%varname(i)
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
 WRITE(*,*)'NLON =   ',myctl%nlon
 WRITE(*,*)'NLAT =   ',myctl%nlat

END SUBROUTINE parse_ctl

SUBROUTINE get_var_index(myctl,varname,varlev,varindex)
IMPLICIT NONE
TYPE(ctl_info),INTENT(IN) :: myctl
CHARACTER(*)  ,INTENT(IN) :: varname
REAL(r_size)  ,INTENT(IN) :: varlev
INTEGER       ,INTENT(OUT):: varindex
integer                   :: i

varindex=0 
DO i=1,myctl%nfields
   IF( myctl%varnameall(i) == varname .and. myctl%levall(i) == varlev )THEN
     varindex=i
   ENDIF
ENDDO

END SUBROUTINE get_var_index

END MODULE common_verification
