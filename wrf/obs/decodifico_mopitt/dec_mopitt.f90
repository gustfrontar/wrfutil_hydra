PROGRAM dec_mopittret
!  USE common
!  USE common_obs_wrf


  IMPLICIT NONE
  include         'hdfeos5.inc'

  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  REAL(r_size),PARAMETER :: pi=3.1415926535d0
  REAL(r_size),PARAMETER :: gg=9.81d0
  REAL(r_size),PARAMETER :: rd=287.0d0
  REAL(r_size),PARAMETER :: cp=1005.7d0
  REAL(r_size),PARAMETER :: re=6371.3d3
  REAL(r_size),PARAMETER :: r_omega=7.292d-5
  REAL(r_size),PARAMETER :: t0c=273.15d0
  REAL(r_sngl),PARAMETER :: undef=-9999.0

  INTEGER,PARAMETER :: nid_obs=7
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
  INTEGER,PARAMETER :: id_tv_obs=3079
  INTEGER,PARAMETER :: id_co_obs=5002
  INTEGER,PARAMETER :: id_totco_obs=5001
!
! surface observations codes > 9999
!
  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_rain_obs=19999
  INTEGER,PARAMETER :: id_tclon_obs=99991
  INTEGER,PARAMETER :: id_tclat_obs=99992
  INTEGER,PARAMETER :: id_tcmip_obs=99993

  LOGICAL           :: file_exists
  LOGICAL,PARAMETER :: verbose = .FALSE.
  INTEGER,PARAMETER :: mopitttype = 31 ! observation type record for MOPITT
  REAL(r_size),PARAMETER :: lon1 =  270.0d0 ! minimum longitude
  REAL(r_size),PARAMETER :: lon2 =  340.0d0 ! maximum longitude
  REAL(r_size),PARAMETER :: lat1 =  -70.0d0 ! minimum latitude
  REAL(r_size),PARAMETER :: lat2 =   20.0d0 ! maximum latitude

  ! use full data in the following domain.
  ! (use thinned data in other domain (iskip)).
  REAL(r_size),PARAMETER :: loni1 =  114.40d0 ! minimum longitude (inner)
  REAL(r_size),PARAMETER :: loni2 =  134.06d0 ! maximum longitude (inner)
  REAL(r_size),PARAMETER :: lati1 =  11.63d0 ! minimum latitude (inner)
  REAL(r_size),PARAMETER :: lati2 =  29.91d0 ! maximum latitude (inner)

  REAL(r_size),PARAMETER :: coef_error = 1.0d0  
  INTEGER :: iskip = 3
  INTEGER :: nunit=70
  INTEGER*8 :: nx
  INTEGER*8 :: nz
  INTEGER*8 :: nvar
  INTEGER*8 :: start(3),stride(3),edge(3)
  INTEGER*8 :: statn
  INTEGER*8 :: fid
  INTEGER*8 :: nswath
  INTEGER*8 :: nchar
  INTEGER*8 :: swid
  INTEGER*8 :: he5_swopen,he5_swinqswath,he5_swattach,he5_swrdattr,he5_swrdfld,he5_swdetach,he5_swclose,he5_swinqdims
  INTEGER*8 :: he5_swfldinfo
  CHARACTER(256) :: swathname
  CHARACTER(100) :: filename
  REAL(r_sngl),ALLOCATABLE :: rlon(:)
  REAL(r_sngl),ALLOCATABLE :: rlat(:)
  REAL(r_size),ALLOCATABLE :: time(:) 
  REAL(r_size) ::  timemin , timemax
  REAL(r_sngl),ALLOCATABLE :: levco(:)
  REAL(r_sngl),ALLOCATABLE :: co(:,:,:),apriorico(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: surfacepressure(:),sdf(:)
  REAL(r_sngl),ALLOCATABLE :: averaging_kernel(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: pressuregrid(:)
  INTEGER*4  ,ALLOCATABLE :: swathindex(:,:),alongtrack(:)
  INTEGER :: i,j,k
  REAL(r_sngl) :: wk(7)
  INTEGER :: iy,im,id,ih,imin
  INTEGER :: iy1,im1,id1,ih1,imin1
  INTEGER :: iy2,im2,id2,ih2,imin2
  REAL(r_sngl) :: sec4
  REAL(r_size) :: sec
  LOGICAL :: ex
  CHARACTER(14) :: outfile='yyyymmddhh.dat'
!                           123456789012345
  CHARACTER(360) :: dimname
  INTEGER,ALLOCATABLE      :: NCO(:)
  INTEGER*8 :: dims(32) , ndims , iunit , ii 
  REAL(r_size) :: tile_lat_mean , tile_lon_mean 
  INTEGER*8      :: datafr(32) 
  INTEGER*8      :: datafdm(32)
  INTEGER*8      :: dataft(1)
  CHARACTER(512) :: datafn , datafmn
  INTEGER :: TIMEDATACOUNT(25)

! FOR DATA THINING
  LOGICAL,PARAMETER :: GOING = .TRUE. , RETURNING=.FALSE.
  INTEGER,PARAMETER :: SKIPACROSS=2 , SKIPALONG=2
 
!=======================================================================
! READING HDF
!=======================================================================
  tile_lat_mean=0.0d0
  tile_lon_mean=0.0d0
  !We will use timedatacount to store the number of profiles that will be
  !asigned to each file.
  TIMEDATACOUNT=0

  !
  ! open
  !
  CALL GETARG( 1, filename )
  !READ(5,'(A)') filename
  PRINT *,filename
  !PRINT *,filename
  fid = he5_swopen(filename,HE5F_ACC_RDONLY )
  IF(fid == -1) THEN
    PRINT *,'IO ERROR: opening file ',filename
    STOP
  END IF
  WRITE(*,*)"FILE ID=",fid
  swid = he5_swattach(fid,'MOP02') 
  !
  ! read
  !

  !Get the dimensions. While some dimensions are fixed (like the number of
  !vertical leveles), the number of valid profiles in a file changes from day to
  !day (since only valid cloud free profiles are incorporated into the files).
  statn = he5_swfldinfo(swid, 'RetrievedCOMixingRatioProfile', datafr,datafdm,dataft,datafn,datafmn)
     nx=datafdm(3)
     nz=datafdm(2)
     nvar=datafdm(1)

  ALLOCATE( rlon(nx) , rlat(nx) , time(nx) , levco(nz) , co(nvar,nz+1,nx) , NCO(nz+1) )
  ALLOCATE( apriorico(nvar,nz+1,nx ) )
  ALLOCATE( surfacepressure(nx) , pressuregrid(nz+1) , sdf(nx) )
  ALLOCATE( averaging_kernel(nz+1,nz+1,nx) )
  ALLOCATE( swathindex(3,nx),alongtrack(nx) )
  NCO=0 
  WRITE(*,*)"NX = ",nx," NZ =",nz," NVAR = ",nvar

  start=0
  stride=1
  edge(1)=nx
  edge(2)=nz
  edge(3)=nvar
  statn = he5_swrdfld(swid,'Longitude',start(1),stride(1),edge(1),rlon)
  statn = he5_swrdfld(swid,'Latitude',start(1),stride(1),edge(1),rlat)
  statn = he5_swrdfld(swid,'Time',start(1),stride(1),edge(1),time)
  edge(1)=nz
  statn = he5_swrdfld(swid,'PressureGrid',start(1),stride(1),edge(1),pressuregrid(2:nz+1))
  edge(1)=nvar
  edge(2)=nz
  edge(3)=nx
  !Read the retrieved and a priori mixing ratio vertical profiles.
  statn =he5_swrdfld(swid,'RetrievedCOMixingRatioProfile',start,stride,edge,co(1:nvar,2:nz+1,nx))
  statn=he5_swrdfld(swid,'APrioriCOMixingRatioProfile',start,stride,edge,apriorico(1:nvar,2:nz+1,nx))

  edge(1)=nvar
  edge(2)=nx
  !Read the surface retrieved and a priori mixing ratio values.
  statn=he5_swrdfld(swid,'RetrievedCOSurfaceMixingRatio',start(1),stride(1),edge(1),co(1:nvar,1,nx))
  statn=he5_swrdfld(swid,'APrioriCOSurfaceMixingRatio',start(1),stride(1),edge(1),apriorico(1:nvar,1,nx))

  edge(1)=nx
  !Get the surface pressure and degrees of freedom.
  statn =he5_swrdfld(swid,'SurfacePressure',start(1),stride(1),edge(1),surfacepressure) 
  statn =he5_swrdfld(swid,'DegreesofFreedomforSignal',start(1),stride(1),edge(1),sdf)

  edge(1)=nz+1
  edge(2)=nz+1
  edge(3)=nx
  statn=he5_swrdfld(swid,'RetrievalAveragingKernelMatrix',start,stride,edge,averaging_kernel)
  edge(1)=3
  edge(2)=nx
  statn=he5_swrdfld(swid,'SwathIndex',start(2),stride(2),edge(2),swathindex)
  alongtrack=swathindex(3,:)*4+swathindex(1,:)  
  !WRITE(*,*)co(1,:,1)
  !WRITE(*,*)co(2,:,1)
  !STOP
  !
  ! close
  !

  timemin = MINVAL(time)
  timemax = MAXVAL(time)

  

  statn = he5_swdetach(swid)
  statn = he5_swclose(fid)
!=======================================================================
! QC -> OUTPUT LETKF OBS FORMAT
!=======================================================================

!  Lets open files corresponding to different times.

!MOPITT data comes in one file per day. We split these data into hourly files
!data in these ourly files are centered at 00 minutes. 
!For example file for 03UTC contains data from 2:30 UTC to 3:30 UTC.

  CALL com_tai2utc(timemin,iy1,im1,id1,ih1,imin1,sec)
  !We open the file corresponding to the previous day
 DO ii=0,24
  !First form file name.
  IF( ii < 24)THEN
    CALL com_tai2utc(timemin,iy1,im1,id1,ih1,imin1,sec)
    WRITE(outfile(1:10),'(I4.4,3I2.2)') iy1,im1,id1,ii
  ELSE
    CALL com_tai2utc(timemax+1800,iy1,im1,id1,ih1,imin1,sec)
    WRITE(outfile(1:10),'(I4.4,3I2.2)') iy1,im1,id1,00
  ENDIF
  !If file exist then open in append position. Else initialize a new file.
  INQUIRE(FILE=outfile,EXIST=file_exists)
  IF( file_exists )THEN
    OPEN(nunit+ii,FILE=outfile,FORM='unformatted',ACCESS='sequential')
    READ(nunit+ii)TIMEDATACOUNT(ii+1)
    CLOSE(nunit+ii)
    WRITE(*,*)"OPENING THE FILE ",outfile
    OPEN(nunit+ii,FILE=outfile,FORM='unformatted',ACCESS='sequential',position='APPEND')
  ELSE
    OPEN(nunit+ii,FILE=outfile,FORM='unformatted',ACCESS='sequential')
    WRITE(nunit+ii)TIMEDATACOUNT(ii+1),id_co_obs,mopitttype
    WRITE(*,*)"CREATING THE FILE ",outfile
  ENDIF

 ENDDO

  !We loop over all the profiles. And write them to the corresponding hourly
  !files.
    DO i=1,nx
   
      !First we check that the profile is within the domain selected by
      !the user.
      IF(rlon(i) < 0) rlon(i) = rlon(i) + 360.0d0
      IF(rlon(i) < lon1) CYCLE
      IF(rlon(i) > lon2) CYCLE
      IF(rlat(i) < lat1) CYCLE
      IF(rlat(i) > lat2) CYCLE

      !Apply thining options.   
      !Thinning strategy is based on MOPITT scan strategy and in the variable
      !SWATHINDEX variable. See
      !http://www.acom.ucar.edu/mopitt/mopitt-qa-plan-v2.shtml for more
      !information on the scan strategy. 
 
      !Filter by crosstrack direction (going or returning)
      IF( .NOT. GOING .AND. swathindex(2,i) < 16)CYCLE
      IF( .NOT. RETURNING .AND. swathindex(2,i) >= 16)CYCLE

      !Apply skiping rules for the along track and cross track
      !directions
      IF( MOD(swathindex(2,i),SKIPACROSS) /= 0.0)CYCLE
      IF( MOD(alongtrack(i),SKIPALONG) /= 0.0)CYCLE

      WRITE(*,*)swathindex(:,i),alongtrack(i)

      !Now we keep only 
      
            
      CALL com_tai2utc(REAL(time(i),r_size),iy,im,id,ih,imin,sec)
           IF( imin > 30 )ih = ih + 1
           iunit=nunit + ih 
      TIMEDATACOUNT(ih+1)=TIMEDATACOUNT(ih+1)+1

      
      WRITE(iunit)rlon(i),rlat(i),time(i)
      WRITE(iunit)nz+1
      WRITE(iunit)co(1,:,i)
      WRITE(iunit)co(2,:,i)
      pressuregrid(1)=surfacepressure(i)
      WRITE(iunit)pressuregrid
      DO ii=1,nz+1
         WRITE(iunit)averaging_kernel(ii,:,i)
      ENDDO
      WRITE(iunit)sdf(i)
      WRITE(iunit)nz+1
      WRITE(iunit)apriorico(1,:,i)
      WRITE(iunit)pressuregrid

  !That's it, lets proceed to the following profile.
  END DO

!We finish processing all the profiles. Before the end we will
!update the number of profiles stored in each binary file.

  !We open the file corresponding to the previous day
  DO ii=0,24
   REWIND(nunit+ii)
   WRITE(nunit+ii)TIMEDATACOUNT(ii+1)
  ENDDO

  !Write a summary...

  WRITE(*,*)"THE NUMBER OF PROFILES IN EACH FILE IS:"
  WRITE(*,*)TIMEDATACOUNT


  STOP
END PROGRAM dec_mopittret

!-----------------------------------------------------------------------
! UTC to TAI93
!-----------------------------------------------------------------------
SUBROUTINE com_utc2tai(iy,im,id,ih,imin,sec,tai93)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: iy,im,id,ih,imin
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  REAL(r_size),INTENT(IN) :: sec
  REAL(r_size),INTENT(OUT) :: tai93
  REAL(r_size),PARAMETER :: mins = 60.0d0
  REAL(r_size),PARAMETER :: hour = 60.0d0*mins
  REAL(r_size),PARAMETER :: day = 24.0d0*hour
  REAL(r_size),PARAMETER :: year = 365.0d0*day
  INTEGER,PARAMETER :: mdays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER :: days,i

  tai93 = REAL(iy-1993,r_size)*year + FLOOR(REAL(iy-1993)/4.0,r_size)*day
  days = id -1
  DO i=1,12
    IF(im > i) days = days + mdays(i)
  END DO
  IF(MOD(iy,4) == 0 .AND. im > 2) days = days + 1 !leap year
  tai93 = tai93 + REAL(days,r_size)*day + REAL(ih,r_size)*hour &
              & + REAL(imin,r_size)*mins + sec
  IF(iy > 1993 .OR. (iy==1993 .AND. im > 6)) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 1994 .OR. (iy==1994 .AND. im > 6)) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 1995) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 1997 .OR. (iy==1997 .AND. im > 6)) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 1998) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 2005) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 2008) tai93 = tai93 + 1.0d0 !leap second

  RETURN
END SUBROUTINE com_utc2tai
!-----------------------------------------------------------------------
! TAI93 to UTC
!-----------------------------------------------------------------------
SUBROUTINE com_tai2utc(tai93,iy,im,id,ih,imin,sec)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER,PARAMETER :: n=7 ! number of leap seconds after Jan. 1, 1993
  INTEGER,PARAMETER :: leapsec(n) = (/  15638399,  47174400,  94608001,&
                                  &    141868802, 189302403, 410227204,&
                                  &    504921605/)
  REAL(r_size),INTENT(IN) :: tai93
  INTEGER,INTENT(OUT) :: iy,im,id,ih,imin
  REAL(r_size),INTENT(OUT) :: sec
  REAL(r_size),PARAMETER :: mins = 60.0d0
  REAL(r_size),PARAMETER :: hour = 60.0d0*mins
  REAL(r_size),PARAMETER :: day = 24.0d0*hour
  REAL(r_size),PARAMETER :: year = 365.0d0*day
  INTEGER,PARAMETER :: mdays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  REAL(r_size) :: wk,tai
  INTEGER :: days,i,leap

  tai = tai93
  sec = 0.0d0
  DO i=1,n
    IF(FLOOR(tai93) == leapsec(i)+1) sec = 60.0d0 + tai93-FLOOR(tai93,r_size)
    IF(FLOOR(tai93) > leapsec(i)) tai = tai -1.0d0
  END DO
  iy = 1993 + FLOOR(tai /year)
  wk = tai - REAL(iy-1993,r_size)*year - FLOOR(REAL(iy-1993)/4.0,r_size)*day
  IF(wk < 0.0d0) THEN
    iy = iy -1
    wk = tai - REAL(iy-1993,r_size)*year - FLOOR(REAL(iy-1993)/4.0,r_size)*day
  END IF
  days = FLOOR(wk/day)
  wk = wk - REAL(days,r_size)*day
  im = 1
  DO i=1,12
    leap = 0
    IF(im == 2 .AND. MOD(iy,4)==0) leap=1
    IF(im == i .AND. days >= mdays(i)+leap) THEN
      im = im + 1
      days = days - mdays(i)-leap
    END IF
  END DO
  id = days +1

  ih = FLOOR(wk/hour)
  wk = wk - REAL(ih,r_size)*hour
  imin = FLOOR(wk/mins)
  IF(sec < 60.0d0) sec = wk - REAL(imin,r_size)*mins

  RETURN
END SUBROUTINE com_tai2utc


