MODULE obsop_tools
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_wrf
  USE common_obs_wrf
  USE common_mpi_wrf
  USE map_utils
  USE common_namelist

  IMPLICIT NONE
  PUBLIC

  INTEGER,SAVE :: nobs , nobsradar
  REAL(r_size),ALLOCATABLE,SAVE :: obselm(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslon(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslev(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obserr(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obstyp(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslot(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsi(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsj(:)
!  REAL(r_size),ALLOCATABLE,SAVE :: obsk(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdep(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obshdxf(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsradar(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsmean(:) !For verification
  REAL(r_size),ALLOCATABLE,SAVE :: obssprd(:) !For verification
  REAL(r_size),ALLOCATABLE,SAVE :: obsaz(:),obsra(:),obsel(:)
  INTEGER,ALLOCATABLE,SAVE :: nobsgrd(:,:)


  !Obsgrid variables
  REAL(r_size) :: minlon,maxlon,minlat,maxlat
  REAL(r_size), ALLOCATABLE :: lonreg(:,:),latreg(:,:),levreg(:),levzreg(:)

  REAL(r_size) :: minlonreg , minlatreg , maxlonreg , maxlatreg
  REAL(r_size) :: minlevreg , maxlevreg 
  REAL(r_size) :: minzlevreg, maxzlevreg

  INTEGER , ALLOCATABLE :: varidall(:)
  REAL(r_size) , ALLOCATABLE :: levall(:)

  INTEGER,SAVE :: nlonreg,nlatreg !Regrid domain size


  integer :: nvreg !Total number of variables in regrid output ( vars times levels)

  !Filenames, variable_initialslot_endslot.grd/txt
  character(20) :: meanfile = 'obsmeanXX_XX.grd' !Ensemble mean
  character(20) :: sprdfile = 'obssprdXX_XX.grd' !Ensemble spread
  character(20) :: merrfile = 'obsmerrXX_XX.grd' !Ensemble mean error
  character(20) :: onumfile = 'obsonumXX_XX.grd' !Number of obs in each grid.

  character(20) :: rmsefile = 'ormseXX_XX.grd'
  character(20) :: biasfile = 'obiasXX_XX.grd'
  character(20) :: absefile = 'oabseXX_XX.grd'
  character(20) :: counfile = 'ocounXX_XX.grd'
  character(20) :: areafile = 'oareaXX_XX.txt'

  logical :: file_exist

  integer :: nv3dver , nv2dver  , nvver
  integer , allocatable :: varidver(:) 

CONTAINS
!-----------------------------------------------------------------------
! Initialize
!-----------------------------------------------------------------------
SUBROUTINE set_letkf_obs
  IMPLICIT NONE
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:,:) 
! <--- QC for tyinfo
  REAL(r_size) :: dz,tg,qg
  REAL(r_size) :: dlon1,dlon2,dlon,dlat
  REAL(r_size),ALLOCATABLE :: tmpelm(:)
  REAL(r_size),ALLOCATABLE :: tmplon(:)
  REAL(r_size),ALLOCATABLE :: tmplat(:)
  REAL(r_size),ALLOCATABLE :: tmplev(:)
  REAL(r_size),ALLOCATABLE :: tmpdat(:)
  REAL(r_size),ALLOCATABLE :: tmperr(:)
  REAL(r_size),ALLOCATABLE :: tmptyp(:)
  REAL(r_size),ALLOCATABLE :: tmplot(:)
  REAL(r_size),ALLOCATABLE :: tmpi(:)
  REAL(r_size),ALLOCATABLE :: tmpj(:)
  REAL(r_size),ALLOCATABLE :: tmpk(:)
  REAL(r_size),ALLOCATABLE :: tmpsprd(:)
  REAL(r_size),ALLOCATABLE :: tmpdep(:)
  REAL(r_size),ALLOCATABLE :: tmpradar(:)
  REAL(r_size),ALLOCATABLE :: tmpaz(:),tmpel(:),tmpra(:)
  REAL(r_size),ALLOCATABLE :: tmphdxf(:,:)
  REAL(r_size),ALLOCATABLE :: tmp2elm(:)
  REAL(r_size),ALLOCATABLE :: tmp2lon(:)
  REAL(r_size),ALLOCATABLE :: tmp2lat(:)
  REAL(r_size),ALLOCATABLE :: tmp2lev(:)
  REAL(r_size),ALLOCATABLE :: tmp2dat(:)
  REAL(r_size),ALLOCATABLE :: tmp2err(:)
  REAL(r_size),ALLOCATABLE :: tmp2typ(:)
  REAL(r_size),ALLOCATABLE :: tmp2lot(:)
  REAL(r_size),ALLOCATABLE :: tmp2i(:)
  REAL(r_size),ALLOCATABLE :: tmp2j(:)
!  REAL(r_size),ALLOCATABLE :: tmp2k(:)
  REAL(r_size),ALLOCATABLE :: tmp2dep(:)
  REAL(r_size),ALLOCATABLE :: tmp2hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: tmp2radar(:)
  REAL(r_size),ALLOCATABLE :: tmp2az(:),tmp2el(:),tmp2ra(:)
  INTEGER , ALLOCATABLE :: nobslotsradar(:,:)
  INTEGER :: nobslots(nslots) !, nobslotsradar(nslots,nradar)
  INTEGER :: n,i,j,ierr,islot,nn,l,im, nini,nend
  INTEGER :: nj(0:nlat-1)
  INTEGER :: njs(1:nlat-1)
  CHARACTER(9) :: obsfile='obsTT.dat'
  CHARACTER(11) :: guesfile='gsTTNNNNN'
  CHARACTER(11) :: radarfile='radTTtt.dat'  !radTimeTimeTypeType.dat
  real(r_size) :: ratio
  LOGICAL      :: ex

  INTEGER      :: ivar
  INTEGER      :: iradar  , current_radar
  REAL(r_size) :: zmodel(nlon,nlat,nlev)

  minref=10.0d0**( minrefdbz / 10.0d0)

  WRITE(6,'(A)') 'Hello from set_letkf_obs'

  nobslots=0
  DO islot=1,nslots
    WRITE(obsfile(4:5),'(I2.2)') islot
    INQUIRE(FILE=obsfile,EXIST=ex)
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',obsfile
     IF( ex )THEN
     CALL get_nobs(obsfile,nobslots(islot))
     ELSE
      WRITE(6,*)'WARNING: FILE ',obsfile,' NOT FOUND'
     ENDIF
  END DO
  nobs = SUM(nobslots)

  !Allocate radar variables
  ALLOCATE( radar_lon(nradar),radar_lat(nradar),radar_z(nradar),radar_mindbz(nradar) )
  ALLOCATE( nobslotsradar(nslots,nradar) )

  !Get the number of radar observations.
  nobsradar=0

  !Take into account derived pseudo observations that 
  !are not in the radar file.
  
  DO islot=1,nslots
   DO iradar=1,nradar
    nobslotsradar(islot,iradar)=0
    WRITE(radarfile(4:5),'(I2.2)')islot
    WRITE(radarfile(6:7),'(I2.2)')iradar
    INQUIRE(FILE=radarfile,EXIST=ex)
     IF( ex )THEN
      !Get radar data size
      CALL get_nobs_radar( radarfile , nobslotsradar(islot,iradar),radar_lon(iradar),radar_lat(iradar),radar_z(iradar) )   
      nobsradar=nobsradar+nobslotsradar(islot,iradar)
     ELSE
      WRITE(6,*)'WARNING: FILE ',radarfile,' NOT FOUND'
     ENDIF
   ENDDO
  ENDDO

  WRITE(6,'(I10,A)') nobs,' TOTAL OBSERVATIONS INPUT'
  WRITE(6,'(I10,A)') nobsradar,' TOTAL RADAR OBSERVATIONS INPUT'

!
! INITIALIZE GLOBAL VARIABLES
!
  ALLOCATE( tmpelm(nobs+nobsradar) )
  ALLOCATE( tmplon(nobs+nobsradar) )
  ALLOCATE( tmplat(nobs+nobsradar) )
  ALLOCATE( tmplev(nobs+nobsradar) )
  ALLOCATE( tmpdat(nobs+nobsradar) )
  ALLOCATE( tmperr(nobs+nobsradar) )
  ALLOCATE( tmptyp(nobs+nobsradar) )
  ALLOCATE( tmplot(nobs+nobsradar) )
  ALLOCATE( tmpi(nobs+nobsradar) )
  ALLOCATE( tmpj(nobs+nobsradar) )
  ALLOCATE( tmpk(nobs+nobsradar) )
  ALLOCATE( tmpdep(nobs+nobsradar) )
  ALLOCATE( tmpsprd(nobs+nobsradar) )
  ALLOCATE( tmphdxf(nobs+nobsradar,nbv) )
  ALLOCATE( tmpradar(nobs+nobsradar) )
  ALLOCATE( tmpaz(nobs+nobsradar) )
  ALLOCATE( tmpra(nobs+nobsradar) )
  ALLOCATE( tmpel(nobs+nobsradar) )

  
  tmpelm=0.0
  tmplon=0.0
  tmplat=0.0
  tmplev=0.0
  tmpdat=0.0
  tmperr=0.0
  tmptyp=0.0
  tmplon=0.0
  tmpi=0.0
  tmpj=0.0
  tmpk=0.0
  tmpdep=0.0
  tmpsprd=0.0
  tmphdxf=0.0
  tmpradar=0.0


!
! LOOP of timeslots 
!
  nn=0

timeslots: DO islot=1,nslots


  nini=nn+1

  !READ CONVENTIONAL DATA.

  IF(nobslots(islot) .GT. 0) THEN 
    WRITE(obsfile(4:5),'(I2.2)') islot
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',obsfile
    CALL read_obs(obsfile,nobslots(islot),&
      & tmpelm(nn+1:nn+nobslots(islot)),tmplon(nn+1:nn+nobslots(islot)),&
      & tmplat(nn+1:nn+nobslots(islot)),tmplev(nn+1:nn+nobslots(islot)),&
      & tmpdat(nn+1:nn+nobslots(islot)),tmperr(nn+1:nn+nobslots(islot)),&
      & tmptyp(nn+1:nn+nobslots(islot)) )
    tmplot(nn+1:nn+nobslots(islot)) = REAL(islot,r_size)
    tmpradar(nn+1:nn+nobslots(islot)) = 0.0d0
  ENDIF
  nn=nn+nobslots(islot)

  !READ RADAR DATA.

     DO iradar=1,nradar
      IF( nobslotsradar(islot,iradar) .GT. 0 )THEN
      WRITE(radarfile(4:5),'(I2.2)')islot
      WRITE(radarfile(6:7),'(I2.2)')iradar
      INQUIRE(FILE=radarfile,EXIST=ex)
       IF( ex )THEN
          !Radar data comes georeferenced in azimuth,elevation and range.
          WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',radarfile
          CALL read_obs_radar(radarfile,nobslotsradar(islot,iradar),&
         & tmpelm(nn+1:nn+nobslotsradar(islot,iradar)),tmpaz(nn+1:nn+nobslotsradar(islot,iradar)),&
         & tmpel(nn+1:nn+nobslotsradar(islot,iradar)),tmpra(nn+1:nn+nobslotsradar(islot,iradar)),&
         & tmpdat(nn+1:nn+nobslotsradar(islot,iradar)),tmperr(nn+1:nn+nobslotsradar(islot,iradar)),&
         & tmptyp(nn+1:nn+nobslotsradar(islot,iradar)) )

         tmplot(nn+1:nn+nobslotsradar(islot,iradar))=REAL(islot,r_size)
         !Add a flag to identify the radar that produced the observation. This is used by the 
         !observation operator in the computation of radial velocity and to compute the lat, lon
         !and height from the azimuth, range and elevation angle data.
         tmpradar(nn+1:nn+nobslotsradar(islot,iradar))=iradar*1.0d0 !Store radar type
        ENDIF !IF over the existance of the file

      ENDIF

      nn=nn+nobslotsradar(islot,iradar)
     ENDDO

    nend=nn
    !nini and nend are the start and end index of all the observations in the current
    !time slot. The observation operator will loop from nini to nend  

    DO n=nini,nend 
      IF( tmpradar(n) .GT. 0.0d0) THEN 
          !This is a radar observation.
          current_radar=NINT(tmpradar(n))
          !Convert from azimuth,elevation,range to lat,lon,z
          CALL aer2llz(radar_lon(current_radar),radar_lat(current_radar),radar_z(current_radar), &
                       tmpaz(n),tmpel(n),tmpra(n),tmplon(n),tmplat(n),tmplev(n))
          !WRITE(6,*)tmpaz(n),tmpel(n),tmpra(n),tmpdat(n)
      ENDIF
      CALL latlon_to_ij(projection,tmplat(n),tmplon(n),tmpi(n),tmpj(n))
      !CALL ll2ij(tmpelm(n),tmplon(n),tmplat(n),tmpi(n),tmpj(n))
    END DO



    l=0
    DO
      im = myrank+1 + nprocs * l
      IF(im > nbv) EXIT
      IF( nini-1 == nend )THEN   !NO OBSERVATIONS FOR THIS TIME SLOT, SKIP OBS. OPERATOR.
           l=l+1
           CYCLE
      ENDIF

      IF(.NOT.ALLOCATED(v3d))ALLOCATE(v3d(nlon,nlat,nlev,nv3d))
      IF(.NOT.ALLOCATED(v2d))ALLOCATE(v2d(nlon,nlat,nv2d))

      WRITE(guesfile(3:9),'(I2.2,I5.5)') islot,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',guesfile
      CALL read_grd(guesfile,v3d,v2d)
      !Keep model height 
      zmodel=v3d(:,:,:,iv3d_ph)/gg

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,dz,tg,qg)
      DO n=nini,nend

        tmphdxf(n,im)=undef !Initialize observation array with undefined values.

        IF(CEILING(tmpi(n)) < 2 .OR. nlon-1 < CEILING(tmpi(n))) THEN
          CYCLE
        END IF
        IF(CEILING(tmpj(n)) < 2 .OR. nlat-1 < CEILING(tmpj(n))) THEN
          CYCLE
        END IF

        IF( tmpradar(n) .GT. 0.0d0 )THEN
         !This is a radar observation. Compute k from z.
         CALL z2k_fast(zmodel,tmpi(n),tmpj(n),tmplev(n),tmpk(n))
        ELSE
         !This is a conventional observation. Compute k from p.
         CALL p2k(v3d(:,:,:,iv3d_p),tmpelm(n),tmpi(n),tmpj(n),tmplev(n),tmpk(n))
        ENDIF

        IF(CEILING(tmpk(n)) > nlev-1) THEN
          CYCLE
        END IF
        IF(CEILING(tmpk(n)) < 2 .AND. NINT(tmpelm(n)) /= id_ps_obs) THEN
          IF(NINT(tmpelm(n)) == id_u_obs .OR.&
           & NINT(tmpelm(n)) == id_v_obs) THEN
            tmpk(n) = 1.00001d0
          ELSE IF(NINT(tmpelm(n)) > 9999) THEN
            tmpk(n) = 0.0d0         
          ELSE
            CYCLE
          END IF
        END IF
        IF(NINT(tmpelm(n)) == id_ps_obs .AND. tmpdat(n) < -100.0d0) CYCLE
        IF(NINT(tmpelm(n)) == id_ps_obs) THEN
          CALL itpl_2d(phi0,tmpi(n),tmpj(n),dz)
          tmpk(n) = tmplev(n) - dz
          IF(ABS(tmpk(n)) > threshold_dz) THEN ! pressure adjustment threshold
            CYCLE
          END IF
        END IF
        !
        ! observational operator
        !
        CALL Trans_XtoY(tmpelm(n),tmptyp(n),tmplon(n),tmplat(n),tmpi(n),tmpj(n),tmpk(n),&
         & tmpaz(n),tmpel(n),v3d,v2d,tmphdxf(n,im))

      END DO
!$OMP END PARALLEL DO

      l = l+1
    END DO

  END DO timeslots
 
  WRITE(6,*)"Finish forward observation operator"

  IF ( ALLOCATED(v3d) ) DEALLOCATE(v3d)
  IF ( ALLOCATED(v2d) ) DEALLOCATE(v2d)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


  CALL allreduce_obs_mpi(nobs+nobsradar,nbv,tmphdxf)


  !Select radar observations to be assimilted and compute pseudo relative humidity
IF ( nobsradar .GT. 0 ) THEN
  CALL process_radar_obs(nobs+nobsradar,nbv,tmphdxf,tmpdat,tmpelm,tmperr)
ENDIF

IF( nobs + nobsradar .GT. 0)THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i)
 DO n=1,nobs+nobsradar
    IF( ANY( tmphdxf(n,:) == undef ) )THEN
      tmpdep(n)=undef
    ENDIF
    tmpdep(n) = 0.0d0
    DO i=1,nbv
      tmpdep(n) = tmpdep(n) + tmphdxf(n,i)
    END DO
    tmpdep(n) = tmpdep(n) / REAL(nbv,r_size)
    DO i=1,nbv
      tmphdxf(n,i) = tmphdxf(n,i) - tmpdep(n) ! Hdx
    END DO
    tmpdep(n) = tmpdat(n) - tmpdep(n) ! y-Hx
    IF(tmpelm(n) == id_tcmip_obs) THEN
!!! (1) gross_error set for tcmip_obs
      IF(ABS(tmpdep(n)) > gross_error_tycmip) THEN !gross error  
        WRITE(6,'(a)') 'xxx TC PC INFO REJECTED IN QC xxx'
        WRITE(6,'(a)') '   dep =', ABS(tmpdep(n))  
        tmpdep(n)=undef
      END IF
    ELSE IF(tmpelm(n) == id_tclon_obs .or. tmpelm(n) == id_tclat_obs) THEN
      IF(ABS(tmpdep(n)) > gross_error_tycll) THEN
        WRITE(6,'(a)') 'xxx TC LL INFO REJECTED IN QC xxx'
        WRITE(6,'(a)') '   dep =', ABS(tmpdep(n))  
        tmpdep(n)=undef
      END IF
    ELSE   
      IF(ABS(tmpdep(n)) > gross_error*tmperr(n)) THEN !gross error
        tmpdep(n)=undef
      END IF
    END IF
! < --- QC for tyinfo
 END DO
!$OMP END PARALLEL DO

ENDIF

!
!PERFORM SUPEROBING ON REQUESTED OBSERVATION TYPES
!

   CALL superobbing()




do n=1,nobs+nobsradar
   CALL com_stdev(nbv,tmphdxf(n,:),tmpsprd(n) )
enddo

CALL monit_dep(nobs+nobsradar,tmpelm,tmpdep,tmpsprd)

  nn = 0
  DO n=1,nobs+nobsradar
    IF(tmpdep(n) == undef ) CYCLE
    nn = nn+1
    tmpelm(nn) = tmpelm(n)
    tmplon(nn) = tmplon(n)
    tmplat(nn) = tmplat(n)
    tmplev(nn) = tmplev(n)
    tmpdat(nn) = tmpdat(n)
    tmperr(nn) = tmperr(n)
    tmptyp(nn) = tmptyp(n)
    tmplot(nn) = tmplot(n)
    tmpi(nn) = tmpi(n)
    tmpj(nn) = tmpj(n)
    tmpdep(nn) = tmpdep(n)
    tmphdxf(nn,:) = tmphdxf(n,:)
    tmpaz(nn) = tmpaz(n)
    tmpel(nn) = tmpel(n)
    tmpra(nn) = tmpra(n)
    tmpradar(nn) = tmpradar(n)
  END DO
  nobs = nn !WARNING: now obs is the total number of observations. Including radar
            !observation 
  WRITE(6,'(I10,A,I3.3)') nobs,' OBSERVATIONS TO BE ASSIMILATED '
!
! SORT
!
  ALLOCATE( tmp2elm(nobs) )
  ALLOCATE( tmp2lon(nobs) )
  ALLOCATE( tmp2lat(nobs) )
  ALLOCATE( tmp2lev(nobs) )
  ALLOCATE( tmp2dat(nobs) )
  ALLOCATE( tmp2err(nobs) )
  ALLOCATE( tmp2typ(nobs) )
  ALLOCATE( tmp2lot(nobs) )
  ALLOCATE( tmp2i(nobs) )
  ALLOCATE( tmp2j(nobs) )
  ALLOCATE( tmp2az(nobs) )
  ALLOCATE( tmp2el(nobs) )
  ALLOCATE( tmp2ra(nobs) )
  ALLOCATE( tmp2radar(nobs) )
  ALLOCATE( tmp2dep(nobs) )
  ALLOCATE( tmp2hdxf(nobs,nbv) )
  ALLOCATE( obselm(nobs) )
  ALLOCATE( obslon(nobs) )
  ALLOCATE( obslat(nobs) )
  ALLOCATE( obslev(nobs) )
  ALLOCATE( obsdat(nobs) )
  ALLOCATE( obserr(nobs) )
  ALLOCATE( obstyp(nobs) )
  ALLOCATE( obslot(nobs) )
  ALLOCATE( obsi(nobs) )
  ALLOCATE( obsj(nobs) )
  ALLOCATE( obsaz(nobs) )
  ALLOCATE( obsra(nobs) )
  ALLOCATE( obsel(nobs) )
  ALLOCATE( obsradar(nobs) )
  ALLOCATE( obsdep(nobs) )
  ALLOCATE( obshdxf(nobs,nbv) )
  ALLOCATE( nobsgrd(nlon,nlat) )
 
  nobsgrd = 0
  nj = 0
!$OMP PARALLEL PRIVATE(i,j,n,nn)
!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    DO n=1,nobs
      IF(tmpj(n) < j .OR. j+1 <= tmpj(n)) CYCLE
      nj(j) = nj(j) + 1
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    njs(j) = SUM(nj(0:j-1))
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    nn = 0
    DO n=1,nobs
      IF(tmpj(n) < j .OR. j+1 <= tmpj(n)) CYCLE
      nn = nn + 1
      tmp2elm(njs(j)+nn) = tmpelm(n)
      tmp2lon(njs(j)+nn) = tmplon(n)
      tmp2lat(njs(j)+nn) = tmplat(n)
      tmp2lev(njs(j)+nn) = tmplev(n)
      tmp2dat(njs(j)+nn) = tmpdat(n)
      tmp2err(njs(j)+nn) = tmperr(n)
      tmp2typ(njs(j)+nn) = tmptyp(n)
      tmp2lot(njs(j)+nn) = tmplot(n)
      tmp2i(njs(j)+nn) = tmpi(n)
      tmp2j(njs(j)+nn) = tmpj(n)
      tmp2az(njs(j)+nn) = tmpaz(n)
      tmp2ra(njs(j)+nn) = tmpra(n)
      tmp2el(njs(j)+nn) = tmpel(n)
      tmp2radar(njs(j)+nn) = tmpradar(n)    
      tmp2dep(njs(j)+nn) = tmpdep(n)
      tmp2hdxf(njs(j)+nn,:) = tmphdxf(n,:)
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    IF(nj(j) == 0) THEN
      nobsgrd(:,j) = njs(j)
      CYCLE
    END IF
    nn = 0
    DO i=1,nlon
      DO n=njs(j)+1,njs(j)+nj(j)
        IF(tmp2i(n) < i .OR. i+1 <= tmp2i(n)) CYCLE
        nn = nn + 1
        obselm(njs(j)+nn) = tmp2elm(n)
        obslon(njs(j)+nn) = tmp2lon(n)
        obslat(njs(j)+nn) = tmp2lat(n)
        obslev(njs(j)+nn) = tmp2lev(n)
        obsdat(njs(j)+nn) = tmp2dat(n)
        obserr(njs(j)+nn) = tmp2err(n)
        obstyp(njs(j)+nn) = tmp2typ(n)
        obslot(njs(j)+nn) = tmp2lot(n)
        obsi(njs(j)+nn) = tmp2i(n)
        obsj(njs(j)+nn) = tmp2j(n)
        obsaz(njs(j)+nn) = tmp2az(n)
        obsra(njs(j)+nn) = tmp2ra(n)
        obsel(njs(j)+nn) = tmp2el(n)
        obsradar(njs(j)+nn) = tmp2radar(n)
        obsdep(njs(j)+nn) = tmp2dep(n)
        obshdxf(njs(j)+nn,:) = tmp2hdxf(n,:)
      END DO
      nobsgrd(i,j) = njs(j) + nn
    END DO
    IF(nn /= nj(j)) THEN
!$OMP CRITICAL
      !WRITE(6,'(A,2I)') 'OBS DATA SORT ERROR: ',nn,nj(j)
      WRITE(6,'(F6.2,A,F6.2)') j,'< J <',j+1
      WRITE(6,'(F6.2,A,F6.2)') MINVAL(tmp2j(njs(j)+1:njs(j)+nj(j))),'< OBSJ <',MAXVAL(tmp2j(njs(j)+1:njs(j)+nj(j)))
!$OMP END CRITICAL
    END IF
  END DO
!$OMP END DO
!$OMP END PARALLEL

  DEALLOCATE( tmp2elm )
  DEALLOCATE( tmp2lon )
  DEALLOCATE( tmp2lat )
  DEALLOCATE( tmp2lev )
  DEALLOCATE( tmp2dat )
  DEALLOCATE( tmp2err )
  DEALLOCATE( tmp2typ )
  DEALLOCATE( tmp2lot )
  DEALLOCATE( tmp2i )
  DEALLOCATE( tmp2j )
  DEALLOCATE( tmp2az , tmp2ra , tmp2el )
  DEALLOCATE( tmp2radar )
!  DEALLOCATE( tmp2k )
  DEALLOCATE( tmp2dep )
  DEALLOCATE( tmp2hdxf )
  DEALLOCATE( tmpelm )
  DEALLOCATE( tmplon )
  DEALLOCATE( tmplat )
  DEALLOCATE( tmplev )
  DEALLOCATE( tmpdat )
  DEALLOCATE( tmperr )
  DEALLOCATE( tmptyp )
  DEALLOCATE( tmplot )
  DEALLOCATE( tmpi )
  DEALLOCATE( tmpj )
  DEALLOCATE( tmpk )
  DEALLOCATE( tmpsprd )
  DEALLOCATE( tmpdep )
  DEALLOCATE( tmphdxf )
  DEALLOCATE( tmpaz , tmpra , tmpel )
  DEALLOCATE( tmpradar )


! WRITE the result of the observation operator (one file per ensemble member)

  call obsope_io(nbv,nobs,2,obselm,obslon,obslat,obslev,           &
                 obsdat,obserr,obstyp,obslot,obsi,obsj,obsdep,obshdxf )


  RETURN
END SUBROUTINE set_letkf_obs

!-----------------------------------------------------------------------
! Process radar obs
! This subroutine performs some additional computations based on input radar
! observations.
! 1) The pseudo rh observations are computed.
! 2) Data is converted into dBz.
! 3) Useless observations are discarded.
! 4) Attenuated observations (that cannot be used as pseudo rh obs) are discarded.
!-----------------------------------------------------------------------
SUBROUTINE process_radar_obs(nobs,nmember,hdxf,value,id,error)
INTEGER, INTENT(IN) :: nobs,nmember
REAL(r_size), INTENT(INOUT) :: hdxf(nobs,nmember),id(nobs),value(nobs),error(nobs)
INTEGER :: iobs , i
REAL(r_size) :: rainratio 
!We will check if we shall use reflectivity or pseudo rh observations.
INTEGER :: ir , irh , irf , ius , imiss , iatt , ilt , izeroref

!minref=10.d0**(minrefdbz/10.0d0)

irh=0
ir=0
ius=0
ilt=0
imiss=0
izeroref=0
iatt=0

  DO iobs=1,nobs

   IF( id(iobs) /= id_reflectivity_obs  )CYCLE
     !First check that there are no missing observations.
     IF( ANY( hdxf(iobs,:) == undef ) )THEN
       imiss=imiss+1
       CYCLE
     ENDIF

     !------COMPUTE REF in DBZ  -------------------
     
     IF( value(iobs) .GT. minref)THEN
        value(iobs) = 10*log10( value(iobs) )
     ELSE
        value(iobs) = minrefdbz
     ENDIF

     rainratio=0.0d0
     DO i=1,nmember
        IF( hdxf(iobs,i) .GT. minref)THEN
           hdxf(iobs,i)=10*log10(hdxf(iobs,i)) 
           rainratio=rainratio+1.0d0
        ELSE
           hdxf(iobs,i)=minrefdbz
        ENDIF
     ENDDO

     rainratio=rainratio/REAL(nmember,r_size)

     IF( rainratio .GE. rainratio_threshold )THEN
       IF( value(iobs) .GT. minrefdbz )THEN
         ir=ir+1
       ELSE
         izeroref=izeroref+1
       ENDIF
     ELSE
       IF( value(iobs) .GT. minrefdbz )THEN
         hdxf(iobs,:)=undef
         ilt=ilt+1
       ELSE
         IF( rainratio .EQ. 0.0d0 )THEN
           hdxf(iobs,:)=undef
           ius=ius+1
         ELSE
           izeroref=izeroref+1
         ENDIF
       ENDIF
     ENDIF


  ENDDO

WRITE(6,*)'========================================='
WRITE(6,*)'DETAILED ROUTINE OUTPUT                  '
WRITE(6,*)'REFLECTIVITY OBS    = ',ir
WRITE(6,*)'MISSING OBS         = ',imiss
WRITE(6,*)'USELESS OBS         = ',ius
WRITE(6,*)'LOW RAIN THRESHOLD  = ',ilt
WRITE(6,*)'ZERO REFLECTIVITY   = ',izeroref
WRITE(6,*)'TOTAL OBS           = ',nobs
WRITE(6,*)'========================================='


END SUBROUTINE process_radar_obs

!Compute regrid grid.
!This grid will be used to generate a grided output
!for the errors computed from the observations.
!------------------------------------------------------
SUBROUTINE get_regrid_grid()
IMPLICIT NONE
INTEGER      :: i,j,ivar,ilev,ilevreg,ilevprev,iregprev
INTEGER      :: counter2d , counter3d , iv , k
REAL(r_size) :: a
REAL(r_size) :: aux
LOGICAL      :: testvar

  nv3dver=0
  nv2dver=0
  !Get the variables that will be verified and whether they are 3D or 2D.
  do ivar=1,maxverifvars
     if( variable_list(ivar) > 0 )then
       call is_3dvar(variable_list(ivar),testvar)
       if ( testvar ) then
          nv3dver=nv3dver+1
       else
          nv2dver=nv2dver+1
       endif
     endif
  enddo


  nvver=nv3dver+nv2dver

  allocate( varidver(nv3dver+nv2dver) )

  !Generate a new variable list with the 3d variables first and then
  !2d variables.
  counter3d=1
  counter2d=1
  do ivar=1,maxverifvars
     if( variable_list(ivar) > 0 )then
       call is_3dvar(variable_list(ivar),testvar)
       if ( testvar ) then
          varidver(counter3d)=variable_list(ivar)
          counter3d=counter3d+1
       else
          varidver(counter2d+nv3dver)=variable_list(ivar)
          counter2d=counter2d+1 
       endif
     endif
  enddo

  minlon=1d30
  maxlon=-1d30
  minlat=1d30
  maxlat=-1d30

 
  DO i=1,nlon-1
   DO j=1,nlat-1
     if( lat(i,j) > maxlat )maxlat=lat(i,j)
     if( lat(i,j) < minlat )minlat=lat(i,j)
     if( lon(i,j) > maxlon )maxlon=lon(i,j)
     if( lon(i,j) < minlon )minlon=lon(i,j)
   ENDDO
  ENDDO

  nlonreg=CEILING( (maxlon-minlon)/regrid_res ) + 1
  nlatreg=CEILING( (maxlat-minlat)/regrid_res ) + 1

  ALLOCATE( lonreg(nlonreg,nlatreg) , latreg(nlonreg,nlatreg) )

  DO i = 1,nlonreg
   DO j = 1,nlatreg

    lonreg(i,j)=minlonreg + REAL((i-1),r_size)*regrid_res
    latreg(i,j)=minlatreg + REAL((j-1),r_size)*regrid_res

   ENDDO
  ENDDO

  minlonreg=minlon
  minlatreg=minlat
  maxlonreg=maxlon
  maxlatreg=maxlat
  maxlonreg=minlon+nlonreg*regrid_res
  maxlatreg=minlat+nlatreg*regrid_res


  maxlevreg=100000.0d0 !Force the lowest grid level to be 1000hPa.
  !minlev=P_TOP    !Model top pressure from input file.

  !nlevreg= CEILING( (maxlev-minlev)/ regrid_vert_res ) + 1

  minlevreg= maxlevreg - REAL(nlevreg,r_size)*regrid_vert_res

  minzlevreg=0.0d0
  maxzlevreg=minzlevreg + REAL(nlevreg,r_size)*regrid_vert_zres

   WRITE(6,*)'=================================================='
   WRITE(6,*)'   Oringial grid limits' 
   WRITE(6,*)'=================================================='
   WRITE(6,*)' MINLON = ',minlon
   WRITE(6,*)' MAXLON = ',maxlon
   WRITE(6,*)' MINLAT = ',minlat
   WRITE(6,*)' MAXLAT = ',maxlat
   WRITE(6,*)'=================================================='
   WRITE(6,*)'    Obsgrid output enambled'
   WRITE(6,*)'=================================================='
   WRITE(6,*)'Regrid output grid                                '
   WRITE(6,*)'INILON= ',minlonreg
   WRITE(6,*)'INILAT= ',minlatreg
   WRITE(6,*)'NLON  = ',nlonreg
   WRITE(6,*)'NLAT  = ',nlatreg
   WRITE(6,*)'HRES  = ',regrid_res,' (degree)'
   WRITE(6,*)'MAXLEV= ',maxlevreg
   WRITE(6,*)'MINLEV= ',minlevreg
   WRITE(6,*)'NLEV  = ',nlevreg
   WRITE(6,*)'LEVRES= ',regrid_vert_res, ' (hPa)'
   WRITE(6,*)'MAXLEVZ= ',maxzlevreg
   WRITE(6,*)'MINLEVZ= ',minzlevreg
   WRITE(6,*)'NLEV  = ',nlevreg
   WRITE(6,*)'LEVRES= ',regrid_vert_zres, ' (m)'
   WRITE(6,*)'3D VARIABLE LIST'
   do ivar=1,nv3dver
   WRITE(6,*)'Var ID ',varidver(ivar)
   enddo
   WRITE(6,*)'2D VARIABLE LIST'
   do ivar=1,nv2dver
   WRITE(6,*)'Var ID ',varidver(nv3dver+ivar)
   enddo
   WRITE(6,*)'=================================================='

  ALLOCATE( levreg(nlevreg) , levzreg(nlevreg) )
  DO i = 1,nlevreg
    levreg(i)=maxlevreg - REAL((i-1),r_size)*regrid_vert_res
  ENDDO
  DO i = 1,nlevreg
    levzreg(i)=minzlevreg + REAL((i-1),r_size)*regrid_vert_zres
  ENDDO
  
  nvreg=nv3dver * nlevreg + nv2dver


  ALLOCATE( levall(nvreg) , varidall(nvreg) )

  DO iv = 1 , nvver

     if( iv > nv3dver )then
       !2D variable
       k = nv3dver * nlevreg + ( iv - nv3dver ) 
       levall( k )=maxlevreg
       varidall( k )=varidver(iv)
     else
      !3D variable
      CALL is_vertp( varidver(iv) , testvar )
      DO ilev=1,nlevreg      
       k= (iv - 1)*nlevreg + ilev
       varidall( k )=varidver(iv)
       
       if( testvar) then
         !P is the vertical coordinate.
         levall(k)=maxlevreg - (ilev - 1)*regrid_vert_res
       else
         !Z is the vertical coordinate.
         levall(k)=minlevreg + (ilev - 1)*regrid_vert_zres
       endif
      ENDDO
     end if

  ENDDO
  WRITE(6,*)' Variable and levels list '
  DO iv= 1 , nvreg
    WRITE(6,*)varidall(iv),levall(iv)
  ENDDO
  WRITE(6,*)'=================================================='


END SUBROUTINE get_regrid_grid

!Regrid var
!This subroutine regrids obs data into a regular grid.
!TODO: This routine might be usefull for data thining.
!-----------------------------------------------------------
SUBROUTINE obs_to_grid(inislot,endslot,obsvar,grdobs,grdnum)
IMPLICIT NONE
INTEGER     , INTENT(IN) :: inislot , endslot
REAL(r_size), INTENT(IN) :: obsvar(nobs)
REAL(r_size), INTENT(OUT):: grdobs(nlonreg,nlatreg,nvreg)
INTEGER                  :: grdnum(nlonreg,nlatreg,nvreg)
INTEGER                  :: i , j , k ,  iobs , iv , jreg , ireg
LOGICAL                  :: testvar

   grdobs=0.0d0
   grdnum=0

   DO i = 1,nobs
      if ( obslot(i) < inislot .OR. obslot(i) > endslot )cycle
      if ( obsvar(i) == undef )cycle

       DO iv=1,nvver
         if( varidver(iv) == obselm(i) ) then
 
           ireg= NINT( ( obslon(i) - minlonreg ) / regrid_res  ) + 1
           jreg= NINT( ( obslat(i) - minlatreg ) / regrid_res  ) + 1
 
           if( iv > nv3dver )then
             !2D variable
             k = nv3dver * nlevreg + iv - nv3dver 
           else
             !3D variable
             CALL is_vertp( varidver(iv) , testvar )
             if( testvar )then  
              !P is the vertical coordinate.
              k= NINT( ( maxlevreg- obslev(i) ) / regrid_vert_res ) + 1
              k=k+nlevreg*(iv-1)
             else
              !Z is the vertical coordinate.
              k= NINT( ( obslev(i) - minzlevreg ) / regrid_vert_zres ) + 1
              k=k+nlevreg*(iv-1)
             endif
           endif
           grdobs(ireg,jreg,k)=grdobs(ireg,jreg,k) + obsvar(i)
           grdnum(ireg,jreg,k)=grdnum(ireg,jreg,k) + 1
           exit
         endif
       ENDDO
   ENDDO

   WHERE( grdnum > 0 )
     grdobs=grdobs/REAL(grdnum,r_size)
   ELSEWHERE
     grdobs=undef
   ENDWHERE

   

END SUBROUTINE obs_to_grid

!-- Write a grid file for obsgrid output
SUBROUTINE write_obsgrid(filename,nx,ny,nz,var)
  IMPLICIT NONE
  INTEGER, INTENT(IN)     :: nx,ny,nz
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: var(nx,ny,nz)
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
   WHERE( var(:,:,n) == undef )
    buf4=undefs
   ENDWHERE
   WRITE(iunit,rec=n) buf4
  ENDDO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obsgrid

!-- Read a grid file ---------------------------------------------------
SUBROUTINE read_obsgrid(filename,nx,ny,nz,var)
  IMPLICIT NONE
  INTEGER,INTENT(IN)      :: nx,ny,nz
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: var(nx,ny,nz)
  REAL(r_sngl) :: buf4(nx,ny)
  INTEGER :: iunit,iolen
  INTEGER :: k,n,reclength,ios

  INQUIRE(IOLENGTH=reclength) reclength
  reclength = reclength * nx * ny


 iunit=33
 INQUIRE(FILE=filename,EXIST=file_exist)
  IF(file_exist) THEN
   OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=reclength)
  ELSE
   var=undef
   RETURN
  ENDIF

  DO n=1,nz
    READ(iunit,rec=n,IOSTAT=ios)buf4
     IF( ios == 0 )THEN
       var(:,:,n) = REAL(buf4,r_size)
       WHERE( buf4 == undefs )
        var(:,:,n)=undef
       ENDWHERE

     ELSE !We reached the end of the file.
       IF( n==1)THEN
         WRITE(6,*)"[Warning]: Error reading, this might be a new file "
         WRITE(6,*)filename
         var=0.0d0
       ELSE
         WRITE(6,*)"[Error]: Unexpectively reached end of file"
         WRITE(6,*)filename
         STOP
       ENDIF
     ENDIF
  END DO

  

  CLOSE(iunit)

  RETURN

END SUBROUTINE read_obsgrid


!Obsgrid is the main subroutine that inputs Yb matrix and the 
!obs metadata and generates the grided output and statistics.
!It also generates the text output.
SUBROUTINE obsgrid()
IMPLICIT NONE
!REAL(r_size)            :: obsmean(nobs) , obssprd(nobs) !Ensemble mean and spread in Ospace
REAL(r_size)            :: mean1(nlonreg,nlatreg,nvreg)   !Ensemble mean. 
REAL(r_size)            :: sprd1(nlonreg,nlatreg,nvreg)  !Ensemble spread.
REAL(r_size)            :: erro1(nlonreg,nlatreg,nvreg)  !Mean error.
INTEGER                 :: onum1(nlonreg,nlatreg,nvreg)
INTEGER                 :: islot , inislot , endslot



ALLOCATE( obsmean(nobs) , obssprd(nobs) )
!Compute the ensemble mean and spread.
obsmean=obsdep+obsdat
obssprd=SQRT( SUM(obshdxf**2,2)/REAL(nbv,r_size) )

WHERE( obsdep == undef )
  obsmean=undef
  obssprd=undef
ENDWHERE

inislot=slotoffset

do
  if( inislot > nslots )exit
  endslot=inislot+slotstep

  !Generate filenames.
  WRITE(meanfile(8:12),'(I2.2,A1,I2.2)')inislot,'_',endslot
  WRITE(sprdfile(8:12),'(I2.2,A1,I2.2)')inislot,'_',endslot
  WRITE(merrfile(8:12),'(I2.2,A1,I2.2)')inislot,'_',endslot
  WRITE(onumfile(8:12),'(I2.2,A1,I2.2)')inislot,'_',endslot


  CALL obs_to_grid(inislot,endslot,obsmean,mean1,onum1)
  WRITE(6,*)"I will write a file ",meanfile
  CALL write_obsgrid(meanfile,nlonreg,nlatreg,nvreg,mean1)

  !Grid and write ensemble spread
  CALL obs_to_grid(inislot,endslot,obssprd,sprd1,onum1)
  WRITE(6,*)"I will write a file ",sprdfile
  CALL write_obsgrid(sprdfile,nlonreg,nlatreg,nvreg,sprd1)

  !Grid and write ensemble mean error
  CALL obs_to_grid(inislot,endslot,obsdep,erro1,onum1)
  WRITE(6,*)"I will write a file ",merrfile
  CALL write_obsgrid(merrfile,nlonreg,nlatreg,nvreg,erro1)
  WRITE(6,*)"I will write a file ",onumfile
  CALL write_obsgrid(onumfile,nlonreg,nlatreg,nvreg,REAL(onum1,r_size))

  !Generate text output
  WRITE(areafile(6:10),'(I2.2,A1,I2.2)')inislot,'_',endslot
  WRITE(6,*)"I will write a file ",areafile
  CALL obsarea(inislot,endslot,areafile)

  inislot=inislot+slotstep+1


end do ![Do over time slots]

DEALLOCATE( obsmean , obssprd )


!To compute text output.

END SUBROUTINE obsgrid

!Obsarea is the main subroutine that inputs Yb matrix and the
!obs metadata and generates a text ouput containing the error statistics
!over selected areas.
!It also generates the text output.

SUBROUTINE obsarea(inislot,endslot,filename)
IMPLICIT NONE
integer , intent(in)     :: inislot,endslot
integer                  :: num_ana(narea,nvreg)  !narea,nv
real(r_size)             :: wei_ana(narea,nvreg)  !narea,nv
real(r_size)             :: bias_ana(narea,nvreg) !narea,nv
real(r_size)             :: abse_ana(narea,nvreg) !narea,nv
real(r_size)             :: rmse_ana(narea,nvreg) !narea,nv
real(r_size)             :: sprd_ana(narea,nvreg) !narea,nv
integer                  :: iarea , i , iv , k , klev , iunit , maxk
character(*), INTENT(IN) :: filename
logical                  :: testvar

iunit=101
  !Compute rmse over selected regions and generate text ouput.

  num_ana = 0
  wei_ana = 0.0d0
  bias_ana = 0.0d0
  abse_ana = 0.0d0
  rmse_ana = 0.0d0
  sprd_ana = 0.0d0

Do iarea = 1,narea

  DO i = 1,nobs
      if( obslot(i) < inislot .OR. obslot(i) > endslot )cycle
      if( obsdep(i) == undef )cycle
      if( obslon(i) < vlon1(iarea) )cycle
      if( obslon(i) > vlon2(iarea) )cycle
      if( obslat(i) < vlat1(iarea) )cycle
      if( obslat(i) > vlat2(iarea) )cycle


       DO iv=1,nvver
         if( varidver(iv) == obselm(i) ) then

           if( iv > nv3dver )then
             !2D variable
             k = nv3dver * nlevreg + iv - nv3dver
           else
             !3D variable
             CALL is_vertp( varidver(iv) , testvar )
             if( testvar) then
               !P is the vertical coordinate.
               k= NINT( ( maxlevreg-obslev(i) ) / regrid_vert_res ) + 1
               k=k+nlevreg*(iv-1)
             else
               !Z is the vertical coordinate.
               k= NINT( ( obslev(i) - minzlevreg ) / regrid_vert_zres ) + 1
               k=k+nlevreg*(iv-1)
             end if
           endif
              num_ana(iarea,k) = num_ana(iarea,k) + 1
              bias_ana(iarea,k) = bias_ana(iarea,k) + obsdep(i)
              abse_ana(iarea,k) = abse_ana(iarea,k) + abs(obsdep(i))
              rmse_ana(iarea,k) = rmse_ana(iarea,k) + obsdep(i)  ** 2
              sprd_ana(iarea,k) = sprd_ana(iarea,k) + obssprd(i) ** 2
           exit
         endif
       ENDDO
   ENDDO

ENDDO



do k = 1, nvreg
   do iarea = 1, narea
      if (num_ana(iarea,k) == 0) then
         bias_ana(iarea,k)=undef
         sprd_ana(iarea,k)=undef
         rmse_ana(iarea,k)=undef
         abse_ana(iarea,k)=undef
      else
         bias_ana(iarea,k)=bias_ana(iarea,k)/real(num_ana(iarea,k),r_size)
         rmse_ana(iarea,k)=SQRT(rmse_ana(iarea,k)/real(num_ana(iarea,k),r_size) )
         abse_ana(iarea,k)=abse_ana(iarea,k)/real(num_ana(iarea,k),r_size)
         sprd_ana(iarea,k)=SQRT(sprd_ana(iarea,k)/real(num_ana(iarea,k),r_size) )
      end if
   end do
end do

    open(iunit,file=filename,form='formatted',recl=1000)
      WRITE(iunit,*)'AREAL ERROR STATISTICS COMPUTED FROM THE ORIGINAL GRID'
      WRITE(iunit,*)'UNDEF ',undef,' NVARS ',nvreg,' N VERIF. VARS ',nvver,' NLEVS ',nlevreg,' NAREA ',narea
      WRITE(iunit,*)'VLON1',vlon1(1:narea)
      WRITE(iunit,*)'VLON2',vlon2(1:narea)
      WRITE(iunit,*)'VLAT1',vlat1(1:narea)
      WRITE(iunit,*)'VLAT2',vlat2(1:narea)

      WRITE(iunit,*)'RMSE SECTION'
      WRITE(iunit,*)'VAR','  LEV  ',(' AREA',i,i=1,narea)
 
     DO iv=1,nvreg
       if( testvar)then
          !P is the vertical coordinate.
          WRITE(iunit,*)varidall(iv),levall(iv),(REAL(rmse_ana(iarea,iv),r_sngl),iarea=1,narea)
       else
          !Z is the vertical coordinate.
          WRITE(iunit,*)varidall(iv),levall(iv),(REAL(rmse_ana(iarea,iv),r_sngl),iarea=1,narea)
       endif
     ENDDO

      WRITE(iunit,*)'SPRD SECTION'
      WRITE(iunit,*)'VAR','  LEV  ',(' AREA',i,i=1,narea)

     DO iv=1,nvreg
       if( testvar)then
          !P is the vertical coordinate.
          WRITE(iunit,*)varidall(iv),levall(iv),(REAL(sprd_ana(iarea,iv),r_sngl),iarea=1,narea)
       else
          !Z is the vertical coordinate.
          WRITE(iunit,*)varidall(iv),levall(iv),(REAL(sprd_ana(iarea,iv),r_sngl),iarea=1,narea)
       endif
     ENDDO

      WRITE(iunit,*)'BIAS SECTION'
      WRITE(iunit,*)'VAR','  LEV  ',(' AREA',i,i=1,narea)

     DO iv=1,nvreg
       if( testvar)then
          !P is the vertical coordinate.
          WRITE(iunit,*)varidall(iv),levall(iv),(REAL(bias_ana(iarea,iv),r_sngl),iarea=1,narea)
       else
          !Z is the vertical coordinate.
          WRITE(iunit,*)varidall(iv),levall(iv),(REAL(bias_ana(iarea,iv),r_sngl),iarea=1,narea)
       endif
     ENDDO

      WRITE(iunit,*)'ABSE SECTION'
      WRITE(iunit,*)'VAR','  LEV  ',(' AREA',i,i=1,narea)

     DO iv=1,nvreg
       if( testvar)then
          !P is the vertical coordinate.
          WRITE(iunit,*)varidall(iv),levall(iv),(REAL(abse_ana(iarea,iv),r_sngl),iarea=1,narea)
       else
          !Z is the vertical coordinate.
          WRITE(iunit,*)varidall(iv),levall(iv),(REAL(abse_ana(iarea,iv),r_sngl),iarea=1,narea)
       endif
     ENDDO

      WRITE(iunit,*)'OBSNUMBER SECTION'
      WRITE(iunit,*)'VAR','  LEV  ',(' AREA',i,i=1,narea)

     DO iv=1,nvreg
       if( testvar)then
          !P is the vertical coordinate.
          WRITE(iunit,*)varidall(iv),levall(iv),(REAL(num_ana(iarea,iv),r_sngl),iarea=1,narea)
       else
          !Z is the vertical coordinate.
          WRITE(iunit,*)varidall(iv),levall(iv),(REAL(num_ana(iarea,iv),r_sngl),iarea=1,narea)
       endif
     ENDDO


    close(iunit)


END SUBROUTINE OBSAREA

!Perform data superobbing for selected observation types.
SUBROUTINE superobbing
IMPLICIT NONE






END SUBROUTINE SUPEROBING


END MODULE obsop_tools
