MODULE letkf_obs
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
  USE common_letkf
  USE common_mpi
  USE common_wrf
  USE common_obs_wrf
  USE common_mpi_wrf
  USE map_utils
  USE common_namelist

  IMPLICIT NONE
  PUBLIC

  INTEGER,SAVE :: nobs , nobsradar
  REAL(r_size),SAVE :: dist_zero
  REAL(r_size),SAVE :: dist_zerov
  REAL(r_size),SAVE :: dist_zeroz
  REAL(r_size),SAVE :: dist_zeroij
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
  REAL(r_size),ALLOCATABLE,SAVE :: obsaz(:),obsra(:),obsel(:)
  INTEGER,ALLOCATABLE,SAVE :: nobsgrd(:,:)


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
  REAL(r_size),ALLOCATABLE :: tmpdep(:)
  REAL(r_size),ALLOCATABLE :: tmpsprd(:)
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

  dist_zero = sigma_obs * SQRT(10.0d0/3.0d0) * 2.0d0
  dist_zerov = sigma_obsv * SQRT(10.0d0/3.0d0) * 2.0d0
  dist_zeroz= sigma_obsz * SQRT(10.0d0/3.0d0) * 2.0d0
  dist_zeroij = dist_zero / dx

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
  ALLOCATE( tmpsprd(nobs+nobsradar) )
  ALLOCATE( tmpdep(nobs+nobsradar) )
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


   !FROM LAT LON TO I J FOR ALL THE OBSERVATIONS.
    DO n=nini,nend 
      IF( tmpradar(n) .GT. 0.0d0) THEN 
          !This is a radar observation.
          current_radar=NINT(tmpradar(n))
          !Convert from azimuth,elevation,range to lat,lon,z
          CALL aer2llz(radar_lon(current_radar),radar_lat(current_radar),radar_z(current_radar),tmpaz(n), &
                       tmpel(n),tmpra(n),tmplon(n),tmplat(n),tmplev(n))
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

      !Check if the observation is within the horizontal domain.
      IF(CEILING(tmpi(n)) < 2 .OR. nlon-1 < CEILING(tmpi(n))) THEN
          WRITE(6,*)'X coordinate out of range'
          WRITE(6,*)tmpelm(n),tmplon(n),tmplat(n),tmplev(n),tmpi(n),tmpj(n),tmpk(n),tmpdat(n),tmperr(n)
        CYCLE
      END IF
      IF(CEILING(tmpj(n)) < 2 .OR. nlat-1 < CEILING(tmpj(n))) THEN
        WRITE(6,*)'Y coordinate out of range'
        WRITE(6,*)tmpelm(n),tmplon(n),tmplat(n),tmplev(n),tmpi(n),tmpj(n),tmpk(n),tmpdat(n),tmperr(n)
        CYCLE
      END IF

      !If it is within the horizontal domain compute the observation level.
      IF( tmpradar(n) .GT. 0.0d0 )THEN
        !This is a radar observation. Compute k from z.
        CALL z2k_fast(zmodel,tmpi(n),tmpj(n),tmplev(n),tmpk(n))
      ELSE
        !This is a conventional observation. Compute k from p.
        CALL p2k(v3d(:,:,:,iv3d_p),tmpelm(n),tmpi(n),tmpj(n),tmplev(n),tmpk(n))
      ENDIF

      !Check if the observation is within the vertical domain.
      IF(CEILING(tmpk(n)) > nlev-1) THEN
        WRITE(6,*)'Z coordinate out of range'
        WRITE(6,*)tmpelm(n),tmplon(n),tmplat(n),tmplev(n),tmpi(n),tmpj(n),tmpk(n),tmpdat(n),tmperr(n)
        CYCLE
      END IF
        IF(CEILING(tmpk(n)) < 2 .AND. NINT(tmpelm(n)) /= id_ps_obs) THEN
          IF(NINT(tmpelm(n)) == id_u_obs .OR.&
           & NINT(tmpelm(n)) == id_v_obs) THEN
            tmpk(n) = 1.00001d0
          ELSE IF(NINT(tmpelm(n)) > 9999) THEN
            tmpk(n) = 0.0d0         
	  ELSE
           WRITE(6,*)'Z below topography'
           WRITE(6,*)tmpelm(n),tmplon(n),tmplat(n),tmplev(n),tmpi(n),tmpj(n),tmpk(n),tmpdat(n),tmperr(n)
           CYCLE
          END IF
        END IF
        IF(NINT(tmpelm(n)) == id_ps_obs .AND. tmpdat(n) < -100.0d0) CYCLE
        IF(NINT(tmpelm(n)) == id_ps_obs) THEN
          CALL itpl_2d(phi0,tmpi(n),tmpj(n),dz)
          tmpk(n) = tmplev(n) - dz
          IF(ABS(tmpk(n)) > threshold_dz) THEN ! pressure adjustment threshold
            WRITE(6,*)'Pressure adjustment fails'
            WRITE(6,*)tmpelm(n),tmplon(n),tmplat(n),tmpi(n),tmpj(n),tmpk(n),tmpdat(n),tmperr(n)
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
IF( nobsradar .GT. 0 )THEN
  CALL process_radar_obs(nobs+nobsradar,nbv,tmphdxf,tmpdat,tmpelm,tmperr)
ENDIF

IF( nobs + nobsradar .GT. 0)THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i)
 DO n=1,nobs+nobsradar


    IF( ANY( tmphdxf(n,:) == undef ) )THEN
      tmpdep(n)=undef
      CYCLE
    ELSE 
      tmpdep(n) = 0.0d0
      DO i=1,nbv
        tmpdep(n) = tmpdep(n) + tmphdxf(n,i)
      END DO
      tmpdep(n) = tmpdep(n) / REAL(nbv,r_size)
      DO i=1,nbv
        tmphdxf(n,i) = tmphdxf(n,i) - tmpdep(n) ! Hdx
      END DO
      tmpdep(n) = tmpdat(n) - tmpdep(n) ! y-Hx
    ENDIF

    IF(tmpelm(n) == id_tcmip_obs) THEN
!!! (1) gross_error set for tcmip_obs
      IF(ABS(tmpdep(n)) > gross_error_tycmip) THEN !gross error  
        tmpdep(n)=undef
      END IF
    ELSE IF(tmpelm(n) == id_tclon_obs .or. tmpelm(n) == id_tclat_obs) THEN
      IF(ABS(tmpdep(n)) > gross_error_tycll) THEN
        tmpdep(n)=undef
      END IF
    ELSE IF(tmpelm(n) == id_reflectivity_obs )THEN
      IF(ABS(tmpdep(n)) > gross_error*tmperr(n) .AND. .NOT. force_norain_assimilation )THEN
        tmpdep(n)=undef
      ENDIF
    ELSE   
      IF(ABS(tmpdep(n)) > gross_error*tmperr(n)) THEN !gross error
        tmpdep(n)=undef
      END IF
    END IF
! < --- QC for tyinfo
 END DO
!$OMP END PARALLEL DO

ENDIF

do n=1,nobs+nobsradar
  CALL com_stdev(nbv,tmphdxf(n,:),tmpsprd(n))
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
            !observation but without rejected observations.
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
INTEGER :: ir , irh , irf , ius , imiss , iatt , ilt , izeroref , iref

!minref=10.d0**(minrefdbz/10.0d0)

irh=0
ir=0
ius=0
ilt=0
imiss=0
izeroref=0
iatt=0
iref=0

  DO iobs=1,nobs

   IF( id(iobs) /= id_reflectivity_obs  )CYCLE
     iref=iref+1
     !First check that there are no missing observations.
     IF( ANY( hdxf(iobs,:) == undef ) )THEN
       hdxf(iobs,:)=undef 
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

     IF( rainratio .GE. rainratio_threshold .AND. value(iobs) .GT. minrefdbz )THEN
         ir=ir+1
     ENDIF
     IF( rainratio .LT. rainratio_threshold .AND. value(iobs) .GT. minrefdbz )THEN
         ilt=ilt+1
         hdxf(iobs,:)=undef  !Observation rejected due to low rain threshold.
     ENDIF
     IF( value(iobs) == minrefdbz .AND. rainratio == 0.0d0 )THEN
         ius=ius+1           !Observation rejected because it is useless
         hdxf(iobs,:)=undef
     ENDIF
     IF( value(iobs) == minrefdbz .AND. rainratio .GT. 0.0d0 )THEN
         izeroref = izeroref + 1
     ENDIF


  ENDDO

WRITE(6,*)'========================================='
WRITE(6,*)'DETAILED ROUTINE OUTPUT                  '
WRITE(6,*)'REFLECTIVITY OBS    = ',ir
WRITE(6,*)'MISSING OBS         = ',imiss
WRITE(6,*)'USELESS OBS         = ',ius
WRITE(6,*)'LOW RAIN THRESHOLD  = ',ilt
WRITE(6,*)'ZERO REFLECTIVITY   = ',izeroref
WRITE(6,*)'TOTAL OBS           = ',iref
WRITE(6,*)'========================================='


END SUBROUTINE process_radar_obs
!-----------------------------------------------------------------------
! Monitor departure from gues/anal mean
!-----------------------------------------------------------------------
SUBROUTINE monit_mean(file,depout)
  IMPLICIT NONE
  CHARACTER(4),INTENT(IN) :: file
  REAL(r_size),INTENT(OUT) :: depout(nobs)
  REAL(r_size) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size) :: v2d(nlon,nlat,nv2d)
  REAL(r_size) :: z3d(nlon,nlat,nlev)
  REAL(r_size) :: elem
  REAL(r_size) :: bias_u,bias_v,bias_t,bias_ps,bias_q,bias_rh,bias_ref,bias_vr
  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_ps,rmse_q,rmse_rh,rmse_ref,rmse_vr
  REAL(r_size) :: hdxf,dep,rk
  INTEGER :: n,iu,iv,it,iq,ips,irh,iref,ivr
  CHARACTER(9) :: filename='filexxxxx'
  

  rmse_u  = 0.0d0
  rmse_v  = 0.0d0
  rmse_t  = 0.0d0
  rmse_q  = 0.0d0
  rmse_ps = 0.0d0
  rmse_rh = 0.0d0
  rmse_ref=0.0d0
  rmse_vr =0.0d0
  bias_u = 0.0d0
  bias_v = 0.0d0
  bias_t = 0.0d0
  bias_q = 0.0d0
  bias_ps = 0.0d0
  bias_rh = 0.0d0
  bias_ref = 0.0d0
  bias_vr  = 0.0d0
  iu  = 0
  iv  = 0
  it  = 0
  iq  = 0
  ips = 0
  irh = 0
  iref= 0
  ivr = 0

  WRITE(filename(1:9),'(A4,I5.5)') file,nbv+1
  WRITE(6,*)"I will read ", filename
  CALL read_grd(filename,v3d,v2d)
  z3d=v3d(:,:,:,iv3d_ph)/gg

  DO n=1,nobs

    IF( obsradar(n) .GT. 0.0d0) THEN
      !This is a radar observation
      CALL z2k_fast(z3d,obsi(n),obsj(n),obslev(n),rk)
    ELSE
       CALL p2k(v3d(:,:,:,iv3d_p),obselm(n),obsi(n),obsj(n),obslev(n),rk)
    ENDIF

    IF(CEILING(rk) > nlev - 1 ) CYCLE
    IF(obslot(n) .NE. nbslot)CYCLE     !Only observations corresponding to the analysis time will be usred for background mean verification.
    IF(CEILING(rk) < 2 .AND. NINT(obselm(n)) /= id_ps_obs) THEN
      IF(NINT(obselm(n)) == id_u_obs .OR. NINT(obselm(n)) == id_v_obs) THEN
        rk = 1.00001d0
      ELSE
        CYCLE
      END IF
    END IF
    IF(NINT(obselm(n)) == id_ps_obs) THEN
      CALL itpl_2d(phi0,obsi(n),obsj(n),rk)
      rk = obslev(n) - rk
    END IF
    CALL Trans_XtoY(obselm(n),obstyp(n),obslon(n),obslat(n),obsi(n),obsj(n),rk,obsaz(n),obsel(n),v3d,v2d,hdxf)

!Transform reflectivity to dBz scale.
if( NINT(obselm(n)) == id_reflectivity_obs )then
   !obs(n)=10*log10(obs(n))
   hdxf=10*log10(hdxf)
   if( obsdat(n) <= minrefdbz )then 
      obsdat(n)=minrefdbz 
   endif
   if( hdxf <= minrefdbz )then
       hdxf=minrefdbz
   endif
endif

    dep = obsdat(n) - hdxf
    depout(n) = dep
    SELECT CASE(NINT(obselm(n)))
    CASE(id_u_obs)
      rmse_u = rmse_u + dep**2
      bias_u = bias_u + dep
      iu = iu + 1
    CASE(id_v_obs)
      rmse_v = rmse_v + dep**2
      bias_v = bias_v + dep
      iv = iv + 1
    CASE(id_t_obs)
      rmse_t = rmse_t + dep**2
      bias_t = bias_t + dep
      it = it + 1
    CASE(id_q_obs)
      rmse_q = rmse_q + dep**2
      bias_q = bias_q + dep
      iq = iq + 1
    CASE(id_ps_obs)
      rmse_ps = rmse_ps + dep**2
      bias_ps = bias_ps + dep
      ips = ips + 1
    CASE(id_rh_obs)
      rmse_rh = rmse_rh + dep**2
      bias_rh = bias_rh + dep
      irh = irh + 1
    CASE(id_reflectivity_obs)
      rmse_ref= rmse_ref + dep**2
      bias_ref= bias_ref + dep
      iref    = iref + 1
    CASE(id_radialwind_obs)
      rmse_vr = rmse_vr  + dep**2
      bias_vr = bias_vr  + dep
      ivr     = ivr + 1
    END SELECT
  END DO

  IF(iu == 0) THEN
    rmse_u = undef
    bias_u = undef
  ELSE
    rmse_u = SQRT(rmse_u / REAL(iu,r_size))
    bias_u = bias_u / REAL(iu,r_size)
  END IF
  IF(iv == 0) THEN
    rmse_v = undef
    bias_v = undef
  ELSE
    rmse_v = SQRT(rmse_v / REAL(iv,r_size))
    bias_v = bias_v / REAL(iv,r_size)
  END IF
  IF(it == 0) THEN
    rmse_t = undef
    bias_t = undef
  ELSE
    rmse_t = SQRT(rmse_t / REAL(it,r_size))
    bias_t = bias_t / REAL(it,r_size)
  END IF
  IF(iq == 0) THEN
    rmse_q = undef
    bias_q = undef
  ELSE
    rmse_q = SQRT(rmse_q / REAL(iq,r_size))
    bias_q = bias_q / REAL(iq,r_size)
  END IF
  IF(ips == 0) THEN
    rmse_ps = undef
    bias_ps = undef
  ELSE
    rmse_ps = SQRT(rmse_ps / REAL(ips,r_size))
    bias_ps = bias_ps / REAL(ips,r_size)
  END IF
  IF(irh == 0) THEN
    rmse_rh = undef
    bias_rh = undef
  ELSE
    rmse_rh = SQRT(rmse_rh / REAL(irh,r_size))
    bias_rh = bias_rh / REAL(irh,r_size)
  END IF
  IF(iref == 0) THEN
    rmse_ref = undef
    bias_ref = undef
  ELSE
    rmse_ref = SQRT(rmse_ref / REAL(iref,r_size))
    bias_ref = bias_ref / REAL(iref,r_size)
  END IF
  IF(ivr == 0) THEN
    rmse_vr = undef
    bias_vr = undef
  ELSE
    rmse_vr = SQRT(rmse_vr / REAL(ivr,r_size))
    bias_vr = bias_vr / REAL(ivr,r_size)
  END IF


  WRITE(6,'(3A)') '== PARTIAL OBSERVATIONAL DEPARTURE (',file,') ========================================'
  WRITE(6,'(8A12)') 'U','V','T','Q','PS','RH','REF','VR'
  WRITE(6,'(8ES12.3)') bias_u,bias_v,bias_t,bias_q,bias_ps,bias_rh,bias_ref,bias_vr
  WRITE(6,'(8ES12.3)') rmse_u,rmse_v,rmse_t,rmse_q,rmse_ps,rmse_rh,rmse_ref,rmse_vr
  WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS ========================================================'
  WRITE(6,'(8A12)') 'U','V','T','Q','PS','RH','REF','VR'
  WRITE(6,'(8I12)') iu,iv,it,iq,ips,irh,iref,ivr
  WRITE(6,'(A)') '=================================================================================='

  RETURN

END SUBROUTINE monit_mean

END MODULE letkf_obs
