MODULE common_wrf_to_bufr
!=======================================================================
!
! [PURPOSE:] Create observations from a "nature" run.
!
! [HISTORY:]
!   01/01/2013 Juan Ruiz  created based on the code provided by Takemasa
!   Miyoshi.
!
!=======================================================================
  USE common
  USE common_wrf
  USE common_obs_wrf
  USE map_utils
  USE common_namelist

  IMPLICIT NONE

  INTEGER :: nobs,nobslev
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:,:)
  REAL(r_size) :: dz,tg,qg
  REAL(r_size) :: dlon1,dlon2,dlon,dlat
  REAL(r_size),ALLOCATABLE :: tmpelm(:)
  REAL(r_size),ALLOCATABLE :: tmplon(:)
  REAL(r_size),ALLOCATABLE :: tmplat(:)
  REAL(r_size),ALLOCATABLE :: tmplev(:)
  REAL(r_size),ALLOCATABLE :: tmpdat(:)
  REAL(r_size),ALLOCATABLE :: tmperr(:)
  REAL(r_size),ALLOCATABLE :: tmptyp(:)
  REAL(r_size),ALLOCATABLE :: tmpi(:)
  REAL(r_size),ALLOCATABLE :: tmpj(:)
  REAL(r_size),ALLOCATABLE :: tmpk(:)
  REAL(r_size),ALLOCATABLE :: tmpqc(:)
  INTEGER :: nobslots
  INTEGER :: n,i,j,ierr,islot,nn,l,im
  INTEGER :: iolen,recloutput
  real(r_size) :: ratio,tmpx,tmpy,tmpz
  INTEGER  :: nslots,count,ii,zz,jj,var
  INTEGER      :: NTOBS,NT2DOBS,NT3DOBS
  INTEGER      :: NU,NV,NT,NQ,NPS

  CHARACTER(9) :: obsfile='obs.dat'
  CHARACTER(9) :: outfile='out.dat'
  CHARACTER(11) :: guesfile='minput'


  CONTAINS


!--------------------------------------------------------------------------
! DETERMINE THE NUMBER OF OBSERVATIONS THAT WILL BE GENERATED

SUBROUTINE get_obs_num()
  IMPLICIT NONE


  NTOBS=0
  NT3DOBS=0
  NT2DOBS=0
  IF(UOBS)THEN
    NT3DOBS=NT3DOBS+1
    NU=NT3DOBS
  ENDIF
  IF(VOBS)THEN
    NT3DOBS=NT3DOBS+1
    NV=NT3DOBS
  ENDIF
  IF(TOBS)THEN
    NT3DOBS=NT3DOBS+1
    NT=NT3DOBS
  ENDIF
  IF(QOBS)THEN
    NT3DOBS=NT3DOBS+1
    NQ=NT3DOBS
  ENDIF
  IF(PSOBS)THEN
    NT2DOBS=NT2DOBS+1
    NPS=NT2DOBS
  ENDIF

  NTOBS=NT3DOBS+NT2DOBS


  IF ( OBSTYPE == 3 )THEN
     !READ OBS LOCATIONS FROM A FILE. GET NUMBER OF OBS FIRST.
     CALL get_nobs(obsfile,nobs)
     WRITE(*,*) nobs,' TOTAL OBSERVATIONS INPUT'
  ELSEIF ( OBSTYPE == 1 )THEN 
     !RANDOM LOCATION FOR THE OBSERVATIONS.
     nobs=(NINT(nij0*OBSDENSITY)*CEILING(REAL(nlev-1)/REAL(SKIPZ))*NT3DOBS) + (NINT(nij0*OBSDENSITY)*NT2DOBS)
     WRITE(*,*) NINT(nij0*OBSDENSITY),'OBSERVATIONS PER LEV'
     WRITE(*,*) nobs,' TOTAL OBSERVATIONS INPUT'
     WRITE(*,*) CEILING(REAL(nlev)/REAL(SKIPZ)),'TOTAL NUMBER OF LEVELS'
  ELSEIF ( OBSTYPE == 2)THEN
     nobs=(CEILING(REAL(nlon)/REAL(SKIP))*CEILING(REAL(nlat)/REAL(SKIP))*CEILING(REAL(nlev-1)/REAL(SKIPZ))*NT3DOBS) + (CEILING(REAL(nlon)/REAL(SKIP))*CEILING(REAL(nlat)/REAL(SKIP)))*NT2DOBS
     WRITE(*,*) CEILING(REAL(nlon)/REAL(SKIP))*CEILING(REAL(nlat)/REAL(SKIP)),'OBSERVATIONS PER LEV'
     WRITE(6,'(I10,A)') nobs,' TOTAL OBSERVATIONS INPUT'
  ENDIF

  !WRITE(*,*)"A total number of ", nobs , " observations will be generated."

END SUBROUTINE get_obs_num 
!

SUBROUTINE get_obs_ijk()
  IMPLICIT NONE
  INTEGER ii,jj,zz,var 
  REAL(r_size) :: randnumber(1)

  ALLOCATE( tmpelm(nobs) )
  ALLOCATE( tmplon(nobs) )
  ALLOCATE( tmplat(nobs) )
  ALLOCATE( tmplev(nobs) )
  ALLOCATE( tmpdat(nobs) )
  ALLOCATE( tmperr(nobs) )
  ALLOCATE( tmptyp(nobs) )
  ALLOCATE( tmpi(nobs) )
  ALLOCATE( tmpj(nobs) )
  ALLOCATE( tmpk(nobs) )
  ALLOCATE( tmpqc(nobs) )

!
! FIND THE VALUE OF i,j,k for EACH OBSERVATION.
!
  IF (OBSTYPE == 1 )THEN
  !RANDOM OBSERVATION LOCATION
     nobslev=NINT(nij0*OBSDENSITY)
     WRITE(*,*)nobs,nobslev
     count=1
     DO ii = 1,nobslev
       CALL com_rand(1,randnumber)
       tmpx=randnumber(1)
       CALL com_rand(1,randnumber)
       tmpy=randnumber(1)
       DO zz=1,nlev-1,SKIPZ
          DO var=1,NT3DOBS 
            tmpi(count)=(nlon-8)*tmpx+4            
            tmpj(count)=(nlat-8)*tmpy+4
            !Get obs lat and lon from i and j.
            CALL itpl_2d(lon,tmpi(n),tmpj(n),tmplon(n))
            CALL itpl_2d(lat,tmpi(n),tmpj(n),tmplat(n))

            tmpk(count)=zz
            tmptyp(count)=1
              IF(var == NU)THEN
                  tmpelm(count)=id_u_obs
                  tmperr(count)=U_ERROR
              ENDIF
              IF(var == NV)THEN 
                  tmpelm(count)=id_v_obs
                  tmperr(count)=V_ERROR
              ENDIF
              IF(var == NT)THEN
                  tmpelm(count)=id_t_obs
                  tmperr(count)=T_ERROR
              ENDIF
              IF(var == NQ)THEN
                  tmpelm(count)=id_q_obs
                  tmperr(count)=Q_ERROR
              ENDIF
              count=count+1
              ENDDO
          ENDDO 
               tmpi(count)=(nlon-8)*tmpx+4
               tmpj(count)=(nlat-8)*tmpy+4
               tmpk(count)=1
               tmptyp(count)=1
              DO var=1,NT2DOBS
                IF(var == NPS)THEN
                  tmpelm(count)=id_ps_obs
                  tmperr(count)=PS_ERROR
                ENDIF
                tmptyp(count)=1
                count=count+1
              ENDDO
       ENDDO

  ELSEIF (OBSTYPE == 2 ) THEN
     !REGULAR OBSERVATION GENERATOR.
     count=1
     DO ii = 1,nlon,SKIP
      DO jj= 1,nlat,SKIP
          DO var=1,NT3DOBS
            DO zz=1,nlev-1,SKIPZ
            tmpi(count)=ii
            tmpj(count)=jj
            tmpk(count)=zz
            tmplon(count)=lon(ii,jj)
            tmplat(count)=lat(ii,jj)
            tmptyp(count)=1
              IF(var == NU)THEN
                  tmpelm(count)=id_u_obs
                  tmperr(count)=U_ERROR
              ENDIF
              IF(var == NV)THEN
                  tmpelm(count)=id_v_obs
                  tmperr(count)=V_ERROR
              ENDIF
              IF(var == NT)THEN
                  tmpelm(count)=id_t_obs
                  tmperr(count)=T_ERROR
              ENDIF
              IF(var == NQ)THEN
                  tmpelm(count)=id_q_obs
                  tmperr(count)=Q_ERROR
              ENDIF
              count=count+1
              ENDDO
            ENDDO

              tmpi(count)=ii
              tmpj(count)=jj
              tmpk(count)=1
              tmplon(count)=lon(ii,jj)
              tmplat(count)=lat(ii,jj)
              tmptyp(count)=1
  
            DO var = 1,NT2DOBS
              IF(var == NPS)THEN
                  tmpelm(count)=id_ps_obs
                  tmperr(count)=PS_ERROR
              ENDIF
              count=count+1
            ENDDO
     ENDDO
    ENDDO

  ELSEIF (OBSTYPE == 3 ) THEN

    CALL read_obs(obsfile,nobs,tmpelm,tmplon,tmplat,tmplev,tmpdat,tmperr,tmptyp )
    DO n=1,nobs
       CALL latlon_to_ij(projection,tmplat(n),tmplon(n),tmpi(n),tmpj(n))
    ENDDO
    tmpk=0.0d0
  ENDIF


  !Get tmplev for random and regular obs and tmpk for input file obs.
  DO n=1,nobs

       IF( CEILING(tmpi(n)) < 2 .OR. nlon-1 < CEILING(tmpi(n)) .OR. CEILING(tmpj(n)) < 2 .OR. nlat-1 < CEILING(tmpj(n)) .OR. CEILING(tmpk(n)) > nlev )THEN
         CYCLE
       ENDIF
       IF( obstype == 3 )THEN
        !Obtain tmpk from tmplev
        CALL p2k(v3d(:,:,:,iv3d_p),tmpelm(n),tmpi(n),tmpj(n),tmplev(n),tmpk(n))
       ELSE
        !Obtain tmplev from tmpk
        CALL itpl_3d(v3d(:,:,:,iv3d_p),tmpi(n),tmpj(n),tmpk(n),tmplev(n))
       ENDIF
  ENDDO 

  !Now independently of the observation distribution type we have tmpi,tmpj,tmpk
  !and tmplon,tmplat,tmplev

END SUBROUTINE get_obs_ijk



!Write the pseudo-observations in the LETKF input format.
SUBROUTINE write_obs(cfile)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  REAL(r_sngl) :: wk(7) 
  INTEGER :: n,iunit,writenobs

  iunit=91
  writenobs=0


  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')

  DO n=1,nobs
    IF(tmpdat(n) == undef )CYCLE

    !IF( tmpelm(n) == id_q_obs )THEN
    !  WRITE(*,*)tmpelm(n),tmpdat(n),tmperr(n)
    !ENDIF

    writenobs=writenobs+1
    wk(1)=REAL(tmpelm(n),r_sngl)
    wk(2)=REAL(tmplon(n),r_sngl)
    wk(3)=REAL(tmplat(n),r_sngl)
    wk(4)=REAL(tmplev(n),r_sngl)
    wk(5)=REAL(tmpdat(n),r_sngl)
    wk(6)=REAL(tmperr(n),r_sngl)
    wk(7)=REAL(tmptyp(n),r_sngl)
    SELECT CASE(NINT(wk(1)))
    CASE(id_u_obs,id_v_obs,id_t_obs,id_tv_obs,id_q_obs)
      wk(4) = wk(4) /  100.0 ! hPa <- Pa
    CASE(id_ps_obs)
      wk(5) = wk(5) / 100.0 ! hPa <- Pa
      wk(6) = wk(6) / 100.0 ! hPa <- Pa
    CASE(id_rh_obs)
      wk(4) = wk(4) / 100.0 ! hPa <- Pa
      wk(5) = wk(5) / 0.01 ! percent input
      wk(6) = wk(6) / 0.01 ! percent input
    CASE(id_tcmip_obs)
      wk(5) = wk(5) / 100.0 ! hPa <- Pa
      wk(6) = wk(6) / 100.0 ! hPa <- Pa
    END SELECT
    WRITE(iunit)wk

    !IF( tmpelm(n) == id_q_obs )THEN
    !  WRITE(*,*)wk
    !ENDIF

  END DO
  CLOSE(iunit)

  WRITE(*,*)"A total number of ",writenobs," observations have been written"
RETURN


END SUBROUTINE write_obs


!========================================================================================
!THIS ROUTINE READS ONE ENSEMBLE MEMBER AND GOES FROM THE MODEL STATE SPACE TO
!THE 
!OBSERVATIONAL SPACE.
!========================================================================================

SUBROUTINE model_to_obs()
IMPLICIT NONE
INTEGER(4)                  :: nxvar,nyvar,nzvar,k,iv,n,ncid,varid,obsid,ii
REAL(r_sngl)                :: fieldg(nlon,nlat),auxfieldg(nlon,nlat)
REAL(r_size)                :: dz,tg,qg,dummy(3)
REAL(r_size)                :: randnumber(nobs)
REAL(r_size)                :: tmp_val


 IF( add_obs_error)THEN
   CALL com_randn(nobs,randnumber)
 ENDIF

 tmpdat= undef

     DO n=1,nobs
 
        IF( CEILING(tmpi(n)) < 2 .OR. nlon-1 < CEILING(tmpi(n)) .OR. CEILING(tmpj(n)) < 2 .OR. nlat-1 < CEILING(tmpj(n)) .OR. CEILING(tmpk(n)) > nlev )THEN
          CYCLE
        ENDIF
        IF( tmpk(n) < 1.0 )THEN
          tmpk(n)=1.0
        ENDIF

        tmp_val=tmpdat(n) !Original obs value.

        CALL Trans_XtoY(tmpelm(n),tmptyp(n),tmplon(n),tmplat(n),tmpi(n),tmpj(n),tmpk(n),0.0d0,0.0d0,v3d,v2d,tmpdat(n))

        
        IF( NINT(tmpelm(n)) == id_ps_obs .AND. tmpdat(n) /= UNDEF )THEN
          CALL itpl_2d(phi0,tmpi(n),tmpj(n),tmplev(n))
        ENDIF

        IF( add_obs_error )THEN
         tmpdat(n)=tmpdat(n)+randnumber(n)*tmperr(n)
        ENDIF
  
        IF( add_obs_error .AND. tmpelm(n) == id_q_obs )THEN
          IF(tmpdat(n) < 0.0d0 )tmpdat(n)=0.0d0
        ENDIF


      END DO !End over the computation of i,j, and k.


END SUBROUTINE model_to_obs

END MODULE common_wrf_to_bufr

