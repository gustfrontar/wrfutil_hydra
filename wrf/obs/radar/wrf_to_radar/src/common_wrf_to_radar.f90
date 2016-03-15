MODULE common_wrf_to_radar
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
  USE common_radar_tools
  USE common_wrf
  USE common_obs_wrf

  IMPLICIT NONE
  PUBLIC

  !Namelist variables
  LOGICAL            :: ADD_OBS_ERROR=.TRUE.
  REAL(r_size)       :: REFLECTIVITY_ERROR=1.0 !Standard deviation of reflectivity error in DBZ
  REAL(r_size)       :: RADIALWIND_ERROR=0.5   !Standard deviation of radialwind error in m/s

  !--------------------------------------------FAKE RADAR PARAMETERS
  LOGICAL   ::   FAKE_RADAR=.TRUE. !We will create a non existant radar.
  REAL(r_size), PARAMETER :: radar_lon = -63.4 , radar_lat = -34.37 !Grid point position of radar.
  REAL(r_size), PARAMETER :: radar_z = 0.0d0
  REAL(r_size), PARAMETER :: radar_az_res = 1.0d0 , radar_r_res = 500.0d0 , radar_el_res = 1.25d0
  REAL(r_size), PARAMETER :: radar_min_az = 0.0d0 , radar_max_az = 360.0d0
  REAL(r_size), PARAMETER :: radar_min_r  = 0.0d0 , radar_max_r  = 240000
  REAL(r_size), PARAMETER :: radar_min_el = 0.1d0 , radar_max_el = 18.0d0

  !Define the type of forward operator.
  LOGICAL :: USE_WT=.TRUE.       !Are we going to use the fall speed of hydrometeors in the computation of RADIAL VELOCITY?

  REAL(r_size) :: radar_lambda=3.0d0  !Wave length of the radar (to be used in the computation of reflectivity.
                                      !3 Xband , 5 Cband , 10 Sband

CONTAINS

!-----------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------

SUBROUTINE read_namelist

NAMELIST / GENERAL / add_obs_error , reflectivity_error , radialwind_error , &
                     fake_radar , radar_lon , radar_lat , radar_z , &
                     radar_min_res , radar_r_res , radar_el_res   , &
                     radar_min_az  , radar_max_az   , &
                     radar_min_el  , radar_max_l    , &
                     radar_min_r   , radar_max_r    , &
                     use_wt , radar_lambda 

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
! Interpolate model data to a radar (or part of a radar) domain
!-----------------------------------------------------------------------
SUBROUTINE model_to_radar( input_radar , v3d , v2d  )
  IMPLICIT NONE
  TYPE(RADAR), INTENT(INOUT) :: input_radar
  REAL(r_size),INTENT(IN)    :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN)    :: v2d(nlon,nlat,nlev,nv2d)
  INTEGER              :: ia , ir , ie , i , j , k
  REAL(r_size)         :: ri , rj , rk , rlev
  REAL(r_size)         :: tmp_z , att , tmp
  REAL(r_size)         :: p , t
  REAL(r_size)         :: qv  , qc , qr
  REAL(r_size)         :: qci , qs , qg 
  REAL(r_size)         :: u   , v  , w
  REAL(r_size)         :: ref , vr ,  cref 
  REAL(r_size)         :: pik !Path integrated attenuation coefficient.
  REAL(r_size)         :: prh !Pseudo relative humidity
  LOGICAL              :: ISALLOC
  REAL(r_size)         :: max_model_z
  REAL(r_size)         :: tmpref,tmperr(1)
  

  !NOTE: input_radar can be the full radar domain or a local radar domain of the same type.
  !This means that this algorithm can be paralelized in AZIMUTH or ELEVATION. Due to the 
  !possibility of radiative transfer extention, paralelization in range is not recomended.

  !Different approaches can be taken. 
  !
  !1) Simple and fast but dangerous: Interpolate condensates / wind to the
  !   position of the center of the beam and then compute radial velocity
  !   and reflectivity.( INTERPOLATION_TECHNIQUE == 1)
  !
  !2) Consider a ray radiative transfer model for several paths at the center and sourronding
  !   the beam. Compute a simple radiative transfer for each path and then perform a weigthed average.
  !   In this case the attenuation can be considered in a more robust way. And the contribution
  !   of each ray to the total radar reflectivity can be weighted according to the distribution
  !   of the power within the beam.


   !Compute maximum model height
   max_model_z=0.0d0
   DO i=1,nlon
     DO j=1,nlat
        IF( v3d(i,j,nlev-1,iv3d_ph).GT. max_model_z )max_model_z=v3d(i,j,nlev-1,iv3d_ph)
     ENDDO
   ENDDO

  IF ( COMPUTE_ATTENUATION .AND. METHOD_REF_CALC .LE. 2)THEN
    WRITE(6,*)'WARNING: Attenuation can not be computed with methods 1 or 2.'
    COMPUTE_ATTENUATION=.false.
  ENDIF

  !Inquire radar band
  IF ( input_radar%lambda == 3.0d0 )THEN
     method_ref_calc=3
  ELSEIF( input_radar%lambda == 5.0d0 )THEN
     method_ref_calc=4
     WRITE(*,*)"[Error]: C-Band reflectivity has not been coded yet"
     STOP
  ELSEIF( input_radar%lambda == 10.0d0 )THEN
     method_ref_calc=2
  ELSE
     WRITE(*,*)"[Error]: Radar lambda not recognized", input_radar%lambda
  ENDIF

   ALLOCATE( input_radar%radarv3d_model(input_radar%na,input_radar%nr,input_radar%ne,input_radar%nv3d) )
   input_radar%radarv3d_model=input_radar%missing

  !Begin with the interpolation. 

  IF   (INTERPOLATION_TECHNIQUE .EQ. 1)THEN
    !SIMPLE AND FAST INTERPOLATION APPROACH.

!$OMP PARALLEL DO DEFAULT(SHARED) FIRSTPRIVATE(ia,ie,ir,pik,tmp_z,ri,rj,rk,qv,qc,qr,qci,qs,qg,t,p,u,v,w,ref,vr,att,cref,prh)

    DO ia=1,input_radar%na 
     DO ie=1,input_radar%ne
      pik=0.0d0 !Path integrated attenuation coefficient.
      DO ir=1,input_radar%nr 
 
        !PAWR can have data up to 60.000 m we speed up the forward operator
        !by checking this simple condition.  
        IF( input_radar%z(ia,ir,ie) .LE. max_model_z )THEN

       !Do not change DO order, radial loop has to be at the end for 
       !Attenuation computation.

       !Find i,j,k for the center of the beam.
       tmp=REAL(id_reflectivity_obs,r_size)


       CALL latlon_to_ij(projection,input_radar%lat(ia,ir,ie),input_radar%lon(ia,ir,ie), ri,rj)

       !CALL  ll2ij(tmp,input_radar%lon(ia,ir,ie),input_radar%lat(ia,ir,ie),ri,rj) !This is the old routine.

       tmp_z=input_radar%z(ia,ir,ie)
       CALL  z2k_fast(v3d(:,:,:,iv3d_ph),ri,rj,tmp_z,rk)
   
        !Interpolate qv,qc,qr,qci,qs,qg,t and p to the beam center.
        IF( ri >= 1 .AND. ri <= nlon .AND. rj >= 1 .AND.    &
            rj <= nlat .AND. rk >= 1 .AND. rk <= nlev )THEN

          CALL itpl_3d(v3d(:,:,:,iv3d_qv),ri,rj,rk,qv)
          CALL itpl_3d(v3d(:,:,:,iv3d_qc),ri,rj,rk,qc)
          CALL itpl_3d(v3d(:,:,:,iv3d_qr),ri,rj,rk,qr)
          CALL itpl_3d(v3d(:,:,:,iv3d_qci),ri,rj,rk,qci)
          CALL itpl_3d(v3d(:,:,:,iv3d_qs),ri,rj,rk,qs)
          CALL itpl_3d(v3d(:,:,:,iv3d_qg),ri,rj,rk,qg)
          CALL itpl_3d(v3d(:,:,:,iv3d_t),ri,rj,rk,t)
          CALL itpl_3d(v3d(:,:,:,iv3d_p),ri,rj,rk,p)
          !Interpolate winds.
          CALL itpl_3d(v3d(:,:,:,iv3d_u),ri+0.5,rj,rk,u)
          CALL itpl_3d(v3d(:,:,:,iv3d_v),ri,rj+0.5,rk,v)
          CALL itpl_3d(v3d(:,:,:,iv3d_w),ri,rj,rk+0.5,w)
          !Rotate ur and vr
          CALL rotwind_letkf(ur,vr,input_radar%lon(ia,ir,ie),1.0d0,projection)

 
          !Compute reflectivity at the beam center.
          CALL calc_ref_vr(qv,qc,qr,qci,qs,qg,u,v,w,t,p,           &
               input_radar%azimuth(ia),input_radar%elevation(ie)   &
                ,method_ref_calc,ref,vr,att)
   
          !ADD ERRORS TO THE OBSERVATIONS
          IF( ADD_OBS_ERROR )THEN
            !Add error to reflectivity.
            IF( ref > minz)THEN
             tmpref=10.0d0*log10(ref)
             CALL com_randn(1,tmperr)
             tmpref=tmpref+tmperr(1)*REFLECTIVITY_ERROR
             ref=10.0d0**(tmpref/10.0d0) 
             IF( ref < minz )ref=minz
            ENDIF
        
            !Add error to doppler velocity.
            CALL com_randn(1,tmperr)
            vr=vr+tmperr(1)*RADIALWIND_ERROR 

          ENDIF

          IF( COMPUTE_ATTENUATION )THEN
          !If compute attenuation then correct the reflectivity using the PIK factor.
            IF( ref .GT. minz )THEN
             cref = ref * EXP( -0.46d0 * pik )
             IF(cref .GT. minz)THEN
               ref=cref
             ELSE
               ref=minz
             ENDIF
            ENDIF
            !Update PIK
            pik = pik + att * input_radar%range_resolution / 1.0d3
          ENDIF
 
          IF( ref .GT. minz )THEN
          input_radar%radarv3d_model(ia,ir,ie,input_radar%iv3d_ref)=ref !refdb
          ENDIF

          !Will generate wind observations only where reflectivity data is good enough (in this case were we have clouds). 
          IF( ref .GT. minz )THEN
          input_radar%radarv3d_model(ia,ir,ie,input_radar%iv3d_wind)=vr    !Radial wind
          ENDIF

        ENDIF  !Endif for domain check

       ENDIF   !Endif for QC_PASS

      ENDDO
     ENDDO
    ENDDO
!$OMP END PARALLEL DO

  ELSEIF(INTERPOLATION_TECHNIQUE .EQ. 2 )THEN
    WRITE(6,*)'OOOPS: Not coded yet'
    STOP


  ELSE

    WRITE(6,*)'ERROR: Not recognized interpolation technique for radar data'
    STOP

  ENDIF
  


  RETURN
END SUBROUTINE model_to_radar


END MODULE common_wrf_to_radar
