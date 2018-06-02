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
  USE common_namelist

  IMPLICIT NONE
  PUBLIC

  TYPE(radar) :: RADAR_1
  INTEGER :: method_ref_calc

  real(r_size),allocatable :: gues3d(:,:,:,:)
  real(r_size),allocatable :: gues2d(:,:,:)
  
  CHARACTER*100  :: model_file_name='input_modelTTTT.nc'  !Model input file. T slots are for different times (or ensemble members)
  CHARACTER*100  :: inputradar='iradarRRR_TTTT.grd'  !Radar data input.
  CHARACTER*100  :: outputradar='oradarRRRR_TTTT.grd'!output files in radar format.(R slots are for radar number, T slots are for time number.)
  CHARACTER*100  :: current_model_file  , current_radar_file                  

CONTAINS

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
  INTEGER              :: INTERPOLATION_TECHNIQUE

  INTERPOLATION_TECHNIQUE=1 

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


   CALL get_method_refcalc( input_radar%lambda , method_ref_calc )

   ALLOCATE( input_radar%radarv3d_model(input_radar%na,input_radar%nr,input_radar%ne,input_radar%nv3d) )
   input_radar%radarv3d_model=input_radar%missing

  !Begin with the interpolation. 

  IF(INTERPOLATION_TECHNIQUE .EQ. 1)THEN
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
          CALL rotwind_letkf(u,v,input_radar%lon(ia,ir,ie),1.0d0,projection)
 
          !Compute reflectivity at the beam center.
          CALL calc_ref_vr(qv,qc,qr,qci,qs,qg,u,v,w,t,p,           &
               input_radar%azimuth(ia),input_radar%elevation(ie)   &
                ,method_ref_calc,ref,vr,att)

           WRITE(*,*)vr
   
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

          !IF( COMPUTE_ATTENUATION )THEN
          !We need to take into account the radar wavelength
          !!If compute attenuation then correct the reflectivity using the PIK factor.
          !  IF( ref .GT. minz )THEN
          !   cref = ref * EXP( -0.46d0 * pik )
          !   IF(cref .GT. minz)THEN
          !     ref=cref
          !   ELSE
          !     ref=minz
          !   ENDIF
          !  ENDIF
          !  !Update PIK
          !  pik = pik + att * input_radar%range_resolution / 1.0d3
          !ENDIF
 
          !IF( ref .GT. minz )THEN
          input_radar%radarv3d_model(ia,ir,ie,input_radar%iv3d_ref)=ref !refdb
          !ENDIF

          !Will generate wind observations only where reflectivity data is good enough (in this case were we have clouds). 
          !IF( ref .GT. minz )THEN
          input_radar%radarv3d_model(ia,ir,ie,input_radar%iv3d_vr)=vr    !Radial wind
          !ENDIF

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

SUBROUTINE interp_to_fake_radar
implicit none
integer       :: im , iradar , i 
integer       :: iyyyy,imm,idd,ihh,imn
real(r_size)  :: ss , tmp

 DO im=1,n_model
  current_model_file=model_file_name
  WRITE(current_model_file(12:15),'(I4.4)')im

  CALL set_common_wrf(current_model_file)
  ALLOCATE( gues3d(nlon,nlat,nlev,nv3d),gues2d(nlon,nlat,nv2d) )

  !Read the model data.
  CALL read_grd(current_model_file,gues3d,gues2d)
  gues3d(:,:,:,iv3d_ph)=gues3d(:,:,:,iv3d_ph)/gg

  CALL read_date_wrf(current_model_file,iyyyy,imm,idd,ihh,imn,ss)
  RADAR_1 %year  =REAL(iyyyy,r_size)
  RADAR_1 %month =REAL(imm,  r_size)
  RADAR_1 %day   =REAL(idd,  r_size)
  RADAR_1 %hour  =REAL(ihh,  r_size)
  RADAR_1 %minute=REAL(imn,  r_size)
  RADAR_1 %second=ss


  DO iradar=1,n_radar

    !Dealing with a fake radar, lets create the radar location and characteristics.
    RADAR_1 % lon0 = fradar_lon(iradar)
    RADAR_1 % lat0 = fradar_lat(iradar)
    RADAR_1 % z0   = fradar_z(iradar)
    RADAR_1 % nv3d = 2
    RADAR_1 % radar_type = 1

    tmp= ( radar_max_az(iradar)-radar_min_az(iradar) ) / radar_az_res(iradar)
    RADAR_1 % na = NINT(tmp)
    tmp= ( radar_max_r(iradar)-radar_min_r(iradar) ) / radar_r_res(iradar)
    RADAR_1 % nr = NINT(tmp)
    tmp= ( radar_max_el(iradar) - radar_min_el(iradar) ) / radar_el_res(iradar)
    RADAR_1 % ne = NINT(tmp)

    ALLOCATE( RADAR_1 % azimuth( RADAR_1 % na) )
    ALLOCATE( RADAR_1 % rrange( RADAR_1 % nr) )
    ALLOCATE( RADAR_1 % elevation( RADAR_1 % ne) )

    DO i=1,RADAR_1 % na
       RADAR_1 % azimuth (i) = (i-1)*radar_az_res(iradar) + radar_min_az(iradar)
    ENDDO
    DO i=1,RADAR_1 % nr
       RADAR_1 % rrange (i) = (i-1)*radar_r_res(iradar) + radar_min_r(iradar)
    ENDDO
    DO i=1,RADAR_1 % ne
      RADAR_1 % elevation (i) = (i-1)*radar_el_res(iradar) + radar_min_el(iradar)
    ENDDO

    ALLOCATE(RADAR_1 % radarv3d(RADAR_1%na,RADAR_1%nr,RADAR_1%ne,RADAR_1%nv3d))
    ALLOCATE(RADAR_1 % qcflag(RADAR_1%na,RADAR_1%nr,RADAR_1%ne))
    ALLOCATE(RADAR_1 % attenuation(RADAR_1%na,RADAR_1%nr,RADAR_1%ne))
 
     RADAR_1 % radarv3d         = 0.0d0
     RADAR_1 % qcflag           = 0
     RADAR_1 % attenuation      = 1.0d0
     RADAR_1 % range_resolution = radar_r_res(iradar)

     ALLOCATE(RADAR_1%beam_wid_h(RADAR_1%na),RADAR_1%beam_wid_v(RADAR_1%ne) )

     RADAR_1%beam_wid_h      = 1.0d0  !Assume this constant.
     RADAR_1%beam_wid_v      = 1.0d0  !Assume this constant.
     RADAR_1%lambda          = radar_lambda(iradar)
     RADAR_1%missing         = undefs
 
     RADAR_1%total_attenuation_factor = 1.0d0
     RADAR_1%iv3d_ref=1    !Reflectivity
     RADAR_1%iv3d_vr =2    !Radial velocity


     WRITE(6,*)' FINISH SETTING THE RADAR '
     WRITE(6,*)' FAKE RADAR SUMMARY FOR RADAR=',iradar
     WRITE(6,*)' NA   ',RADAR_1%na
     WRITE(6,*)' NR   ',RADAR_1%nr
     WRITE(6,*)' NE   ',RADAR_1%ne
     WRITE(6,*)' RADAR LON ',radar_1%lon0
     WRITE(6,*)' RADAR LAT ',radar_1%lat0
     WRITE(6,*)' RADAR HEIGHT ',radar_1%z0
     WRITE(6,*)' INTERPOLATE MODEL OUTPUT TO THE CURRENT RADAR'

     CALL radar_georeference( RADAR_1 )

     CALL model_to_radar( RADAR_1 , gues3d , gues2d )


     current_radar_file=outputradar
     WRITE(current_radar_file(7:10),'(I4.4)')iradar
     WRITE(current_radar_file(12:15),'(I4.4)')im

     CALL radar_write_file( radar_1 , radar_1%radarv3d_model(:,:,:,radar_1%iv3d_ref), radar_1%radarv3d_model(:,:,:,radar_1%iv3d_vr), &
                                      radar_1%qcflag , radar_1%attenuation ,  current_radar_file , endian )

  ENDDO ![Endo over radars]
  
 ENDDO ![Endo over model files]

END SUBROUTINE interp_to_fake_radar


SUBROUTINE interp_to_real_radar
implicit none
integer       :: im , iradar
integer       :: iyyyy,imm,idd,ihh,imn
real(r_size)  :: ss

 DO im=1,n_model
  current_model_file=model_file_name
  WRITE(current_model_file(12:15),'(I4.4)')im

  CALL read_date_wrf(current_model_file,iyyyy,imm,idd,ihh,imn,ss)

  CALL set_common_wrf(current_model_file)
  ALLOCATE( gues3d(nlon,nlat,nlev,nv3d),gues2d(nlon,nlat,nv2d) )

  !Read the model data.
  CALL read_grd(current_model_file,gues3d,gues2d)
  gues3d(:,:,:,iv3d_ph)=gues3d(:,:,:,iv3d_ph)/gg


  DO iradar=1,n_radar

    current_radar_file=inputradar
    WRITE(current_radar_file(7:10),'(I4.4)')iradar
    WRITE(current_radar_file(12:15),'(I4.4)')im

    RADAR_1 % radar_type =1
    CALL radar_read_data( RADAR_1 , current_radar_file , endian )
    
    radar_1%qcflag      = 0
    radar_1%attenuation = 1.0d0


     WRITE(6,*)' FINISH READING THE RADAR '
     WRITE(6,*)' REAL RADAR SUMMARY FOR RADAR=',iradar
     WRITE(6,*)' NA   ',RADAR_1%na
     WRITE(6,*)' NR   ',RADAR_1%nr
     WRITE(6,*)' NE   ',RADAR_1%ne
     WRITE(6,*)' RADAR LON ',radar_1%lon0
     WRITE(6,*)' RADAR LAT ',radar_1%lat0
     WRITE(6,*)' RADAR HEIGHT ',radar_1%z0
     WRITE(6,*)' INTERPOLATE MODEL OUTPUT TO THE CURRENT RADAR'

     CALL radar_georeference( RADAR_1 )

     CALL model_to_radar( RADAR_1 , gues3d , gues2d )

     RADAR_1 %year  =REAL(iyyyy,r_size)
     RADAR_1 %month =REAL(imm,  r_size)
     RADAR_1 %day   =REAL(idd,  r_size)
     RADAR_1 %hour  =REAL(ihh,  r_size)
     RADAR_1 %minute=REAL(imn,  r_size)
     RADAR_1 %second=ss
   
     current_radar_file=outputradar
     WRITE(current_radar_file(7:10),'(I4.4)')iradar
     WRITE(current_radar_file(12:15),'(I4.4)')im

     CALL radar_write_file( radar_1 , radar_1%radarv3d_model(:,:,:,radar_1%iv3d_ref), radar_1%radarv3d_model(:,:,:,radar_1%iv3d_vr), &
                                      radar_1%qcflag , radar_1%attenuation ,  current_radar_file , endian )

  ENDDO ![Endo over radars]
  
 ENDDO ![Endo over model files]

END SUBROUTINE interp_to_real_radar




END MODULE common_wrf_to_radar
