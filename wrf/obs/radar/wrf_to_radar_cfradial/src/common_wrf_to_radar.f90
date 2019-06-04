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
  !INTEGER :: method_ref_calc

  real(r_size),allocatable :: gues3d(:,:,:,:)
  real(r_size),allocatable :: gues2d(:,:,:)

  real(r_size)             :: min_ref_dbz
  
  CHARACTER*100  :: model_file_name='modelTTTT.nc'  !Model input file. T slots are for different times (or ensemble members)
  CHARACTER*100  :: inputradar='radarRRR_TTTT.nc'  !Radar data input.
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
  REAL(r_size)         :: ref , rv ,  cref  , refdb
  REAL(r_size)         :: pik !Path integrated attenuation coefficient.
  REAL(r_size)         :: prh !Pseudo relative humidity
  LOGICAL              :: ISALLOC
  REAL(r_size)         :: max_model_z
  REAL(r_size)         :: tmpref,tmperr(1)
  INTEGER              :: INTERPOLATION_TECHNIQUE

  write(*,*)"Hello from model_to_radar"

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
  !   of the power within the beam. (this one is not implemented yet)
  !Compute maximum model height

   min_ref_dbz=10.0*log10(min_ref)

   max_model_z=0.0d0
   DO i=1,nlon
     DO j=1,nlat
        IF( v3d(i,j,nlev-1,iv3d_ph).GT. max_model_z )max_model_z=v3d(i,j,nlev-1,iv3d_ph)
     ENDDO
   ENDDO


   !CALL get_method_refcalc( input_radar%lambda , method_ref_calc )

   ALLOCATE( input_radar%radarv3d_model(input_radar%nr,input_radar%na,input_radar%nv3d_model) )
   input_radar%radarv3d_model=input_radar%missing

   write(*,*)input_radar%nr , input_radar%na 

  !Begin with the interpolation. 

  IF(INTERPOLATION_TECHNIQUE .EQ. 1)THEN
  !SIMPLE AND FAST INTERPOLATION APPROACH.


!$OMP PARALLEL DO DEFAULT(SHARED) FIRSTPRIVATE(ia,ie,ir,pik,tmp_z,ri,rj,rk,qv,qc,qr,qci,qs,qg,t,p,u,v,w,ref,rv,att,cref,prh)

    DO ia=1,input_radar%na 
      pik=0.0d0 !Path integrated attenuation coefficient.
      DO ir=1,input_radar%nr 
 
        !PAWR can have data up to 60.000 m we speed up the forward operator
        !by checking this simple condition.  
        IF( input_radar%z(ia,ir) .LE. max_model_z )THEN

       !Do not change DO order, radial loop has to be at the end for 
       !Attenuation computation.

       !Find i,j,k for the center of the beam.
       tmp=REAL(id_reflectivity_obs,r_size)

       CALL latlon_to_ij(projection,input_radar%lat(ia,ir),input_radar%lon(ia,ir), ri,rj)

       !CALL  ll2ij(tmp,input_radar%lon(ia,ir,ie),input_radar%lat(ia,ir,ie),ri,rj) !This is the old routine.

       tmp_z=input_radar%z(ia,ir)
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
          !Rotate ur and rv
          CALL rotwind_letkf(u,v,input_radar%lon(ia,ir),1.0d0,projection)
 
          !Compute reflectivity at the beam center.
          CALL calc_ref_vr(qv,qc,qr,qci,qs,qg,u,v,w,t,p,           &
               input_radar%azimuth(ia),input_radar%elevation(ie)   &
                ,method_ref_calc,ref,rv,att)

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
            rv=rv+tmperr(1)*RADIALWIND_ERROR 

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

 

          IF( ref .GT. min_ref )THEN
            input_radar%radarv3d_model(ir,ia,input_radar%iv3d_ref)=10.0d0*log10(ref) !refdb
          ELSE
            input_radar%radarv3d_model(ir,ia,input_radar%iv3d_ref)=min_ref_dbz
          ENDIF

          input_radar%radarv3d_model(ir,ia,input_radar%iv3d_rv)=rv


        ENDIF  !Endif for domain check

       ENDIF   !Endif for QC_PASS

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

  write(*,*)"Good bye from model_to_radar"
  


  RETURN
END SUBROUTINE model_to_radar

SUBROUTINE interp_to_real_radar
implicit none
integer       :: im , iradar
integer       :: iyyyy,imm,idd,ihh,imn
real(r_size)  :: ss
logical       :: model_split_output

  !Get grid properties.
  write(*,*)"Getting model info"
  current_model_file=model_file_name
  im=1
  WRITE(current_model_file(6:9),'(I4.4)')im
  CALL set_common_wrf( current_model_file )
  ALLOCATE( gues3d(nlon,nlat,nlev,nv3d),gues2d(nlon,nlat,nv2d) )


  if( ntime .gt. 1 )then
    !The number of times in the file is greather than one.
    model_split_output = .false.
    !Check is the number of requested times is equal to the number 
    !of times inside the file.
    if( ntime /= n_times)then
      n_times=min( ntime , n_times)
      write(*,*)"Warning! number of times in file is different from n_times."
      write(*,*)"Reseting n_times to: ",n_times
    endif
  else
    model_split_output = .true.
  endif



 DO im=1,n_times
  !Por el momento no uso la fecha, pero mas adelante la idea es usarla para
  !hacer interpolacion 4D.
  !CALL read_date_wrf(current_model_file,iyyyy,imm,idd,ihh,imn,ss)

  !Read the model data.
   write(*,*)"Reading model data"
   if( .not. model_split_output ) then

     CALL read_grd(current_model_file,im,gues3d,gues2d)
     gues3d(:,:,:,iv3d_ph)=gues3d(:,:,:,iv3d_ph)/gg

   else

     WRITE(model_file_name(6:9),'(I4.4)')im
     CALL read_grd(current_model_file,1,gues3d,gues2d)
     gues3d(:,:,:,iv3d_ph)=gues3d(:,:,:,iv3d_ph)/gg

   endif


  DO iradar=1,n_radar

    current_radar_file = inputradar
    WRITE(current_radar_file(6:8),'(I3.3)')iradar
    WRITE(current_radar_file(10:13),'(I4.4)')im

    write(*,*)"Reading file ",current_radar_file

    CALL radar_set_common( RADAR_1 , current_radar_file  )

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

     write(*,*)"Writing output data"
     CALL radar_write_model( RADAR_1 , current_radar_file )

  ENDDO ![Endo over radars]
  
 ENDDO ![Endo over model files]

END SUBROUTINE interp_to_real_radar

SUBROUTINE interp_to_real_radar_member(my_model_file,my_radar_file,my_radar,my_time)
implicit none
character(100)  , intent(IN)  :: my_model_file , my_radar_file
integer         , intent(IN)  :: my_time
type(RADAR)     , intent(OUT) :: my_radar
real(r_size)                  :: my_gues3d(nlon,nlat,nlev,nv3d) , my_gues2d(nlon,nlat,nv2d)
integer                       :: iyyyy,imm,idd,ihh,imn
real(r_size)                  :: ss


  CALL read_date_wrf(my_model_file,iyyyy,imm,idd,ihh,imn,ss)

  CALL set_common_wrf(my_model_file)

  !Read the model data.
  CALL read_grd(current_model_file,my_time,my_gues3d,my_gues2d)
  my_gues3d(:,:,:,iv3d_ph)=my_gues3d(:,:,:,iv3d_ph)/gg


  CALL radar_read_data( my_radar , my_radar_file )

  CALL radar_georeference( my_radar )

  CALL model_to_radar( my_radar , my_gues3d , my_gues2d )


END SUBROUTINE interp_to_real_radar_member


END MODULE common_wrf_to_radar
