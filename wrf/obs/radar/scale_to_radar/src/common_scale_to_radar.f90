MODULE common_scale_to_radar
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
  USE common_scale
  USE common_obs_scale
  USE common_namelist

  IMPLICIT NONE
  PUBLIC

  TYPE(radar) :: RADAR_1
  !INTEGER :: method_ref_calc

  real(r_size),allocatable :: gues3d(:,:,:,:)
  real(r_size),allocatable :: gues2d(:,:,:)
  
  CHARACTER*100  :: model_file_prefix = 'modelTTTT'  !Model input file. T slots are for different times (or ensemble members)
  CHARACTER*100  :: radar_file = 'radarRRR_TTTT.nc' !Radar data input in cfradial format
  CHARACTER*100  :: topo_file_prefix = 'topo'       !Model topography input file prefix.
  CHARACTER*100  :: current_model_prefix  , current_radar_file                  

  real(r_size)  :: min_ref_dbz

CONTAINS

!-----------------------------------------------------------------------
! Interpolate model data to a radar (or part of a radar) domain
!-----------------------------------------------------------------------
SUBROUTINE model_to_radar( input_radar , v3d , v2d  )
  IMPLICIT NONE
  TYPE(RADAR), INTENT(INOUT) :: input_radar
  REAL(r_size),INTENT(IN)    :: v3d(nlev,nlonh,nlath,nv3d)
  REAL(r_size),INTENT(IN)    :: v2d(nlonh,nlath,nv2d)
  INTEGER              :: ia , ir , ie , i , j , k
  REAL(r_size)         :: ri , rj , rk , rlev
  REAL(r_size)         :: tmp_z ,  tmp
  REAL(r_size)         :: p , t
  REAL(r_size)         :: qv  , qc , qr
  REAL(r_size)         :: qi  , qs , qg 
  REAL(r_size)         :: u   , v  , w
  REAL(r_size)         :: ref , rv ,  cref , refdb , minrefdb
  REAL(r_size)         :: pik !Path integrated attenuation coefficient.
  REAL(r_size)         :: prh !Pseudo relative humidity
  LOGICAL              :: ISALLOC
  REAL(r_size)         :: max_model_z
  REAL(r_size)         :: tmpref,tmperr(1)
  INTEGER              :: INTERPOLATION_TECHNIQUE

  REAL(r_size)         :: test_lon , test_lat

  INTERPOLATION_TECHNIQUE=1 

  minrefdb=10.0*log10(min_ref)

  ALLOCATE( input_radar%radarv3d_model(input_radar%nr,input_radar%nt,input_radar%nv3d) )
  input_radar%radarv3d_model=undef

  !Begin with the interpolation. 

  IF(INTERPOLATION_TECHNIQUE .EQ. 1)THEN
  !SIMPLE AND FAST INTERPOLATION APPROACH.

  max_model_z = 0.0d0
  do i = 1 , nlonh
    do j = 1 , nlath 
       if( v3d(nlev,i,j,iv3d_hgt) > max_model_z )then
         max_model_z = v3d(nlev,i,j,iv3d_hgt)
       endif
    enddo
  enddo


!$OMP PARALLEL DO DEFAULT(SHARED) FIRSTPRIVATE(ia,ir,pik,tmp_z,ri,rj,rk,qv,qc,qr,qci,qs,qg,t,p,u,v,w,ref,rv,cref,prh)

    DO ia=1,input_radar%na 
      pik=0.0d0 !Path integrated attenuation coefficient.
      DO ir=1,input_radar%nr 
 
        !PAWR can have data up to 60.000 m we speed up the forward operator
        !by checking this simple condition.  
        IF( input_radar%z(ia,ir) .LE. max_model_z )THEN

        CALL  phys2ij(input_radar%lon(ia,ir),input_radar%lat(ia,ir),ri,rj)
        CALL  phys2ijkz(v3d(:,:,:,iv3d_hgt),ri,rj,input_radar%z(ia,ir),rk)
 
        !For debug observation location.
        !call itpl_2d(v2d(:,:,iv2d_lon),ri,rj,test_lon)
        !call itpl_2d(v2d(:,:,iv2d_lat),ri,rj,test_lat)
        !WRITE(*,*)input_radar%lon(ia,ir),input_radar%lat(ia,ir),test_lon-360.0d0,test_lat
        !STOP

        if( rk /= undef )THEN
   
           !Interpolate qv,qc,qr,qci,qs,qg,t and p to the beam center.

           CALL itpl_3d(v3d(:,:,:,iv3d_u),rk,ri-0.5d0,rj,u)  
           CALL itpl_3d(v3d(:,:,:,iv3d_v),rk,ri,rj-0.5d0,v) 
           CALL itpl_3d(v3d(:,:,:,iv3d_w),rk-0.5d0,ri,rj,w) 
           CALL itpl_3d(v3d(:,:,:,iv3d_t),rk,ri,rj,t)
           CALL itpl_3d(v3d(:,:,:,iv3d_p),rk,ri,rj,p)
           CALL itpl_3d(v3d(:,:,:,iv3d_q),rk,ri,rj,qv)
           CALL itpl_3d(v3d(:,:,:,iv3d_qc),rk,ri,rj,qc)
           CALL itpl_3d(v3d(:,:,:,iv3d_qr),rk,ri,rj,qr)
           CALL itpl_3d(v3d(:,:,:,iv3d_qi),rk,ri,rj,qi)
           CALL itpl_3d(v3d(:,:,:,iv3d_qs),rk,ri,rj,qs)
           CALL itpl_3d(v3d(:,:,:,iv3d_qg),rk,ri,rj,qg)
 
            !Compute reflectivity at the beam center.
           CALL calc_ref_rv(qv,qc,qr,qi,qs,qg,u,v,w,t,p,            &
               input_radar%azimuth(ia),input_radar%elevation(ia)    &
                ,ref,rv)
   
           !ADD ERRORS TO THE OBSERVATIONS
           IF( ADD_OBS_ERROR )THEN
             !Add error to reflectivity.
             IF( ref > min_ref)THEN
              tmpref=10.0d0*log10(ref)
              CALL com_randn(1,tmperr)
              tmpref=tmpref+tmperr(1)*REFLECTIVITY_ERROR
              ref=10.0d0**(tmpref/10.0d0) 
              IF( ref < min_ref )ref=min_ref
             ENDIF
        
             !Add error to doppler velocity.
             CALL com_randn(1,tmperr)
             rv=rv+tmperr(1)*RADIALWIND_ERROR 

           ENDIF

          refdb=10.0d0*log10(ref)

          IF( ref .GT. min_ref )THEN
              input_radar%radarv3d_model(ir,ia,input_radar%iv3d_ref)=refdb !refdb
          ELSE
              input_radar%radarv3d_model(ir,ia,input_radar%iv3d_ref)=minrefdb
          ENDIF

          input_radar%radarv3d_model(ir,ia,input_radar%iv3d_rv)=rv    !Radial wind


        ENDIF  !Domain check

       ENDIF   !Max height check

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


SUBROUTINE interp_to_real_radar
implicit none
integer       :: im , iradar
integer       :: iyyyy,imm,idd,ihh,imn
real(r_size)  :: ss
logical       :: model_split_output

  !Get grid properties.
  write(*,*)"Getting model info"
  current_model_prefix=model_file_prefix
  WRITE(current_model_prefix(6:9),'(I4.4)')1
  CALL set_common_scale( current_model_prefix , topo_file_prefix , input_type )

  if( ntime .gt. 1 )then
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

  ALLOCATE( gues3d(nlev,nlonh,nlath,nv3d),gues2d(nlonh,nlath,nv2d) )

!We will process the first n_times as requested from namelist.
DO im=1,n_times

   write(*,*)"Reading model data"
   if( .not. model_split_output ) then
 
     CALL read_file(current_model_prefix,gues3d,gues2d,im,input_type)

   else

     WRITE(current_model_prefix(6:9),'(I4.4)')im
     CALL read_file(current_model_prefix,gues3d,gues2d,1,input_type)

   endif

   write(*,*)"Reading radar data"
   DO iradar=1,n_radar

    current_radar_file = radar_file 
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
  
 ENDDO ![Endo over model times]

END SUBROUTINE interp_to_real_radar




END MODULE common_scale_to_radar
