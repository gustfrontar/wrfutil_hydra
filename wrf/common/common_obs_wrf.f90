MODULE common_obs_wrf
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!
!=======================================================================
  USE common
  USE common_wrf
  USE common_namelist
  USE map_utils

  IMPLICIT NONE
  PUBLIC

  REAL(r_size) , ALLOCATABLE :: radar_lon(:),radar_lat(:),radar_z(:),radar_mindbz(:)


  !3d observations codes 
  INTEGER , PARAMETER :: id_u_obs=2819
  INTEGER , PARAMETER :: id_v_obs=2820
  INTEGER , PARAMETER :: id_t_obs=3073
  INTEGER , PARAMETER :: id_q_obs=3330
  INTEGER , PARAMETER :: id_rh_obs=3331
  INTEGER , PARAMETER :: id_tv_obs=3079
  !surface observations codes > 9999
  INTEGER , PARAMETER :: id_ps_obs=14593
  INTEGER , PARAMETER :: id_us_obs=82819
  INTEGER , PARAMETER :: id_vs_obs=82820
  INTEGER , PARAMETER :: id_ts_obs=83073
  INTEGER , PARAMETER :: id_qs_obs=83330
  INTEGER , PARAMETER :: id_rhs_obs=83331
  INTEGER , PARAMETER :: id_pwv_obs=83344
  INTEGER , PARAMETER :: id_rain_obs=19999
  INTEGER , PARAMETER :: id_tclon_obs=99991
  INTEGER , PARAMETER :: id_tclat_obs=99992
  INTEGER , PARAMETER :: id_tcmip_obs=99993
  INTEGER , PARAMETER :: id_tcr15_obs=99994
  INTEGER , PARAMETER :: id_tcr25_obs=99995
  !radar observation codes
  INTEGER , PARAMETER :: id_reflectivity_obs=4001
  INTEGER , PARAMETER :: id_radialwind_obs  =4002
  INTEGER , PARAMETER :: id_pseudorh_obs    =4003 
  !chem observation codes
  INTEGER, PARAMETER :: id_totco_obs=5001
  INTEGER, PARAMETER :: id_co_obs = 5002

 
  REAL(r_size) :: minref !=10.0d0**( minrefdbz / 10.0d0)

  !INTEGER, PARAMETER :: type_sband_radar = 1 
  !INTEGER, PARAMETER :: type_xband_radar = 2
  !INTEGER, PARAMETER :: type_cband_radar = 3

CONTAINS

!-----------------------------------------------------------------------
! Define observations groups for variable localization
!-----------------------------------------------------------------------

SUBROUTINE is_vertp( obstype_in , output )
IMPLICIT NONE
integer , intent(in) :: obstype_in
logical , intent(out) :: output
!Classify observation variables as p vertical coordinate (output is true)
!or z vertical coordinate (ouput is false)
  output=.true.
  SELECT CASE ( obstype_in )
    CASE(  id_u_obs , id_v_obs , id_t_obs , id_q_obs , id_rh_obs , id_tv_obs ,      &
           id_pwv_obs, id_rain_obs , id_totco_obs )
         output=.true.
    CASE( id_reflectivity_obs , id_radialwind_obs,         &
          id_pseudorh_obs , id_us_obs , id_vs_obs ,        & 
          id_ts_obs , id_qs_obs , id_rhs_obs, id_ps_obs )
         output=.false.
    CASE DEFAULT
         WRITE(6,*)"[Warning]: Not recognized obs type in is_vertp ",obstype_in
  END SELECT
END SUBROUTINE is_vertp


SUBROUTINE is_3dvar( obstype_in , output )
IMPLICIT NONE
integer , intent(in) :: obstype_in
logical , intent(out) :: output
!Clasify observation variables as 3D variables or 2D (eg surface) variables.
 output=.true.
 SELECT CASE ( obstype_in )
   CASE( id_u_obs , id_v_obs , id_t_obs , id_q_obs , id_rh_obs , id_tv_obs , &
         id_reflectivity_obs , id_radialwind_obs, &
         id_pseudorh_obs )
     output=.true.
   CASE( id_ps_obs , id_us_obs , id_vs_obs , id_ts_obs , id_qs_obs , id_rhs_obs , &
         id_pwv_obs, id_rain_obs , id_totco_obs )
     output=.false.
   CASE default
       WRITE(6,*)"[Warning]: Not recognized obs type in is_3dvar ",obstype_in
 END SELECT

END SUBROUTINE is_3dvar


SUBROUTINE get_iobs( obstype_in , iobs_out  )
IMPLICIT NONE
INTEGER, INTENT(IN) :: obstype_in
INTEGER, INTENT(OUT) :: iobs_out


 SELECT CASE(obstype_in)

      CASE(id_u_obs,id_v_obs,id_us_obs,id_vs_obs)
        iobs_out=1
      CASE(id_t_obs,id_ts_obs)
        iobs_out=2
      CASE(id_tv_obs)
        iobs_out=3
      CASE(id_q_obs,id_rh_obs,id_pseudorh_obs,id_pwv_obs,id_rhs_obs,id_qs_obs)
        iobs_out=4
      CASE(id_ps_obs)
        iobs_out=5
      CASE(id_rain_obs)
        iobs_out=6
      CASE(id_tclon_obs,id_tclat_obs,id_tcmip_obs,id_tcr15_obs,id_tcr25_obs)
        iobs_out=7
      CASE(id_reflectivity_obs)
        iobs_out=8
      CASE(id_radialwind_obs)
        iobs_out=9
      CASE(id_co_obs,id_totco_obs)
        iobs_out=10
      END SELECT

END SUBROUTINE get_iobs

SUBROUTINE get_method_refcalc( lambda , method_ref_calc )
REAL(r_size),INTENT(IN)  :: lambda 
INTEGER     ,INTENT(OUT) :: method_ref_calc


     !Check radar band (radar band is in cm and is stored in typ.
     IF( lambda > 3.75 .AND. lambda <= 7.5 )THEN
       write(6,*)"[Error]: Cband reflectivity calculation is not coded yet"
       method_ref_calc=-9
     ELSEIF( lambda > 7.5 .AND. lambda <= 15)THEN
       method_ref_calc=2
     ELSEIF( lambda > 2.5 .AND. lambda >= 3.75)THEN
       method_ref_calc=3
     ELSE
       write(6,*)"[Error]: Not recognized radar wave-length"
       method_ref_calc=-9
     ENDIF

END SUBROUTINE get_method_refcalc

!-----------------------------------------------------------------------
! Transformation from model variables to an observation
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY(elm,typ,olon,olat,odat,ri,rj,rk,raz,rel,v3d,v2d,yobs)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elm , olon , olat
  REAL(r_size),INTENT(IN) :: ri,rj,rk
  REAL(r_size),INTENT(IN) :: raz,rel    !For radial velocity computation.
  INTEGER                 :: method_ref_calc
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: yobs
  REAL(r_size),INTENT(INOUT)  :: odat  !Odat will be modified in case of simulated obs.
  REAL(r_size),INTENT(INOUT)  :: typ !Observation type, currently used to identify radar.
  REAL(r_size) :: rh(nlon,nlat,nlev)
  REAL(r_size) :: sh(nlon,nlat,nlev)
  REAL(r_size) :: tv(nlon,nlat,nlev)    
  REAL(r_size) :: tg,qg,pg,dz
  REAL(r_size) :: dummy(3)
  INTEGER :: i,j,k
  INTEGER :: is,ie,js,je,ks,ke,ityp
  REAL(r_size) :: qvr,qcr,qrr,qcir,qsr,qgr,ur,vr,wr,tr,pr,rhr
  REAL(r_size) :: dist , dlon , dlat , az , elev , ref , radialv 
  REAL(r_size) :: ti 
  REAL(r_size) :: att_coef !Attenuation factor.


  ie = CEILING( ri )
  is = ie-1
  je = CEILING( rj )
  js = je-1
  ke = CEILING( rk )
  ks = ke-1

  SELECT CASE (NINT(elm))
  CASE(id_u_obs )  ! U
    !Grid can be rotated
    ti=ri+0.5d0
    CALL itpl_3d(v3d(:,:,:,iv3d_u),ti,rj,rk,ur)
    ti=rj+0.5d0
    CALL itpl_3d(v3d(:,:,:,iv3d_v),ri,ti,rk,vr)
    CALL rotwind_letkf(ur,vr,olon,1.0d0,projection)
    yobs=ur

  CASE(id_v_obs)  ! V
    !Grid can be rotated
    ti=ri+0.5d0
    CALL itpl_3d(v3d(:,:,:,iv3d_u),ti,rj,rk,ur)
    ti=rj+0.5d0
    CALL itpl_3d(v3d(:,:,:,iv3d_v),ri,ti,rk,vr)
    CALL rotwind_letkf(ur,vr,olon,1.0d0,projection) 
    yobs=vr
 
  CASE(id_t_obs)  ! T
    CALL itpl_3d(v3d(:,:,:,iv3d_t),ri,rj,rk,yobs)
  CASE(id_q_obs)  ! Q
    sh(is:ie,js:je,ks:ke) = v3d(is:ie,js:je,ks:ke,iv3d_qv) / (1.0d0 + v3d(is:ie,js:je,ks:ke,iv3d_qv))
    CALL itpl_3d(sh,ri,rj,rk,yobs)
  CASE(id_ps_obs) ! PS
    CALL itpl_2d(v2d(:,:,iv2d_t2),ri,rj,tg)
    CALL itpl_2d(v2d(:,:,iv2d_q2),ri,rj,qg)
    CALL itpl_2d(v2d(:,:,iv2d_ps),ri,rj,yobs)
    !write(6,*)'Before press adjustment ',rk,yobs
    CALL prsadj(yobs,rk,tg,qg) !Recall that for Surface Pressure rk is equal to topo - station height
    !write(6,*)'After press adjustment ',rk,yobs
  CASE( id_us_obs ) ! U10
    !Grid can be rotated
    ti=ri+0.5d0
    CALL itpl_2d(v2d(:,:,iv2d_u10),ti,rj,ur)
    ti=rj+0.5d0
    CALL itpl_2d(v2d(:,:,iv2d_v10),ri,ti,vr)
    CALL rotwind_letkf(ur,vr,olon,1.0d0,projection)
    yobs=ur
  CASE(id_vs_obs)  ! V10
    !Grid can be rotated
    ti=ri+0.5d0
    CALL itpl_2d(v2d(:,:,iv2d_u10),ti,rj,ur)
    ti=rj+0.5d0
    CALL itpl_2d(v2d(:,:,iv2d_v10),ri,ti,vr)
    CALL rotwind_letkf(ur,vr,olon,1.0d0,projection)
    yobs=vr
  CASE( id_ts_obs )
    CALL itpl_2d(v2d(:,:,iv2d_t2),ri,rj,yobs)
  CASE( id_rhs_obs ) 
    CALL itpl_2d(v2d(:,:,iv2d_t2),ri,rj,tg)
    CALL itpl_2d(v2d(:,:,iv2d_q2),ri,rj,qg)
    CALL itpl_2d(v2d(:,:,iv2d_ps),ri,rj,pg)
    CALL calc_rh(tg,qg,pg,yobs)
  CASE( id_qs_obs )
    CALL itpl_2d(v2d(:,:,iv2d_q2),ri,rj,yobs)
  CASE(id_rh_obs) ! RH
    DO k=ks,ke
      DO j=js,je
        IF(ie <= nlon ) THEN
          CALL calc_rh(v3d(is,j,k,iv3d_t),v3d(is,j,k,iv3d_qv),&
            & v3d(is,j,k,iv3d_p),rh(is,j,k))
          CALL calc_rh(v3d(ie,j,k,iv3d_t),v3d(ie,j,k,iv3d_qv),&
            & v3d(ie,j,k,iv3d_p),rh(ie,j,k))
        ELSE
          CALL calc_rh(v3d(is,j,k,iv3d_t),v3d(is,j,k,iv3d_qv),&
            & v3d(is,j,k,iv3d_p),rh(is,j,k))
          CALL calc_rh(v3d( 1,j,k,iv3d_t),v3d( 1,j,k,iv3d_qv),&
            & v3d( 1,j,k,iv3d_p),rh( 1,j,k))
        END IF
      END DO
    END DO
    CALL itpl_3d(rh,ri,rj,rk,yobs)
  CASE(id_tv_obs) ! TV
    DO k=ks,ke
      DO j=js,je
        IF(ie <= nlon ) THEN
          CALL calc_tv(v3d(is,j,k,iv3d_t),v3d(is,j,k,iv3d_qv),&
            & tv(is,j,k))
          CALL calc_tv(v3d(ie,j,k,iv3d_t),v3d(ie,j,k,iv3d_qv),&
            & tv(ie,j,k))
        ELSE
          CALL calc_tv(v3d(is,j,k,iv3d_t),v3d(is,j,k,iv3d_qv),&
            & tv(is,j,k))
          CALL calc_tv(v3d( 1,j,k,iv3d_t),v3d( 1,j,k,iv3d_qv),&
            & tv( 1,j,k))
        END IF
      END DO
    END DO
    CALL itpl_3d(tv,ri,rj,rk,yobs)      
  CASE(id_tclon_obs)
    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
    yobs = dummy(1)
  CASE(id_tclat_obs)
    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
    yobs = dummy(2)
  CASE(id_tcmip_obs)
    CALL tctrk(v2d(:,:,iv2d_ps),v2d(:,:,iv2d_t2),ri,rj,dummy)
    yobs = dummy(3)

  CASE(id_reflectivity_obs,id_radialwind_obs)

     CALL itpl_3d(v3d(:,:,:,iv3d_qv),ri,rj,rk,qvr)
     CALL itpl_3d(v3d(:,:,:,iv3d_qc),ri,rj,rk,qcr)
     CALL itpl_3d(v3d(:,:,:,iv3d_qr),ri,rj,rk,qrr)
     CALL itpl_3d(v3d(:,:,:,iv3d_qci),ri,rj,rk,qcir)
     CALL itpl_3d(v3d(:,:,:,iv3d_qs),ri,rj,rk,qsr)
     CALL itpl_3d(v3d(:,:,:,iv3d_qg),ri,rj,rk,qgr)
     ti=ri+0.5d0
     CALL itpl_3d(v3d(:,:,:,iv3d_u),ti,rj,rk,ur)
     ti=rj+0.5d0
     CALL itpl_3d(v3d(:,:,:,iv3d_v),ri,ti,rk,vr)
     CALL rotwind_letkf(ur,vr,olon,1.0d0,projection) !Rotate winds
     ti=rk+0.5d0
     CALL itpl_3d(v3d(:,:,:,iv3d_w),ri,rj,ti,wr)
     CALL itpl_3d(v3d(:,:,:,iv3d_t),ri,rj,rk,tr)
     CALL itpl_3d(v3d(:,:,:,iv3d_p),ri,rj,rk,pr)


     !Get the right reflectivity computation method based on radar
     !wavelength.
     CALL get_method_refcalc(typ , method_ref_calc )

     CALL calc_ref_vr(qvr,qcr,qrr,qcir,qsr,qgr,ur,vr,wr,tr,pr,raz,rel,method_ref_calc,ref,radialv,att_coef)

     IF( NINT(elm) .EQ. id_reflectivity_obs )THEN
         yobs=ref
     ELSEIF( NINT(elm) .EQ. id_radialwind_obs)THEN
         yobs=radialv
     ENDIF

     IF( use_pseudorh .and. NINT(elm) == id_reflectivity_obs .and. yobs == minref )THEN

       DO k=ks,ke
         DO j=js,je
           IF(ie <= nlon ) THEN
             CALL calc_rh(v3d(is,j,k,iv3d_t),v3d(is,j,k,iv3d_qv),&
               & v3d(is,j,k,iv3d_p),rh(is,j,k))
             CALL calc_rh(v3d(ie,j,k,iv3d_t),v3d(ie,j,k,iv3d_qv),&
               & v3d(ie,j,k,iv3d_p),rh(ie,j,k))
            ELSE
             CALL calc_rh(v3d(is,j,k,iv3d_t),v3d(is,j,k,iv3d_qv),&
               & v3d(is,j,k,iv3d_p),rh(is,j,k))
             CALL calc_rh(v3d( 1,j,k,iv3d_t),v3d( 1,j,k,iv3d_qv),&
               & v3d( 1,j,k,iv3d_p),rh( 1,j,k))
           END IF
         END DO
       END DO
       CALL itpl_3d(rh,ri,rj,rk,yobs)

       !Set yobs as a large negative value + rh

       yobs=  yobs - 1000.0d0

     ENDIF

  END SELECT




  RETURN
END SUBROUTINE Trans_XtoY
!-----------------------------------------------------------------------
! TC center search
!  [AUTHORS:] T. Miyoshi and M. Kunii
!-----------------------------------------------------------------------
SUBROUTINE tctrk(ps,t2,ri,rj,trk)
  IMPLICIT NONE
  INTEGER,PARAMETER :: isearch = 10 !search radius [grid points]
  REAL(r_size),INTENT(IN) :: ps(nlon,nlat)
  REAL(r_size),INTENT(IN) :: t2(nlon,nlat)
  REAL(r_size),INTENT(IN) :: ri,rj
  REAL(r_size),INTENT(OUT) :: trk(3) !1:lon, 2:lat, 3:minp
  REAL(r_size) :: wk(100,3)
  REAL(r_size) :: slp(nlon,nlat)
  REAL(r_size) :: p1,p2,p3,p4,p5,a,c,d,xx,yy
  INTEGER :: i,j,i0,i1,j0,j1,n

  i0 = MAX(1,FLOOR(ri)-isearch)
  i1 = MIN(nlon,CEILING(ri)+isearch)
  j0 = MAX(1,FLOOR(rj)-isearch)
  j1 = MIN(nlat,CEILING(rj)+isearch)
  trk = undef

  DO j=j0,j1
    DO i=i0,i1
      slp(i,j) = ps(i,j) * (1.0d0 - 0.0065d0 * phi0(i,j) / &
        & (t2(i,j) + 0.0065d0 * phi0(i,j))) ** (-5.257d0)
    END DO
  END DO

  n=0
  DO j=j0+1,j1-1
    DO i=i0+1,i1-1
      IF(slp(i,j) > slp(i  ,j-1)) CYCLE
      IF(slp(i,j) > slp(i  ,j+1)) CYCLE
      IF(slp(i,j) > slp(i-1,j  )) CYCLE
      IF(slp(i,j) > slp(i+1,j  )) CYCLE
      IF(slp(i,j) > slp(i-1,j-1)) CYCLE
      IF(slp(i,j) > slp(i+1,j-1)) CYCLE
      IF(slp(i,j) > slp(i-1,j+1)) CYCLE
      IF(slp(i,j) > slp(i+1,j+1)) CYCLE
      p1 = slp(i,j)
      p2 = slp(i-1,j)
      p3 = slp(i+1,j)
      p4 = slp(i,j-1)
      p5 = slp(i,j+1)
      c = (p3-p2)*0.5d0
      d = (p5-p4)*0.5d0
      a = (p2+p3+p4+p5)*0.25d0 - p1
      IF(a == 0.0d0) CYCLE
      xx = -0.5d0 * c / a
      yy = -0.5d0 * d / a
      n = n+1
      wk(n,3) = p1 - a*(xx*xx + yy*yy)
      wk(n,2) = lat(i,j) * (1.0d0 - yy) + lat(i,j+1) * yy
      wk(n,1) = lon(i,j) * (1.0d0 - xx) + lon(i+1,j) * xx
    END DO
  END DO

  j=1
  IF(n > 1) THEN
    a = wk(1,3)
    DO i=2,n
      IF(wk(i,3) < a) THEN
        a = wk(i,3)
        j = i
      END IF
    END DO
  END IF
  trk = wk(j,:)

END SUBROUTINE tctrk
!-----------------------------------------------------------------------
! Compute relative humidity (RH)
!-----------------------------------------------------------------------
SUBROUTINE calc_rh(t,q,p,rh)
  IMPLICIT NONE
  REAL(r_size),PARAMETER :: t0=273.15d0
  REAL(r_size),PARAMETER :: e0c=6.11d0
  REAL(r_size),PARAMETER :: al=17.3d0
  REAL(r_size),PARAMETER :: bl=237.3d0
  REAL(r_size),PARAMETER :: e0i=6.1121d0
  REAL(r_size),PARAMETER :: ai=22.587d0
  REAL(r_size),PARAMETER :: bi=273.86d0
  REAL(r_size),INTENT(IN) :: t,q,p
  REAL(r_size),INTENT(OUT) :: rh
  REAL(r_size) :: e,es,tc

  e = q * p * 0.01d0 / (0.378d0 * q + 0.622d0)

  tc = t-t0
  IF(tc >= 0.0d0) THEN
    es = e0c * exp(al*tc/(bl+tc))
  ELSE IF(tc <= -15.d0) THEN
    es = e0i * exp(ai*tc/(bi+tc))
  ELSE
    es = e0c * exp(al*tc/(bl+tc)) * (15.0d0+tc)/15.0d0 &
       + e0i * exp(ai*tc/(bi+tc)) * (-tc) / 15.0d0
  END IF

  rh = e/es

  RETURN
END SUBROUTINE calc_rh
!-----------------------------------------------------------------------
! Compute virtual temperature (TV)
!-----------------------------------------------------------------------
SUBROUTINE calc_tv(t,q,tv)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: t,q
  REAL(r_size),INTENT(OUT) :: tv
  REAL(r_size) :: sh
  
  sh = q / (1.0d0 + q)
  tv = t * (1.0d0 + 0.608d0 * sh)

  RETURN
END SUBROUTINE calc_tv

!-----------------------------------------------------------------------
! Compute radar reflectivity and radial wind.
! Radial wind computations for certain methods depend on model reflectivity
! so both functions has been merged into a single one.
! First reflectivity is computed, and the the radial velocity is computed.
!-----------------------------------------------------------------------
SUBROUTINE calc_ref_vr(qv,qc,qr,qci,qs,qg,u,v,w,t,p,az,elev,method,ref,vr,k)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: qv !Water vapor
  REAL(r_size), INTENT(IN) :: qc,qr  !Cloud and rain water
  REAL(r_size), INTENT(IN) :: qci,qs,qg !Cloud ice, snow and graupel
  REAL(r_size), INTENT(IN) :: t,p    !Temperature and pressure.
  REAL(r_size), INTENT(IN) :: u,v,w  !velocities with respecto to earth.
  REAL(r_size), INTENT(OUT) :: ref   !Reflectivity
  REAL(r_size), INTENT(OUT) :: k            !Attenuation factor.
  REAL(r_size)              :: ro 
  REAL(r_size), INTENT(IN) :: az     !Azimuth respect to the radar.
  REAL(r_size), INTENT(IN) :: elev   !Elevation angle respect to the surface.
  INTEGER     , INTENT(IN) :: method !Method used to compute radial velocity.
  REAL(r_size), INTENT(OUT) :: vr    !Radial velocity.
  REAL(r_size)  :: qms , qmg !Melting species concentration (method 3)
  REAL(r_size)  :: qt        !Total condensate mixing ratio (method 1)
  REAL(r_size)  :: zr , zs , zg !Rain, snow and graupel's reflectivities.
  REAL(r_size)  :: wr , ws , wg !Rain, snow and graupel's mean terminal velocities.
  REAL(r_size)  :: wt           !Total mean terminal velocity.
  REAL(r_size)  :: nor, nos, nog !Rain, snow and graupel's intercepting parameters.
  REAL(r_size)  :: ror, ros, rog , roi !Rain, snow and graupel, ice densities.
  REAL(r_size)  :: a,b,c,d,Cd    !Constant for fall speed computations.
  REAL(r_size),SAVE  :: cf,cf2,cf3,cf4, pip , roo
  REAL(r_size)  :: ki2 , kr2
  REAL(r_size)  :: lr , ls , lg
  REAL(r_size)  :: tmp_factor , rofactor
  REAL(r_size)  :: p0
  REAL(r_size)  :: Fs, Fg , zms , zmg , fws , fwg !Method 3
  REAL(r_size)  :: qrp , qsp , qgp
  REAL(r_size)  :: kr , ks , kg , kms , kmg !Partial contribution of different species to
                                            !the attenuation of the radar beam.
  REAL(r_size)  :: maxf                     !Maximum mixture relative concentration. (method 3)

  LOGICAL :: INITIALIZED=.FALSE.
  !Note: While equivalent reflectivity is assumed to be independent of the radar, in 
  !practice short wavelengths as those associated with K band radars migh frequently
  !experience Mie scattering. In that case, the equivalent reflectivity is not longer
  !radar independent and an appropiate relationship between the forecasted concentrations
  !and the reflectivity should be used.
  
  
  !REAL(r_size)  :: trqr,trqs,trqg
  ![P] Pa
  ![T] K
  ![q...] dimensionless
  ![ref] mm^6/m^3
  ![U,V,W] m/s
  ![az, elev] degree
  ![vr] m/s
  
  !This model reflectivity won't be lower than this value.

  !Some options do not include a model for beam attenuation.
  !In that case a 0 value for the attenuation coefficient is returned.
  k=0.0d0

  !Initialize reflectivities
  zr=0.0d0
  zs=0.0d0
  zg=0.0d0
  zms=0.0d0
  zmg=0.0d0
  ref=0.0d0
  k=0.0d0
  kr=0.0d0
  ks=0.0d0
  kg=0.0d0
  kmg=0.0d0
  kms=0.0d0

  !Compute air density (all methods use this)

  ro =  p / (rd * t)

  !Begin computation of reflectivity and vr

  IF( method .EQ. 1 )THEN

  !WRF method: See for example Sugimoto et al. Evaluation of the Performance of Ra
  !dial Velocity Assimilation with WRF-3DVAR System and Simulated Multiple-Doppler
  !Radar Data
  !Only rain is used to estimate the terminal velocity of hidrometeors.
  !Only rain is used to compute equivalent reflectivity.
  !Marshall-Palmer distributions are assumed to find the relationship bestween
  !concentration and equivalent reflectivity.
  ! Sun and Crook 1997 , 1998.
  !Derived for C-band radars.

  !Reflectivity
  nor=8.0d6      ![m^-4]
  ror=1000.0d0   ![Kg/m3]
  pip=pi ** 1.75 !factor
  cf =1.0d18 * 720 !factor
  p0=1.0d5            !Reference pressure.

  qt=qr + qs + qg  !Assume that the total condensate is rain water
                   !But ignore cloud ice and cloud water

  IF( qt .GT. 0.0d0 )THEN
  ref = cf * ( ( ro * qt )**1.75 )
  ref = ref / ( pip * ( nor ** 0.75 ) * ( ror ** 1.75 ) )
  !ref= 2.04d4 *( ( ro * qt * 1.0d3 ) ** 1.75 ) !Original Sun and Crook expresion.
  ELSE
  ref=minref
  ENDIF

  !Radial wind

  IF ( qt .GT. 0.0d0 )THEN
  a=(p0/p)**0.4
  wt = 5.40d0 * a * ( qt ** 0.125 )
  ELSE
  wt=0d0
  ENDIF
  !WRITE(6,*)qr,qs,qg
  

  ELSEIF( method .EQ. 2)THEN
  !Observation operator from Tong and Xue 2006, 2008 a and b.
  !Based on Smith et al 1975.
  !Ensemble Kalman Filter Assimilation of Doppler Radar Data with a Compressible
  !Nonhydrostatic Model: OSS Experiments. MWR. 133, 1789-187.
  !It includes reflectivity contribution by all the microphisical species.
  !is assumes Marshall and Palmer distributions.
  !Based on S band radars.
    nor=8.0d6      ![m^-4]
    nos=2.0d6      ![m^-4] This value has been modified according to WRF WSM6
    nog=4.0d6      ![m^-4] This value has been modified according to WRF WSM6
    ror=1000.0d0   ![Kg/m3]
    ros=100.0d0    ![Kg/m3]
    rog=913.0d0    ![Kg/m3] 
    roi=917.0d0    ![Kg/m3]
    roo=1.0d0      ![Kg/m3] Surface air density.
    ki2=0.176d0    !Dielectric factor for ice.
    kr2=0.930d0    !Dielectric factor for water.
   IF( .NOT. INITIALIZED)THEN
    pip=pi ** 1.75 !factor
    cf=1.0d18*720/(pip*(nor**0.75d0)*(ror**1.75d0))
    cf2=1.0d18*720*ki2*( ros ** 0.25 )/(pip*kr2*(nos ** 0.75)*( roi ** 2 ) )
    cf3=1.0d18*720/( pip * ( nos ** 0.75 ) * ( ros ** 1.75 ) )  
    cf4=(1.0d18*720/( pip * ( nog ** 0.75) * ( rog ** 1.75 ) ) ) ** 0.95
   ENDIF

    zr=0.0d0
    zs=0.0d0
    zg=0.0d0
    IF( qr .GT. 0.0d0 )THEN
    !rain contribution
    zr= cf * ( ( ro * qr )**1.75 )
    ENDIF
    IF( qs .GT. 0.0d0 )THEN
    IF ( t <= 273.16 )THEN
     !Dry snow
     zs = cf2*( ( ro * qs ) ** 1.75 )
    ELSE
     !Wet snow
     zs = cf3 * ( ( ro * qs ) ** 1.75 )
    ENDIF
    ENDIF

    !Only wet graupel contribution is ussed.
    IF( qg .GT. 0.0d0 )THEN
    zg= cf4 * ( ( ro * qg ) ** 1.6625 )
    
    ENDIF

    ref = zr + zs + zg  

    !Compute reflectivity weigthed terminal velocity.
    !Lin et al 1983.
    IF( ref > minref )THEN
    !There are hidrometeors, compute their terminal velocity.
    !Change units to be consistent with Lin et al 1983 and
    !to obtain wt in m/s
    nor=nor*1e-3      ![cm^-4]
    nos=nos*1e-3      ![cm^-4]
    nog=nog*1e-3      ![cm^-4]
    ror=ror*1e-3        ![g/cm3]
    ros=ros*1e-3        ![g/cm3]
    rog=rog*1e-3      ![g/cm3] 
    roo=roo*1e-3      ![g/cm3] Surface air density.
    ro= ro*1e-3

    a=2115d0   ![cm**1-b / s]
    b=0.8d0
    c=152.93d0 ![cm**1-b / s]
    d=0.25d0
    Cd=0.6d0

    lr= ( pi * ror * nor / ( ro * qr ) ) ** 0.25
    ls= ( pi * ros * nos / ( ro * qs ) ) ** 0.25
    lg= ( pi * rog * nog / ( ro * qg ) ) ** 0.25
   
    rofactor= ( roo / ro  ) ** 0.25 
    wr=0.0d0
    ws=0.0d0
    wg=0.0d0

    IF( qr .GT. 0.0d0)THEN
    CALL com_gamma( 4.0d0 + b , tmp_factor )
    wr= a * tmp_factor / ( 6.0d0 * ( lr ** b ) )
    wr= 1.0d-2*wr * rofactor
    ENDIF
    IF( qs .GT. 0.0d0)THEN
    CALL com_gamma( 4.0d0 + d , tmp_factor ) 
    ws= c * tmp_factor / ( 6.0d0 * ( ls ** d ) )
    ws= 1.0d-2*ws * rofactor
    ENDIF
    IF( qg .GT. 0.0d0 )THEN
    CALL com_gamma( 4.5d0 , tmp_factor )
    wg= tmp_factor * ( ( ( 4.0d0 * gg * rog )/( 3.0d0 * Cd * ro ) ) ** 0.5 )
    wg= 1.0d-2*wg / ( 6.0d0 * ( lg ** 0.5 ) )
    ENDIF
    !Reflectivity weighted terminal velocity. 
    wt = ( wr * zr + ws * zs + wg * zg )/ ( zr + zs + zg )
    ELSE

    wt=0.0d0

    ENDIF   

 ELSEIF( method .EQ. 3) THEN
 !Observation operator from Xue et al 2007
 !Asumes power law between Z and mass concentration of different 
 !hydrometeor categories. 
 !Derived for X-Band radars tacking into account Mie scattering.
 !Includes a computation of the attenuation. (k)

 MAXF=0.5d0

  !First we need to compute the mixtures between rain, snow and hail
  !Following Jung et al 2007 eq (2) and (3)
  Fg=0.0d0
  Fs=0.0d0
  fwg=0.0d0
  fws=0.0d0
  IF( qr .GT. 0.0d0 .AND. qg .GT. 0.0d0)THEN
     Fg=MAXF * ( min( qr/qg , qg/qr ) )**(1.0d0/3.0d0)
     fwg= qr / ( qr + qg )
  ENDIF
  IF( qr .GT. 0.0d0 .AND. qs .GT. 0.0d0)THEN
     Fs=MAXF * ( min( qr/qs , qs/qr ) )**(1.0d0/3.0d0)
     fws= qr / ( qr + qs )
  ENDIF


  !Correct the rain, snow and hail mixing ratios assuming
  !that we have a mixture due to melting.

  qrp=(1.0d0-Fs-Fg)*qr
  
  qsp=(1.0d0-Fs)*qs

  qgp=(1.0d0-Fg)*qg

  !Compute the new species resulting from melting.

  qms=Fs * (qr + qs) !Melting snow concentration.

  qmg=Fg * (qr + qg) !Melting hail concentration.

  !Compute reflectivities for each species including the melting species.

  IF( qrp .GT. 0.0d0)THEN
  zr= 2.53d4 * ( ro * qrp * 1.0d3 )**1.84
  kr= 0.319 * ( ro * qrp * 1.0d3 )**1.38
  ENDIF
  IF( qsp .GT. 0.0d0)THEN
  zs= 3.48d3 * ( ro * qsp * 1.0d3 )**1.66
  ks=0.00483 * ( ro * qsp * 1.0d3 )**1.28
  ENDIF
  IF( qgp .GT. 0.0d0)THEN
  zg= 8.18d4 * ( ro * qgp * 1.0d3 )**1.50
  kg=0.159 * ( ro * qgp * 1.0d3 )**1.64
  ENDIF
  IF( qms .GT. 0.0d0 )THEN
  zms=( 0.00491 + 5.75*fws - 5.588*(fws**2) )*1.0d5
  zms= zms * ( ro * qms * 1.0d3 )**( 1.67 - 0.202*fws + 0.398*(fws**2) )
  kms= 0.0413 + 22.7*fws -50.5*(fws**2) + 28.6*(fws**3)
  kms= kms * ( ro * qms * 1.0d3 )**( 1.06-0.579*fws+2.03*(fws**2) -1.24*(fws**3) )

  ENDIF
  IF( qmg .GT. 0.0d0 )THEN
  zmg=( 0.809 + 10.13*fwg -5.98*(fwg**2) )*1.0d5
  zmg= zmg * ( ro * qmg * 1.0d3 )**( 1.48 + 0.0448*fwg - 0.0313*(fwg**2) )
  kmg= 0.256 + 6.28*fwg -11.36*(fwg**2) +6.01*(fwg**3)
  kmg= kmg * ( ro * qmg * 1.0d3 )**(1.26 -0.659*fwg +1.44*(fwg**2)-0.817*(fwg**3) )
  ENDIF
 
  ref = zr +  zg  + zs + zms + zmg

  k   = kr +  kg  + ks + kms + kmg

  !Compute reflectivity weigthed terminal velocity.
  !Lin et al 1983. (The distribution parameters are 
  !consistent with the work of Jung et al 2007)

  IF( ref > minref )THEN
    !There are hidrometeors, compute their terminal velocity.
    !Units according to Lin et al 1983.
    nor=8.0d-2      ![cm^-4]
    nos=3.0d-2      ![cm^-4]
    nog=4.0d-4      ![cm^-4]
    ror=1.0d0        ![g/cm3]
    ros=0.1d0        ![g/cm3]
    rog=0.917d0      ![g/cm3] 
    roo=0.001d0      ![g/cm3] Surface air density.
    ro=1.0d-3 * ro
    a=2115d0   ![cm**1-b / s]
    b=0.8d0
    c=152.93d0 ![cm**1-b / s]
    d=0.25d0
    Cd=0.6d0
  
    rofactor= ( roo / ro  ) ** 0.5

    IF ( qr .GT. 0.0d0 )THEN 
    CALL com_gamma( 4.0d0 + b , tmp_factor )
    lr= ( pi * ror * nor / ( ro * qr ) ) ** 0.25
    wr= a * tmp_factor / ( 6.0d0 * ( lr ** b ) )
    wr= 1.0d-2 * wr * rofactor
    ELSE
    wr=0.0d0
    ENDIF

    IF( qs .GT. 0.0d0 )THEN
    ls= ( pi * ros * nos / ( ro * qs ) ) ** 0.25
    CALL com_gamma( 4.0d0 + d , tmp_factor )
    ws= c * tmp_factor / ( 6.0d0 * ( ls ** d ) )
    ws= 1.0d-2 * ws * rofactor
    ELSE
    ws=0.0d0
    ENDIF

    IF ( qg .GT. 0.0d0 )THEN
    lg= ( pi * rog * nog / ( ro * qg ) ) ** 0.25
    CALL com_gamma( 4.5d0 , tmp_factor )
    wg= tmp_factor * ( ( ( 4.0d0 * gg * rog )/( 3.0d0 * Cd * ro ) ) ** 0.5 )
    wg= 1.0d-2 * wg / ( 6.0d0 * ( lg ** 0.5 ) )
    ELSE
    wg=0.0d0
    ENDIF

    !Reflectivity weighted terminal velocity. 
    !The melting species are assumed to fail as their non-melting counterparts.
    !however this might not be the case for melting snow.
    wt = ( wr * zr + ws *  zs +  ws * zms + wg *  zg + wg * zmg ) / ( zr + zs + zg + zms + zmg )
  
 ELSE

    wt=0.0d0

  ENDIF


 ELSE  !IF OVER DIFFERENT OPTIONS

    WRITE(6,*)'ERROR: Not recognized method for radar reflectivity and wind computation'
    
    ref=undef
    vr =undef

 ENDIF !END IF OVER DIFFERENT COMPUTATION OPTIONS.

 IF( ref /= undef )THEN
 !Compute radial velocity
    vr = u * cos(elev*deg2rad) * sin(az*deg2rad)
    vr = vr + v * cos(elev*deg2rad) * cos(az*deg2rad)
    IF( USE_WT )THEN
      vr = vr + (w - wt)*sin(elev*deg2rad)
    ELSE
      vr = vr + w*sin(elev*deg2rad)  
    ENDIF
 ELSE
    vr = undef 
 ENDIF


 INITIALIZED=.TRUE.
RETURN
END SUBROUTINE calc_ref_vr



SUBROUTINE aer2llz(radarlon,radarlat,radarz,az,el,ra,lon,lat,z)
IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: radarlon,radarlat,radarz,az,ra
  REAL(r_size),INTENT(OUT) :: lon,lat,z
  REAL(r_size),INTENT(INOUT) :: el  !Elevation will be corrected in this routine.
  REAL(r_size)             :: tmp1 , tmp2 , ke
  REAL(r_size)             :: distance_to_radar
  !This routine computes lon,lat,z of a radar beam at a certain
  !azimuth,elevation and range with respect to the radar location.
  !This version uses beam propagation assuming the temperature and
  !moisture profile of a standard atmosphere.
  !Future versions can use the background vertical profile of t and moisture
  !to produce a better representation of the beam propagation.

  ke=(4d0/3d0)

  !Perform standard height beam heigth computation.

      tmp1= ra**2 + (ke*Re)**2 + 2*ra*ke*Re*sin( el*deg2rad )
      z=radarz + SQRT( tmp1 )-ke*Re;
      !Compute the local elevation angle tacking into account the efective earth radius.
      tmp1= ra * cos( el * deg2rad )
      tmp2= ra * sin( el * deg2rad ) + ke * Re
      !Correct local elevation angle with respect to the sourface tacking into account
      !beam propagation.
      el = el + atan(tmp1/tmp2)

  !Compute the latitude and longitude corresponding to each radar point.
  distance_to_radar=ke*Re*asin( ra * cos(el * deg2rad) / (ke * Re) )

  CALL com_ll_arc_distance(radarlon,radarlat,distance_to_radar,az,lon,lat)

END SUBROUTINE aer2llz

!-----------------------------------------------------------------------
! Pressure adjustment for a different height level
!-----------------------------------------------------------------------
SUBROUTINE prsadj(p,dz,t,q)
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: p
  REAL(r_size),INTENT(IN) :: dz ! height difference (target - original) [m]
  REAL(r_size),INTENT(IN) :: t  ! temperature [K] at original level
  REAL(r_size),INTENT(IN) :: q  ! humidity [kg/kg] at original level
  REAL(r_size),PARAMETER :: gamma=5.0d-3 ! lapse rate [K/m]
  REAL(r_size) :: tv

  IF(dz /= 0) THEN
    tv = t * (1.0d0 + 0.608d0 * q)
    p = p * ((-gamma*dz+tv)/tv)**(gg/(gamma*rd)) !tv is at original level
!    p = p * (tv/(tv+gamma*dz))**(gg/(gamma*rd)) !tv is at target level

  END IF

  RETURN
END SUBROUTINE prsadj
!-----------------------------------------------------------------------
! Horizontal coordinate conversion
!-----------------------------------------------------------------------
!SUBROUTINE ll2ij(elem,rlon,rlat,ri,rj)
!  IMPLICIT NONE
!  REAL(r_size),INTENT(IN) :: elem
!  REAL(r_size),INTENT(IN) :: rlon(1)
!  REAL(r_size),INTENT(IN) :: rlat(1)
!  REAL(r_size),INTENT(OUT) :: ri(1)
!  REAL(r_size),INTENT(OUT) :: rj(1)
!  INTEGER(4),PARAMETER :: isw = 3  
!  INTEGER :: i,j
  !
  ! rlon, rlat --> ri, rj
  !
  !SELECT CASE(NINT(elem))
  !CASE(id_u_obs)
  !  i=nlon
  !  j=nlat-1
  !  CALL com_pos2ij(isw,i,j,lonu(1:i,1:j),latu(1:i,1:j),1,rlon,rlat,ri,rj)
  !CASE(id_v_obs)
  !  i=nlon-1
  !  j=nlat
  !  CALL com_pos2ij(isw,i,j,lonv(1:i,1:j),latv(1:i,1:j),1,rlon,rlat,ri,rj)
  !CASE DEFAULT
!    i=nlon-1
!    j=nlat-1
!    CALL com_pos2ij(isw,i,j,lon(1:i,1:j),lat(1:i,1:j),1,rlon,rlat,ri,rj)
!  !END SELECT
!
!  RETURN
!END SUBROUTINE ll2ij
!-----------------------------------------------------------------------
! Vertical coordinate conversion
!-----------------------------------------------------------------------
SUBROUTINE p2k(p_full,elem,ri,rj,rlev,rk)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: p_full(nlon,nlat,nlev)
  REAL(r_size),INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rlev ! pressure levels
  REAL(r_size),INTENT(OUT) :: rk
  REAL(r_size) :: ak
  REAL(r_size) :: lnps(nlon,nlat)
  REAL(r_size) :: plev(nlev)
  INTEGER :: i,j,k,n
  !
  ! rlev --> rk
  !
  IF(NINT(elem) > 9999) THEN ! surface observation
    rk = 0.0d0
  ELSE
    !
    ! horizontal interpolation
    !
    i = CEILING(ri)
    j = CEILING(rj)
    IF(p_full(i,j,1) < 0.1) THEN
      rk = -1.0d0
      RETURN
    END IF
    DO k=1,nlev-1
      lnps(i-1:i,j-1:j) = LOG(p_full(i-1:i,j-1:j,k))
      CALL itpl_2d(lnps,ri,rj,plev(k))
    END DO
    !
    ! Log pressure
    !
    rk = LOG(rlev)
    !
    ! find rk
    !
    DO k=2,nlev-1
      IF(plev(k) < rk) EXIT ! assuming descending order of plev
    END DO
    ak = (rk - plev(k-1)) / (plev(k) - plev(k-1))
    rk = REAL(k-1,r_size) + ak
    !WRITE(*,*)EXP(plev(k-1)),EXP(plev(k))
  END IF

  RETURN
END SUBROUTINE p2k

!-----------------------------------------------------------------------
! Fast vertical coordinate conversion (Z-level)
! Instead of perform horizontal interpolation take the 
! nearest model grid profile of z. (Takes advantange of 
! the high resolution of the model)
!-----------------------------------------------------------------------
SUBROUTINE z2k_fast(z_full,ri,rj,rlev,rk)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: z_full(nlon,nlat,nlev)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rlev 
  REAL(r_size),INTENT(OUT):: rk
  REAL(r_size) :: ak
  !REAL(r_size) :: z2d(nlon,nlat)
  REAL(r_size) :: zlev(nlev-1)
  INTEGER :: i,j,k,n,zindex
  rk=0.0d0
  !
  ! rlev --> rk
  !
    !
    ! horizontal interpolation
    !
    !i = CEILING(ri)
    !j = CEILING(rj)
    !IF(z_full(i,j,2) < phi0(i,j)) THEN
    !  rk = -1.0d0
    !  RETURN
    !END IF
    !DO k=1,nlev
    !  z2d(i-1:i,j-1:j) = z_full(i-1:i,j-1:j,k)
    !  CALL itpl_2d(z2d,ri,rj,zlev(k))
    !END DO
    i=INT(ri)
    j=INT(rj)
    zlev(1:nlev-1)=z_full(i,j,1:nlev-1)
    !Find rk
    DO k=1,nlev-2
      zindex=k
      IF(zlev(k) > rlev) EXIT ! assuming increasing order of zlev
    END DO

    IF( zindex > 1)THEN
      ak = (rlev - zlev(zindex-1)) / (zlev(zindex) - zlev(zindex-1))
      rk = REAL(zindex-1,r_size) + ak
    ELSE
      rk = -1.0d0 !Data is below topography.
    ENDIF

  RETURN
END SUBROUTINE z2k_fast

!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
SUBROUTINE itpl_2d(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlon,nlat)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  IF(i <= nlon) THEN
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(i  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(i  ,j  ) *    ai  *    aj
  ELSE
    var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
       & + var(1  ,j-1) *    ai  * (1-aj) &
       & + var(i-1,j  ) * (1-ai) *    aj  &
       & + var(1  ,j  ) *    ai  *    aj
  END IF

  RETURN
END SUBROUTINE itpl_2d

SUBROUTINE itpl_3d(var,ri,rj,rk,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlon,nlat,nlev)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rk
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj,ak
  INTEGER :: i,j,k

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)
  k = CEILING(rk)
  ak = rk - REAL(k-1,r_size)

  IF(i <= nlon) THEN
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(i  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(i  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(i  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(i  ,j  ,k  ) *    ai  *    aj  *    ak
  ELSE
    var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
       & + var(1  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
       & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
       & + var(1  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
       & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
       & + var(1  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
       & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
       & + var(1  ,j  ,k  ) *    ai  *    aj  *    ak
  END IF

  RETURN
END SUBROUTINE itpl_3d

!-----------------------------------------------------------------------
! Monitor departure
!-----------------------------------------------------------------------
SUBROUTINE monit_dep(nn,elm,dep,sprd)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elm(nn)
  REAL(r_size),INTENT(IN) :: dep(nn)
  REAL(r_size),INTENT(IN) :: sprd(nn)
  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_q,rmse_ps,rmse_rh,rmse_pwv,rmse_ref,rmse_vr
  REAL(r_size) :: sprd_u,sprd_v,sprd_t,sprd_q,sprd_ps,sprd_rh,sprd_pwv,sprd_ref,sprd_vr
  REAL(r_size) :: bias_u,bias_v,bias_t,bias_q,bias_ps,bias_rh,bias_pwv,bias_ref,bias_vr
  INTEGER :: n,iu,iv,it,iq,ips,irh,ipwv,iref,ivr

  rmse_u = 0.0d0
  rmse_v = 0.0d0
  rmse_t = 0.0d0
  rmse_q = 0.0d0
  rmse_ps = 0.0d0
  rmse_rh = 0.0d0
  rmse_pwv = 0.0d0 
  rmse_ref = 0.0d0
  rmse_vr  = 0.0d0

  sprd_u = 0.0d0
  sprd_v = 0.0d0
  sprd_t = 0.0d0
  sprd_q = 0.0d0
  sprd_ps = 0.0d0
  sprd_rh = 0.0d0
  sprd_pwv = 0.0d0 
  sprd_ref = 0.0d0
  sprd_vr  = 0.0d0

  bias_u = 0.0d0
  bias_v = 0.0d0
  bias_t = 0.0d0
  bias_q = 0.0d0
  bias_ps = 0.0d0
  bias_rh = 0.0d0
  bias_pwv = 0.0d0   
  bias_ref = 0.0d0
  bias_vr  = 0.0d0

  iu = 0
  iv = 0
  it = 0
  iq = 0
  ips = 0
  irh = 0
  ipwv = 0  
  iref = 0
  ivr  = 0
  
  DO n=1,nn
    IF( dep(n) == undef ) CYCLE
    SELECT CASE(NINT(elm(n)))
    CASE(id_u_obs, id_us_obs)
      rmse_u = rmse_u + dep(n)**2
      bias_u = bias_u + dep(n)
      sprd_u = sprd_u + sprd(n)**2
      iu = iu + 1
    CASE(id_v_obs, id_vs_obs)
      rmse_v = rmse_v + dep(n)**2
      bias_v = bias_v + dep(n)
      sprd_v = sprd_v + sprd(n)**2
      iv = iv + 1
    CASE(id_t_obs, id_ts_obs , id_tv_obs )
      rmse_t = rmse_t + dep(n)**2
      bias_t = bias_t + dep(n)
      sprd_t = sprd_t + sprd(n)**2
      it = it + 1
    CASE(id_q_obs, id_qs_obs)
      rmse_q = rmse_q + dep(n)**2
      bias_q = bias_q + dep(n)
      sprd_q = sprd_q + sprd(n)**2
      iq = iq + 1
    CASE(id_ps_obs)
      rmse_ps = rmse_ps + dep(n)**2
      bias_ps = bias_ps + dep(n)
      sprd_ps = sprd_ps + sprd(n)**2
      ips = ips + 1
    CASE(id_rh_obs, id_rhs_obs , id_pseudorh_obs )
      rmse_rh = rmse_rh + dep(n)**2
      bias_rh = bias_rh + dep(n)
      sprd_rh = sprd_rh + sprd(n)**2
      irh = irh + 1
    CASE(id_pwv_obs)
      rmse_pwv = rmse_pwv + dep(n)**2
      bias_pwv = bias_pwv + dep(n)
      sprd_pwv = sprd_pwv + sprd(n)**2
      ipwv = ipwv + 1    
     CASE(id_reflectivity_obs )
      rmse_ref = rmse_ref + dep(n)**2 
      bias_ref = bias_ref + dep(n)
      sprd_ref = sprd_ref + sprd(n)**2
      iref = iref + 1
     CASE(id_radialwind_obs)
      rmse_vr = rmse_vr + dep(n)**2
      bias_vr = bias_vr + dep(n)
      sprd_vr = sprd_vr + sprd(n)**2
      ivr = ivr + 1
     CASE DEFAULT
      WRITE(6,*) 'Observation type is not in the list'
      WRITE(6,*)NINT(elm(n)),dep(n),n  
    END SELECT
  END DO
  IF(iu == 0) THEN
    rmse_u = undef
    bias_u = undef
    sprd_u = undef
  ELSE
    rmse_u = SQRT(rmse_u / REAL(iu,r_size))
    bias_u = bias_u / REAL(iu,r_size)
    sprd_u = sqrt(sprd_u / REAL(iu,r_size))
  END IF
  IF(iv == 0) THEN
    rmse_v = undef
    bias_v = undef
    sprd_v = undef
  ELSE
    rmse_v = SQRT(rmse_v / REAL(iv,r_size))
    bias_v = bias_v / REAL(iv,r_size)
    sprd_v = sqrt(sprd_v / REAL(iv,r_size))
  END IF
  IF(it == 0) THEN
    rmse_t = undef
    bias_t = undef
    sprd_v = undef
  ELSE
    rmse_t = SQRT(rmse_t / REAL(it,r_size))
    bias_t = bias_t / REAL(it,r_size)
    sprd_t = sqrt(sprd_t / REAL(it,r_size))
  END IF
  IF(iq == 0) THEN
    rmse_q = undef
    bias_q = undef
    sprd_q = undef 
  ELSE
    rmse_q = SQRT(rmse_q / REAL(iq,r_size))
    bias_q = bias_q / REAL(iq,r_size)
    sprd_q = sqrt(sprd_q / REAL(iq,r_size))
  END IF
  IF(ips == 0) THEN
    rmse_ps = undef
    bias_ps = undef
    sprd_ps = undef 
  ELSE
    rmse_ps = SQRT(rmse_ps / REAL(ips,r_size))
    bias_ps = bias_ps / REAL(ips,r_size)
    sprd_ps = sqrt(sprd_ps / REAL(ips,r_size))
  END IF
  IF(irh == 0) THEN
    rmse_rh = undef
    bias_rh = undef
    sprd_rh = undef 
  ELSE
    rmse_rh = SQRT(rmse_rh / REAL(irh,r_size))
    bias_rh = bias_rh / REAL(irh,r_size)
    sprd_rh = sqrt(sprd_rh / REAL(irh,r_size))
  END IF
  IF(ipwv == 0) THEN
    rmse_pwv = undef
    bias_pwv = undef
    sprd_pwv = undef 
  ELSE
    rmse_pwv = SQRT(rmse_pwv / REAL(ipwv,r_size))
    bias_pwv = bias_pwv / REAL(ipwv,r_size)
    sprd_pwv = sqrt(sprd_pwv / REAL(ipwv,r_size))
  END IF
  IF(iref == 0) THEN
    rmse_ref = undef
    bias_ref = undef
    sprd_ref = undef 
  ELSE
    rmse_ref = SQRT(rmse_ref / REAL(iref,r_size))
    bias_ref = bias_ref / REAL(iref,r_size)
    sprd_ref = sqrt(sprd_ref / REAL(iref,r_size))
  END IF
  IF(ivr == 0) THEN
    rmse_vr = undef
    bias_vr = undef
    sprd_vr = undef 
  ELSE
    rmse_vr = SQRT(rmse_vr / REAL(ivr,r_size))
    bias_vr = bias_vr / REAL(ivr,r_size)
    sprd_vr = sqrt(sprd_vr / REAL(ivr,r_size))
  END IF


  WRITE(6,'(A)') '== OBSERVATIONAL DEPARTURE ========================================================='
  WRITE(6,'(9A12)') 'U','V','T','Q','PS','RH','PWV','Ze','Vr'
  WRITE(6,'(9ES12.3)') bias_u,bias_v,bias_t,bias_q,bias_ps,bias_rh,bias_pwv,bias_ref,bias_vr
  WRITE(6,'(9ES12.3)') rmse_u,rmse_v,rmse_t,rmse_q,rmse_ps,rmse_rh,rmse_pwv,rmse_ref,rmse_vr
  WRITE(6,'(9ES12.3)') sprd_u,sprd_v,sprd_t,sprd_q,sprd_ps,sprd_rh,sprd_pwv,sprd_ref,sprd_vr
  WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS TO BE ASSIMILATED ========================================'
  WRITE(6,'(9A12)') 'U','V','T','Q','PS','RH','PWV','Ze','Vr'
  WRITE(6,'(9I12)') iu,iv,it,iq,ips,irh,ipwv,iref,ivr
  WRITE(6,'(A)') '===================================================================================='

  RETURN
END SUBROUTINE monit_dep

!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs(cfile,nn)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl) :: wk(7)
  INTEGER :: ios
  INTEGER :: iu,iv,it,iq,irh,itv,ips,itc
  INTEGER :: iunit
  LOGICAL :: ex

  nn = 0
  iu = 0
  iv = 0
  it = 0
  iq = 0
  irh = 0
  itv = 0
  ips = 0
  itc = 0
  iunit=91
  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
    DO
      READ(iunit,IOSTAT=ios) wk
      IF(ios /= 0) EXIT
      SELECT CASE(NINT(wk(1)))
      CASE(id_u_obs,id_us_obs)
        iu = iu + 1
      CASE(id_v_obs,id_vs_obs)
        iv = iv + 1
      CASE(id_t_obs,id_ts_obs)
        it = it + 1
      CASE(id_q_obs,id_qs_obs)
        iq = iq + 1
      CASE(id_rh_obs,id_rhs_obs)
        irh = irh + 1
      CASE(id_tv_obs)
        itv = itv + 1      
      CASE(id_ps_obs)
        ips = ips + 1
      CASE(id_tclon_obs)
        itc = itc + 1
      END SELECT
      nn = nn + 1
    END DO
    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
    WRITE(6,'(A12,I10)') '          U:',iu
    WRITE(6,'(A12,I10)') '          V:',iv
    WRITE(6,'(A12,I10)') '          T:',it
    WRITE(6,'(A12,I10)') '          Q:',iq
    WRITE(6,'(A12,I10)') '         RH:',irh
    WRITE(6,'(A12,I10)') '         TV:',itv    
    WRITE(6,'(A12,I10)') '         Ps:',ips
    WRITE(6,'(A12,I10)') '   TC TRACK:',itc
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF

  RETURN
END SUBROUTINE get_nobs

SUBROUTINE read_obs(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,otyp)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(OUT) :: elem(nn) ! element number
  REAL(r_size),INTENT(OUT) :: rlon(nn)
  REAL(r_size),INTENT(OUT) :: rlat(nn)
  REAL(r_size),INTENT(OUT) :: rlev(nn)
  REAL(r_size),INTENT(OUT) :: odat(nn)
  REAL(r_size),INTENT(OUT) :: oerr(nn)
  REAL(r_size),INTENT(OUT) :: otyp(nn)
  REAL(r_sngl) :: wk(7)
  INTEGER :: n,iunit

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,nn
    READ(iunit) wk
    SELECT CASE(NINT(wk(1)))
    CASE(id_u_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_v_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_t_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_q_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_tv_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
    CASE(id_ps_obs)
      wk(5) = wk(5) * 100.0 ! hPa -> Pa
      wk(6) = wk(6) * 100.0 ! hPa -> Pa
    CASE(id_rh_obs)
      wk(4) = wk(4) * 100.0 ! hPa -> Pa
      wk(5) = wk(5) * 0.01 ! percent input
      wk(6) = wk(6) * 0.01 ! percent input
    CASE(id_tcmip_obs)
      wk(5) = wk(5) * 100.0 ! hPa -> Pa
      wk(6) = wk(6) * 100.0 ! hPa -> Pa
    END SELECT
    elem(n) = REAL(wk(1),r_size)
    rlon(n) = REAL(wk(2),r_size)
    rlat(n) = REAL(wk(3),r_size)
    rlev(n) = REAL(wk(4),r_size)
    odat(n) = REAL(wk(5),r_size)
    oerr(n) = REAL(wk(6),r_size)
    otyp(n) = REAL(wk(7),r_size)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs



!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs_radar(cfile,nn,radarlon,radarlat,radarz)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl) :: wk(7),tmp
  INTEGER :: ios
  INTEGER :: ir,iv
  INTEGER :: iunit
  LOGICAL :: ex
  REAL(r_size) radarlon,radarlat,radarz

  nn = 0
  iv = 0
  ir = 0
  iunit=91
  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
    READ(iunit)tmp
    radarlon=REAL(tmp,r_size)
    READ(iunit)tmp
    radarlat=REAL(tmp,r_size)
    READ(iunit)tmp
    radarz=REAL(tmp,r_size)
    DO
      READ(iunit,IOSTAT=ios) wk
      IF(ios /= 0) EXIT
      !WRITE(6,*)wk
      SELECT CASE(NINT(wk(1)))
      CASE(id_reflectivity_obs)
        ir = ir + 1
      CASE(id_radialwind_obs)
        iv = iv + 1
      END SELECT
      nn = nn + 1
    END DO
    WRITE(6,*)' RADAR FILE ', cfile
    WRITE(6,*)' RADAR LON = ',radarlon
    WRITE(6,*)' RADAR LAT = ',radarlat
    WRITE(6,*)' RADAR Z   = ',radarz
    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
    WRITE(6,'(A12,I10)') '   REFLECTIVITY:',ir
    WRITE(6,'(A12,I10)') 'RADIAL VELOCITY:',iv
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF

  RETURN
END SUBROUTINE get_nobs_radar

SUBROUTINE read_obs_radar(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,otyp)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(OUT) :: elem(nn) ! element number
  REAL(r_size),INTENT(OUT) :: rlon(nn)
  REAL(r_size),INTENT(OUT) :: rlat(nn)
  REAL(r_size),INTENT(OUT) :: rlev(nn)
  REAL(r_size),INTENT(OUT) :: odat(nn)
  REAL(r_size),INTENT(OUT) :: oerr(nn)
  REAL(r_size),INTENT(OUT) :: otyp(nn)
  REAL(r_sngl) :: wk(7)
  INTEGER :: n,iunit
  REAL(r_sngl) :: tmp

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  READ(iunit)tmp
  READ(iunit)tmp
  READ(iunit)tmp
  DO n=1,nn
    READ(iunit) wk
    elem(n) = REAL(wk(1),r_size)
    rlon(n) = REAL(wk(2),r_size)
    rlat(n) = REAL(wk(3),r_size)
    rlev(n) = REAL(wk(4),r_size)
    odat(n) = REAL(wk(5),r_size)
    oerr(n) = REAL(wk(6),r_size)
    otyp(n) = REAL(wk(7),r_size)
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs_radar

SUBROUTINE monit_obs(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,otyp,omb,oma,slot)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile
  INTEGER,INTENT(IN) :: nn
  REAL(r_size),INTENT(IN) :: elem(nn) ! element number
  REAL(r_size),INTENT(IN) :: rlon(nn)
  REAL(r_size),INTENT(IN) :: rlat(nn)
  REAL(r_size),INTENT(IN) :: rlev(nn)
  REAL(r_size),INTENT(IN) :: odat(nn)
  REAL(r_size),INTENT(IN) :: oerr(nn)
  REAL(r_size),INTENT(IN) :: otyp(nn)
  REAL(r_size),INTENT(IN) :: omb(nn)
  REAL(r_size),INTENT(IN) :: oma(nn)
  REAL(r_size),INTENT(IN) :: slot(nn)
  REAL(r_sngl) :: wk(10)
  INTEGER :: n,iunit

  iunit=91
  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential')
  DO n=1,nn
    wk(1) = REAL(elem(n),r_sngl)
    wk(2) = REAL(rlon(n),r_sngl)
    wk(3) = REAL(rlat(n),r_sngl)
    wk(4) = REAL(rlev(n),r_sngl)
    wk(5) = REAL(odat(n),r_sngl)
    wk(6) = REAL(oerr(n),r_sngl)
    wk(7) = REAL(otyp(n),r_sngl)
    wk(8) = REAL(omb(n),r_sngl)
    wk(9) = REAL(oma(n),r_sngl)
    wk(10) = REAL(slot(n),r_sngl)
    SELECT CASE(NINT(wk(1)))
    CASE(id_u_obs)
      wk(4) = wk(4) / 100.0 ! Pa -> hPa
    CASE(id_v_obs)
      wk(4) = wk(4) / 100.0 ! Pa -> hPa
    CASE(id_t_obs)
      wk(4) = wk(4) / 100.0 ! Pa -> hPa
    CASE(id_q_obs)
      wk(4) = wk(4) / 100.0 ! Pa -> hPa
    CASE(id_tv_obs)
      wk(4) = wk(4) / 100.0 ! Pa -> hPa
    CASE(id_ps_obs)
      wk(5) = wk(5) / 100.0 ! Pa -> hPa
      wk(6) = wk(6) / 100.0 ! Pa -> hPa
      wk(8) = wk(8) / 100.0 ! Pa -> hPa
      wk(9) = wk(9) / 100.0 ! Pa -> hPa
    CASE(id_rh_obs)
      wk(4) = wk(4) / 100.0 ! Pa -> hPa
      wk(5) = wk(5) / 0.01 ! percent output
      wk(6) = wk(6) / 0.01 ! percent output
      wk(8) = wk(8) / 0.01 ! percent output
      wk(9) = wk(9) / 0.01 ! percent output
    CASE(id_tcmip_obs)
      wk(5) = wk(5) / 100.0 ! Pa -> hPa
      wk(6) = wk(6) / 100.0 ! Pa -> hPa
      wk(8) = wk(8) / 100.0 ! Pa -> hPa
      wk(9) = wk(9) / 100.0 ! Pa -> hPa
    END SELECT
    WRITE(iunit) wk
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE monit_obs





END MODULE common_obs_wrf
