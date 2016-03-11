PROGRAM WRF_TO_RADAR
!=======================================================================
!
! [PURPOSE:] Main program of WRF_TO_RADAR
! This program creates radial velocity and reflectivyt observations at
! each model grid point.
!
! [HISTORY:]
!   02/19/2014 Juan Ruiz created
!=======================================================================

USE common
USE common_wrf
!USE common_wrf_to_radar
!USE common_radar_tools

 IMPLICIT NONE

 REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:)
 REAL(r_size),ALLOCATABLE :: gues2d(:,:,:)

! REAL(r_size),ALLOCATABLE :: radar3d(:,:,:,:)
 REAL(r_size) , ALLOCATABLE :: z3d(:,:,:) , dv3d(:,:,:)
 REAL(r_size)             :: rtimer00 , rtimer
 INTEGER                  :: i , ii , jj , kk

 CHARACTER*20  :: model_file_name='input_model.nc'

  REAL(r_size)         :: ri , rj , rk , rlev
  REAL(r_size)         :: tmp_z , att , tmp
  REAL(r_size)         :: p , t
  REAL(r_size)         :: qv  , qc , qr
  REAL(r_size)         :: qci , qs , qg
  REAL(r_size)         :: u   , v  , w
  REAL(r_size)         :: az , elev , dlon , dlat , dist
  INTEGER              :: NOBS

  INTEGER     , PARAMETER :: id_z_obs   =4001
  INTEGER     , PARAMETER :: id_dv_obs  =4002

  REAL(r_sngl) :: wk(7)

  REAL(r_size) , PARAMETER :: ERROR_Z = 5.0d0 , ERROR_DV= 1.0d0

 !--------------------------------------------FAKE RADAR PARAMETERS
 INTEGER , PARAMETER  :: radar_cen_x = 65 , radar_cen_y = 65 !Grid point position of radar.
 REAL(r_size), PARAMETER :: radar_z = 0.0d0
 REAL(r_size) :: radar_lat , radar_lon
 !Define the type of radar.
 INTEGER            :: METHOD_REF_CALC=3

 !------------------------------------------------------------------

CALL CPU_TIME(rtimer00)


CALL set_common_wrf(model_file_name)

ALLOCATE( gues3d(nlon,nlat,nlev,nv3d),gues2d(nlon,nlat,nv2d) )
ALLOCATE( z3d(nlon,nlat,nlev) , dv3d(nlon,nlat,nlev) )


CALL CPU_TIME(rtimer)
WRITE(6,'(A,2F10.2)') '### TIMER(SET MODEL):',rtimer,rtimer-rtimer00
rtimer00=rtimer


CALL read_grd(model_file_name,gues3d,gues2d)

gues3d(:,:,:,iv3d_ph)=gues3d(:,:,:,iv3d_ph)/gg

CALL CPU_TIME(rtimer)
WRITE(6,'(A,2F10.2)') '### TIMER(READ MODEL):',rtimer,rtimer-rtimer00
rtimer00=rtimer


 radar_lon = lon(radar_cen_x,radar_cen_y)
 radar_lat = lat(radar_cen_x,radar_cen_y)

!Compute dv and z at each model grid point.

 do ii=1,nlon-1
  do jj=1,nlat-1
   do kk=1,nlev-1
      

     dlon=lon(ii,jj)-radar_lon
     dlat=lat(ii,jj)-radar_lat

     rlev=gues3d(ii,jj,kk,iv3d_ph) 
     ri=real(ii,r_size)
     rj=real(jj,r_size)

     az=rad2deg*atan2(dlon*cosd(radar_lat),dlat)
     IF( az .LT. 0)az=360.0d0 + az
     !elevation (note that we do not correct elevation according to
     !beam propagation)
     CALL com_distll_1(lon(ii,jj),lat(ii,jj),radar_lon,radar_lat,dist)
     elev=rad2deg*atan2(rlev-radar_z,dist)

     qv=gues3d(ii,jj,kk,iv3d_qv)
     qc=gues3d(ii,jj,kk,iv3d_qc)
     qr=gues3d(ii,jj,kk,iv3d_qr)
     qci=gues3d(ii,jj,kk,iv3d_qci)
     qs= gues3d(ii,jj,kk,iv3d_qs)
     qg= gues3d(ii,jj,kk,iv3d_qg)
     t = gues3d(ii,jj,kk,iv3d_t)
     p = gues3d(ii,jj,kk,iv3d_p)

     
     u = 0.5*(gues3d(ii,jj,kk,iv3d_u)+gues3d(ii+1,jj,kk,iv3d_u))
     v = 0.5*(gues3d(ii,jj,kk,iv3d_v)+gues3d(ii,jj+1,kk,iv3d_v))
     w = 0.5*(gues3d(ii,jj,kk,iv3d_w)+gues3d(ii,jj,kk+1,iv3d_w))

     !Compute z and rv at this particular grid point.
     CALL calc_ref_vr(qv,qc,qr,qci,qs,qg,u,v,w,t,p,           &
          az,elev   &
         ,METHOD_REF_CALC,z3d(ii,jj,kk),dv3d(ii,jj,kk),att)

 
   enddo
  enddo
 enddo


CALL CPU_TIME(rtimer)
WRITE(6,'(A,2F10.2)') '### TIMER(FORWARD OPERATOR):',rtimer,rtimer-rtimer00
rtimer00=rtimer


!Write the results in LETKF-OBS format.
   !WRITE THE DATA
   !WRITE FILE HEADER.
   OPEN(UNIT=99,FILE='radarobs.dat',STATUS='unknown',FORM='unformatted')
   !Introduce a small header with the radar possition and two values that might be useful.
   write(99)REAL(radar_lon,r_sngl)
   write(99)REAL(radar_lat,r_sngl)
   write(99)REAL(radar_z,r_sngl)

    nobs=0
    DO ii=1,nlon-1
     DO jj=1,nlat-1
      DO kk=1,nlev-1

       !correspond to the location where the stronger echoes are located.
          !IF( grid_count_ref(ii,jj,kk) .GT. 0.0d0 )THEN
           wk(1)=REAL(id_z_obs,r_sngl)
           wk(6)=REAL(error_z,r_sngl)
           wk(2)=REAL(lon(ii,jj),r_sngl)
           wk(3)=REAL(lat(ii,jj),r_sngl)
           wk(4)=REAL(gues3d(ii,jj,kk,iv3d_ph),r_sngl)
           wk(5)=REAL(z3d(ii,jj,kk),r_sngl)
           wk(7)=REAL(1.0,r_sngl)
           WRITE(99)wk
           nobs = nobs + 1
          !ENDIF
          !IF( grid_count_vr(ii,jj,kk) .GT. 0.0d0 )THEN
         !DV is only available where reflectivity is over 5DBZ
         IF( 10.0d0*log( z3d(ii,jj,kk) ) .GT. 5.0d0 )THEN
           wk(1)=REAL(id_dv_obs,r_sngl)
           wk(6)=REAL(error_dv,r_sngl)
           wk(2)=REAL(lon(ii,jj),r_sngl)
           wk(3)=REAL(lat(ii,jj),r_sngl)
           wk(4)=REAL(gues3d(ii,jj,kk,iv3d_ph),r_sngl)
           wk(5)=REAL(dv3d(ii,jj,kk),r_sngl)
           wk(7)=REAL(1.0,r_sngl)
           WRITE(99)wk
           nobs=nobs +1
          ENDIF
      ENDDO
     ENDDO
    ENDDO

WRITE(*,*)'A TOTAL NUMBER OF ', nobs , ' HAS BEEN WRITTEN TO THE OBSERVATION FILE'


CALL CPU_TIME(rtimer)
WRITE(6,'(A,2F10.2)') '### TIMER(WRITE DATA):',rtimer,rtimer-rtimer00
rtimer00=rtimer

END PROGRAM WRF_TO_RADAR

!-----------------------------------------------------------------------
! Compute radar reflectivity and radial wind.
! Radial wind computations for certain methods depend on model reflectivity
! so both functions has been merged into a single one.
! First reflectivity is computed, and the the radial velocity is computed.
!-----------------------------------------------------------------------
SUBROUTINE calc_ref_vr(qv,qc,qr,qci,qs,qg,u,v,w,t,p,az,elev,method,ref,vr,k)
  USE common
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
  REAL(r_size)  :: cf, pip , roo
  REAL(r_size)  :: ki2 , kr2
  REAL(r_size)  :: lr , ls , lg
  REAL(r_size)  :: tmp_factor , rofactor
  REAL(r_size)  :: p0
  REAL(r_size)  :: Fs, Fg , zms , zmg , fws , fwg !Method 3
  REAL(r_size)  :: qrp , qsp , qgp
  REAL(r_size)  :: kr , ks , kg , kms , kmg !Partial contribution of different species to
                                            !the attenuation of the radar beam.
  REAL(r_size)  :: maxf                     !Maximum mixture relative concentration. (method 3)

  REAL(r_size) , PARAMETER  :: minz=0.01d0
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
  !Attenuation is not computed.

  !Reflectivity
  nor=8.0d6      ![m^-4]
  ror=1000.0d0   ![Kg/m3]
  pip=pi ** 1.75 !factor
  cf =10.0d18 * 72 !factor
  p0=1.0d5            !Reference pressure.

  qt=qr + qs + qg  !Assume that the total condensate is rain water
                   !But ignore cloud ice and cloud water

  IF( qt .GT. 0.0d0 )THEN
  ref = cf * ( ( ro * qt )**1.75 )
  ref = ref / ( pip * ( nor ** 0.75 ) * ( ror ** 1.75 ) )
  !ref= 2.04d4 *( ( ro * qt * 1.0d3 ) ** 1.75 ) !Original Sun and Crook expresion.
  ELSE
  ref=minz
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
  !It includes reflectivity contribution by all the microphisical species.
  !is assumes Marshall and Palmer distributions.
  !Based on C band radars.
  !Attenuation is not computed.
    nor=8.0d6      ![m^-4]
    nos=3.0d6      ![m^-4]
    nog=4.0d4      ![m^-4]
    ror=1000.0d0   ![Kg/m3]
    ros=100.0d0    ![Kg/m3]
    rog=913.0d0    ![Kg/m3]
    roi=917.0d0    ![Kg/m3]
    roo=1.0d0      ![Kg/m3] Surface air density.
    ki2=0.176d0    !Dielectric factor for ice.
    kr2=0.930d0    !Dielectric factor for water.
    pip=pi ** 1.75 !factor
    cf =1.0d18 * 720 !factor

    IF( qr .GT. 0.0d0 )THEN
    zr= cf * ( ( ro * qr )**1.75 )
    zr= zr / ( pip * ( nor ** 0.75 ) * ( ror ** 1.75 ) )
    ENDIF
    !The contribution of snow depends on temperature (bright band effect)
    IF( qs .GT. 0.0d0 )THEN
    IF ( t <= 273.16 )THEN
     zs = cf * ki2 * ( ros ** 0.25 ) * ( ( ro * qs ) ** 1.75 )
     zs = zs / ( pip * kr2 * ( nos ** 0.75  ) * ( roi ** 2 ) )
    ELSE
     !WARNING: This formulation has to be checked the paper says that
     !ros instead of roi should be used in thes expresion, but that
     !leads to unrealistic jumps in the model derived reflectivity.
     zs = cf * ( ( ro * qs ) ** 1.75 )
     zs = zs / ( pip * ( nos ** 0.75 ) * ( roi ** 1.75 ) )
    ENDIF
    ENDIF

    !Only dry graupel contribution is ussed.
    IF( qg .GT. 0.0d0 )THEN
    zg= ( cf / ( pip * ( nog ** 0.75) * ( rog ** 1.75 ) ) ) ** 0.95
    zg= zg * ( ( ro * qg ) ** 1.6625 )

    !zg=(ki2/kr2)*( cf / ( pip * ( nog ** 0.75) * ( rog ** 1.75 ) ) )
    !zg= zg * ( ( ro * qg ) ** 1.75 )

    !if( qg * ro * 1000 .GT. 1.0d0 )THEN
    !  WRITE(6,*)qg, zg, ro, rd , t , p, 10*log10(zg)
    !endif
    ENDIF

    ref = zr + zs + zg

   !Compute reflectivity weigthed terminal velocity.
    !Lin et al 1983.
    IF( ref > minz )THEN
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
    CALL com_gamma( 4.0d0 + b , tmp_factor )
    wr= a * tmp_factor / ( 6.0d0 * ( lr ** b ) )
    wr= 1.0d-2*wr * rofactor
    CALL com_gamma( 4.0d0 + d , tmp_factor )
    ws= c * tmp_factor / ( 6.0d0 * ( ls ** d ) )
    ws= 1.0d-2*ws * rofactor
    CALL com_gamma( 4.5d0 , tmp_factor )
    wg= tmp_factor * ( ( ( 4.0d0 * gg * rog )/( 3.0d0 * Cd * ro ) ) ** 0.5 )
    wg= 1.0d-2*wg / ( 6.0d0 * ( lg ** 0.5 ) )

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

 IF( ref > minz )THEN
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
    STOP

  ENDIF !END IF OVER DIFFERENT COMPUTATION OPTIONS.


  !Compute radial velocity
    vr = u * cos(elev*deg2rad) * sin(az*deg2rad)
    vr = vr + v * cos(elev*deg2rad) * cos(az*deg2rad)
    vr = vr + w*sin(elev*deg2rad)  !The contribution of the terminal velocity is not taken into account
    !vr = vr + (w - wt)*sin(elev*deg2rad)

    !WRITE(*,*)wt,p

    !WRITE(6,*) u , v , w
    !WRITE(6,*) wt , vr
    !WRITE(6,*) elev , az , deg2rad
    !STOP


  RETURN
END SUBROUTINE calc_ref_vr

SUBROUTINE itpl_3d(var,ri,rj,rk,var5,nlon,nlat,nlev)
  USE common
  IMPLICIT NONE
  INTEGER , INTENT(IN) :: nlon , nlat , nlev
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

