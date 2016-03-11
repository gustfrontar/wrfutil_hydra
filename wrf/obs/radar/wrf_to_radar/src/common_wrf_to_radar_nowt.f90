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

  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: nid_obs=7
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073  
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
  INTEGER,PARAMETER :: id_tv_obs=3079
!
! surface observations codes > 9999
!
  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_rain_obs=19999
  INTEGER,PARAMETER :: id_tclon_obs=99991
  INTEGER,PARAMETER :: id_tclat_obs=99992
  INTEGER,PARAMETER :: id_tcmip_obs=99993

  !RADAR DATA ASSIMILATION PARAMETERS!!!!!!!!!!!!!!!!!!!!!
  !Values assigned to model derived reflectivity
  !in areas with no hydrometeors.
!  REAL(r_size), PARAMETER :: minz=0.01d0 , 
  !Ids for radar obserbations and pseudo observations.
  INTEGER     , PARAMETER :: id_reflectivity_obs=4001
  INTEGER     , PARAMETER :: id_radialwind_obs  =4002
  INTEGER     , PARAMETER :: id_pseudorh_obs    =4003

  !These 2 flags affects the computation of model reflectivity and
  !radial velocity.
  INTEGER            :: INTERPOLATION_TECHNIQUE=1
  INTEGER            :: METHOD_REF_CALC=3
  !This flag is to include or not the effect of attenuation in the
  !model derived reflectivities (this only makes sense if reflectivity
  !is going to be directly assimilated into the model)
  LOGICAL            :: COMPUTE_ATTENUATION=.FALSE.



CONTAINS

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


!-----------------------------------------------------------------------
! Horizontal coordinate conversion
!-----------------------------------------------------------------------
SUBROUTINE ll2ij(elem,rlon,rlat,ri,rj)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: elem
  REAL(r_size),INTENT(IN) :: rlon(1)
  REAL(r_size),INTENT(IN) :: rlat(1)
  REAL(r_size),INTENT(OUT) :: ri(1)
  REAL(r_size),INTENT(OUT) :: rj(1)
  INTEGER(4),PARAMETER :: isw = 3
  INTEGER :: i,j
  !
  ! rlon, rlat --> ri, rj
  !
  SELECT CASE(NINT(elem))
  CASE(id_u_obs)
    i=nlon
    j=nlat-1
    CALL com_pos2ij(isw,i,j,lonu(1:i,1:j),latu(1:i,1:j),1,rlon,rlat,ri,rj)
  CASE(id_v_obs)
    i=nlon-1
    j=nlat
    CALL com_pos2ij(isw,i,j,lonv(1:i,1:j),latv(1:i,1:j),1,rlon,rlat,ri,rj)
  CASE DEFAULT
    i=nlon-1
    j=nlat-1
    CALL com_pos2ij(isw,i,j,lon(1:i,1:j),lat(1:i,1:j),1,rlon,rlat,ri,rj)
  END SELECT

  RETURN
END SUBROUTINE ll2ij

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
  REAL(r_size) :: zlev(nlev)
  INTEGER :: i,j,k,n,zindex
  rk=0.0d0
  !
  ! rlev --> rk
  !
    i=INT(ri)
    j=INT(rj)
    IF( i < 1 )i=1
    IF( i > nlon )i=nlon
    IF( j < 1 )j=1
    IF( j > nlat )j=nlat
    zlev(1:nlev)=z_full(i,j,1:nlev)
    !Find rk
    DO k=FLOOR(rk),nlev
      zindex=k
      IF(zlev(k) > rlev) EXIT ! assuming increasing order of zlev
    END DO
    ak = (rlev - zlev(zindex-1)) / (zlev(zindex) - zlev(zindex-1))
    rk = REAL(zindex-1,r_size) + ak

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
! Interpolate model data to a radar (or part of a radar) domain
!-----------------------------------------------------------------------
SUBROUTINE model_to_radar( input_radar , v3d , v2d  )
  IMPLICIT NONE
  TYPE(RADAR), INTENT(INOUT) :: input_radar
  REAL(r_size),INTENT(IN)    :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN)    :: v2d(nlon,nlat,nlev,nv2d)
  INTEGER              :: ia , ir , ie , i , j , k
  REAL(r_size)         :: ri(1) , rj(1) , rk(1) , rlev(1)
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


  ! First allocate arrays where the model derived reflectivity and radial velocity will be stored.
   !ALLOCATE( input_radar%i( input_radar%na , input_radar%nr , input_radar%ne ) )
   !ALLOCATE( input_radar%j( input_radar%na , input_radar%nr , input_radar%ne ) )
   !ALLOCATE( input_radar%k( input_radar%na , input_radar%nr , input_radar%ne ) )
   !ALLOCATE( input_radar%p( input_radar%na , input_radar%nr , input_radar%ne ) )

   !input_radar%i = 0.0d0
   !input_radar%j = 0.0d0
   !input_radar%k = 0.0d0
   !input_radar%p = 0.0d0
  
   !Compute maximum model height
   !WRITE(*,*)input_radar%lat(1,input_radar%nr,1)
   !STOP
   max_model_z=0.0d0
   DO i=1,nlon
     DO j=1,nlat
        IF( v3d(i,j,nlev-1,iv3d_ph).GT. max_model_z )max_model_z=v3d(i,j,nlev-1,iv3d_ph)
     ENDDO
   ENDDO




  IF ( COMPUTE_ATTENUATION .AND. METHOD_REF_CALC .LE. 2)THEN
    WRITE(6,*)'WARNING: Attenuation can not be computed with methods 1 or 2.'
  !  COMPUTE_ATTENUATION=.FALSE. 
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
 
       CALL  ll2ij(tmp,input_radar%lon(ia,ir,ie),input_radar%lat(ia,ir,ie),ri,rj)

       tmp_z=input_radar%z(ia,ir,ie)
       CALL  z2k_fast(v3d(:,:,:,iv3d_ph),ri(1),rj(1),tmp_z,rk(1))
        !Interpolate qv,qc,qr,qci,qs,qg,t and p to the beam center.
        IF( ri(1) >= 1 .AND. ri(1) <= nlon .AND. rj(1) >= 1 .AND.    &
            rj(1) <= nlat .AND. rk(1) >= 1 .AND. rk(1) <= nlev )THEN


          !input_radar%i(ia,ir,ie) = ri
          !input_radar%j(ia,ir,ie) = rj
          !input_radar%k(ia,ir,ie) = rk

          !CALL k2p(v3d(:,:,:,iv3d_p),id_reflectivity_obs,ri,rj,rk,rlev)
          !input_radar%p = rlev

          CALL itpl_3d(v3d(:,:,:,iv3d_qv),ri(1),rj(1),rk(1),qv)
          CALL itpl_3d(v3d(:,:,:,iv3d_qc),ri(1),rj(1),rk(1),qc)
          CALL itpl_3d(v3d(:,:,:,iv3d_qr),ri(1),rj(1),rk(1),qr)
          CALL itpl_3d(v3d(:,:,:,iv3d_qci),ri(1),rj(1),rk(1),qci)
          CALL itpl_3d(v3d(:,:,:,iv3d_qs),ri(1),rj(1),rk(1),qs)
          CALL itpl_3d(v3d(:,:,:,iv3d_qg),ri(1),rj(1),rk(1),qg)
          CALL itpl_3d(v3d(:,:,:,iv3d_t),ri(1),rj(1),rk(1),t)
          CALL itpl_3d(v3d(:,:,:,iv3d_p),ri(1),rj(1),rk(1),p)
          !Interpolate winds.
          CALL itpl_3d(v3d(:,:,:,iv3d_u),ri(1),rj(1),rk(1),u)
          CALL itpl_3d(v3d(:,:,:,iv3d_v),ri(1),rj(1),rk(1),v)
          CALL itpl_3d(v3d(:,:,:,iv3d_w),ri(1),rj(1),rk(1),w)

          !input_radar%p(ia,ir,ie) = p

          !Compute reflectivity at the beam center.
          !a=input_radar%azimuth(ia)
          !e=input_radar%elevation(ie) 
          CALL calc_ref_vr(qv,qc,qr,qci,qs,qg,u,v,w,t,p,           &
               input_radar%azimuth(ia),input_radar%elevation(ie)   &
                ,METHOD_REF_CALC,ref,vr,att)
    
          !IF(qr .gt. 1e-5)WRITE(*,*)qv,qc,qr,qci,qs,qg,u,v,w,t,p,vr
          !Compute radial velocity at the beam center. 

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
            !WRITE(*,*) ia, ie, ir , pik , att 
          ENDIF
 
          IF( input_radar%iv3d_ref .GT. 0 )THEN
          input_radar%radarv3d_model(ia,ir,ie,input_radar%iv3d_ref) =ref !refdb
          ENDIF

          !Will generate wind observations only where reflectivity data is good enough (in this case were we have clouds). 
          IF( input_radar%iv3d_wind .GT. 0  .AND. input_radar%radarv3d_model(ia,ir,ie,input_radar%iv3d_ref) .GT. minz )THEN
          input_radar%radarv3d_model(ia,ir,ie,input_radar%iv3d_wind)=vr    !Radial wind
          ENDIF

          !IF( input_radar%iv3d_prh .GT. 0 )THEN
          !Compute model relative humidity to be compared with radar pseudo relative humidity
          !CALL calc_rh(t,qv,p,prh)
          !input_radar%radarv3d_model(ia,ir,ie,input_radar%iv3d_prh)=prh    !Pseudo rh
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


END MODULE common_wrf_to_radar
