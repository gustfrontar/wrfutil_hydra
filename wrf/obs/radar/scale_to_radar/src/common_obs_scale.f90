MODULE common_obs_scale
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   10/04/2012 Guo-Yuan Lien     modified for GFS model
!   07/25/2014 Guo-Yuan Lien     modified for SCALE model
!   .......... See git history for the following revisions
!   14/07/2017 Juan Ruiz         modified for SCALTE to RADAR
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE netcdf
  USE common_namelist
  USE common_scale

  IMPLICIT NONE
  PUBLIC

!
! radar observations
!
  INTEGER,PARAMETER :: id_radar_ref_obs=4001
  INTEGER,PARAMETER :: id_radar_vr_obs=4002

!  INTEGER,PARAMETER :: iqc_good=0
!  INTEGER,PARAMETER :: iqc_gross_err=5
!  INTEGER,PARAMETER :: iqc_ps_ter=10
!  INTEGER,PARAMETER :: iqc_ref_low=11
!  INTEGER,PARAMETER :: iqc_ref_mem=12
!  INTEGER,PARAMETER :: iqc_radar_vhi=19
!  INTEGER,PARAMETER :: iqc_out_vhi=20
!  INTEGER,PARAMETER :: iqc_out_vlo=21
!  INTEGER,PARAMETER :: iqc_obs_bad=50
!  INTEGER,PARAMETER :: iqc_otype=90
!  INTEGER,PARAMETER :: iqc_time=97
!  INTEGER,PARAMETER :: iqc_out_h=98
!  INTEGER,PARAMETER :: iqc_undef=99

  

CONTAINS

!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
SUBROUTINE Trans_XtoY_radar(elm,radar_lon,radar_lat,radar_z,ri,rj,rk,lon,lat,lev,v3d,v2d,yobs,stggrd)
!  USE common_mpi
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: elm
  REAL(r_size),INTENT(IN) :: ri,rj,rk,radar_lon,radar_lat,radar_z !!!!! Use only, ri, rj, rk eventually... (radar_lon,lat,z in ri,rj,rk)
  REAL(r_size),INTENT(IN) :: lon,lat,lev
  REAL(r_size),INTENT(IN) :: v3d(nlev,nlonh,nlath,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlonh,nlath,nv2d)
  REAL(r_size),INTENT(OUT) :: yobs
  INTEGER,INTENT(IN),OPTIONAL :: stggrd
  INTEGER :: stggrd_ = 0

  REAL(r_size) :: qvr,qcr,qrr,qir,qsr,qgr,ur,vr,wr,tr,pr,rhr
  REAL(r_size) :: dist , dlon , dlat , az , elev , radar_ref,radar_rv

!  integer :: ierr
!  REAL(r_dble) :: rrtimer00,rrtimer
!  rrtimer00 = MPI_WTIME()


  if (present(stggrd)) stggrd_ = stggrd


  yobs = undef
  !qc = iqc_good

  if (stggrd_ == 1) then
    CALL itpl_3d(v3d(:,:,:,iv3d_u),rk,ri-0.5d0,rj,ur)  !###### should modity itpl_3d to prevent '1.0' problem....??
    CALL itpl_3d(v3d(:,:,:,iv3d_v),rk,ri,rj-0.5d0,vr)  !######
    CALL itpl_3d(v3d(:,:,:,iv3d_w),rk-0.5d0,ri,rj,wr)  !######
  else
    CALL itpl_3d(v3d(:,:,:,iv3d_u),rk,ri,rj,ur)
    CALL itpl_3d(v3d(:,:,:,iv3d_v),rk,ri,rj,vr)
    CALL itpl_3d(v3d(:,:,:,iv3d_w),rk,ri,rj,wr)
  end if
  CALL itpl_3d(v3d(:,:,:,iv3d_t),rk,ri,rj,tr)
  CALL itpl_3d(v3d(:,:,:,iv3d_p),rk,ri,rj,pr)
  CALL itpl_3d(v3d(:,:,:,iv3d_q),rk,ri,rj,qvr)
  CALL itpl_3d(v3d(:,:,:,iv3d_qc),rk,ri,rj,qcr)
  CALL itpl_3d(v3d(:,:,:,iv3d_qr),rk,ri,rj,qrr)
  CALL itpl_3d(v3d(:,:,:,iv3d_qi),rk,ri,rj,qir)
  CALL itpl_3d(v3d(:,:,:,iv3d_qs),rk,ri,rj,qsr)
  CALL itpl_3d(v3d(:,:,:,iv3d_qg),rk,ri,rj,qgr)


  !Compute az and elevation for the current observation.
  !Simple approach (TODO: implement a more robust computation)

  !Azimuth
  dlon=lon-radar_lon
  dlat=lat-radar_lat
  IF ( dlon == 0.0d0 .and. dlat == 0.0d0  )THEN
!      WRITE(6,*)'OOPS',dlon,dlat,lon,lat,radar_lon,radar_lat
!    qc = iqc_out_h
    RETURN
  ELSE
    az = rad2deg*atan2(dlon*cos(radar_lat*deg2rad),dlat)
  ENDIF
  !WRITE(6,*)dlon,dlat,lon,lat,radar_lon(ityp),radar_lat(ityp),dlon*cos(radar_lat(ityp))
  !IF( abs(dlon) > maxdlon )maxdlon=abs(dlon)
  !IF( abs(dlat) > maxdlat )maxdlat=abs(dlat)
  !WRITE(6,*)maxdlon,maxdlat
  IF( az < 0) az = 360.0d0 + az
  !elevation
  CALL com_distll_1(lon,lat,radar_lon,radar_lat,dist)
  elev=rad2deg*atan2(lev-radar_z,dist)

  CALL calc_ref_rv(qvr,qcr,qrr,qir,qsr,qgr,ur,vr,wr,tr,pr,az,elev,radar_ref,radar_rv)



!  SELECT CASE (elm)
!  CASE(id_radar_ref_obs)
!    if (radar_ref < MIN_RADAR_REF) then
!      qc = iqc_ref_low
!      yobs = MIN_RADAR_REF_DBZ + LOW_REF_SHIFT  !!! even if the above qc is bad, still return the value
!    else
!      yobs = 10.0d0 * log10(radar_ref)
!    end if
!  CASE(id_radar_vr_obs)
!    if (radar_ref < MIN_RADAR_REF) then
!      qc = iqc_ref_low
!    end if
!    yobs = radar_rv  !!! even if the above qc is bad, still return the value
!  CASE DEFAULT
!    qc = iqc_otype
!  END SELECT



  RETURN
END SUBROUTINE Trans_XtoY_radar

!-----------------------------------------------------------------------
! Compute radar reflectivity and radial wind.
! Radial wind computations for certain methods depend on model reflectivity
! so both functions has been merged into a single one.
! First reflectivity is computed, and the the radial velocity is computed.
!-----------------------------------------------------------------------
SUBROUTINE calc_ref_rv(qv,qc,qr,qci,qs,qg,u,v,w,t,p,az,elev,ref,vr)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: qv        !Water vapor
  REAL(r_size), INTENT(IN) :: qc,qr     !Cloud and rain water
  REAL(r_size), INTENT(IN) :: qci,qs,qg !Cloud ice, snow and graupel
  REAL(r_size), INTENT(IN) :: t,p       !Temperature and pressure.
  REAL(r_size), INTENT(IN) :: u,v,w     !velocities with respecto to earth.
  REAL(r_size), INTENT(INOUT) :: ref      !Reflectivity.
  REAL(r_size)              :: ro
  REAL(r_size), INTENT(IN) :: az     !Azimuth respect to the radar.
  REAL(r_size), INTENT(IN) :: elev   !Elevation angle respect to the surface.
  REAL(r_size), INTENT(INOUT) :: vr    !Radial velocity.
  REAL(r_size)  :: qms , qmg !Melting species concentration (METHOD_REF_CALC 3)
  REAL(r_size)  :: qt        !Total condensate mixing ratio (METHOD_REF_CALC 1)
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
  REAL(r_size)  :: Fs, Fg , zms , zmg , fws , fwg !METHOD_REF_CALC 3
  REAL(r_size)  :: qrp , qsp , qgp
  REAL(r_size)  :: maxf                     !Maximum mixture relative concentration. (METHOD_REF_CALC 3)

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
  ![refdb] refdb
  ![U,V,W] m/s
  ![az, elev] degree
  ![vr] m/s

  !This model reflectivity won't be lower than this value.

  !Initialize reflectivities
  zr=0.0d0
  zs=0.0d0
  zg=0.0d0
  zms=0.0d0
  zmg=0.0d0
  ref=0.0d0

  !Compute air density (all methods use this)

  ro =  p / (rd * t)

  !Begin computation of reflectivity and vr

  if (METHOD_REF_CALC == 1) then

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
    ref=0.0d0
    ENDIF

    !Radial wind

    IF ( qt .GT. 0.0d0 )THEN
    a=(p0/p)**0.4
    wt = 5.40d0 * a * ( qt ** 0.125 )
    ELSE
    wt=0d0
    ENDIF
   !WRITE(6,*)qr,qs,qg


  else if (METHOD_REF_CALC == 2) then
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
     !ros instead of roi should be used in this expresion, but that 
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
    IF( ref > 0.0d0 )THEN
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


    rofactor= ( roo / ro  ) ** 0.25
    IF( qr .gt. 0.0d0 )THEN
     lr= ( pi * ror * nor / ( ro * qr ) ) ** 0.25
     CALL com_gamma( 4.0_r_size + b , tmp_factor )
     wr= a * tmp_factor / ( 6.0d0 * ( lr ** b ) )
     wr= 1.0d-2*wr * rofactor
    ELSE
     wr=0.0d0
    ENDIF
    IF( qs .gt. 0.0d0 )THEN 
     ls= ( pi * ros * nos / ( ro * qs ) ) ** 0.25
     CALL com_gamma( 4.0_r_size + d , tmp_factor )
     ws= c * tmp_factor / ( 6.0d0 * ( ls ** d ) )
     ws= 1.0d-2*ws * rofactor
    ELSE
     ws=0.0d0
    ENDIF
    IF( qg .gt. 0.0d0 )THEN
     lg= ( pi * rog * nog / ( ro * qg ) ) ** 0.25
     CALL com_gamma( 4.5_r_size , tmp_factor )
     wg= tmp_factor * ( ( ( 4.0d0 * gg * 100.0d0 * rog )/( 3.0d0 * Cd * ro ) ) ** 0.5 )
     wg= 1.0d-2*wg / ( 6.0d0 * ( lg ** 0.5 ) )
    ELSE
     wg=0.0d0
    ENDIF
    !Reflectivity weighted terminal velocity. 
    wt = ( wr * zr + ws * zs + wg * zg )/ ( zr + zs + zg )

    ELSE

    wt=0.0d0

    ENDIF

  else if (METHOD_REF_CALC == 3) then
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
    ENDIF
    IF( qsp .GT. 0.0d0)THEN
    zs= 3.48d3 * ( ro * qsp * 1.0d3 )**1.66
    ENDIF
    IF( qgp .GT. 0.0d0)THEN
    zg= 8.18d4 * ( ro * qgp * 1.0d3 )**1.50
    ENDIF
    IF( qms .GT. 0.0d0 )THEN
    zms=( 0.00491 + 5.75*fws - 5.588*(fws**2) )*1.0d5
    zms= zms * ( ro * qms * 1.0d3 )**( 1.67 - 0.202*fws + 0.398*(fws**2) )

    ENDIF
    IF( qmg .GT. 0.0d0 )THEN
    zmg=( 0.809 + 10.13*fwg -5.98*(fwg**2) )*1.0d5
    zmg= zmg * ( ro * qmg * 1.0d3 )**( 1.48 + 0.0448*fwg - 0.0313*(fwg**2) )
    ENDIF

    ref = zr +  zg  + zs + zms + zmg

    !Compute reflectivity weigthed terminal velocity.
    !Lin et al 1983. (The distribution parameters are 
    !consistent with the work of Jung et al 2007)

    IF( ref > 0.0d0 )THEN
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
      CALL com_gamma( 4.0_r_size + b , tmp_factor )
      lr= ( pi * ror * nor / ( ro * qr ) ) ** 0.25
      wr= a * tmp_factor / ( 6.0d0 * ( lr ** b ) )
      wr= 1.0d-2 * wr * rofactor
      ELSE
      wr=0.0d0
      ENDIF

      IF( qs .GT. 0.0d0 )THEN
      ls= ( pi * ros * nos / ( ro * qs ) ) ** 0.25
      CALL com_gamma( 4.0_r_size + d , tmp_factor )
      ws= c * tmp_factor / ( 6.0d0 * ( ls ** d ) )
      ws= 1.0d-2 * ws * rofactor
      ELSE
      ws=0.0d0
      ENDIF

      IF ( qg .GT. 0.0d0 )THEN
      lg= ( pi * rog * nog / ( ro * qg ) ) ** 0.25
      CALL com_gamma( 4.5_r_size , tmp_factor )
      wg= tmp_factor * ( ( ( 4.0d0 * gg * 100.0d0 * rog )/( 3.0d0 * Cd * ro ) ) ** 0.5 )
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

  else  !IF OVER DIFFERENT OPTIONS

    WRITE(6,*)'ERROR: Not recognized method for radar reflectivity and wind computation'
    STOP

  end if ! [METHOD_REF_CALC == ?]


  !Compute radial velocity
  !WRITE(6,*)'ICRV',u,v,w,wt,az,elev
  vr = u * cos(elev*deg2rad) * sin(az*deg2rad)
  vr = vr + v * cos(elev*deg2rad) * cos(az*deg2rad)
  IF( USE_WT )THEN
    vr = vr + (w - wt)*sin(elev*deg2rad)
  ELSE
    vr = vr + (w)*sin(elev*deg2rad)
  ENDIF

  !WRITE(6,*) u , v , w
  !WRITE(6,*) wt , vr
  !WRITE(6,*) elev , az , deg2rad


  RETURN
END SUBROUTINE calc_ref_rv

!-----------------------------------------------------------------------
! Coordinate conversion (find rk in height levels)
!
! rk = 0.0d0  : surface observation
!-----------------------------------------------------------------------
SUBROUTINE phys2ijkz(z_full,ri,rj,rlev,rk)
!  use scale_grid_index, only: &
!      KHALO
!  use common_mpi
  IMPLICIT NONE

  REAL(r_size),INTENT(IN) :: z_full(nlev,nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(IN) :: rlev ! height levels
  REAL(r_size),INTENT(OUT) :: rk
  REAL(r_size) :: ak
  REAL(r_size) :: zlev(nlev)
  REAL(r_size) :: ztmp
  INTEGER :: i,j,k, ii, jj, ks

  !qc = iqc_good
!
! rlev -> rk
!
  if (ri < 1.0d0 + IHALO .or. ri > nlonh - IHALO .or. rj < 1.0d0 + JHALO .or. rj > nlath - JHALO ) then
    write (6,'(A)') 'warning: observation is outside of the horizontal domain'
    rk = undef
    !qc = iqc_out_h
    return
  end if

  !
  ! horizontal interpolation
  !
  i = CEILING(ri)
  j = CEILING(rj)
  !
  ! Find the lowest valid level
  !
  !ks = 1+KHALO
  !do jj = j-1, j
  !  do ii = i-1, i
  !    DO k=1+KHALO,nlev+KHALO
  !      if (z_full(k,ii,jj) > -300.0d0 .and. z_full(k,ii,jj) < 10000.0d0) exit
  !    END DO
  !    if (k > ks) ks = k
  !  end do
  !end do

  call itpl_2d_column(z_full,ri,rj,zlev)

  !
  ! determine if rlev is within bound.
  !
  IF(rlev > zlev(nlev)) THEN
    call itpl_2d(z_full(nlev,:,:),ri,rj,ztmp)
    write(6,'(A,F8.1,A,F8.1)') 'warning: observation is too high: ztop=', ztmp, ', lev=', rlev
    rk = undef
    !qc = iqc_out_vhi
    RETURN
  END IF
  IF(rlev < zlev(1)) THEN
    call itpl_2d(z_full(1,:,:),ri,rj,ztmp)
    write(6,'(A,F8.1,A,F8.1)') 'warning: observation is too low: zbottom=', ztmp, ', lev=', rlev
    rk = undef
    !qc = iqc_out_vlo
    RETURN
  END IF
  !
  ! find rk
  !
  DO k=2,nlev
    IF(zlev(k) > rlev) EXIT ! assuming ascending order of zlev
  END DO

  ak = (rlev - zlev(k-1)) / (zlev(k) - zlev(k-1))
  rk = REAL(k-1,r_size) + ak


  RETURN
END SUBROUTINE phys2ijkz
!-----------------------------------------------------------------------
! Coordinate conversion
!-----------------------------------------------------------------------
SUBROUTINE phys2ij(rlon,rlat,rig,rjg)
!  use scale_grid, only: &
!      GRID_CXG, &
!      GRID_CYG, &
!      DX, &
!      DY
  use scale_mapproj, only: &
      MPRJ_lonlat2xy
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(OUT) :: rig
  REAL(r_size),INTENT(OUT) :: rjg
!
! rlon,rlat -> ri,rj
!
  call MPRJ_lonlat2xy(rlon*deg2rad,rlat*deg2rad,rig,rjg)
 
  rig = (rig - CXG(IHALO)) / DX + 1.0d0 
  rjg = (rjg - CYG(JHALO)) / DY + 1.0d0 

  RETURN
END SUBROUTINE phys2ij

SUBROUTINE ij2phys(rig,rjg,rlon,rlat)
  !use scale_grid, only: &
  !    GRID_CXG, &
  !    GRID_CYG, &
  !    DX, &
  !    DY
  use scale_mapproj, only: &
      MPRJ_xy2lonlat
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rig
  REAL(r_size),INTENT(IN) :: rjg
  REAL(r_size),INTENT(OUT) :: rlon ! (deg)
  REAL(r_size),INTENT(OUT) :: rlat ! (deg)
  REAL(r_size) :: x, y ! (m)
!
! ri,rj -> rlon,rlat
!
  x = (rig - 1.0d0) * DX + CXG(1) 
  y = (rjg - 1.0d0) * DY + CYG(1) 

  call MPRJ_xy2lonlat(x,y,rlon,rlat)

  rlon = rlon * rad2deg
  rlat = rlat * rad2deg

  RETURN
END SUBROUTINE ij2phys
!
!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------
SUBROUTINE itpl_2d(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
     & + var(i  ,j-1) *    ai  * (1-aj) &
     & + var(i-1,j  ) * (1-ai) *    aj  &
     & + var(i  ,j  ) *    ai  *    aj

  RETURN
END SUBROUTINE itpl_2d

SUBROUTINE itpl_2d_column(var,ri,rj,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlev,nlonh,nlath)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5(nlev)
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  var5(:) = var(:,i-1,j-1) * (1-ai) * (1-aj) &
        & + var(:,i  ,j-1) *    ai  * (1-aj) &
        & + var(:,i-1,j  ) * (1-ai) *    aj  &
        & + var(:,i  ,j  ) *    ai  *    aj

  RETURN
END SUBROUTINE itpl_2d_column

SUBROUTINE itpl_3d(var,ri,rj,rk,var5)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: var(nlev,nlonh,nlath)
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

  var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
     & + var(i  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
     & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
     & + var(i  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
     & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
     & + var(i  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
     & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
     & + var(i  ,j  ,k  ) *    ai  *    aj  *    ak

  RETURN
END SUBROUTINE itpl_3d




END MODULE common_obs_scale
