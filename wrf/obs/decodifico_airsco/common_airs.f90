MODULE common_airs
USE common
USE common_mtx

IMPLICIT NONE
PUBLIC

INTEGER,PARAMETER :: id_co_obs=5002
INTEGER,PARAMETER :: id_totco_obs=5001
!
! surface observations codes > 9999
!
LOGICAL,PARAMETER :: verbose = .FALSE.
INTEGER,PARAMETER :: airstype = 21 ! observation type record for AIRS
REAL(r_size),PARAMETER :: lon1 =  270.0d0 ! minimum longitude
REAL(r_size),PARAMETER :: lon2 =  340.0d0 ! maximum longitude
REAL(r_size),PARAMETER :: lat1 =  -70.0d0 ! minimum latitude
REAL(r_size),PARAMETER :: lat2 =   20.0d0 ! maximum latitude

REAL(r_size),PARAMETER :: coef_error = 1.0d0  
INTEGER,PARAMETER :: iskip = 3
INTEGER,PARAMETER :: nunit=70
INTEGER,PARAMETER :: nx=30
INTEGER,PARAMETER :: ny=45
INTEGER,PARAMETER :: nzt=28
INTEGER,PARAMETER :: nzco=100
INTEGER,PARAMETER :: n_trapezoids=9 !Number of CO Layers in AIRS Retrievals.
 
REAL(r_dble) :: rlon(nx,ny)
REAL(r_dble) :: rlat(nx,ny)
REAL(r_dble) :: time(nx,ny),time1,time2
REAL(r_sngl) :: covmr(nzco,nx,ny)
REAL(r_sngl) :: covmre(nzco,nx,ny)
REAL(r_sngl) :: covmrqc(nzco,nx,ny)
REAL(r_sngl) :: co_ave_kernel(n_trapezoids,n_trapezoids,nx,ny)
REAL(r_size) :: co_prior(12,nzco,2) 
REAL(r_size) :: smooth_profile(n_trapezoids)
REAL(r_size) :: pressure_trapezoids(n_trapezoids+1)   !Pressure limits of AIRS CO Trapezoid function.
REAL(r_size) :: pressure_levels(nzco)                 !Pressure levels of AIRS CO Suplementary retrievals.
INTEGER      :: level_trapezoids(n_trapezoids+1)

REAL(r_sngl) :: psurf(nx,ny)
INTEGER(2) :: qq(nx,ny)

CONTAINS
!=======================================================================
! READING HDF
!=======================================================================
SUBROUTINE read_hdf(filename)
IMPLICIT NONE
CHARACTER(*) :: filename
INTEGER :: statn
INTEGER :: fid
INTEGER :: nswath
INTEGER :: nchar
INTEGER :: swid
INTEGER :: swopen,swinqswath,swattach,swrdattr,swrdfld,swdetach,swclose
CHARACTER(256) :: swathname
INTEGER :: i,j,k
  !
  ! open
  !
  !READ(5,'(A)') filename
  PRINT *,filename
  !PRINT *,filename
  fid = swopen(filename,1)
  IF(fid == -1) THEN
    PRINT *,'IO ERROR: opening file ',filename
    STOP
  END IF
  nswath = swinqswath(filename,swathname,nchar)
  IF(nswath /= 1) THEN
    PRINT *,'FILE ERROR: bad nswath ',nswath
    STOP
  END IF
  IF(swathname /= 'L2_Standard_atmospheric&surface_product') THEN
    PRINT *,'FILE ERROR: bad swath name ',swathname
    STOP
  END IF
  swid = swattach(fid, swathname)
  IF(swid == -1) THEN
    PRINT *,'FILE ERROR: failed to attach to swath ',swathname
    STOP
  END IF
  !
  ! read
  !
!  sec = sec4
  statn = swrdattr(swid,'start_Time',time1)
  statn = swrdattr(swid,'end_Time',time2)
  statn = swrdfld(swid,'Longitude',(/0,0/),(/1,1/),(/nx,ny/),rlon)
  statn = swrdfld(swid,'Latitude',(/0,0/),(/1,1/),(/nx,ny/),rlat)
  statn = swrdfld(swid,'Time',(/0,0/),(/1,1/),(/nx,ny/),time)

  statn = swrdfld(swid,'COVMRLevSup',(/0,0,0/),(/1,1,1/),(/nzco,nx,ny/),covmr)
  statn = swrdfld(swid,'COVMRLevSupErr',(/0,0,0/),(/1,1,1/),(/nzco,nx,ny/),covmre)
  statn = swrdfld(swid,'COVMRLevSup_QC',(/0,0,0/),(/1,1,1/),(/nzco,nx,ny/),covmrqc)
  statn = swrdfld(swid,'CO_ave_kern',(/0,0,0/),(/1,1,1/),(/n_trapezoids,n_trapezoids,nx,ny/),co_ave_kernel)
  statn = swrdfld(swid,'PSurfStd',(/0,0/),(/1,1/),(/nx,ny/),psurf)
  !
  ! close
  !
  statn = swdetach(swid)
  statn = swclose(fid)

END SUBROUTINE read_hdf


SUBROUTINE read_first_guess()
IMPLICIT NONE
INTEGER :: ii , jj 
  !
  !  READ first guess data
  !
  !Read Northern Hemisphere prior.
  OPEN(44,FILE='airs_co_prior_nh.tbl',FORM='FORMATTED')
   READ(44,*)
   READ(44,*)pressure_levels
   DO ii=1,12
     READ(44,*)co_prior(ii,:,1)
   ENDDO
  CLOSE(44)
  !Read Southern Hemisphere prior.
  OPEN(44,FILE='airs_co_prior_sh.tbl',FORM='FORMATTED')
   READ(44,*)
   READ(44,*)
   DO ii=1,12
     READ(44,*)co_prior(ii,:,2)
   ENDDO
  CLOSE(44)

  !Invert original data order
  DO jj=1,2
   DO ii=1,12
      CALL invert_order(co_prior(ii,:,jj),nzco)
   ENDDO
  ENDDO
  CALL invert_order(pressure_levels,nzco)
  

END SUBROUTINE read_first_guess

!Compute first guess for a particular latitutude and month.
SUBROUTINE get_first_guess(latitude,month,first_guess)
IMPLICIT NONE
REAL(r_size),INTENT(IN) :: latitude
INTEGER     ,INTENT(IN) :: month
REAL(r_size),INTENT(OUT):: first_guess
REAL(r_size)            :: alphash,alphanh

IF( latitude > 15.0d0)THEN
  alphash=0.0d0
  alphanh=1.0d0
ELSEIF( latitude < -15.0d0)THEN
  alphash=1.0d0
  alphanh=0.0d0
ELSE
  alphash=(15.0d0-latitude)/30.0d0
  alphanh=1.0d0-alphash
ENDIF

 first_guess=co_prior(month,:,1)*alphanh+co_prior(month,:,2)*alphash

END SUBROUTINE get_first_guess

SUBROUTINE read_trapezoid_limits()
IMPLICIT NONE
INTEGER :: ii
  ! 
  !  READ trapezoid limits in pressure (hPa)
  !
  !Read in the levels corresponding to trapezoid limits in AIRS data.
   OPEN(33,FILE='co_trapezoid_layers.tbl',FORM='FORMATTED')
    DO ii=n_trapezoids+1,1,-1
      READ(33,*)level_trapezoids(ii),pressure_trapezoids(ii)
    ENDDO
    
    CALL invert_order(pressure_trapezoids,n_trapezoids)
    CALL invert_order(level_trapezoids,n_trapezoids)

END SUBROUTINE read_trapezoid_limits

!-----------------------------------------------------------------------
! From vertical profile to trapezoid layers. 
!-----------------------------------------------------------------------
!This subroutine generates the transformation matrix between a vertical profile
!and AIRS trapezoid layers.
!Pressure has to be in hPa.

SUBROUTINE model_to_airs(co_profile,co_prior_profile,co_profile_airs,surface_pressure,surface_co,surface_co_prior)
IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: co_profile(nzco) , co_prior_profile(nzco)
  REAL(r_size), INTENT(IN) :: surface_pressure , surface_co  , surface_co_prior
  REAL(r_size),ALLOCATABLE :: F(:,:),Finv(:,:)
  INTEGER              :: ii  , jj , kk
  REAL(r_size)         :: tmp_trapezoids_limit(4),tmp_w_value(4)
  REAL(r_size)         :: pmin , pmax
  REAL(r_size),INTENT(OUT) :: co_profile_airs(n_trapezoids)
  REAL(r_size),INTENT(IN)  :: co_ave_kernel(n_trapezoids,n_trapezoids)
  REAL(r_size),ALLOCATABLE :: work1(:,:),work2(:,:),work3(:),eivec(:,:),eival(:),tmp_co_profile(:)
  INTEGER              :: include_level(n_trapezoids), min_trap , max_trap , n_eff_trap
  INTEGER              :: n_pos_eival
  REAL(r_size),ALLOCATABLE :: co_profile_eff(:) , co_profile_airs_eff(:) , co_prior_profile_eff(:)
  REAL(r_size),ALLOCATABLE :: pressure_trapezoids_eff(:) ,pressure_levels_eff(:)
  REAL(r_size),ALLOCATABLE :: log_pressure_trapezoids_eff(:) , log_pressure_levels_eff(:)
  REAL(r_size),ALLOCATABLE :: log_co_profile_eff(:) , log_co_profile_airs_eff(:)
  REAL(r_size),ALLOCATABLE :: log_co_prior_profile_eff(:) 
  REAL(r_size),ALLOCATABLE :: co_ave_kernel_eff(:,:)
  INTEGER                  :: ini_p_lev , ini_trap , nzco_eff , n_trapezoids_eff
  
  !First we need to detect the efective number of levels and the efective number
  !of trapezoids.
  !We assume that trapezoids and co_profile levels are ordered in decreasing P.
  !First element corresponds to the level closer to the surface.

  !Generate the effective trapezoid layers.
  DO ii=1,ntrapezoid
     IF( pressure_trapezoid(ii) >= surface_pressure .AND. pressure_trapezoid(ii) <= surface_pressure )THEN
     ini_trap=ii
     ENDIF
  ENDDO
  n_trapezoids_eff=n_trapezodis-ini_trap+1
  ALLOCATE( pressure_trapezoids_eff(n_trapezoids_eff+1))
  ALLOCATE( log_pressure_trapezoids_eff(n_trapezoids_eff+1))
  ALLOCATE( co_profile_airs_eff(n_trapezoids_eff))
  ALLOCATE( log_co_profile_airs_eff(n_trapezoids_eff))
  ALLOCATE( co_ave_kernel_eff(n_trapezoids_eff,n_trapezoids_eff) )
  pressure_trapezoids_eff(1)=surface_pressure
  DO ii=ini_trap+1,n_trapezoid+1
    pressure_trapezoids_eff(ii-ini_trap+1)=pressure_trapezoids(ii)
  ENDDO
  co_profile_airs_eff=undef
  log_co_profile_airs_eff=undef
  co_ave_kernel_eff=co_ave_kernel(ini_trap:n_trapezoids,ini_trap:n_trapezoids)

  !Generate the effective pressure levels
  DO ii=1,nzco
     IF( pressure_levels(ii) >= surface_pressure .AND. pressure_levels(ii+1) <= surface_pressure )THEN
     ini_p_lev=ii
     ENDIF
  ENDDO
  n_pressure_levels_eff=nzco - ini_p_lev + 1
  ALLOCATE( pressure_levels_eff(n_pressure_levels_eff))
  ALLOCATE( log_pressure_levels_eff(n_pressure_levels_eff))
  ALLOCATE( co_profile_eff(n_pressure_levels_eff))
  ALLOCATE( co_prior_profile_eff(n_pressure_levels_eff))
  ALLOCATE( log_co_prior_profile_eff(n_pressure_levels_eff))
  ALLOCATE( log_co_profile_eff(n_pressure_levels_eff))

  pressure_levels_eff(1)=surface_pressure
  co_profile_eff(1)=surface_co
  co_prior_profile_eff(1)=surface_co
  
  DO ii=ini_p_lev+1,nzco
    pressure_levels_eff(ii-ini_p_lev+1)=pressure_levels(ii)
    co_profile_eff(ii-ini_p_lev+1)=co_profile(ii)
    co_prior_profile_eff(ii-ini_p_lev+1)=co_prior_profile(ii)
  ENDDO
  log_co_prior_profile_eff=log10(co_prior_profile_eff)
  log_co_profile_eff=log10(co_profile_eff)-log_co_prior_profile_eff


  !Compute pressure in logaritmic scale.
  log_pressure_levels_eff=log10(pressure_levels_eff)
  log_pressure_trapezoids_eff=log10(pressure_trapezoids_eff)


  tmp_w_value=0.0d0
  !Define F
  ALLOCATE( F(nzco_eff,n_trapezoids_eff) )

  DO ii =1,n_trapezoids_eff
    IF( ii == 1)THEN
      tmp_trapezoids_limit(1)=log_pressure_trapezoids_eff(1)-1.0e-10
      tmp_trapezoids_limit(2:4)=log_pressure_trapezoids_eff(ii:ii+2)
      tmp_w_value(3)=0.5d0
      tmp_w_value(2)=1.0d0
    ELSEIF( ii == n_trapezoids )THEN
      tmp_trapezoids_limit(1:3)=log_pressure_trapezoids_eff(ii-1:ii+1)
      tmp_trapezoids_limit(4)=log_pressure_trapezoids_eff(n_trapezoids_eff+1)+1.0e-10
      tmp_w_value(3)=1.0d0
      tmp_w_value(2)=0.5d0
    ELSE
      tmp_trapezoids_limit(1:4)=log_pressure_trapezoids_eff(ii-1:ii+2)
      tmp_w_value(2:3)=0.5d0
    ENDIF

    CALL linear_interp(n_trapezoids_eff,nzco_eff, &
  & log_pressure_trapezoids_eff,log_pressure_levels_eff ,      &
  & tmp_w_value,F(:,ii))

  ENDDO

  WHERE( F == undef )
     F=0.0d0
  ENDWHERE


  !Compute vertical smoothing of the model profile based on ARIS products
  !documentation V6.
  ALLOCATE(Finv(n_trapezoids_eff,nzco_eff) )
  ALLOCATE(work1(n_trapezoids_eff,n_trapezoids_eff),eivec(n_trapezoids_eff,n_trapezoids_eff),eival(n_trapezoids_eff) )
  ALLOCATE(work2(n_trapezoids_eff,n_trapezoids_eff) )
  ALLOCATE(work3(n_trapezoids_eff) )
  !Compute F^T * F
  CALL dgemm('t','n',n_trapezoids_eff,n_trapezoids_eff,nzco_eff,1.0d0,F,nzco_eff,F,&
    & nzco_eff,0.0d0,work1,n_trapezoids_eff)
  !Compute inv( F^T*F )
  CALL mtx_eigen(1,n_trapezoids_eff,work1,eival,eivec,n_pos_eival)
  DO jj=1,n_trapezoids_eff
    DO ii=1,n_trapezoids_eff
      work1(ii,jj) = eivec(ii,jj) / eival(jj)
    END DO
  END DO
  !Compute psedo inverse of F inv(F^T*F)*F^T
  CALL dgemm('n','t',n_trapezoids_eff,n_trapezoids_eff,n_trapezoids_eff,1.0d0,work1,n_trapezoids_eff,eivec,&
    & n_trapezoids_eff,0.0d0,work2,n_trapezoids_eff)
  CALL dgemm('n','t',n_trapezoids_eff,nzco_eff,n_trapezoids_eff,1.0d0,work2,n_trapezoids_eff,F,&
    & nzco_eff,0.0d0,Finv,nzco_eff)
  !Finv * log_co_profile
  CALL dgemm('n','n',n_trapezoids_eff,1,nzco_eff,1.0d0,Finv,n_trapezoids_eff,log_co_profile,&
    & nzco_eff,0.0d0,log_co_profile_eff,n_trapezoids_eff)

  ! A * Finv * log_co_profile
   CALL dgemm('n','n',n_trapezoids_eff,1,n_trapezoids_eff,1.0d0,co_ave_kernel_eff,n_trapezoids_eff,log_co_profile_eff,&
    & n_trapezoids_eff,0.0d0,work3,n_trapezoids_eff)
    log_co_profile_eff = work3

  ! F * A * Finv * log_co_profile 
     CALL dgemm('n','n',nzco_eff,1,n_trapezoids_eff,1.0d0,F,n_trapezoids_eff,log_co_profile_eff,&
    & nzco_eff,0.0d0,log_co_profile_eff,nzco_eff)

  DO ii=1,nzco
    IF( ii >= ini_p_lev )THEN
      co_profile(ii)=10**(log_co_prior_profile(ii-ini_p_lev+1) + log_co_profile_eff(ii-ini_p_lev+1))
    ELSE
      co_profile(ii)=undef
    ENDIF 
 
  ENDDO
  

  


  DEALLOCATE(F,Finv)
  DEALLOCATE(work1,work2,work3,eivec,eival)
  !Faltan deallocates!

END SUBROUTINE model_to_airs

!SUBROUTINE moleccmsq_to_mixingratio(moleccmsq,deltaP,mixingratio)
!Converts from molecules per cm2 to mixing ratio.
!Pressure input in Pa.
!  REAL(r_size), INTENT(IN) :: moleccsmsq , deltaP
!  REAL(r_size), INTENT(OUT):: mixingratio
!  REAL(r_size), PARAMETER :: Av=6.0221409d+23

  !mixingratio=moleccmsq*R*g/(Av*10-3*Rd*deltaP)
  !The factor 10-4 transforms from Pa=N/m2 to N/cm2
  !R is the ideal gas constant.
  !Rd is the ideal gas constant for dry air.
  !g is the acceleration of gravity.

!  mixingratio=moleccmsq*8.31*gg/(Av*10d-4*Rd*deltaP)

!END SUBROUTINE moleccmsq_to_mixingratio

Subroutine linear_interp(    &
     & klevi ,&  ! in
     & klevf ,&  ! in
     & xi ,&  ! in
     & xf ,&  ! in
     & yi  ,&  ! in
     & yf)  ! out  

  ! Description:
  ! To interpolate the array vec from the presi levels to presf levels
  !
  ! Copyright:
  !
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2002, EUMETSAT, All Rights Reserved.
  !
  !
  ! Method:
  ! Linear interpolation 
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date        Comment
  ! -------   ----        -------
  ! 1         09/2002     ECMWF

  Implicit None

  !
  ! Subroutine arguments
  !
  Integer, Intent(in) :: klevi      ! number of levels of the initial grid
  Integer, Intent(in) :: klevf      ! number of levels of the final grid
  !
  Real(r_size), Intent(in), Dimension(klevi)  :: xi ! initial grid
  Real(r_size), Intent(in), Dimension(klevf)  :: xf ! final grid
  !
  Real(r_size), Intent(in), Dimension(klevi)  :: yi  ! initial vec array
  Real(r_size), Intent(out), Dimension(klevf) :: yf  ! final vec array
  !
  ! Local scalars :
  !
  Integer :: jki, jkf
  Real(r_size)    :: slope, t1, t2, p1, p2, lp1, lp2
  !
  !- End of header --------------------------------------------------------

  yf(:) = undef

  Do jkf = 1,klevf
     Do jki = 1,klevi-1
        p1 = xi(jki)
        p2 = xi(jki+1)
        If (xf(jkf) >= p1 .And. xf(jkf) < p2) Then
           t1 = yi(jki)
           t2 = yi(jki+1)
           slope = (t1-t2)/(p1-p2)
           If (t2 == 0.) slope = 0.
           yf(jkf) = t1 + slope*(yf(jkf)-p1)
           !
        Else If (jki == 1 .And. xf(jkf) < p1) Then
           yf(jkf) = undef
        Else If (jki == (klevi-1) .And. yf(jkf) == undef ) Then
           yf(jkf) = undef
        End If
     End Do
  End Do

  Return
End Subroutine linear_interp


END MODULE common_airs

