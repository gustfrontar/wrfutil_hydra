MODULE common_breeding
!=======================================================================
!
! [PURPOSE:] Common Breeding routines for SCALE
!
! [HISTORY:]
!   10/02/2017 Juan Ruiz  created  - Based on rescaling norms by Shigenori Otsuka
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_scale
  USE common_namelist
  USE common_smooth2d

  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE compute_norm(v3d_p,v2d_p,v3d_n,v2d_n,norm)

IMPLICIT NONE
!character(len=*)  :: norm_type   !Norm to be used: UV , T ,  UVT
real(r_size), intent(in) :: v3d_p(nlev,nlong,nlatg,nv3d) , v3d_n(nlev,nlong,nlatg,nv3d)
real(r_size), intent(in) :: v2d_p(nlong,nlatg,nv2d) , v2d_n(nlong,nlatg,nv2d)
!integer     , intent(in) :: is,ie,js,je,zs,ze  !Subdomain for norm computation (only for global norms)

real(r_size) , intent(out) :: norm(nlev,nlong,nlatg) !Norm is always 3D but in the case of global norms
                                                   !all its values will be the same.
real(r_size) :: normg
integer      :: maskh(nlong,nlatg),maskv(nlev) , filter_size_x , filter_size_z 
integer      :: ii,jj,kk

character(256) :: filter_type='Lanczos' 

norm=0.0d0 !Initialize norm

!First compute local norm (this will be averaged later for global norms)

!write(*,*)v3d_p(20,20,20,iv3d_u),v3d_n(20,20,20,iv3d_u)

if ( TRIM( norm_type ) == 'UV'  )then

   write(*,*) "We will use UV norm"

   norm = ( ( v3d_p(:,:,:,iv3d_u) - v3d_n(:,:,:,iv3d_u) ) / wind_scale )**2  &
        + ( ( v3d_p(:,:,:,iv3d_v) - v3d_n(:,:,:,iv3d_v) ) / wind_scale )**2

   norm = norm / 2.0d0 


elseif( TRIM( norm_type ) == 'T' )then

   write(*,*) "We will use T norm"

   norm = ( ( v3d_p(:,:,:,iv3d_t) - v3d_n(:,:,:,iv3d_t) ) / t_scale )**2

elseif( TRIM( norm_type ) == 'UVT'  )then

   write(*,*) "We will use UVT norm"

   norm = ( ( v3d_p(:,:,:,iv3d_u) - v3d_n(:,:,:,iv3d_u) ) / wind_scale )**2  &
        + ( ( v3d_p(:,:,:,iv3d_v) - v3d_n(:,:,:,iv3d_v) ) / wind_scale )**2  & 
        + ( ( v3d_p(:,:,:,iv3d_t) - v3d_n(:,:,:,iv3d_t) ) / t_scale )**2

   norm = norm / 3.0d0

endif

   call norm_check( norm )


!If ave_norm_xy option is enabled, then the norm will be horizontally averaged.

if( ave_norm_xy )then

 if( is == 0 )is = 1
 if( ie == 0 )ie = nlong
 if( js == 0 )js = 1
 if( je == 0 )je = nlatg

 DO kk = 1 , nlev 
  normg=0.0d0
   DO ii = is , ie
    DO jj = js , je
         
       normg = normg + norm(kk,ii,jj)

    ENDDO
   ENDDO

   norm(kk,:,:) = normg / REAL( (ie-is+1)*(je-js+1) ,r_size )  

 ENDDO

endif

!If ave_norm_z option is enabled, then the norm will be vertically averaged.

if( ave_norm_z )then

  if( ks == 0 )ks = 1
  if( ke == 0 )ke = nlev

  DO ii = 1 , nlong
   DO jj = 1 , nlatg

    normg=0
    
    DO kk = ks , ke
       
       normg=normg + norm(kk,ii,jj)

    ENDDO

    norm(:,ii,jj)=normg / REAL( (ke-ks+1) , r_size ) 

   ENDDO
  ENDDO


endif

!If smooth_xy is enabled, perform horizontall smoothing of the norm.

if( smooth_norm_xy ) then

!Use a simple lanczos filter to smooth the norm.
 maskh       = 1
 maskv       = 1

 filter_size_x = NINT( smooth_sx / dx )

 DO kk = 1 , nlev
  CALL filter_2d(norm(kk,:,:),norm(kk,:,:),maskh,filter_size_x,filter_type,nlong,nlatg)
 ENDDO

endif 

!If smooth_z is enabled, perform vertical smoothing of the norm.

if( smooth_norm_z ) then
    
 filter_size_z = NINT( smooth_sz / dz )

 DO ii  = 1 , nlatg
  DO jj = 1 , nlong

   CALL filter_2d(norm(:,ii,jj),norm(:,ii,jj),maskv,filter_size_z,filter_type,nlev,1)

  ENDDO
 ENDDO
 

endif  

  norm = sqrt( norm )

END SUBROUTINE compute_norm

SUBROUTINE norm_check(norm)

IMPLICIT NONE
real(r_size) , intent(in) :: norm(nlev,nlong,nlatg) !Norm is always 3D but in the case of global norms
real(r_size) :: norm_max , norm_min , norm_mean
integer      :: ii,jj,kk

 norm_min=9.99d9
 norm_max=0.0d0
 norm_mean=0.0d0

  DO ii = 1 , nlong
   DO jj = 1 , nlatg
    DO kk = 1 , nlev

       if( norm(kk,ii,jj) > norm_max) norm_max=norm(kk,ii,jj)
       if( norm(kk,ii,jj) < norm_min) norm_min=norm(kk,ii,jj)
       norm_mean = norm_mean + norm(kk,ii,jj)

    ENDDO
   ENDDO
  ENDDO

  norm_mean=norm_mean / REAL( nlong*nlatg*nlev, r_size)

  WRITE(*,*) "--------------- NORM CHECK ------------------"
  WRITE(*,*) "norm_max = ",norm_max
  WRITE(*,*) "norm_min = ",norm_min
  WRITE(*,*) "norm_mean= ",norm_mean
  WRITE(*,*) "norm_threshold = ",norm_threshold
  WRITE(*,*) "--------------- NORM CHECK ------------------"


END SUBROUTINE norm_check

SUBROUTINE rescale( v3d_c , v2d_c , v3d_p , v2d_p , v3d_n , v2d_n , norm )
!The main goal of this routin is to rescale the perturbations.
IMPLICIT NONE
real(r_size), intent(in)    :: v3d_c(nlev,nlong,nlatg,nv3d)    !Ctrl 3D variables.
real(r_size), intent(inout) :: v3d_p(nlev,nlong,nlatg,nv3d)    !Positive perturbation 3D variables.
real(r_size), intent(inout) :: v3d_n(nlev,nlong,nlatg,nv3d)    !Negative perturbation 3D variables.
real(r_size), intent(in)    :: v2d_c(nlong,nlatg,nv2d)         !Ctrl 2D variables.
real(r_size), intent(inout) :: v2d_p(nlong,nlatg,nv2d)         !Positive perturbation 2D variables.
real(r_size), intent(inout) :: v2d_n(nlong,nlatg,nv2d)         !Negative perturbation 2D variables.

real(r_size), intent(inout) :: norm(nlev,nlong,nlatg)       !Perturbation norm.

real(r_size)  :: v3d_pert(nlev,nlong,nlatg,nv3d) , v2d_pert(nlong,nlatg,nv2d)  !Auxiliary perturbation array.

integer :: ii , jj , kk , iv , pcount


!If the bred vector is not initialize then initialize it.
!Rescaling won't be applied to the initial perturbation.

if ( init_bv ) then

   write(*,*)"We are going to initialize the bred vector"

!Initialize perturbation according to namelist parameters.

   call  get_random_pert(v3d_pert,v2d_pert)

   !Damp perturbation amplitude near the boundaries.
   call damp_latbnd(v3d_pert,v2d_pert)
 
   v3d_p = v3d_c + v3d_pert / 2.0d0
   v2d_p = v2d_c + v2d_pert / 2.0d0
   v3d_n = v3d_c - v3d_pert / 2.0d0
   v2d_n = v2d_c - v2d_pert / 2.0d0

   !Get perturbations in rho, rhou, rhov, rhow consistent with the
   !initial random perturbations introduced in u, v, w and t.
   call state_trans_inv(v3d_p)
   call state_trans_inv(v3d_n)

   !Compute the initial perturbation norm. 
   call compute_norm( v3d_p , v2d_p , v3d_n , v2d_n , norm )

else

   !Compute nomr according to namelist parameters.
   call compute_norm( v3d_p , v2d_p , v3d_n , v2d_n , norm )

endif

   !Compute perturbation (if bv is initialized recompute perturbation
   !to include perturbations in the model native variables)
   v3d_pert = v3d_p - v3d_n
   v2d_pert = v2d_p - v2d_n



!Rescale perturbation

 pcount = 0

  do ii = 1 , nlong
   do jj = 1 , nlatg

     if( norm(1,ii,jj) >= norm_threshold ) then

       v2d_pert(ii,jj,:) = v2d_pert(ii,jj,:) * ( norm_threshold / norm(1,ii,jj) )
  
     endif

     do kk = 1, nlev
      if (  norm(kk,ii,jj) >= norm_threshold )then

         pcount = pcount + 1

         v3d_pert(kk,ii,jj,:) = v3d_pert(kk,ii,jj,:) * ( norm_threshold / norm(kk,ii,jj) )

      endif
     enddo 
  enddo
 enddo

 !Damp perturbation boundary near the boundaries.
 call damp_latbnd(v3d_pert,v2d_pert)

 !Now generate the new perturbation pair centered at the control state.

 v3d_p = v3d_c + v3d_pert / 2.0d0
 v3d_n = v3d_c - v3d_pert / 2.0d0

 v2d_p = v2d_c + v2d_pert / 2.0d0
 v2d_n = v2d_c - v2d_pert / 2.0d0
  
 write(*,*)"---------------------------------------------------------------"
 write(*,*)" Perturbation has been rescaled at ",pcount," grid points"
 write(*,*)" That is ",(100*real(pcount,r_size)/real(nlong*nlatg*nlev,r_size))," % of the total grid points"
 write(*,*)"---------------------------------------------------------------"


END SUBROUTINE rescale

SUBROUTINE write_norm( norm )
IMPLICIT NONE
character(len=8) :: filename='norm.grd'
real(r_size) , intent(in) :: norm(nlev,nlong,nlatg)
integer :: ii , jj , kk , iunit , reclength
integer :: irec , iv 
!Write the norm to a direct access unfromatted binary file for ploting.

 iunit=44
 INQUIRE(IOLENGTH=reclength) reclength
 reclength = 4* reclength * nlong * nlatg

 irec=1

 OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=reclength)
 

   DO kk = 1 , nlev
  
     WRITE(iunit,rec=irec)real(norm(kk,:,:) ,r_sngl)

     irec = irec + 1

  ENDDO

  write(*,*)"NX = ",nlong," NY = ",nlatg," NLEV = ",nlev


END SUBROUTINE write_norm

SUBROUTINE  get_random_pert(pert3d,pert2d)
  IMPLICIT NONE
  REAL(r_size) , INTENT(OUT) :: pert3d(nlev,nlong,nlatg,nv3d) , pert2d(nlong,nlatg,nv2d)
  REAL(r_size) :: tmpfield(nlev,nlong,nlatg) , tmpval(1)
  INTEGER      :: i , j , k , iv
  INTEGER      :: filter_size_x , filter_size_z 
  REAL(r_size) :: AMP , SCLH , SCLV , tmp_amp_factor 
  REAL(r_size) :: dz
  REAL(r_size) :: max_pert , min_pert 
  INTEGER      :: maskh(nlong,nlatg),maskv(nlev)
  CHARACTER(256) :: filter_type 
  INTEGER      :: ii,jj,kk
  REAL(r_size) :: factor , pert_var , pert_mean
  !In this case on output we get the perturbation (not the full field)

  pert3d=0.0d0
  pert2d=0.0d0 

  filter_type = 'Lanczos'
  maskh       = 1
  maskv       = 1


     DO iv = 1 , nv3d
          
           IF( iv == iv3d_t  .AND.  PERTURB_T )THEN
             AMP=PERTURB_T_AMP
             SCLH=PERTURB_T_SCLH
             SCLV=PERTURB_T_SCLV
           ELSEIF( ( iv == iv3d_u .OR. iv == iv3d_v ) .AND. PERTURB_WIND )THEN
             AMP=PERTURB_WIND_AMP
             SCLH=PERTURB_WIND_SCLH
             SCLV=PERTURB_WIND_SCLV
           ELSE
             CYCLE
           ENDIF

           WRITE(*,*)'Adding random perturbation to field ',v3d_name(iv),' iv= ',iv

           !PERTURB VARIBLE IV
           filter_size_x = NINT( SCLH / dx )     
           !dz= var3d(2,2,nlev-1,iv3d_ght)/gg / REAL(nlev,r_size)  !Approximate dz asuming equispaced levels.      
           filter_size_z = NINT( SCLV / dz )

           DO k=1,nlev
             DO j=1,nlatg
                CALL com_randn2(nlong,tmpfield(k,:,j),0)
             ENDDO
             CALL filter_2d(tmpfield(k,:,:),tmpfield(k,:,:),maskh,filter_size_x,filter_type,nlong,nlatg)
           ENDDO 
           DO i =1,nlatg
            DO j =1,nlong
             CALL filter_2d(tmpfield(:,i,j),tmpfield(:,i,j),maskv,filter_size_z,filter_type,nlev,1)
            ENDDO
           ENDDO

         !Adjust the amplitud of the perturbation.
         pert_var=0.0d0
         pert_mean=0.0d0
         DO i = 1 , nlong
          DO j = 1 , nlatg
           DO k = 1, nlev
           pert_var=pert_var+tmpfield(k,i,j)**2
           pert_mean=pert_mean+tmpfield(k,i,j)
           ENDDO
          ENDDO
         ENDDO
         pert_mean=pert_mean/REAL(nlong*nlatg*nlev,r_size)
         pert_var =sqrt(pert_var/REAL(nlong*nlatg*nlev,r_size) - pert_mean**2)

         !Match the perturbation standard deviation with the requested standard
         !deviation.
         tmp_amp_factor = AMP / ( pert_var )
         tmpfield= tmpfield * tmp_amp_factor

         pert3d(:,:,:,iv)=pert3d(:,:,:,iv) + tmpfield

     ENDDO

  CALL pert_max_min(pert3d,pert2d)

  RETURN
END SUBROUTINE get_random_pert

SUBROUTINE pert_max_min(pert3d,pert2d)
IMPLICIT NONE
REAL(r_size) , INTENT(IN) :: pert3d(nlev,nlong,nlatg,nv3d) , pert2d(nlong,nlatg,nv2d)
INTEGER                   :: iv , i , j , k
REAL(r_size)              :: pert_max , pert_min

  DO iv=1,nv3d
     DO k=1,nlev
      pert_max=-9.0d9
      pert_min=9.0d9

      DO i=1,nlong-1
       DO j=1,nlatg-1
   
          IF( pert3d(k,i,j,iv) < pert_min)pert_min=pert3d(k,i,j,iv)
          IF( pert3d(k,i,j,iv) > pert_max)pert_max=pert3d(k,i,j,iv)

       ENDDO
      ENDDO

     IF( iv == iv3d_t  )THEN
       WRITE(*,*)"Temperature pert. range at level ",k," is ",pert_min,pert_max
     ELSEIF( iv == iv3d_u )THEN
       WRITE(*,*)"U-wind  pert. range at level ",k," is ",pert_min,pert_max
     ELSEIF( iv == iv3d_v )THEN
       WRITE(*,*)"V-wind  pert. range at level ",k," is ",pert_min,pert_max
     !ELSEIF( iv == iv3d_rh )THEN
     !  WRITE(*,*)"RH  pert. range at level ",k," is ",pert_min,pert_max
     ENDIF
 
     ENDDO

  ENDDO

RETURN
END SUBROUTINE pert_max_min


!-----------------------------------------------------------------------
! Damp lateral boundary
!-----------------------------------------------------------------------
SUBROUTINE damp_latbnd(v3d,v2d)
  IMPLICIT NONE
  REAL(r_size),INTENT(inout) :: v3d(nlev,nlong,nlatg,nv3d),v2d(nlong,nlatg,nv2d)
  INTEGER :: i,j,k,iv
  INTEGER :: ist,ied,jst,jed

  ist = 1
  ied = nlong - 1
  jst = 1
  jed = nlatg - 1

 DO iv = 1 , nv3d
  DO k = 1, nlev
    DO i = ist, bdy_width
      v3d(k,i,1:nlatg,iv) = v3d(k,i,1:nlatg,iv) * real(i-1) / real(bdy_width)
    END DO
    DO i = ied-bdy_width+1, ied
      v3d(k,i,1:nlatg,iv) = v3d(k,i,1:nlatg,iv) * real(ied-i) / real(bdy_width)
    END DO
    DO j = jst, bdy_width
      v3d(k,1:nlong,j,iv) = v3d(k,1:nlong,j,iv) * real(j-1) / real(bdy_width)
    END DO
    DO j = jed-bdy_width+1, jed
      v3d(k,1:nlong,j,iv) = v3d(k,1:nlong,j,iv) * real(jed-j) / real(bdy_width)
    END DO
  END DO
 END DO

 DO iv = 1 , nv2d
    DO i = ist, bdy_width
      v2d(i,1:nlatg,iv) = v2d(i,1:nlatg,iv) * real(i-1) / real(bdy_width)
    END DO
    DO i = ied-bdy_width+1, ied
      v2d(i,1:nlatg,iv) = v2d(i,1:nlatg,iv) * real(ied-i) / real(bdy_width)
    END DO
    DO j = jst, bdy_width
      v2d(1:nlong,j,iv) = v2d(1:nlong,j,iv) * real(j-1) / real(bdy_width)
    END DO
    DO j = jed-bdy_width+1, jed
      v2d(1:nlong,j,iv) = v2d(1:nlong,j,iv) * real(jed-j) / real(bdy_width)
    END DO
 END DO


  RETURN
END SUBROUTINE




END MODULE common_breeding
