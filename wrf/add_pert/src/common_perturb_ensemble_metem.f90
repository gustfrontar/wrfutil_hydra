MODULE common_perturb_ensemble
!=======================================================================
!
! [PURPOSE:] Common Information for WRF
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   07/14/2010 Takemasa Miyoshi  modified for WRF
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_wrf
  USE common_namelist
  USE common_smooth2d
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------

!  LOGICAL :: PERTURB_T , PERTURB_RH , PERTURB_WIND   !Whether this variables will be perturbed.
!  REAL(r_size) :: PERTURB_T_AMP , PERTURB_RH_AMP , PERTURB_WIND_AMP   !Amplitude of the perturbations.
!  REAL(r_size) :: PERTURB_T_SCLH , PERTURB_RH_SCLH , PERTURB_WIND_SCLH  !Horizontal scale of the perturbations.
!  REAL(r_size) :: PERTURB_T_SCLV , PERTURB_RH_SCLV , PERTURB_WIND_SCLV  !Horizontal scale of the perturbations.

!  LOGICAL :: PERTURB_ATMOSPHERE= .TRUE.  !If atmospheric variables will be perturbed
!  LOGICAL :: PERTURB_SST       = .TRUE.  !If sst will be perturbed
!  LOGICAL :: PERTURB_SOIL      = .TRUE.  !If soiltemp and soilmoisture will be perturbed


!  REAL(r_size) :: amp_factor , random_amp_factor

CONTAINS

!-----------------------------------------------------------------------
! Get time interpolation factor.
!-----------------------------------------------------------------------
SUBROUTINE time_interpol_factor(date1,date2,dateint,interpolation_factor)
character(19),INTENT(IN) :: date1,date2,dateint
integer                  :: iy,im,id,ih,imin
real(r_size)             :: tai1,tai2,taiint , sec
real(r_size)             :: interpolation_factor
!expected date format "2008-08-07_00:00:00"
   read(date1(1:4) ,*)iy
   read(date1(6:7) ,*)im
   read(date1(9:10),*)id
   read(date1(12:13),*)ih
   read(date1(15:16),*)imin
   read(date1(18:19),*)sec
   call com_utc2tai(iy,im,id,ih,imin,sec,tai1)

   read(date2(1:4) ,*)iy
   read(date2(6:7) ,*)im
   read(date2(9:10),*)id
   read(date2(12:13),*)ih
   read(date2(15:16),*)imin
   read(date2(18:19),*)sec
   call com_utc2tai(iy,im,id,ih,imin,sec,tai2)

   read(dateint(1:4) ,*)iy
   read(dateint(6:7) ,*)im
   read(dateint(9:10),*)id
   read(dateint(12:13),*)ih
   read(dateint(15:16),*)imin
   read(dateint(18:19),*)sec
   call com_utc2tai(iy,im,id,ih,imin,sec,taiint)

  if ( taiint - tai1 /= 0 )then
   interpolation_factor= REAL( taiint - tai1, r_size) / REAL( tai2 - tai1,r_size)
  else
   interpolation_factor= 0.0d0
  endif

  if ( interpolation_factor < 0.0d0 )then
    write(*,*)'[Error] : interpolation factor is less than one '
    write(*,*)'Date 2 has to be later than date 1 '
    write(*,*)'Date2: ',date2
    write(*,*)'Date1: ',date1
    stop
  endif
  if ( interpolation_factor > 1.0d0 )then
    write(*,*)'[Error] : can not extrapolate perturbations in time '
    write(*,*)'Date2: ',date2
    write(*,*)'Date1: ',date1
    write(*,*)'Date int: ',dateint
    stop
  endif

END SUBROUTINE time_interpol_factor

!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE  PERTURB_MEMBER_BALANCED(var3d1,var2d1,var3d2,var2d2,pert3d,pert2d)
IMPLICIT NONE
REAL(r_size), INTENT(IN) :: var3d1(nlon,nlat,nlev,nv3d),var3d2(nlon,nlat,nlev,nv3d)
REAL(r_size), INTENT(IN) :: var2d1(nlon,nlat,nv2d)     ,var2d2(nlon,nlat,nv2d)
REAL(r_size), INTENT(OUT):: pert3d(nlon,nlat,nlev,nv3d),pert2d(nlon,nlat,nv2d)
INTEGER :: iv

pert3d=0.0d0
pert2d=0.0d0

DO iv=1,nv3d
   if( PERTURB_ATMOSPHERE )then
      WRITE(*,*)'Adding balanced perturbation to field ',element(iv),' iv=',iv," amp_factor = ",amp_factor
      pert3d=amp_factor*( var3d1 - var3d2)
         
   endif
ENDDO

DO iv=1,nv2d
   SELECT CASE(iv)
    CASE( iv2d_psfc , iv2d_pmsl )
     if( PERTURB_ATMOSPHERE )then
       WRITE(*,*)'Adding balanced perturbation to field ',element(nv3d+iv),' iv= ',iv
       pert2d(:,:,iv)=amp_factor * ( var2d1(:,:,iv)-var2d2(:,:,iv) )
     endif

    CASE( iv2d_skintemp)
     if( PERTURB_SST)then
       WRITE(*,*)'Adding balanced perturbation to field ',element(nv3d+iv),' iv= ',iv
       pert2d(:,:,iv)=amp_factor * ( var2d1(:,:,iv)-var2d2(:,:,iv) )
     endif
    CASE( iv2d_st100200 , iv2d_st040100 , iv2d_st010040 , iv2d_st000010 , &
          iv2d_sm100200 , iv2d_sm040100 , iv2d_sm010040 , iv2d_sm000010  )
     WRITE(*,*)'Adding balanced perturbation to field ',element(nv3d+iv),' iv= ',iv
     if( PERTURB_SOIL)then
       pert2d(:,:,iv)=amp_factor * ( var2d1(:,:,iv)-var2d2(:,:,iv) )
     endif
   END SELECT

ENDDO


!Get perturbation maximum and minimum

CALL pert_max_min(pert3d,pert2d)


END SUBROUTINE PERTURB_MEMBER_BALANCED

SUBROUTINE  PERTURB_MEMBER_RANDOM(var3d,var2d,pert3d,pert2d)
  IMPLICIT NONE
  REAL(r_size) :: var3d(nlon,nlat,nlev,nv3d) , var2d(nlon,nlat,nv2d)
  REAL(r_size) :: pert3d(nlon,nlat,nlev,nv3d) , pert2d(nlon,nlat,nv2d)
  REAL(r_size) :: tmpfield(nlon,nlat,nlev) , tmpval(1)
  INTEGER      :: i , j , k , iv
  INTEGER      :: filter_size_x , filter_size_z 
  REAL(r_size) :: AMP , SCLH , SCLV , tmp_amp_factor 
  REAL(r_size) :: dz
  REAL(r_size) :: max_pert , min_pert 
  INTEGER      :: maskh(nlon,nlat),maskv(nlev)
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
           ELSEIF( iv == iv3d_rh .AND.  PERTURB_RH )THEN
             AMP=PERTURB_RH_AMP
             SCLH=PERTURB_RH_SCLH
             SCLV=PERTURB_RH_SCLV
           ELSEIF( ( iv == iv3d_u .OR. iv == iv3d_v ) .AND. PERTURB_WIND )THEN
             AMP=PERTURB_WIND_AMP
             SCLH=PERTURB_WIND_SCLH
             SCLV=PERTURB_WIND_SCLV
           ELSE
             CYCLE
           ENDIF
           WRITE(*,*)'Adding random perturbation to field ',element(iv),' iv= ',iv

           !PERTURB VARIBLE IV
           filter_size_x = NINT( SCLH / dx )     
           dz= var3d(2,2,nlev-1,iv3d_ght)/gg / REAL(nlev,r_size)  !Approximate dz asuming equispaced levels.      
           filter_size_z = NINT( SCLV / dz )

           DO k=1,nlev
             DO j=1,nlat
                CALL com_randn2(nlon,tmpfield(:,j,k),0)
             ENDDO
             CALL filter_2d(tmpfield(:,:,k),tmpfield(:,:,k),maskh,filter_size_x,filter_type,nlon,nlat)
           ENDDO 
           DO i =1,nlat
            DO j =1,nlon
             CALL filter_2d(tmpfield(i,j,:),tmpfield(i,j,:),maskv,filter_size_z,filter_type,nlev,1)
            ENDDO
           ENDDO

         !Adjust the amplitud of the perturbation.
         pert_var=0.0d0
         pert_mean=0.0d0
         DO i = 1 , nlon
          DO j = 1 , nlat
           DO k = 1, nlev
           pert_var=pert_var+tmpfield(i,j,k)**2
           pert_mean=pert_mean+tmpfield(i,j,k)
           ENDDO
          ENDDO
         ENDDO
         pert_mean=pert_mean/REAL(nlon*nlat*nlev,r_size)
         pert_var =sqrt(pert_var/REAL(nlon*nlat*nlev,r_size) - pert_mean**2)

         !Match the perturbation standard deviation with the requested standard
         !deviation.
         tmp_amp_factor =  AMP / ( pert_var )
         tmpfield= tmpfield * tmp_amp_factor



         pert3d(:,:,:,iv)=pert3d(:,:,:,iv)+tmpfield

     ENDDO

  CALL pert_max_min(pert3d,pert2d)

  RETURN
END SUBROUTINE PERTURB_MEMBER_RANDOM

SUBROUTINE PERT_MAX_MIN(pert3d,pert2d)
IMPLICIT NONE
REAL(r_size) , INTENT(IN) :: pert3d(nlon,nlat,nlev,nv3d) , pert2d(nlon,nlat,nv2d)
INTEGER                   :: iv , i , j , k
REAL(r_size)              :: pert_max , pert_min

  DO iv=1,nv3d
     DO k=1,nlev
      pert_max=-9.0d9
      pert_min=9.0d9

      DO i=1,nlon-1
       DO j=1,nlat-1
   
          IF( pert3d(i,j,k,iv) < pert_min)pert_min=pert3d(i,j,k,iv)
          IF( pert3d(i,j,k,iv) > pert_max)pert_max=pert3d(i,j,k,iv)

       ENDDO
      ENDDO

     IF( iv == iv3d_t  )THEN
       WRITE(*,*)"Temperature pert. range at level ",k," is ",pert_min,pert_max
     ELSEIF( iv == iv3d_u )THEN
       WRITE(*,*)"U-wind  pert. range at level ",k," is ",pert_min,pert_max
     ELSEIF( iv == iv3d_v )THEN
       WRITE(*,*)"V-wind  pert. range at level ",k," is ",pert_min,pert_max
     ELSEIF( iv == iv3d_rh )THEN
       WRITE(*,*)"RH  pert. range at level ",k," is ",pert_min,pert_max
     ENDIF
 
     ENDDO

  ENDDO

RETURN
END SUBROUTINE PERT_MAX_MIN

SUBROUTINE CHECK_PERT(var3d,var2d)
IMPLICIT NONE
REAL(r_size) :: var3d(nlon,nlat,nlev,nv3d)
REAL(r_size) :: var2d(nlon,nlat,nv2d)
INTEGER      :: iv

DO iv=1,nv3d
   if(  iv == iv3d_rh )then
      where( var3d(:,:,:,iv) <= 1.0d-2 )
          var3d(:,:,:,iv)=1.0d-2
      endwhere
   endif
ENDDO

DO iv=1,nv2d
   if(  iv == iv2d_st100200 .or. iv == iv2d_st040100  &
   .or. iv == iv2d_st010040 .or. iv == iv2d_st000010 )then
      where( var2d(:,:,iv) < 0.0d0 )
         var2d(:,:,iv)=0.0d0
      endwhere
   endif
ENDDO


END SUBROUTINE CHECK_PERT


END MODULE common_perturb_ensemble
