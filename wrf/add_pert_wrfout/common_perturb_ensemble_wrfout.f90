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
  USE common_smooth2d
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------

  LOGICAL :: PERTURB_T = .true. , PERTURB_Q = .true.  , PERTURB_WIND = .true.   !Whether this variables will be perturbed.
  REAL(r_size) :: PERTURB_T_AMP = 1.0d0 , PERTURB_Q_AMP = 1.0d-3 , PERTURB_WIND_AMP = 1.0d0   !Amplitude of the perturbations.
  REAL(r_size) :: PERTURB_T_SCLH = 1.0d4 , PERTURB_Q_SCLH = 1.0d4 , PERTURB_WIND_SCLH = 1.0d4  !Horizontal scale of the perturbations.
  REAL(r_size) :: PERTURB_T_SCLV = 0.5d4 , PERTURB_Q_SCLV = 0.5d4 , PERTURB_WIND_SCLV = 0.5d4  !Horizontal scale of the perturbations.
  
  REAL(r_size) :: dt = 1
  REAL(r_size) :: SCLT  = 10
  INTEGER :: ntimes = 1

  CHARACTER(100) :: namelist_file='perturb.namelist'

  LOGICAL :: file_exist
  INTEGER :: ierr

CONTAINS


!-----------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------

SUBROUTINE read_namelist

NAMELIST / GENERAL / perturb_t , perturb_q , perturb_wind , &
                     perturb_t_amp , perturb_q_amp , perturb_wind_amp , &
                     perturb_t_sclh, perturb_q_sclh, perturb_wind_sclh, &
                     perturb_t_sclv, perturb_q_sclv, perturb_wind_sclv, &
                     dt , sclt , ntimes

INQUIRE(FILE=NAMELIST_FILE, EXIST=file_exist)

IF( .NOT. file_exist )THEN
  WRITE(*,*)"ERROR: Could not find namelist"
ENDIF

OPEN(54,FILE=NAMELIST_FILE)
READ(54,NML=GENERAL,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at GENERAL section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)

END SUBROUTINE read_namelist



SUBROUTINE  PERTURB_MEMBER_RANDOM(var3d,var2d)
  IMPLICIT NONE
  REAL(r_size) :: var3d(nlon,nlat,nlev,nv3d,ntimes),var2d(nlon,nlat,nv2d,ntimes)
  REAL(r_size) :: tmpfield(nlon,nlat,nlev,ntimes) , tmpval(1)
  INTEGER      :: i , j , k , iv
  INTEGER      :: filter_size_x , filter_size_z  , filter_size_t
  REAL(r_size) :: AMP , SCLH , SCLV  
  REAL(r_size) :: dz 
  REAL(r_size) :: amp_factor
  INTEGER      :: maskh(nlon,nlat),maskv(nlev)
  CHARACTER(256) :: filter_type 
  INTEGER      :: ii,jj,kk,it
  REAL(r_size) :: pert_var , pert_mean
  !In this case on output we get the perturbation (not the full field)

  filter_type = 'Lanczos'
  maskh       = 1
  maskv       = 1

     DO iv = 1 , nv3d
          
           IF( iv == iv3d_t  .AND.  PERTURB_T )THEN
             AMP=PERTURB_T_AMP
             SCLH=PERTURB_T_SCLH
             SCLV=PERTURB_T_SCLV
           ELSEIF( iv == iv3d_qv .AND.  PERTURB_Q )THEN
             AMP=PERTURB_Q_AMP
             SCLH=PERTURB_Q_SCLH
             SCLV=PERTURB_Q_SCLV
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
         dz= ( var3d(2,2,nlev-1,iv3d_ph,1)+var3d(2,2,nlev-1,iv3d_phb,1) ) /gg / REAL(nlev,r_size)  !Approximatedzasuming equispaced levels.
         filter_size_z = NINT( SCLV / dz )
         filter_size_t = NINT( SCLT / dt ) 

         DO it=1,ntimes
           DO k=1,nlev
             DO j=1,nlat
                CALL com_randn2(nlon,tmpfield(:,j,k,it),0)
             ENDDO
             CALL filter_2d(tmpfield(:,:,k,it),tmpfield(:,:,k,it),maskh,filter_size_x,filter_type,nlon,nlat)
           ENDDO

           DO i =1,nlat
            DO j =1,nlon
             CALL filter_2d(tmpfield(i,j,:,it),tmpfield(i,j,:,it),maskv,filter_size_z,filter_type,nlev,1)
            ENDDO
           ENDDO
         ENDDO
  
          DO i =1,nlat
            DO j =1,nlon
              DO k = 1,nlev
               CALL filter_2d(tmpfield(i,j,k,:),tmpfield(i,j,k,:),maskv,filter_size_t,filter_type,ntimes,1)
              ENDDO
            ENDDO
          ENDDO



     ! Adjust the amplitud of the perturbation.

     pert_var=0.0d0
     pert_mean=0.0d0
     DO i = 1 , nlon
      DO j = 1 , nlat
       DO k = 1, nlev
        DO it = 1 , ntimes
           pert_var=pert_var+tmpfield(i,j,k,it)**2
           pert_mean=pert_mean+tmpfield(i,j,k,it)
        ENDDO
       ENDDO
      ENDDO
     ENDDO
     pert_mean=pert_mean/REAL(nlon*nlat*nlev*ntimes,r_size)
     pert_var =sqrt(pert_var/REAL(nlon*nlat*nlev*ntimes,r_size) - pert_mean**2)

     !Match the perturbation standard deviation with the requested standard
     !deviation.
     amp_factor =  AMP / ( pert_var )
     tmpfield= tmpfield * amp_factor

     !Damp perturbation near the boundaries.
     !CALL damp_latbnd_double(spec_bdy_width,nlon,nlat,nlev,tmpfield)
     !Apply perturbation

     var3d(:,:,:,iv,:)=var3d(:,:,:,iv,:)+tmpfield

     ENDDO ![End do over variables]



  RETURN
END SUBROUTINE PERTURB_MEMBER_RANDOM

SUBROUTINE CHECK_PERT(var3d,var2d)
IMPLICIT NONE
REAL(r_size) :: var3d(nlon,nlat,nlev,nv3d)
REAL(r_size) :: var2d(nlon,nlat,nv2d)
INTEGER      :: iv

DO iv=1,nv3d
   if(  iv == iv3d_qv )then
      where( var3d(:,:,:,iv) < 0.0d0 )
          var3d(:,:,:,iv)=0.0d0
      endwhere
   endif
ENDDO


END SUBROUTINE CHECK_PERT


END MODULE common_perturb_ensemble
