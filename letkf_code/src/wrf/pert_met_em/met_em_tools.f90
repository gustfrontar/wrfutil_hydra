MODULE met_em_tools
!=======================================================================
!
! [PURPOSE:] Module for LETKF with WRF
!
! [HISTORY:]
!   01/26/2009 Takemasa Miyoshi  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_met_em
  USE common_mpi_met_em

  IMPLICIT NONE

  PUBLIC 

  REAL(r_sngl) , ALLOCATABLE :: eori3d(:,:,:,:) , eori2d(:,:,:) , eoris3d(:,:,:,:)   !Original ensemble
  REAL(r_sngl) , ALLOCATABLE :: eper3d(:,:,:,:) , eper2d(:,:,:) , epers3d(:,:,:,:)   !Ensemble with additional perturbations

  REAL(r_sngl) , ALLOCATABLE :: mean3dg(:,:,:) , mean2dg(:,:) , means3dg(:,:,:)


  INTEGER :: iworkp0d,np0d

CONTAINS

!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
SUBROUTINE transform_pert( wmatrix_iter , nbv_iter )
  IMPLICIT NONE
  INTEGER :: m,n,iv
  INTEGER     ,INTENT(IN)  :: nbv_iter
  REAL(r_sngl),INTENT(IN)  :: wmatrix_iter(nbv_ori,nbv_iter)
  REAL(r_sngl),ALLOCATABLE :: permean3dg(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: permean2dg(:,:)
  REAL(r_sngl),ALLOCATABLE :: permeans3dg(:,:,:)

  WRITE(6,'(A)') 'Hello pert_met_em'

  !
  ! MAIN TRNSFORMATION
  !
  eper3d=0.0e0
  eper2d=0.0e0
  epers3d=0.0e0
  DO n=1,nbv_iter
    DO m=1,nbv_ori
       eper3d(:,:,n,:)  = eper3d(:,:,n,:)  + eori3d(:,:,m,:)  * wmatrix_iter(m,n)  
       eper2d(:,n,:)    = eper2d(:,n,:)    + eori2d(:,m,:)    * wmatrix_iter(m,n)
       epers3d(:,:,n,:) = epers3d(:,:,n,:) + eoris3d(:,:,m,:) * wmatrix_iter(m,n)
    END DO
  END DO 

  !
  !  RECENTER PERTURBATIONS SO THEY HAVE 0-MEAN AND ADD THE ORIGINAL ENSEMBLE MEAN
  !
  ALLOCATE(permean3dg(nij1,nlev,nv3d),permean2dg(nij1,nv2d),permeans3dg(nij1,nlev_soil,ns3d))
  CALL ensmean_grd(nbv_iter,nij1,eper3d,eper2d,epers3d,permean3dg,permean2dg,permeans3dg)
  DO m=1,nbv_iter
     eper3d(:,:,m,:) = eper3d(:,:,m,:) - permean3dg + mean3dg
  END DO
  DO m=1,nbv_iter
     eper2d(:,m,:) = eper2d(:,m,:) - permean2dg + mean2dg
  END DO
  DO m=1,nbv_iter
     epers3d(:,:,m,:) = epers3d(:,:,m,:) - permeans3dg + means3dg
  END DO

  !CONDENSATES AND QVAPOR CANNOT BE 0 AFTER ASSIMILATION.
  DO iv = 1 , nv3d
   IF( VAR_NAME_3D(iv)=='SPECHUMD'   .or. VAR_NAME_3D(iv)=='QR' .or.  &
    &    VAR_NAME_3D(iv)=='QC'       .or. VAR_NAME_3D(iv)=='QS' .or.  &
    &    VAR_NAME_3D(iv)=='QI'       .or. VAR_NAME_3D(iv)=='QG' .or.  &
    &    VAR_NAME_3D(iv)=='QH' )THEN
     WHERE( eper3d(:,:,:,iv) <= 0.0e0 )eper3d(:,:,:,iv) = 0.0e0 
   ENDIF
   IF( VAR_NAME_3D(iv)=='RH' ) THEN
     WHERE( eper3d(:,:,:,iv) <= 0.0e0 )eper3d(:,:,:,iv) = 0.0e0
     WHERE( eper3d(:,:,:,iv) >= 100.0e0 )eper3d(:,:,:,iv) = 100.0e0
   ENDIF          
  END DO

  DEALLOCATE(permean3dg,permean2dg,permeans3dg)

  RETURN
END SUBROUTINE transform_pert


END MODULE met_em_tools
