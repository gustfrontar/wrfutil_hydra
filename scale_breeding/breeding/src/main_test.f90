PROGRAM test
!=======================================================================
!
! [PURPOSE:] Generate an ensemble of perturbed WRF runs.
!
! [HISTORY:]
!   12/02/2017 Juan Ruiz created
!=======================================================================
USE common
USE common_scale

 IMPLICIT NONE
 REAL(r_size),ALLOCATABLE :: v3d_p(:,:,:,:)  , v2d_p(:,:,:)


!Get model domain from the first ensemble member.
CALL set_common_scale('ctrl')

ALLOCATE( v3d_p(nlev,nlong,nlatg,nv3d)  )
ALLOCATE( v2d_p(nlong,nlatg,nv2d)            )


!Reading postive perturvation

 WRITE(*,*)"Reading positive perturbation"
 call read_restart('pp',v3d_p,v2d_p)
 call state_trans(v3d_p)

 write(*,*)v3d_p(20,20,20,iv3d_u),v3d_p(20,20,20,iv3d_v)

 WRITE(*,*)"Writing positive perturbation"
 call write_restart('pp',v3d_p,v2d_p)

 write(*,*)v3d_p(20,20,20,iv3d_u),v3d_p(20,20,20,iv3d_v)

END PROGRAM test

