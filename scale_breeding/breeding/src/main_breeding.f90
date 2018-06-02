PROGRAM BREEDING_RESCALE
!=======================================================================
!
! [PURPOSE:] Generate an ensemble of perturbed WRF runs.
!
! [HISTORY:]
!   12/02/2017 Juan Ruiz created
!=======================================================================
USE common
USE common_scale
USE common_namelist
USE common_breeding

 IMPLICIT NONE
 REAL(r_size),ALLOCATABLE :: v3d_c(:,:,:,:) , v3d_p(:,:,:,:) , v3d_n(:,:,:,:) 
 REAL(r_size),ALLOCATABLE :: v2d_c(:,:,:)   , v2d_p(:,:,:)   , v2d_n(:,:,:)

 REAL(r_size),ALLOCATABLE :: v3d_p_old(:,:,:,:) , v3d_n_old(:,:,:,:) !For breeding_in_place
 REAL(r_size),ALLOCATABLE :: v2d_p_old(:,:,:)   , v2d_n_old(:,:,:)   !For breeding_in_place

 REAL(r_size),ALLOCATABLE :: v3d_aux(:,:,:,:) , v2d_aux(:,:,:)

 REAL(r_size),ALLOCATABLE :: norm(:,:,:)

!Get configuration from namelist.
CALL read_namelist()



!Get model domain from the first ensemble member.
CALL set_common_scale(input_file_prefix_ctrl)

ALLOCATE( v3d_c(nlev,nlong,nlatg,nv3d) , v3d_p(nlev,nlong,nlatg,nv3d) , v3d_n(nlev,nlong,nlatg,nv3d) )
ALLOCATE( v2d_c(nlong,nlatg,nv2d)      , v2d_p(nlong,nlatg,nv2d)      , v2d_n(nlong,nlatg,nv2d)      )


ALLOCATE( norm(nlev,nlong,nlatg) )

!Reading control file

 WRITE(*,*)"Reading control file"
 call read_restart(input_file_prefix_ctrl,v3d_c,v2d_c)
 call state_trans(v3d_c)


CALL set_common_scale(input_file_prefix_pp)  !We can set a different number of processes for the ctrl data and the perturbation data.
!Only the number of subdomains can change between the control run and the perturbations. The size (vertical and horizontal) and location of the domain is assumed to 
!be the same. 

!Reading postive perturvation

 WRITE(*,*)"Reading positive perturbation"
 call read_restart(input_file_prefix_pp,v3d_p,v2d_p)
 call state_trans(v3d_p)

!Reading negative perturvation

 WRITE(*,*)"Reading negative perturbation"
 call read_restart(input_file_prefix_pn,v3d_n,v2d_n)
 call state_trans(v3d_n)

!Rescaling perturbation
 
 WRITE(*,*)"Rescaling perturbation"
 call rescale( v3d_c , v2d_c , v3d_p , v2d_p , v3d_n , v2d_n , norm )

!Writing positive perturbation

 WRITE(*,*)"Writing positive perturbation"
 call write_restart(input_file_prefix_pp,v3d_p,v2d_p)

!Writing negative perturbation

 WRITE(*,*)"Writing negative perturbation"
 call write_restart(input_file_prefix_pn,v3d_n,v2d_n)

!Write norm

 WRITE(*,*)"Writing the perturbation norm."
 call write_norm( norm )


 if ( breeding_in_place ) then
  write(*,*)"Breeding in place is active!"
 !Allocate memory space to store the perturbation at the previous time.

  ALLOCATE( v3d_p_old(nlev,nlong,nlatg,nv3d) , v3d_n_old(nlev,nlong,nlatg,nv3d)  )
  ALLOCATE( v2d_p_old(     nlong,nlatg,nv2d) , v2d_n_old(     nlong,nlatg,nv2d)  )

  ALLOCATE( v3d_aux(nlev,nlong,nlatg,nv3d) , v2d_aux(nlong,nlatg,nv2d) )
  !Reading postive perturvation

  WRITE(*,*)"Reading old positive perturbation"
  call read_restart(input_file_prefix_pp_old,v3d_p_old,v2d_p_old)

  !Reading negative perturvation

  WRITE(*,*)"Reading old negative perturbation"
  call read_restart(input_file_prefix_pn_old,v3d_n_old,v2d_n_old)

  !Remove the perturbation from the fields corresponding to the begining of the cycle 
  !and add the new rescaled perturbation obtained at the end of the cycle.
  v3d_aux = ( v3d_p_old + v3d_n_old ) * 0.5d0
  v2d_aux = ( v2d_p_old + v2d_n_old ) * 0.5d0

  v3d_p_old = v3d_aux + ( v3d_p - v3d_n ) * 0.5d0
  v3d_n_old = v3d_aux - ( v3d_p - v3d_n ) * 0.5d0

  v2d_p_old = v2d_aux + ( v2d_p - v2d_n ) * 0.5d0
  v2d_n_old = v2d_aux - ( v2d_p - v2d_n ) * 0.5d0

  !Writing positive perturbation

  WRITE(*,*)"Writing old positive perturbation"
  call write_restart(input_file_prefix_pp_old,v3d_p_old,v2d_p_old)

  !Writing negative perturbation

  WRITE(*,*)"Writing old negative perturbation"
  call write_restart(input_file_prefix_pn_old,v3d_n_old,v2d_n_old)


 endif


END PROGRAM BREEDING_RESCALE

