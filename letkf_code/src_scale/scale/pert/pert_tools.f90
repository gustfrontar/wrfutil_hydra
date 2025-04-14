MODULE pert_tools
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
  USE common_scale
  USE common_mpi_scale

  IMPLICIT NONE

  PUBLIC 

  REAL(r_size) , ALLOCATABLE :: eori3d(:,:,:,:) , eori2d(:,:,:) , eoris3d(:,:,:,:)   !Original ensemble
  REAL(r_size) , ALLOCATABLE :: eper3d(:,:,:,:) , eper2d(:,:,:) , epers3d(:,:,:,:)   !Ensemble with additional perturbations

  REAL(r_sngl) , ALLOCATABLE :: mean3dg(:,:,:) , mean2dg(:,:) , means3dg(:,:,:)

  REAL(r_sngl) , ALLOCATABLE :: wmatrix(:,:)

  INTEGER :: iworkp0d,np0d

  CHARACTER(vname_max),PARAMETER :: v3db_name(nv3d) = &
     (/'DENS      ', 'VELZ      ', 'VELX      ', 'VELY      ', 'PT        ', &
       'QV        ', 'QC        ', 'QR        ', 'QI        ', 'QS        ', 'QG        '/)


CONTAINS

!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
SUBROUTINE transform_pert( wmatrix_iter , nbv_iter )
  IMPLICIT NONE
  INTEGER :: m,n,iv
  INTEGER     ,INTENT(IN)  :: nbv_iter
  REAL(r_sngl),INTENT(IN)  :: wmatrix_iter(MEMBER,nbv_iter)
  REAL(r_sngl),ALLOCATABLE :: permean3dg(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: permean2dg(:,:)
  REAL(r_sngl),ALLOCATABLE :: permeans3dg(:,:,:)

  !
  ! MAIN TRNSFORMATION
  !
  eper3d=0.0e0
  eper2d=0.0e0
!!!!  epers3d=0.0e0
  DO n=1,nbv_iter
    DO m=1,MEMBER
       eper3d(:,:,n,:)  = eper3d(:,:,n,:)  + eori3d(:,:,m,:)  * wmatrix_iter(m,n)  
       eper2d(:,n,:)    = eper2d(:,n,:)    + eori2d(:,m,:)    * wmatrix_iter(m,n)
!!!!       epers3d(:,:,n,:) = epers3d(:,:,n,:) + eoris3d(:,:,m,:) * wmatrix_iter(m,n)
    END DO
  END DO 

  !
  !  RECENTER PERTURBATIONS SO THEY HAVE 0-MEAN AND ADD THE ORIGINAL ENSEMBLE MEAN
  !
!  CALL ensmean_grd(nbv_iter,nij1,eper3d,eper2d,epers3d,permean3dg,permean2dg,permeans3dg)
  CALL ensmean_grd(nbv_iter,nbv_iter+1,nij1,eper3d,eper2d)
  DO m=1,nbv_iter
     eper3d(:,:,m,:) = eper3d(:,:,m,:) - eper3d(:,:,mmean,:) + eori3d(:,:,mmean,:)
  END DO
  DO m=1,nbv_iter
     eper2d(:,m,:) = eper2d(:,m,:) - eper2d(:,mmean,:) + eori2d(:,mmean,:)
  END DO
!!  DO m=1,nbv_iter
!     epers3d(:,:,m,:) = epers3d(:,:,m,:) - permeans3dg + means3dg
!  END DO

  !CONDENSATES AND QVAPOR CANNOT BE 0 AFTER ASSIMILATION.
  DO iv = 1 , nv3d
   IF( any(trim(V3D_NAME(iv))==(/'QV','QR','QC','QS','QI','QG'/)) ) then
     WHERE( eper3d(:,:,:,iv) <= 0.0e0 ) eper3d(:,:,:,iv) = 0.0e0 
   ENDIF
  END DO

  RETURN
END SUBROUTINE transform_pert
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------

SUBROUTINE get_weights 
IMPLICIT NONE
REAL(r_size)               :: bufr( MEMBER )
INTEGER                    :: ii , jj , im

!This routine generates a transformation matrix that produce the desired perturbation.

allocate(wmatrix(MEMBER,MEMBER_TAR))

IF ( PERT_method == 1 .or. PERT_method == 2 ) THEN  !Particle rejuvenation.
   DO jj = 1 , member_tar
     IF ( PERT_method == 1 ) THEN 
        CALL com_randn( MEMBER , bufr )  !With Gaussian distribution
     ELSE
        CALL com_rand(  MEMBER , bufr )  !With Uniform distribution
     ENDIF
     wmatrix(:,jj) = REAL( bufr , r_sngl )
   ENDDO 
   wmatrix = wmatrix * ( PERT_SIGMA / REAL( MEMBER ) ) ** 0.5 
   im = 1
   DO jj = 1 , member_tar 
      wmatrix(im,jj) = wmatrix(im,jj) + 1.0e0
      im = im + 1
      IF ( im > MEMBER ) THEN
         im = 1
      ENDIF
      WRITE(6,*)'Weights for member ',jj 
      WRITE(6,*)wmatrix(:,jj) 
   ENDDO
ENDIF

IF ( PERT_method == 3 ) THEN   !Specular 

   IF ( member_tar /= 2 * MEMBER ) THEN
      WRITE(*,*)'Error: In specular method the target ensemble size should be '
      WRITE(*,*)'       twice the original ensmble size'
      WRITE(*,*)'       Original ensemble size ',member
      WRITE(*,*)'       Target ensemble size   ',member_tar
      STOP 1
   ENDIF
   wmatrix = 0.0e0 
   im = 0
   DO ii = 1 , MEMBER
      wmatrix(ii,ii)=1.0e0
      wmatrix(ii,ii+MEMBER)=-1.0e0
   ENDDO
ENDIF

END SUBROUTINE get_weights

!-------------------------------------------------------------------------------
! Write ensemble analysis data after collecting from processes
! Support extra output members 
!-------------------------------------------------------------------------------
subroutine write_ens_mpi_ext(v3d, v2d, ishift, filebase, opt_bdy)
  implicit none

  real(r_size), intent(in) :: v3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(in) :: v2d(nij1,nens,nv2d)
  integer, intent(in) :: ishift
  character(len=filelenmax),intent(in) :: filebase
  logical,intent(in) :: opt_bdy
  real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg(nlon,nlat,nv2d)
  character(len=filelenmax) :: filename
  integer :: it, im, mstart, mend

  integer :: ierr
  integer :: im_ext

  do it = 1, nitmax
    call mpi_timer('', 2, barrier=MPI_COMM_e)

    im = myrank_to_mem(it)

    mstart = 1 + (it-1)*nprocs_e
    mend = min(it*nprocs_e, nens)
    if (mstart <= mend) then
      call gather_grd_mpi_alltoall(mstart, mend, v3d, v2d, v3dg, v2dg)
    end if

    call mpi_timer('write_ens_mpi:gather_grd_mpi_alltoall:', 2)

    ! Note: write all members + mean + mdet
    ! 

    im_ext=im+ishift*MEMBER
    if (im >= 1 .and. im <= MEMBER .and. im_ext <= MEMBER_TAR) then !!! Do not overwrite mean and mdet 
      filename = filebase
      call filename_replace_mem(filename, im_ext)

!      write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',myrank_d,'.nc'
      if (.not. opt_bdy) call state_trans_inv(v3dg)

      call mpi_timer('write_ens_mpi:state_trans_inv:', 2)

      if (opt_bdy) then
        call write_boundary(filename, v3dg, v2dg)
      else
        call write_restart(filename, v3dg, v2dg)
      end if

      call mpi_timer('write_ens_mpi:write_restart:', 2)
    end if
  end do ! [ it = 1, nitmax ]

  return
end subroutine write_ens_mpi_ext
!-------------------------------------------------------------------------------
! [File I/O] Write SCALE boundary files
!-------------------------------------------------------------------------------

SUBROUTINE write_boundary(filename,v3dg,v2dg)
  use netcdf
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_HAS_W,  &
    PRC_HAS_S
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX
  use common_mpi, only: myrank
  use common_ncio
  implicit none

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(IN) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: iv3d, iv2d, ncid, varid
  integer :: is, js, ks

  integer :: istat

  is = 1
  js = 1
  if (.not. PRC_HAS_W) then
    is = is + IHALO
  end if
  if (.not. PRC_HAS_S) then
    js = js + JHALO
  end if

!  write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is writing a file ',filename,'.pe',PRC_myrank,'.nc'

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  if ( LOG_OUT) write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is writing a file ',trim(filename) // filesuffix
  call ncio_open(trim(filename) // filesuffix, NF90_WRITE, ncid)

  do iv3d = 1, nv3d
    if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
      write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3db_name(iv3d))
    end if

    if (iv3d == iv3d_rhou) then !!! VELZ is 2nd position
      ks = 2 ! ignore zh=1 (the surface where VELZ = 0.0)
    else
      ks = 1
    endif

    istat = nf90_inq_varid(ncid, trim(v3db_name(iv3d)), varid)
    if ( istat == 0 ) then
      call ncio_check(nf90_put_var(ncid, varid, v3dg(:,:,:,iv3d), &
                                   start = (/ ks, is, js, 1 /),    &
                                   count = (/ KMAX, IMAX, JMAX, 1 /)))
    else
      write(6,'(A,A15,A)') " 3D var ", trim(v3db_name(iv3d))," not found. skip!"
    end if
!    call ncio_check(nf90_inq_varid(ncid, trim(v3d_name(iv3d)), varid))
!    call ncio_check(nf90_put_var(ncid, varid, v3dg(:,:,:,iv3d), &
!                                 start = (/ ks, is, js, 1 /),    &
!                                 count = (/ KMAX, IMAX, JMAX, 1 /)))
  end do

  do iv2d = 1, nv2d !!! currently not used 
    if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
      write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
    end if

    
    istat = nf90_inq_varid(ncid, trim(v2d_name(iv2d)), varid)
    if ( istat == 0 ) then
      call ncio_check(nf90_put_var(ncid, varid, v2dg(:,:,iv2d), &
                                   start = (/ is, js, 1 /),     &
                                   count = (/ IMAX, JMAX, 1 /)))
    else
      write(6,'(A,A15,A)') " 2D var ", trim(v2d_name(iv2d))," not found. skip!"
    end if
!    call ncio_check(nf90_inq_varid(ncid, trim(v2d_name(iv2d)), varid))
!    call ncio_check(nf90_put_var(ncid, varid, v2dg(:,:,iv2d), &
!                                 start = (/ is, js, 1 /),     &
!                                 count = (/ IMAX, JMAX, 1 /)))
  end do

  call ncio_close(ncid)

  RETURN
END SUBROUTINE write_boundary

!
!-------------------------------------------------------------------------------
! Read ensemble boundary data and distribute to processes
!-------------------------------------------------------------------------------
subroutine read_ens_bdy_mpi(v3d, v2d, filebase)
  implicit none

  real(r_size), intent(out) :: v3d(nij1,nlev,nens,nv3d)
  real(r_size), intent(out) :: v2d(nij1,nens,nv2d)
  character(len=filelenmax),intent(in) :: filebase
 
  real(RP) :: v3dg(nlev,nlon,nlat,nv3d)
  real(RP) :: v2dg(nlon,nlat,nv2d)
  character(len=filelenmax) :: filename
  integer :: it, im, mstart, mend

  call mpi_timer('', 2)

  do it = 1, nitmax
    im = myrank_to_mem(it)

    ! Note: read all members + mdetin
    ! 
    if ( ( im >= 1 .and. im <= MEMBER ) .or. im == mmdetin ) then
      filename = filebase
      call filename_replace_mem(filename, im)

      call read_boundary(filename, v3dg, v2dg)

      call mpi_timer('read_ens_mpi:read_restart:', 2)

    else if (im <= nens) then ! This is to avoid the undefined value problem;
      v3dg = undef            ! it has no impact to the results
      v2dg = undef            ! 
    end if

    call mpi_timer('', 2, barrier=MPI_COMM_e)

    mstart = 1 + (it-1)*nprocs_e
    mend = min(it*nprocs_e, nens)
    if (mstart <= mend) then
      call scatter_grd_mpi_alltoall(mstart, mend, v3dg, v2dg, v3d, v2d)
    end if

    call mpi_timer('read_ens_bdy_mpi:scatter_grd_mpi_alltoall:', 2)
  end do ! [ it = 1, nitmax ]

  return
end subroutine read_ens_bdy_mpi

SUBROUTINE read_boundary(filename,v3dg,v2dg)
  use netcdf
  use scale_prc, only: &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_HAS_W,  &
    PRC_HAS_S
  use scale_atmos_grid_cartesC_index, only: &
    IHALO, JHALO, &
    IMAX, JMAX, KMAX
  use common_mpi, only: myrank
  use common_ncio
  IMPLICIT NONE

  CHARACTER(*),INTENT(IN) :: filename
  REAL(RP),INTENT(OUT) :: v3dg(nlev,nlon,nlat,nv3d)
  REAL(RP),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  character(len=12) :: filesuffix = '.pe000000.nc'
  integer :: iv3d, iv2d, ncid, varid
  integer :: is, js, ks

  integer :: istat

  is = 1
  js = 1
  if (.not. PRC_HAS_W) then
    is = is + IHALO
  end if
  if (.not. PRC_HAS_S) then
    js = js + JHALO
  end if

!  write (6,'(A,I6.6,3A,I6.6,A)') 'MYRANK ',myrank,' is reading a file ',filename,'.pe',PRC_myrank,'.nc'

  write (filesuffix(4:9),'(I6.6)') PRC_myrank
  if ( LOG_OUT ) write (6,'(A,I6.6,2A)') 'MYRANK ',myrank,' is reading a file ',trim(filename) // filesuffix
  call ncio_open(trim(filename) // filesuffix, NF90_NOWRITE, ncid)

  do iv3d = 1, nv3d
    if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
      write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3db_name(iv3d))
    end if

    if (iv3d == iv3d_rhou) then !!!! VELZ is 2nd position
      ks = 2 ! ignore zh=1 (the surface where VELZ = 0.0)
    else
      ks = 1
    endif

    istat = nf90_inq_varid(ncid, trim(v3db_name(iv3d)), varid)
    if ( istat == 0 ) then
      call ncio_check(nf90_get_var(ncid, varid, v3dg(:,:,:,iv3d), &
                                   start = (/ ks, is, js, 1 /),    &
                                   count = (/ KMAX, IMAX, JMAX, 1 /)))
    else
      write(6,'(A,A15,A)') " 3D var ", trim(v3db_name(iv3d))," not found."
!      if ( FILL_BY_ZERO_MISSING_VARIABLES ) then
!        v3dg(:,:,:,iv3d) = 0.0_RP
!      else
        stop
!      endif
    end if
!    call ncio_check(nf90_inq_varid(ncid, trim(v3d_name(iv3d)), varid))
!    call ncio_check(nf90_get_var(ncid, varid, v3dg(:,:,:,iv3d), &
!                                 start = (/ ks, is, js, 1 /),    &
!                                 count = (/ KMAX, IMAX, JMAX, 1 /)))
  end do
!!! not considered (nv2d=0)
  do iv2d = 1, nv2d
    if ( LOG_LEVEL >= 1 .and. LOG_OUT ) then
      write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
    end if
    istat = nf90_inq_varid(ncid, trim(v2d_name(iv2d)), varid)
    if ( istat == 0 ) then
    call ncio_check(nf90_get_var(ncid, varid, v2dg(:,:,iv2d), &
                                 start = (/ is, js, 1 /),     &
                                 count = (/ IMAX, JMAX, 1 /)))
    else
      write(6,'(A,A15,A)') " 2D var ", trim(v3d_name(iv3d))," not found."
!      if ( FILL_BY_ZERO_MISSING_VARIABLES ) then
!        v2dg(:,:,iv2d) = 0.0_RP
!      else
        stop
!      endif
    end if
!    call ncio_check(nf90_inq_varid(ncid, trim(v2d_name(iv2d)), varid))
!    call ncio_check(nf90_get_var(ncid, varid, v2dg(:,:,iv2d), &
!                                 start = (/ is, js, 1 /),     &
!                                 count = (/ IMAX, JMAX, 1 /)))
  end do

  call ncio_close(ncid)

  RETURN
END SUBROUTINE read_boundary


END MODULE pert_tools
