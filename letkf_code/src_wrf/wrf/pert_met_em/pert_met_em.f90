PROGRAM pert_met_em
!=======================================================================
!
! [PURPOSE:] Main program to compute me_em ensemble perturbations.
!
! [HISTORY:]
!   01/16/2009 Takemasa Miyoshi  created
!   10/16/2023 Juan Ruiz adapted to met em perturbation
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_met_em
  USE common_mpi_met_em
  USE met_em_tools
  USE common_namelist_met_em

  IMPLICIT NONE
  REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr , itime  , iter , ini_mem , end_mem , mem_iter , nbv_iter
  INTEGER :: m
  CHARACTER(8) :: stdoutf='NOUT-000'
  CHARACTER(4) :: prefix='____'
  CHARACTER(20):: argument
!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
  CALL CPU_TIME(rtimer00)
  CALL initialize_mpi

  !Get configuration from namelist
  CALL read_namelist()
!
  WRITE(stdoutf(6:8), '(I3.3)') myrank
!  WRITE(6,'(3A,I3.3)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I3.3,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf
!
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '  MET EM PERTURBATION TOOL                   '
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '              PARAMETERS                     '
  WRITE(6,'(A)') ' ------------------------------------------- '
  WRITE(6,'(A,I15)')   '  nbv_ori      :',nbv_ori
  WRITE(6,'(A,I15)')   '  nbv_tar      :',nbv_tar
  WRITE(6,'(A,I15)')   '  nbv_ntimes   :',ntimes  
  WRITE(6,'(A,I15)')   '  method       :',method
  WRITE(6,'(A,I15)')   '  niter        :',niter
  WRITE(6,'(A)') '============================================='
  CALL set_common_met_em('eo0100001')
  CALL set_common_mpi_met_em

  mem_iter = INT( CEILING( REAL(nbv_tar,r_sngl) / REAL(niter,r_sngl) ) )
  DO itime = 1 , ntimes 
!
    CALL CPU_TIME(rtimer)
    WRITE(6,'(A,2F10.2)') '### TIMER(INITIALIZE):',rtimer,rtimer-rtimer00
    rtimer00=rtimer

!-----------------------------------------------------------------------
! Original ensemble
!-----------------------------------------------------------------------
  !
  ! READ ORIGINAL MET EM FILES
  !

    ALLOCATE(eori3d(nij1,nlev,nbv_ori,nv3d))
    ALLOCATE(eori2d(nij1,nbv_ori,nv2d))
    ALLOCATE(eoris3d(nij1,nlev_soil,nbv_ori,ns3d))

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    WRITE(prefix(1:4),'(A,I2.2)')'eo',itime
    CALL read_ens_mpi_alltoall(prefix,nbv_ori,eori3d,eori2d,eoris3d)

    CALL CPU_TIME(rtimer)
    WRITE(6,'(A,2F10.2)') '### TIMER(READ_GUES):',rtimer,rtimer-rtimer00
    rtimer00=rtimer


  !
  ! COMPUTE ORIGINAL ENSEMBLE MEAN AND PERTURBATIONS
  !
 
    ALLOCATE(mean3dg(nij1,nlev,nv3d),mean2dg(nij1,nv2d),means3dg(nij1,nlev_soil,ns3d))
    CALL ensmean_grd(nbv_ori,nij1,eori3d,eori2d,eoris3d,mean3dg,mean2dg,means3dg)

    DO m=1,nbv_ori
     eori3d(:,:,m,:) = eori3d(:,:,m,:) - mean3dg
    END DO
    DO m=1,nbv_ori
     eori2d(:,m,:) = eori2d(:,m,:) - mean2dg
    END DO
    DO m=1,nbv_ori
     eoris3d(:,:,m,:) = eoris3d(:,:,m,:) - means3dg
    END DO

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    DO iter = 1 , niter !Loop over the iterations to compute ensemble members
       ini_mem = (iter-1) * mem_iter + 1
       end_mem = ini_mem + mem_iter -1 
       WRITE(6,*) '### INITIAL / FINAL MEMBER FOR THIS ITERATION ARE:',ini_mem,end_mem
       IF ( end_mem > nbv_tar ) THEN
          end_mem = nbv_tar
       ENDIF
       nbv_iter = end_mem - ini_mem + 1


   !-----------------------------------------------------------------------
   ! Perturbation computation
   !-----------------------------------------------------------------------
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

       ALLOCATE(eper3d(nij1,nlev,nbv_iter,nv3d))
       ALLOCATE(eper2d(nij1,nbv_iter,nv2d))
       ALLOCATE(epers3d(nij1,nlev_soil,nbv_iter,ns3d))

       CALL transform_pert(wmatrix(:,ini_mem:end_mem),nbv_iter)

       CALL CPU_TIME(rtimer)
       WRITE(6,'(A,2F10.2)') '### TIMER(PERTURBATION TRANSFORM):',rtimer,rtimer-rtimer00
       rtimer00=rtimer

  !
  ! WRITE PERTURBED METEM FILES
  !
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      WRITE(prefix(1:4),'(A,I2.2)')'ep',itime
      CALL write_ens_mpi_alltoall(prefix,ini_mem,end_mem,eper3d,eper2d,epers3d)

      CALL CPU_TIME(rtimer)
      WRITE(6,'(A,2F10.2)') '### TIMER(WRITE_ANAL_ens):',rtimer,rtimer-rtimer00
      rtimer00=rtimer

      DEALLOCATE(eper3d,eper2d,epers3d)
    END DO !End do over iterations   

    DEALLOCATE(eori3d,eori2d,eoris3d)
    DEALLOCATE(mean3dg,mean2dg,means3dg)


  END DO ! End do over times  
!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  

  WRITE(6,'(A)') ' MET EM PERTURBATION COMPLETED SUCCESSFULLY !!!!'
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi


STOP
END PROGRAM pert_met_em
