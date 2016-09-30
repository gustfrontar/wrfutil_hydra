MODULE letkf_tools
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
  USE common_wrf
  USE common_mpi_wrf
  USE common_letkf
  USE letkf_obs

  IMPLICIT NONE

  PRIVATE
  PUBLIC ::  das_letkf,gues3d,gues2d,anal3d,anal2d,guesp2d,analp2d

  INTEGER,SAVE :: nobstotal
  INTEGER,SAVE :: var_local_n2n(nv3d+nv2d)
  INTEGER,SAVE :: var_local_par_n2n(np2d)
  INTEGER,SAVE :: var_local_par_n2n_0d(np2d)
  REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:) ! background ensemble
  REAL(r_size),ALLOCATABLE :: gues2d(:,:,:)   !  output: destroyed
  REAL(r_size),ALLOCATABLE :: guesp2d(:,:,:)  !Parameters
  REAL(r_size),ALLOCATABLE :: anal3d(:,:,:,:) ! analysis ensemble
  REAL(r_size),ALLOCATABLE :: anal2d(:,:,:)   !
  REAL(r_size),ALLOCATABLE :: analp2d(:,:,:)  !Parameters
  INTEGER      , ALLOCATABLE :: np2dtop0d(:)
  INTEGER :: iworkp0d,np0d

CONTAINS

!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
SUBROUTINE das_letkf
  IMPLICIT NONE
  CHARACTER(12) :: inflfile='infl_mul.grd'
  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:)
  REAL(r_size),ALLOCATABLE :: rloc(:)
  REAL(r_size),ALLOCATABLE :: dep(:)
  REAL(r_size),ALLOCATABLE :: work3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: work2d(:,:)
  REAL(r_size),ALLOCATABLE :: workp2d(:,:)
  REAL(r_sngl),ALLOCATABLE :: workg(:,:,:)
  REAL(r_size),ALLOCATABLE :: work3da(:,:,:)     !GYL
  REAL(r_size),ALLOCATABLE :: work2da(:,:)       !GYL
  REAL(r_size),ALLOCATABLE :: logpfm(:,:),zfm(:,:)
  REAL(r_size) :: parm
  REAL(r_size) :: trans(nbv,nbv,nv3d+nv2d)
  REAL(r_size) :: trans_p(nbv,nbv,np2d)
  REAL(r_size) :: transm(nbv,nv3d+nv2d)    !GYL
  REAL(r_size) :: transrlx(nbv,nbv)        !GYL
  REAL(r_size) :: pa(nbv,nbv,nv3d+nv2d)    !GYL
  REAL(r_size) :: q_mean,q_sprd            !GYL
  REAL(r_size) :: q_anal(nbv)              !GYL

  LOGICAL :: ex
  INTEGER :: ij,ilev,n,m,i,j,k,nobsl,ierr,iv,iolen,irec
  INTEGER :: counterp
  INTEGER :: npoints

  REAL(r_size),ALLOCATABLE :: mean3dg(:,:,:),mean3da(:,:,:)
  REAL(r_size),ALLOCATABLE :: mean2dg(:,:),mean2da(:,:)
  REAL(r_size),ALLOCATABLE :: meanp2d(:,:)
  REAL(r_size),ALLOCATABLE :: sprdp2d(:,:)

  REAL(r_size)             :: tmpsprd

  WRITE(6,'(A)') 'Hello from das_letkf'
  WRITE(6,'(A,F15.2)') '  cov_infl_mul = ',cov_infl_mul, '  relax_alpha = ' &
                 ,relax_alpha , ' relax_alpha_spread =', relax_alpha_spread
  nobstotal = nobs !+ ntvs
  WRITE(6,'(A,I8)') 'Target observation numbers : NOBS=',nobs!,', NTVS=',ntvs
  !
  ! In case of no obs
  !
  !IF(nobstotal == 0) THEN
  !  WRITE(6,'(A)') 'No observation assimilated'
    !anal3d = gues3d
    !anal2d = gues2d
    !IF(ESTPAR)analp2d= guesp2d
    !RETURN
  !END IF

  !
  ! INITIALIZE PARAMETERS (FIRST CYCLE ONLY)
  !

  IF(ESTPAR)THEN
    CALL init_par_mpi(guesp2d,nbv,update_parameter_2d,update_parameter_0d,param_default_value,param_sprd_init)
  ENDIF


  !
  ! Variable localization
  !
  var_local_n2n(1) = 1
  DO n=2,nv3d+nv2d
    DO i=1,n
      var_local_n2n(n) = i
      IF(MAXVAL(ABS(var_local(i,:)-var_local(n,:))) < TINY(var_local)) EXIT
    END DO
  END DO
  IF(ESTPAR)THEN
  var_local_par_n2n(1) = 1
  DO n=2,np2d
    DO i=1,n
      var_local_par_n2n(n) = i
      IF(MAXVAL(ABS(var_local_par(i,:)-var_local_par(n,:))) < TINY(var_local_par) &
         .AND. update_parameter_2d(i) == 1 ) EXIT
    END DO
  END DO
  ENDIF
!print *,var_local_n2n
!print *,var_local_par_n2n
  !
  ! FCST PERTURBATIONS
  !
  ALLOCATE(mean3dg(nij1,nlev,nv3d),mean2dg(nij1,nv2d))
  ALLOCATE(mean3da(nij1,nlev,nv3d),mean2da(nij1,nv2d))
  IF(ESTPAR)ALLOCATE(meanp2d(nij1,np2d))
  CALL ensmean_grd(nbv,nij1,gues3d,gues2d,mean3dg,mean2dg)
  IF(ESTPAR)CALL ensmean_ngrd(nbv,nij1,guesp2d,meanp2d,np2d)
  DO n=1,nv3d
    DO m=1,nbv
      DO k=1,nlev
        DO i=1,nij1
          gues3d(i,k,m,n) = gues3d(i,k,m,n) - mean3dg(i,k,n)
        END DO
      END DO
    END DO
  END DO
  DO n=1,nv2d
    DO m=1,nbv
      DO i=1,nij1
        gues2d(i,m,n) = gues2d(i,m,n) - mean2dg(i,n)
      END DO
    END DO
  END DO


  !
  ! multiplicative inflation
  !
  IF(cov_infl_mul > 0.0d0) THEN ! fixed multiplicative inflation parameter
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    work3d = cov_infl_mul
    work2d = cov_infl_mul
  END IF
  IF(cov_infl_mul <= 0.0d0) THEN ! 3D parameter values are read-in
    ALLOCATE( workg(nlon,nlat,nlev) )
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    IF(myrank == 0) THEN
      INQUIRE(FILE=inflfile,EXIST=ex)
      IF(ex) THEN
        WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',inflfile
        INQUIRE(IOLENGTH=iolen) iolen
        OPEN(55,FILE=inflfile,FORM='unformatted',ACCESS='direct',RECL=nlon*nlat*iolen)
        !OPEN(55,FILE=inflfile,FORM='unformatted',ACCESS='sequential')
      END IF
    END IF
    irec=1
    DO iv=1,nv3d
      DO ilev=1,nlev
        IF(myrank ==0)THEN
          IF(ex)THEN
            READ(55,rec=irec)workg(:,:,ilev)
            irec=irec+1
          ELSE
            workg(:,:,ilev)=-1.0d0 * cov_infl_mul
          ENDIF
        ENDIF
       ENDDO

       CALL scatter_ngrd_mpi(0,workg,work3d(:,:,iv),nlev)
    ENDDO
    DO iv=1,nv2d
        IF(myrank ==0)THEN
          IF(ex)THEN
            READ(55,rec=irec)workg(:,:,1)
            irec=irec+1
          ELSE
            workg(:,:,1)=-1.0d0 * cov_infl_mul
          ENDIF
        ENDIF

       CALL scatter_ngrd_mpi(0,workg(:,:,1),work2d(:,iv),1)
    ENDDO
    IF( myrank == 0 )CLOSE(55)
  END IF

  !
  ! RTPS relaxation
  !
  IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN
    ALLOCATE( work3da(nij1,nlev,nv3d) )
    ALLOCATE( work2da(nij1,nv2d) )
    work3da = 1.0d0
    work2da = 1.0d0
  END IF

  !WRITE(*,*)work3d(10,10,2),work2d(10,1)
  !
  ! p_full for background ensemble mean
  !
  ALLOCATE(logpfm(nij1,nlev),zfm(nij1,nlev))
  logpfm = DLOG(mean3dg(:,:,iv3d_p))
  zfm    = mean3dg(:,:,iv3d_ph)/gg  !Compute z in meters.


  !
  ! ESTIMATE 0D PARAMETERS
  !

  !Get number of 0d estimated parameters
  np0d=SUM(update_parameter_0d)
  IF(np0d .GT. 0 .AND. ESTPAR )THEN
    CALL estimate_0d_par
  ENDIF

  !
  ! MAIN ASSIMILATION LOOP
  !
  ALLOCATE( hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
  DO ilev=1,nlev
    WRITE(6,'(A,I3)') 'ilev = ',ilev
 

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ij,n,hdxf,rdiag,rloc,dep,nobsl,parm,trans_p,trans,transm,transrlx,pa,m,k,q_mean,q_sprd,q_anal)
    DO ij=1,nij1
      DO n=1,nv3d
        IF(var_local_n2n(n) < n) THEN
          trans(:,:,n) = trans(:,:,var_local_n2n(n))
          transm(:,n) = transm(:,var_local_n2n(n))                                     !GYL
          IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                         !GYL
            pa(:,:,n) = pa(:,:,var_local_n2n(n))                                       !GYL
          END IF                                                                       !GYL
          work3d(ij,ilev,n) = work3d(ij,ilev,var_local_n2n(n))
        ELSE
          CALL obs_local(ij,ilev,n,hdxf,rdiag,rloc,dep,nobsl,logpfm,zfm)
          parm = work3d(ij,ilev,n)
          IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                         !GYL
            CALL letkf_core(nbv,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &            !GYL
                            trans(:,:,n),transm=transm(:,n),pao=pa(:,:,n),minfl=MIN_INFL_MUL) !GYL
          ELSE                                                                         !GYL
            CALL letkf_core(nbv,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &            !GYL
                            trans(:,:,n),transm=transm(:,n),minfl=MIN_INFL_MUL)        !GYL
          END IF                                                                       !GYL
          work3d(ij,ilev,n) = parm
        END IF
        IF((n == iv3d_qv .OR. n == iv3d_qc .OR. n == iv3d_qr .OR. n == iv3d_qci .OR. n == iv3d_qs .OR. n == iv3d_qg) &
           .AND. ilev > LEV_UPDATE_Q) THEN !GYL, do not update upper-level q,qc
          anal3d(ij,ilev,:,n) = mean3dg(ij,ilev,n) + gues3d(ij,ilev,:,n)                !GYL
        ELSE                                                                           !GYL
          IF(RELAX_ALPHA /= 0.0d0) THEN                                                !GYL - RTPP method (Zhang et al. 2005)
            CALL weight_RTPP(trans(:,:,n),transrlx)                                    !GYL
          ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                    !GYL - RTPS method (Whitaker and Hamill 2012)
            CALL weight_RTPS(trans(:,:,n),pa(:,:,n),gues3d(ij,ilev,:,n),transrlx,work3da(ij,ilev,n)) !GYL
          ELSE                                                                         !GYL
            transrlx = trans(:,:,n)                                                    !GYL
          END IF                                                                       !GYL
          DO m=1,nbv
            anal3d(ij,ilev,m,n) = mean3dg(ij,ilev,n)
            DO k=1,nbv
              anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &                              !GYL - sum trans and transm here
                & + gues3d(ij,ilev,k,n) * (transrlx(k,m) + transm(k,n))                !GYL
            END DO
          END DO
          END IF                                                                         !GYL
      END DO ! [ n=1,nv3d ]


      IF(ilev == 1) THEN !update 2d variable at ilev=1
        DO n=1,nv2d
          IF(var_local_n2n(nv3d+n) < nv3d+n) THEN
            trans(:,:,nv3d+n) = trans(:,:,var_local_n2n(nv3d+n))
            transm(:,nv3d+n) = transm(:,var_local_n2n(nv3d+n))                         !GYL
            IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                       !GYL
              pa(:,:,nv3d+n) = pa(:,:,var_local_n2n(nv3d+n))                           !GYL
            END IF                                                                     !GYL
            IF(var_local_n2n(nv3d+n) <= nv3d) THEN                                     !GYL - correct the bug of the 2d variable update
              work2d(ij,n) = work3d(ij,ilev,var_local_n2n(nv3d+n))                     !GYL
            ELSE                                                                       !GYL
              work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n)-nv3d)                     !GYL
            END IF                                                                     !GYL
          ELSE
            CALL obs_local(ij,ilev,n,hdxf,rdiag,rloc,dep,nobsl,logpfm,zfm)
            parm = work2d(ij,n)
            IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                       !GYL
              CALL letkf_core(nbv,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &       !GYL
                              trans(:,:,nv3d+n),transm=transm(:,nv3d+n),pao=pa(:,:,nv3d+n),minfl=MIN_INFL_MUL) !GYL
            ELSE                                                                       !GYL
              CALL letkf_core(nbv,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm, &          !GYL
                              trans(:,:,nv3d+n),transm=transm(:,nv3d+n),minfl=MIN_INFL_MUL) !GYL
            END IF                                                                     !GYL
            work2d(ij,n) = parm
          END IF
          IF(RELAX_ALPHA /= 0.0d0) THEN                                                !GYL - RTPP method (Zhang et al. 2005)
            CALL weight_RTPP(trans(:,:,nv3d+n),transrlx)                               !GYL
          ELSE IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN                                    !GYL - RTPS method (Whitaker and Hamill 2012)
            CALL weight_RTPS(trans(:,:,nv3d+n),pa(:,:,nv3d+n),gues2d(ij,:,n),transrlx,work2da(ij,n)) !GYL
          ELSE                                                                         !GYL
            transrlx = trans(:,:,nv3d+n)                                               !GYL
          END IF                                                                       !GYL
          DO m=1,nbv
            anal2d(ij,m,n) = mean2dg(ij,n)
            DO k=1,nbv
              anal2d(ij,m,n) = anal2d(ij,m,n) &                                        !GYL - sum trans and transm here
                & + gues2d(ij,k,n) * (transrlx(k,m) + transm(k,nv3d+n))                !GYL
            END DO
          END DO
        END DO ! [ n=1,nv2d ]
      END IF ! [ ilev == 1 ]


      IF(ilev == 1 .AND. ESTPAR) THEN !update 2d parameters.
        DO n=1,np2d
          !Check that we have to estimate this parameter and if the parameter is not undef.
          IF( update_parameter_2d(n) == 1 .AND. guesp2d(ij,1,n) .NE. REAL(REAL(UNDEF,r_sngl),r_size) )THEN
          !Weigths will be recomputed for the first parameter.

           IF(TRANSPAR)THEN !Transform parameter if necessary.
            CALL PARAM_TRANSFORM(guesp2d(ij,:,n),param_max_value(n),param_min_value(n),1)
           ENDIF

           CALL com_mean(nbv,guesp2d(ij,:,n),meanp2d(ij,n))  !Compute ensemble mean.
           DO m=1,nbv
            guesp2d(ij,m,n) = guesp2d(ij,m,n) - meanp2d(ij,n) !Compute ensemble perturbations.
           ENDDO


           IF(var_local_par_n2n(n) < n ) THEN
            trans_p(:,:,n) = trans_p(:,:,var_local_par_n2n(n))
           ELSE
            CALL obs_local_p2d(parameter_localization_type_2d(n),ij,ilev,n,hdxf,rdiag,rloc,dep,nobsl,logpfm,zfm)

            CALL letkf_core(nbv,nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans_p(:,:,n))
           END IF
           DO m=1,nbv
            analp2d(ij,m,n)  = meanp2d(ij,n)
            DO k=1,nbv
              analp2d(ij,m,n) = analp2d(ij,m,n) + guesp2d(ij,k,n) * trans_p(k,m,n)
            END DO
           END DO
          
           !Reconstruct parameter gues.
           DO m=1,nbv
            guesp2d(ij,m,n) = meanp2d(ij,n) + guesp2d(ij,m,n)
           END DO

           !Apply multiplicative inflation
           IF( SUM(guesp2d(ij,:,n)) .GT. 0)THEN
            CALL PARAMETER_INFLATION(analp2d(ij,:,n),guesp2d(ij,:,n),trans_p(:,:,n),parameter_fixinflation(n), &
                                     param_sprd_init(n),parameter_inflation_type)
           ENDIF
           !Apply additive inflation
           IF(ADDINFPAR)THEN
             CALL PARAMETER_INFLATION_ADDITIVE(analp2d(ij,:,n),additive_inflation_factor)
           ENDIF
           IF(TRANSPAR)THEN
            CALL PARAM_TRANSFORM(analp2d(ij,:,n),param_max_value(n),param_min_value(n),2)   !Anti transform parameter analysis.
            CALL PARAM_TRANSFORM(guesp2d(ij,:,n),param_max_value(n),param_min_value(n),2)   !Anti transform parameter gues.
           ELSE 
            !If the ensemble mean is outside the meaningful range bring it back.
            CALL com_mean(nbv,analp2d(ij,:,n),meanp2d(ij,n))  !Compute ensemble mean.
            IF(meanp2d(ij,n) > param_max_value(n))THEN
             analp2d(ij,:,n)=analp2d(ij,:,n)+(param_max_value(n)-meanp2d(ij,n))
             CALL com_mean(nbv,analp2d(ij,:,n),meanp2d(ij,n))  !Re compute ensemble mean.
            ENDIF
            IF(meanp2d(ij,n) < param_min_value(n))THEN
             analp2d(ij,:,n)=analp2d(ij,:,n)+(param_min_value(n)-meanp2d(ij,n))
             CALL com_mean(nbv,analp2d(ij,:,n),meanp2d(ij,n))  !Re compute ensemble mean.
            ENDIF
           ENDIF
          ENDIF
       ENDDO  ! [ n=1,np2d ]
      ENDIF  ! [ ilev == 1 , estpar=.true. ]

     END DO  !En do over ij
!$OMP END PARALLEL DO
  ENDDO !En do over ilev


  !CONDENSATES AND QVAPOR CANNOT BE 0 AFTER ASSIMILATION.
  WHERE( anal3d(:,:,:,iv3d_qr) <= 0.0d0)anal3d(:,:,:,iv3d_qr)=0.0d0
  WHERE( anal3d(:,:,:,iv3d_qg) <= 0.0d0)anal3d(:,:,:,iv3d_qg)=0.0d0
  WHERE( anal3d(:,:,:,iv3d_qs) <= 0.0d0)anal3d(:,:,:,iv3d_qs)=0.0d0
  WHERE( anal3d(:,:,:,iv3d_qc) <= 0.0d0)anal3d(:,:,:,iv3d_qc)=0.0d0
  WHERE( anal3d(:,:,:,iv3d_qci) <= 0.0d0)anal3d(:,:,:,iv3d_qci)=0.0d0
  WHERE( anal3d(:,:,:,iv3d_qv) <= 0.0d0)anal3d(:,:,:,iv3d_qv)=0.0d0

  DEALLOCATE(hdxf,rdiag,rloc,dep)

  IF(RELAX_ALPHA_SPREAD /= 0.0d0) THEN
    DEALLOCATE(work3da,work2da)
  END IF

  !!!!!
  !SMOOTH 2D PARAMETER UPDATE TO AVOID NOISE ACCUMULATION IN SMALL SCALES.
  !!!!!
  IF(ESTPAR)THEN
  IF(smooth_par_update_flag)THEN
  WRITE(6,*)"BEGIN PARAMETER SMOOTHING"
  !Smoothing can affect the mean or the spread of the ensemble. 
  !After smoothing the spread should be restored to their previous value.
  CALL smooth_par_update(analp2d,guesp2d,nbv,update_parameter_2d)

  WRITE(6,*)"END PARAMETER SMOOTHING"
  ENDIF 
  ENDIF

  !WRITE MULTIPLICATIVE INFLATION
  IF(cov_infl_mul < 0.0d0) THEN
    IF(myrank == 0)THEN
    INQUIRE(IOLENGTH=iolen) iolen
    OPEN(55,FILE=inflfile,FORM='unformatted',ACCESS='direct',RECL=nlon*nlat*iolen)
    !OPEN(55,FILE=inflfile,FORM='unformatted',ACCESS='sequential')
    ENDIF
    irec=1
    DO iv=1,nv3d
     !DO ilev=1,nlev
       CALL gather_ngrd_mpi(0,work3d(:,:,iv),workg,nlev)
       IF(myrank == 0) THEN
        DO ilev=1,nlev
         WRITE(55,rec=irec)workg(:,:,ilev)
         irec=irec+1
        ENDDO
       ENDIF
    ENDDO

    DO iv=1,nv2d
       CALL gather_ngrd_mpi(0,work2d(:,iv),workg(:,:,1),1)
       IF(myrank == 0) THEN
        WRITE(55,rec=irec)workg(:,:,1)
        irec=irec+1
       END IF
    ENDDO

    IF(myrank==0)CLOSE(55)
    DEALLOCATE(workg,work3d,work2d)
  END IF


  DEALLOCATE(logpfm,zfm,mean3da,mean2da,mean3dg,mean2dg)
  IF(ESTPAR)DEALLOCATE(meanp2d)
  RETURN
END SUBROUTINE das_letkf

!-----------------------------------------------------------------------
! Trnasform parameters to keep them within their physical meaningful range.
!-----------------------------------------------------------------------
SUBROUTINE PARAM_TRANSFORM(par_ens,pmax,pmin,ts)
REAL(r_size),INTENT(INOUT) ::   par_ens(nbv)
REAL(r_size),INTENT(IN)    ::   pmax,pmin
INTEGER      ::   i
REAL(r_size) ::   a,b,x
INTEGER      ::   ts  !If ts=1 transform, if ts=2 anti transform

 a=0.5d0*(pmax-pmin)
 b=0.5d0*(pmax+pmin)

 IF(ts==1)THEN
   DO i=1,nbv
    x=(par_ens(i)-b)/a
    IF(abs(x) >= 1)THEN
      WRITE(*,*)"WARNING, INVERSE HYPERBOLIC TANGENT CANNOT BE COMPUTED"
      WRITE(*,*)"PARAMETER= ",par_ens(i)," X = ",x," a = ",a, " b= ", b
      !This usually happens for boundary points.
      CYCLE
    ENDIF
      par_ens(i)=0.5*log( (x+1)/(1-x) )
   ENDDO
 ELSEIF(ts==2)THEN
   DO i=1,nbv
    x=par_ens(i)
    x=( exp(2*x) - 1 )/( exp(2*x) + 1)
    par_ens(i)=a*x+b
   ENDDO
 ELSE

 ENDIF

END SUBROUTINE PARAM_TRANSFORM

!------------------------------------------------------------------------
! PARAMETER_INFLATION: Apply different types of inflation to the parameter
! ensemble.
!------------------------------------------------------------------------

SUBROUTINE PARAMETER_INFLATION(anal_par_ens,gues_par_ens,w,fixedlambda,fix_sprd,pinf_type)

IMPLICIT NONE
INTEGER            :: n,m,k
REAL(r_size), INTENT(INOUT)       :: anal_par_ens(nbv)
REAL(r_size), INTENT(IN)          :: gues_par_ens(nbv)
REAL(r_size), INTENT(IN)          :: fixedlambda
REAL(r_size), INTENT(IN)          :: fix_sprd
INTEGER, INTENT(IN)               :: pinf_type
REAL(r_size)                      :: tmp_siga,tmp_sigg,tmp_meana,tmp_meang,lambda,fixstd
REAL(r_size), INTENT(IN)          :: w(nbv,nbv)
REAL(r_size)                      :: w_par(nbv,nbv)
REAL(r_size)                      :: w_var(nbv),w_mean(nbv)


IF(pinf_type==1)THEN
!PARAMETER INFLATION TYPE: AKSOY FIXED PARAMETER SPREAD.    
      CALL com_mean(nbv,anal_par_ens,tmp_meana)
      !CALL com_stdev(nbv,gues_par_ens,tmp_sigg)
      CALL com_stdev(nbv,anal_par_ens,tmp_siga) 
      !Only apply this inflation is the background spread is greather than 0.     
      IF(tmp_siga .GT. 0)THEN
      lambda=fix_sprd/tmp_siga
      DO m=1,nbv
        anal_par_ens(m)=(anal_par_ens(m)-tmp_meana)*lambda+tmp_meana
      END DO
      ENDIF

!FIXED PARAMETER INFLATION AS IN KW2010 
ELSEIF(pinf_type==2 )THEN
     CALL com_mean(nbv,anal_par_ens,tmp_meana)
     DO m=1,nbv
        anal_par_ens(m)=(anal_par_ens(m)-tmp_meana)*(fixedlambda)+tmp_meana
     END DO

ELSEIF(pinf_type==3)THEN

!W INFLATION TECHNIQUE FOR THE PARAMETER
! Inflate variance of W columns so the equal the varians of the 
! W corresponding to the background (i.e. the identity)
      w_par=w    !To preserve the original w value.
      w_mean=0.0d0
      DO m=1,nbv
         CALL com_mean(nbv,w_par(m,:),w_mean(m))
       DO k=1,nbv
            w_par(m,k)=w_par(m,k)-w_mean(m)
       ENDDO
      ENDDO

      !Compute W column variance
      w_var=0.0d0
      DO m=1,nbv
          CALL com_covar(nbv,w_par(:,m),w_par(:,m),w_var(m))
      ENDDO


      !Defino un coeficiente de inflacion tal que la std de las columnas de W sean iguales a las que tenia originalmente.   
      lambda=SQRT( 1 / SUM(w_var) )
      CALL com_mean(nbv,anal_par_ens,tmp_meana)
      DO m=1,nbv
        anal_par_ens(m)=(anal_par_ens(m)-tmp_meana)*lambda+tmp_meana
      ENDDO
ELSE
     WRITE(*,*) "ERROR: UNRECOGNIZED PARAMETER INFLATION TYPE!!: ",pinf_type
     STOP 99
ENDIF

RETURN
END SUBROUTINE PARAMETER_INFLATION

!------------------------------------------------------------------------
! PARAMETER ADDITIVE INFLATION
!------------------------------------------------------------------------

SUBROUTINE PARAMETER_INFLATION_ADDITIVE(anal_par_ens,additive_inflation_factor)

IMPLICIT NONE
INTEGER            :: n,m,k
REAL(r_size), INTENT(INOUT)       :: anal_par_ens(nbv)
REAL(r_size), INTENT(IN)          :: additive_inflation_factor
REAL(r_size)                      :: tmp_meana,tmp_siga
REAL(r_size)                      :: random_noise(nbv)

      CALL com_stdev(nbv,anal_par_ens,tmp_siga)  
      CALL com_randn(nbv,random_noise)     
      !APPLY ADDITIVE INFLATION.
        
      DO m=1,nbv
        anal_par_ens(m)=anal_par_ens(m)+random_noise(m)*additive_inflation_factor*tmp_siga
      END DO

END SUBROUTINE PARAMETER_INFLATION_ADDITIVE

!-----------------------------------------------------------------------
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
!-----------------------------------------------------------------------
SUBROUTINE obs_local(ij,ilev,nvar,hdxf,rdiag,rloc,dep,nobsl,logpfm,zfm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ij,ilev,nvar
  REAL(r_size),INTENT(IN) :: logpfm(nij1,nlev),zfm(nij1,nlev)
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  REAL(r_size)  ::  sigmah , sigmav
  INTEGER,INTENT(OUT) :: nobsl
  REAL(r_size) :: dist,dlev
  INTEGER,ALLOCATABLE:: nobs_use(:)
  INTEGER :: imin,imax,jmin,jmax,im,ichan,ielm
  INTEGER :: n,nn,iobs
!
! INITIALIZE
!
  
  IF( nobs > 0 ) THEN
    ALLOCATE(nobs_use(nobs))
  END IF
!
! data search
!
  imin = MAX(NINT(ri1(ij) - dist_zeroij),1)
  imax = MIN(NINT(ri1(ij) + dist_zeroij),nlon)
  jmin = MAX(NINT(rj1(ij) - dist_zeroij),1)
  jmax = MIN(NINT(rj1(ij) + dist_zeroij),nlat)

  nn = 1
  IF( nobs > 0 ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
  nn = nn-1
  IF(nn < 1) THEN
    nobsl = 0
    RETURN
  END IF
!
! CONVENTIONAL
!
  nobsl = 0
  IF(nn > 0) THEN
    DO n=1,nn
      !
      ! vertical localization
      !
      ielm=NINT(obselm(nobs_use(n)))
      IF( ielm >= id_tclon_obs) THEN !TC track obs
        dlev = 0.0d0	
      ELSE IF( ielm == id_ps_obs .AND. ilev > 1) THEN
        dlev = ABS(LOG(obsdat(nobs_use(n))) - logpfm(ij,ilev))
        IF(dlev > dist_zerov) CYCLE
      ELSE IF( ielm == id_u_obs .OR. ielm == id_v_obs .OR. ielm == id_t_obs .OR. ielm == id_q_obs & 
      .OR. ielm == id_tv_obs .OR. ielm == id_rh_obs .OR. ielm==id_co_obs .OR. ielm==id_totco_obs ) THEN
        dlev = ABS(LOG(obslev(nobs_use(n))) - logpfm(ij,ilev))
        IF(dlev > dist_zerov) CYCLE
      ELSE IF( ielm == id_reflectivity_obs .OR.   &
        ielm == id_radialwind_obs .OR. ielm == id_pseudorh_obs ) THEN
        dlev = ABS( obslev(nobs_use(n)) - zfm(ij,ilev) ) !Vertical localization in Z.
        IF(dlev > dist_zeroz) CYCLE
      ELSE
        dlev = 0.0d0
      END IF

      SELECT CASE(NINT(obselm(nobs_use(n))))
      CASE(id_u_obs,id_v_obs,id_t_obs,id_tv_obs,id_q_obs,id_rh_obs)
        sigmah=sigma_obs
        sigmav=sigma_obsv
      CASE(id_ps_obs)
        sigmah=sigma_obs
        sigmav=sigma_obsv
      CASE(id_reflectivity_obs,id_radialwind_obs,id_pseudorh_obs)
        sigmah=sigma_obs
        sigmav=sigma_obsz
      END SELECT

      

      !
      ! horizontal localization
      !
      CALL com_distll_1(obslon(nobs_use(n)),obslat(nobs_use(n)),&
        & lon1(ij),lat1(ij),dist)
      IF(dist > dist_zero) CYCLE
   
      !
      ! variable localization
      !
      CALL get_iobs( NINT(obselm(nobs_use(n))) , iobs)

      IF(var_local(nvar,iobs) < TINY(var_local)) CYCLE

      nobsl = nobsl + 1
      hdxf(nobsl,:) = obshdxf(nobs_use(n),:)
      dep(nobsl)    = obsdep(nobs_use(n))
      !
      ! Observational localization
      !

      rdiag(nobsl) = obserr(nobs_use(n))**2
      rloc(nobsl) =EXP(-0.5d0 * ((dist/sigmah)**2 + (dlev/sigmav)**2)) &
                  & * var_local(nvar,iobs)

    END DO
  END IF

  IF( nobsl > nobstotal ) THEN
    WRITE(6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL=',nobsl,' > NOBSTOTAL=',nobstotal
    WRITE(6,*) 'IJ,NN=', ij, nn
    STOP 99
  END IF
!
  IF( nobs > 0 ) THEN
    DEALLOCATE(nobs_use)
  END IF
!
  RETURN
END SUBROUTINE obs_local

!-----------------------------------------------------------------------
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
!-----------------------------------------------------------------------
SUBROUTINE obs_local_p2d(loctype,ij,ilev,nvar,hdxf,rdiag,rloc,dep,nobsl,logpfm,zfm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ij,ilev,nvar,loctype
  REAL(r_size),INTENT(IN) :: logpfm(nij1,nlev),zfm(nij1,nlev)
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  REAL(r_size)             :: sigmah , sigmav
  INTEGER,INTENT(OUT) :: nobsl
  REAL(r_size) :: dist,dlev
  INTEGER,ALLOCATABLE:: nobs_use(:)
  INTEGER :: imin,imax,jmin,jmax,im,ichan
  INTEGER :: n,nn,iobs,ielm
!
! INITIALIZE
!
  IF( nobs > 0 ) THEN
    ALLOCATE(nobs_use(nobs))
  END IF
!
! data search
!
  imin = MAX(NINT(ri1(ij) - dist_zeroij),1)
  imax = MIN(NINT(ri1(ij) + dist_zeroij),nlon)
  jmin = MAX(NINT(rj1(ij) - dist_zeroij),1)
  jmax = MIN(NINT(rj1(ij) + dist_zeroij),nlat)

  nn = 1
  IF( nobs > 0 ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
  nn = nn-1
  IF(nn < 1) THEN
    nobsl = 0
    RETURN
  END IF
!
! CONVENTIONAL
!
  
  nobsl = 0
  IF(nn > 0) THEN
    DO n=1,nn
      !
      ! vertical localization
      !
     IF ( loctype == 1 )THEN !Vertical localization will be used.

      !
      ! vertical localization
      !
      ielm=NINT(obselm(nobs_use(n)))
      IF( ielm >= id_tclon_obs) THEN !TC track obs
        dlev = 0.0d0    
      ELSE IF( ielm == id_ps_obs .AND. ilev > 1) THEN
        dlev = ABS(LOG(obsdat(nobs_use(n))) - logpfm(ij,ilev))
        IF(dist > dist_zerov) CYCLE
      ELSE IF( ielm == id_u_obs .OR. ielm == id_v_obs .OR. ielm == id_t_obs .OR. ielm == id_q_obs .OR. & 
      ielm == id_tv_obs .OR. ielm == id_rh_obs .OR. ielm==id_co_obs .OR. ielm==id_totco_obs ) THEN
        dlev = ABS(LOG(obslev(nobs_use(n))) - logpfm(ij,ilev))
        IF(dist > dist_zerov) CYCLE
      ELSE IF( ielm == id_reflectivity_obs  .OR. & 
        ielm == id_radialwind_obs .OR. ielm == id_pseudorh_obs ) THEN
        dlev = ABS( obslev(nobs_use(n)) - zfm(ij,ilev) ) !Vertical localization in Z.
        IF(dist > dist_zeroz) CYCLE
      ELSE
        dlev = 0.0d0
      END IF

     ELSEIF( loctype == 2)THEN !No vertical localization.
       dlev = 0.0d0

     ELSE
       WRITE(6,*)'WARNING: Not recognized parameter localization type.',loctype
       STOP
     ENDIF

      !
      ! horizontal localization
      !
      CALL com_distll_1(obslon(nobs_use(n)),obslat(nobs_use(n)),&
        & lon1(ij),lat1(ij),dist)


      IF(dist > dist_zero) CYCLE

      !
      ! variable localization
      !
      !Localization length can be set independently for each observation type.
      CALL get_iobs( NINT(obselm(nobs_use(n))) , iobs)
    

      IF(var_local_par(nvar,iobs) < TINY(var_local)) CYCLE

      nobsl = nobsl + 1
      hdxf(nobsl,:) = obshdxf(nobs_use(n),:)
      dep(nobsl)    = obsdep(nobs_use(n))
      !
      ! Observational localization
      !
      sigmah=sigma_obs
      sigmav=sigma_obsv

      rdiag(nobsl) = obserr(nobs_use(n))**2
      rloc(nobsl) =EXP(-0.5d0 * ((dist/sigmah)**2 + (dlev/sigmav)**2)) &
                  & * var_local(nvar,iobs)

    END DO
  END IF
!
! DEBUG
! IF( ILEV == 1 .AND. ILON == 1 ) &
! & WRITE(6,*) 'ILEV,ILON,ILAT,NN,TVNN,NOBSL=',ilev,ij,nn,tvnn,nobsl
!
  IF( nobsl > nobstotal ) THEN
    WRITE(6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL=',nobsl,' > NOBSTOTAL=',nobstotal
    WRITE(6,*) 'IJ,NN=', ij, nn
    STOP 99
  END IF
!
  IF( nobs > 0 ) THEN
    DEALLOCATE(nobs_use)
  END IF
!
  RETURN
END SUBROUTINE obs_local_p2d
!-----------------------------------------------------------------------
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
!-----------------------------------------------------------------------
SUBROUTINE obs_local_p0d(loctype,nvar,hdxf,rdiag,rloc,dep,nobsl)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nvar,loctype
  !REAL(r_size),INTENT(IN) :: logpfm(nij1,nlev)
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT) :: nobsl
  REAL(r_size) :: dist,dlev
  INTEGER,ALLOCATABLE:: nobs_use(:)
  INTEGER :: imin,imax,jmin,jmax,im,ichan
  INTEGER :: n,nn,iobs
  INTEGER :: tmpi,tmpj

nobsl = 0
  DO n=1,nobs


IF ( loctype == 1 )   THEN
      !
      ! vertical localization (localize as a surface variable)
      !
      IF(NINT(obselm(n)) >= id_tclon_obs) THEN !TC track obs
        dlev = 0.0d0
      ELSE IF(NINT(obselm(n)) == id_ps_obs ) THEN
        dlev = 0.0d0
      ELSE IF(NINT(obselm(n)) /= id_ps_obs) THEN
        dlev = ABS(LOG(obslev(n)) - LOG(100000.0e0)) 

      ELSE
        dlev = 0.0d0
      END IF
      IF(dlev > dist_zerov) CYCLE


END IF  !End for localization type 1

IF  ( loctype == 2)  THEN
            !Do not apply vertical localization

      dlev = 0.0d0

END IF


      !
      ! horizontal localization
      !
      !Discard land observations
      tmpi=INT(obsi(n))
      tmpj=INT(obsj(n))
      IF ( landmask(tmpi,tmpj) /= 0 )CYCLE
      !
      ! variable localization
      !
      CALL get_iobs( NINT(obselm(n)) , iobs)
    
      IF(var_local_par(nvar,iobs) < TINY(var_local_par)) CYCLE

      nobsl = nobsl + 1
      hdxf(nobsl,:) = obshdxf(n,:)
      dep(nobsl)    = obsdep(n)
      !
      ! Observational localization
      !
      rdiag(nobsl) = obserr(n)**2
      rloc(nobsl) =EXP( (dlev/sigma_obsv)**2) &
                  & * var_local(nvar,iobs)
    END DO
!
!
!
!
  RETURN
END SUBROUTINE obs_local_p0d

SUBROUTINE obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
  INTEGER,INTENT(INOUT) :: nn, nobs_use(nobs)
  INTEGER :: j,n,ib,ie,ip

  DO j=jmin,jmax
    IF(imin > 1) THEN
      ib = nobsgrd(imin-1,j)+1
    ELSE
      IF(j > 1) THEN
        ib = nobsgrd(nlon,j-1)+1
      ELSE
        ib = 1
      END IF
    END IF
    ie = nobsgrd(imax,j)
    n = ie - ib + 1
    IF(n == 0) CYCLE
    DO ip=ib,ie
      IF(nn > nobs) THEN
        WRITE(6,*) 'FATALERROR, NN > NOBS', NN, NOBS
      END IF
      nobs_use(nn) = ip
      nn = nn + 1
    END DO
  END DO

  RETURN

END SUBROUTINE obs_local_sub

SUBROUTINE estimate_0d_par()
IMPLICIT NONE
INTEGER   n,i,npoints,m,counterp,ierr,nobsl,k
REAL(r_size) , ALLOCATABLE :: tmpmeanp0d(:),tmpsprdp0d(:),tmpguesp0d(:,:),tmpanalp0d(:,:),workp0d(:,:),tmptrans_p(:,:,:)
REAL(r_size),ALLOCATABLE :: hdxf0d(:,:)
REAL(r_size),ALLOCATABLE :: rdiag0d(:)
REAL(r_size),ALLOCATABLE :: rloc0d(:)
REAL(r_size),ALLOCATABLE :: dep0d(:)
REAL(r_size)             :: parm


ALLOCATE(tmpmeanp0d(np0d))
ALLOCATE(tmpsprdp0d(np0d))
ALLOCATE(tmpguesp0d(nbv,np0d))
ALLOCATE(tmpanalp0d(nbv,np0d))
ALLOCATE(workp0d(nbv,np0d))
ALLOCATE(np2dtop0d(np0d))
ALLOCATE(tmptrans_p(nbv,nbv,np0d))

  !First spatially average 2d parameters.
  counterp=0
  tmpguesp0d=0.0d0
  DO n=1,np2d
   IF(update_parameter_0d(n)==1)THEN
     counterp=counterp+1
     npoints=0
     np2dtop0d(counterp)=n  !Correspondance between p0d parameters and total parameters array.
     DO i=1,nij1
      IF(landmask1(i) == 0 .AND. guesp2d(i,1,n) /= REAL(REAL(undef,r_sngl),r_size) )THEN
        DO m=1,nbv
         tmpguesp0d(m,counterp)=tmpguesp0d(m,counterp)+guesp2d(i,m,n)
        ENDDO
       npoints=npoints+1
       ENDIF
     ENDDO
   ENDIF
  ENDDO


  !All the processes will have the same spatially averaged parameters.
  workp0d=0.0d0
  CALL MPI_ALLREDUCE(tmpguesp0d,workp0d,np0d*nbv,MPI_DOUBLE_PRECISION,MPI_SUM,&
          & MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(npoints,iworkp0d,1,MPI_INTEGER,MPI_SUM,&
          & MPI_COMM_WORLD,ierr)
  !Get the parameter mean.

   !WRITE(6,*)"workp0d,iworkp0d : ",workp0d,iworkp0d
   tmpguesp0d=workp0d/iworkp0d


IF(myrank == 0 )THEN
   WRITE(6,*)"UPDATING 0D PARAMETERS"

   !Transform 0d parameters if requested
   IF(TRANSPAR)THEN
    WRITE(6,*)'Transforming parameters'
     DO n=1,np0d
     !Transform parameters using tanh function.
       CALL PARAM_TRANSFORM(tmpguesp0d(:,n),param_max_value(np2dtop0d(n)),param_min_value(np2dtop0d(n)),1)
     ENDDO
   ENDIF




  !Compute mean and perturbations.
    tmpmeanp0d=0.0d0
    DO n=1,np0d
      DO m=1,nbv
         tmpmeanp0d(n)=tmpmeanp0d(n)+tmpguesp0d(m,n)
     ENDDO
       tmpmeanp0d(n)=tmpmeanp0d(n)/REAL(nbv,r_size)
      DO m=1,nbv
       tmpguesp0d(m,n)=tmpguesp0d(m,n)-tmpmeanp0d(n)
      ENDDO
    ENDDO
   ALLOCATE(hdxf0d(nobstotal,nbv),rdiag0d(nobstotal),rloc0d(nobstotal),dep0d(nobstotal))
   DO n=1,np0d
    WRITE(6,*)"PARAMETER ENSEMBLE MEAN (GUES): ",tmpmeanp0d(n)
    WRITE(6,*)"PARAMETER ENSEMBLE (GUES):",tmpguesp0d(:,n)+tmpmeanp0d(n)
   ENDDO
      !Main estimation loop over 0d parameters.
      DO n=1,np0d
          !Observation localization for 0d parameter estimation.
          hdxf0d=0.0d0
          rdiag0d=0.0d0
          rloc0d=0.0d0
          dep0d=0.0d0
          CALL obs_local_p0d(parameter_localization_type_0d(np2dtop0d(n)),np2dtop0d(n),hdxf0d,rdiag0d,rloc0d,dep0d,nobsl)

          WRITE(6,*)nobsl,' observations are used to estimate 0d parameters'
          parm = 1.0d0
          CALL letkf_core(nbv,nobstotal,nobsl,hdxf0d,rdiag0d,rloc0d,dep0d,parm,tmptrans_p(:,:,n))


           !Compute 0d parameter analysis
           DO m=1,nbv
             tmpanalp0d(m,n)  = tmpmeanp0d(n)
              DO k=1,nbv
               tmpanalp0d(m,n) = tmpanalp0d(m,n) + tmpguesp0d(k,n) * tmptrans_p(k,m,n)
              END DO
            END DO
           !Reconstruct parameter gues.


            DO m=1,nbv
              tmpguesp0d(m,n) = tmpmeanp0d(n) + tmpguesp0d(m,n)
            END DO
            tmpmeanp0d=0.0d0
            CALL com_mean(nbv,tmpanalp0d(:,n),tmpmeanp0d(n))
           !Apply multiplicative inflation to the parameter.
           WRITE(6,*)"PARAMETER ENSEMBLE MEAN BEFORE INFLATION (ANAL): ",tmpmeanp0d(n)
           WRITE(6,*)"PARAMETER ENSEMBLE BEFORE INFLATION (ANAL): ",tmpanalp0d(:,n)
           WRITE(6,*)"PARAMETER FIX SPREAD",param_sprd_init(np2dtop0d(n))
           CALL PARAMETER_INFLATION(tmpanalp0d(:,n),tmpguesp0d(:,n),tmptrans_p(:,:,n),  &
           & parameter_fixinflation(np2dtop0d(n)),param_sprd_init(np2dtop0d(n)),parameter_inflation_type)
!
           !Apply additive inflation to the parameter.
           IF(ADDINFPAR)THEN
             CALL PARAMETER_INFLATION_ADDITIVE(tmpanalp0d(:,n),additive_inflation_factor)
           ENDIF
      ENDDO

    tmpmeanp0d=0.0d0
    DO n=1,np0d
      DO m=1,nbv
         tmpmeanp0d(n)=tmpmeanp0d(n)+tmpanalp0d(m,n)
      ENDDO
      tmpmeanp0d(n)=tmpmeanp0d(n)/REAL(nbv,r_size)
    ENDDO

  DEALLOCATE( hdxf0d,rdiag0d,rloc0d,dep0d )

    DO n=1,np0d
    WRITE(6,*)"PARAMETER ENSEMBLE MEAN (ANAL): ",tmpmeanp0d(n)
    WRITE(6,*)"PARAMETER ENSEMBLE (ANAL): ",tmpanalp0d(:,n)
    ENDDO

ENDIF  !END FOR THE IF MY_RANK == 0


  !broadcast the estimated parameters to all process.
  CALL MPI_BCAST( tmpanalp0d, nbv*np0d, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  !Convert 0d parameters to 2d parameters (only parameter value over the sea
  !will be modified.
  !over the land the default parameter value is assumed. 
  counterp=0
  DO n=1,np2d
   IF(update_parameter_0d(n)==1)THEN
     counterp=counterp+1
     DO i=1,nij1
      IF(landmask1(i) == 0)THEN
        DO m=1,nbv
         analp2d(i,m,n)=tmpanalp0d(m,counterp)
        ENDDO
       ELSE
       !This is to fix the parameter value over land for 0d estimated
       !parameters.
        analp2d(i,:,n)=param_default_value(n)
       ENDIF
     ENDDO
   ENDIF
  ENDDO


RETURN
END SUBROUTINE estimate_0d_par




!-----------------------------------------------------------------------
! Relaxation via LETKF weight - RTPP method
!-----------------------------------------------------------------------
subroutine weight_RTPP(w, wrlx)
  implicit none
  real(r_size), intent(in) :: w(nbv,nbv)
  real(r_size), intent(out) :: wrlx(nbv,nbv)
  integer :: m

  wrlx = (1.0d0 - RELAX_ALPHA) * w
  do m = 1, nbv
    wrlx(m,m) = wrlx(m,m) + RELAX_ALPHA
  end do

  return
end subroutine weight_RTPP


!-----------------------------------------------------------------------
! Relaxation via LETKF weight - RTPS method
!-----------------------------------------------------------------------
subroutine weight_RTPS(w, pa, xb, wrlx, infl)
  implicit none
  real(r_size), intent(in) :: w(nbv,nbv)
  real(r_size), intent(in) :: pa(nbv,nbv)
  real(r_size), intent(in) :: xb(nbv)
  real(r_size), intent(out) :: wrlx(nbv,nbv)
  real(r_size), intent(out) :: infl
  real(r_size) :: var_g, var_a
  integer :: m, k

  var_g = 0.0d0
  var_a = 0.0d0
  do m = 1, nbv
    var_g = var_g + xb(m) * xb(m)
    do k = 1, nbv
      var_a = var_a + xb(k) * pa(k,m) * xb(m)
    end do
  end do
  if (var_g > 0.0d0 .and. var_a > 0.0d0) then
    infl = RELAX_ALPHA_SPREAD * sqrt(var_g / (var_a * real(nbv-1,r_size))) - RELAX_ALPHA_SPREAD + 1.0d0   ! Whitaker and Hamill 2012
!    infl = sqrt(RELAX_ALPHA_SPREAD * (var_g / (var_a * real(nbv-1,r_size))) - RELAX_ALPHA_SPREAD + 1.0d0) ! Hamrud et al. 2015 (slightly modified)
    wrlx = w * infl
  else
    wrlx = w
    infl = 1.0d0
  end if

  return
end subroutine weight_RTPS





END MODULE letkf_tools
