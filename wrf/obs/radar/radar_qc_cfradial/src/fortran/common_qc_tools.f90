MODULE QC_CONST
!=======================================================================
!
! [PURPOSE:]  Constants for QC
!
! [HISTORY:]
!   09/01/2014 Juan Ruiz created
!
!=======================================================================
!-----------------------------------------------------------------------
! Variable size definitions
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
!-----------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------
  REAL(r_size),PARAMETER :: pi=3.141592653589793d0
  REAL(r_size),PARAMETER :: gg=9.81d0
  REAL(r_size),PARAMETER :: rd=287.0d0
  REAL(r_size),PARAMETER :: cp=7.0d0 / 2.0d0 * rd
  REAL(r_size),PARAMETER :: re=6371.3d3
  REAL(r_size),PARAMETER :: r_omega=7.292d-5
  REAL(r_size),PARAMETER :: t0c=273.15d0
  REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
  REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0
  REAL(r_size),PARAMETER :: clight=299792458.0d0 !Speed of light

  INTEGER      , PARAMETER :: NPAR_ECHO_TOP_3D=6 , NPAR_ECHO_TOP_2D=7 !Number of parameters in output arrays.
  REAL(r_size) , PARAMETER :: MAX_Z_ECHO_TOP=20.0d4 , MAX_R_ECHO_TOP=240.0d03
  REAL(r_size) , PARAMETER :: DZ_ECHO_TOP = 500.0d0 , DX_ECHO_TOP = 500.0d0

  INTEGER      , PARAMETER :: MAX_ECHO_TOP_LEVS=5 
  REAL(r_size) , PARAMETER :: DBZ_THRESHOLD_ECHO_TOP=5.0d0  !Echo top detection value.

  REAL(r_size)             :: undef 


!  REAL(r_size) , ALLOCATABLE  :: QCARRAY(:,:,:)  !Array to store qccodes

END MODULE QC_CONST

MODULE QC
!=======================================================================
!
! [PURPOSE:] Common quality control parameters computation
!
! [HISTORY:]
!   09/01/2014 Juan Ruiz created
!
!=======================================================================
!$USE OMP_LIB
  USE qc_const
  IMPLICIT NONE
  PUBLIC

 CONTAINS


SUBROUTINE SPECKLE_FILTER(var,na,nr,ne,nx,ny,nz,threshold,speckle)
IMPLICIT NONE
INTEGER, INTENT(IN)        :: na,nr,ne,nx,ny,nz        !Var dims and box dims
REAL(r_size),INTENT(IN)    :: threshold                !Threshold 
REAL(r_size),INTENT(IN)    :: var(na,nr,ne)            !Input variable
REAL(r_size),INTENT(OUT)   :: speckle(na,nr,ne)        !Temporal array
INTEGER                    :: ia,ir,ie 

  CALL BOX_FUNCTIONS_2D(var,na,nr,ne,nx,ny,nz,'COUN',threshold,speckle)

RETURN
END SUBROUTINE SPECKLE_FILTER

SUBROUTINE RHO_FILTER(var,na,nr,ne,nx,ny,nz,rho_smooth)
IMPLICIT NONE
INTEGER, INTENT(IN)        :: na,nr,ne,nx,ny,nz        !Var dims and box dims
!REAL(r_size),INTENT(IN)    :: threshold                !Threshold 
REAL(r_size),INTENT(INOUT) :: var(na,nr,ne)            !Input variable
REAL(r_size),INTENT(OUT)   :: rho_smooth(na,nr,ne)     !Temporal array
INTEGER                    :: ia,ir,ie

  CALL BOX_FUNCTIONS_2D(var,na,nr,ne,nx,ny,nz,'MEAN',0.0d0,rho_smooth)

!  where( rho_smooth < threshold )
!    var=undef 
!    qcarray=QCCODE_RHOFILTER
!  endwhere

RETURN

END SUBROUTINE RHO_FILTER


SUBROUTINE GET_ATTENUATION(var,na,nr,ne,beaml,cal_error,attenuation)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  This function estimates the attenuation percentaje due to metereological
!  echoes and also computes a correction.
!  Input:
!  radar an structure containing radar information.
!  reflectivity: a 3D array (possibly 4D) containing reflectivity data with
!  most of the ground clutter already removed (we will assume that the
!  reflectivity is associatedi with weather echoes only.
!  Output:
!  attenuation which is the Path Integrated Attenuation (PIA) A. Berne and R.Uijlenhoet
!  2006.
!  correction is the correction (in dbz) that has to be added to each grid
!  point. 
!  To avoid the very well known instability of the forward attenuation
!  computation algorithm the algorithm is stopped when the attenuation
!  factor reaches a certain threshold. 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IMPLICIT NONE
INTEGER     ,INTENT(IN)    :: na,nr,ne
REAL(r_size),INTENT(IN)    :: var(na,nr,ne) !Input reflectivity
REAL(r_size),INTENT(OUT)   :: attenuation(na,nr,ne) !Attenuation factor.
REAL(r_size),INTENT(IN)    :: beaml  !Beam length (m)
REAL(r_size)               :: a_coef , b_coef , c_coef , d_coef , alfa , beta  !Attenuation parameters
REAL(r_size),INTENT(IN)    :: cal_error !Calibration erro (use 1.0 if we dont know it)
REAL(r_size)               :: tmp_data_3d(na,nr,ne)
INTEGER                    :: ia,ir,ie
REAL(r_size)               :: power(2) , mean_k


a_coef=543;
b_coef=1.36;
c_coef=1.55e-3;
d_coef=1.30;

alfa=(c_coef**d_coef)/(a_coef**(d_coef/b_coef))
beta=(d_coef/b_coef)


tmp_data_3d=10.0d0**(var/10);
where( var == undef )
     tmp_data_3d=0.0d0
endwhere

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Iterative algorithm to compute attenuation. Based on the forward
! attenuation estimation of HB.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!We iterate forward in the range direction.
attenuation(:,1,:)=1.0d0

DO ia=1,na
 DO ie=1,ne
  DO ir=1,nr-1

    !First: correct Z using the current value for the attenuation.

    mean_k=1.0d-3*(alfa*( (0.5)*(tmp_data_3d(ia,ir,ie)+tmp_data_3d(ia,ir+1,ie)) )**beta ) 
    !Compute mean k between ir and ir+1 (k is dbz/m);
    attenuation(ia,ir+1,ie)=attenuation(ia,ir,ie)*exp(-0.46d0*mean_k*beaml)*cal_error

   ENDDO
 ENDDO
ENDDO

attenuation=-10*log(attenuation)  !Compute PIA

END SUBROUTINE GET_ATTENUATION


SUBROUTINE COMPUTE_TDBZ(var,na,nr,ne,nx,ny,nz,texture)
!This routine performs the radar QC computing the requested fields.
IMPLICIT NONE
INTEGER     ,INTENT(IN) :: na , nr , ne    !Grid dimension
INTEGER     ,INTENT(IN) :: nx , ny , nz  !Box dimension
REAL(r_size),INTENT(INOUT)  :: var(na,nr,ne) 
!REAL(r_size),INTENT(IN)     :: threshold
REAL(r_size),INTENT(OUT)    :: texture(na,nr,ne)
REAL(r_size)             :: tmp_data_3d(na,nr,ne) 
INTEGER                  :: ii , jj , kk

!Compute the difference along the radial direction.
tmp_data_3d=undef

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ii,jj,kk)
 DO ii = 1,na
   DO jj = 1 ,nr-1
     DO kk = 1, ne
        IF( var(ii,jj,kk) /= undef .AND. var(ii,jj+1,kk) /= undef)THEN
          tmp_data_3d(ii,jj,kk) = ( var(ii,jj+1,kk)-var(ii,jj,kk) )**2
        ENDIF
     ENDDO
   ENDDO
 ENDDO 
!$OMP END PARALLEL DO

 !Average the squared radial differences.
 CALL BOX_FUNCTIONS_2D(tmp_data_3d,na,nr,ne,nx,ny,nz,'MEAN',0.0d0,texture)

! where( texture > threshold .or. texture == undef )
!      var=undef
!      qcarray=QCCODE_TEXTURE
! endwhere

RETURN
END SUBROUTINE COMPUTE_TDBZ

SUBROUTINE COMPUTE_SIGN(var,na,nr,ne,nx,ny,nz,varsign)
!This routine computes the sign parameter
!Kessinger et al 2003
IMPLICIT NONE
INTEGER     ,INTENT(IN)     :: na , nr , ne    !Grid dimension
INTEGER     ,INTENT(IN)     :: nx , ny , nz  !Box dimension
REAL(r_size),INTENT(INOUT)  :: var(na,nr,ne)
!REAL(r_size),INTENT(IN)     :: threshold 
REAL(r_size)                :: tmp_data_3d(na,nr,ne) , diff
INTEGER                     :: ii , jj , kk
REAL(r_size),INTENT(OUT)    :: varsign(na,nr,ne)

!Compute the difference along the radial direction.
tmp_data_3d=undef

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ii,jj,kk,diff)
 DO ii = 1,na
   DO jj = 1 ,nr-1
     DO kk = 1, ne
        IF( var(ii,jj,kk) /= undef .AND. var(ii,jj+1,kk) /= undef )THEN
           diff= var(ii,jj,kk) - var(ii,jj+1,kk) 
           IF( ABS(diff) > 0 )THEN
             tmp_data_3d(ii,jj,kk) = diff / ABS(diff)
           ELSE
             tmp_data_3d(ii,jj,kk) = 0.0d0
           ENDIF
        ENDIF
     ENDDO
   ENDDO
 ENDDO
!$OMP END PARALLEL DO

 !Average the squared radial differences.
 CALL BOX_FUNCTIONS_2D(tmp_data_3d,na,nr,ne,nx,ny,nz,'MEAN',0.0d0,varsign)

! where( varsign > threshold .or. varsign == undef )
!    var=undef
!    qcarray=QCCODE_SIGN
! endwhere 


RETURN
END SUBROUTINE COMPUTE_SIGN

SUBROUTINE BOX_FUNCTIONS_2D(datain,na,nr,ne,boxx,boxy,boxz,operation,threshold,dataout)

IMPLICIT NONE
INTEGER     ,INTENT(IN) :: na , nr , ne    !Grid dimension
INTEGER     ,INTENT(IN) :: boxx,boxy,boxz  !Box dimension
REAL(r_size),INTENT(IN) :: datain(na,nr,ne)
CHARACTER(4),INTENT(IN) :: operation     
REAL(r_size),INTENT(IN) :: threshold
REAL(r_size),INTENT(OUT) :: dataout(na,nr,ne) !Result
REAL(r_size),ALLOCATABLE :: tmp_field(:) 
REAL(r_size)             :: tmp_mean , tmp_var
INTEGER                  :: NITEMS
INTEGER                  :: ii , jj , kk , bii , bjj , bkk , box_size , iin ,ii_index , data_count
dataout=UNDEF

box_size=(2*boxx+1)*(2*boxy+1)*(2*boxz+1);

ALLOCATE( tmp_field(box_size) )

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(kk,ii,jj,bkk,bii,bjj,ii_index,NITEMS,tmp_field,tmp_mean,tmp_var,data_count,iin)
DO kk=1,ne
 DO ii=1,na
  DO jj=1,nr
     NITEMS=0
     tmp_field=0.0d0
     DO bkk=kk-boxz,kk+boxz
       DO bii=ii-boxx,ii+boxx
         DO bjj=jj-boxy,jj+boxy
            IF( bkk >= 1 .AND. bkk <= ne .AND. bjj >= 1 .AND. bjj <= nr )THEN
              !Boundary condition in X
              ii_index=bii
              IF( bii < 1 )ii_index=bii+na
              IF( bii > na)ii_index=bii-na 
              IF( OPERATION == 'MEA2' .AND. bii== ii .AND. bjj==jj .AND. bkk==kk)CYCLE !We will not count the center of the box.
                NITEMS=NITEMS+1
                tmp_field(NITEMS)=datain(ii_index,bjj,bkk)
            ENDIF
         ENDDO
       ENDDO
      ENDDO 
      !Perform the operation and save the result in dataout
      IF( OPERATION .EQ. 'MEAN' .OR. OPERATION .EQ. 'MEA2' )THEN
        !Undef values won't be considered.
        data_count=0
        tmp_mean=0.0d0
        DO iin=1,NITEMS
          IF( tmp_field(iin) /= UNDEF )THEN
            data_count=data_count+1
            tmp_mean=tmp_mean+tmp_field(iin)
          ENDIF
        ENDDO
        IF( data_count .GT. 0)THEN
          dataout(ii,jj,kk)=tmp_mean/REAL(data_count,r_size)
        ENDIF
      ELSEIF( OPERATION == 'SIGM')THEN
        !Undef values won't be considered.
        tmp_mean=0.0d0
        tmp_var=0.0d0
        data_count=0
        DO iin=1,NITEMS
          IF( tmp_field(iin) .ne. UNDEF )THEN
            data_count=data_count+1
            tmp_mean=tmp_mean+tmp_field(iin)
            tmp_var=tmp_var+tmp_field(iin) ** 2
          ENDIF
        ENDDO
        IF( data_count .GT. 0)THEN
         tmp_mean=tmp_mean/REAL(data_count,r_size)
         tmp_var=tmp_var/REAL(data_count,r_size)
         dataout(ii,jj,kk)=SQRT(tmp_var - tmp_mean**2 )
        ENDIF

      ELSEIF( OPERATION == 'COUN')THEN
        !Count values over a certain threshold (note that undef values will be 
        !always below the threshold.
        tmp_mean=0.0d0
        DO iin=1,NITEMS
          IF( tmp_field(iin)  >= threshold .AND. tmp_field(iin) /= UNDEF )THEN
            tmp_mean=tmp_mean+1.0d0
          ENDIF
        ENDDO
        IF( NITEMS > 0)dataout(ii,jj,kk)=tmp_mean/REAL(NITEMS,r_size)

      ELSEIF( OPERATION == 'MAXN')THEN
        !Local maximum tacking care of undef values.
        tmp_mean=UNDEF
        DO iin=1,NITEMS
          IF( tmp_mean == UNDEF .AND. tmp_field(iin) /= UNDEF )tmp_mean=tmp_field(iin)
          IF( tmp_field(iin)  > tmp_mean .AND. tmp_mean /= UNDEF .AND. tmp_field(iin) /= UNDEF )THEN
            tmp_mean = tmp_field(iin)
          ENDIF
        ENDDO
        dataout(ii,jj,kk)=tmp_mean
      ELSEIF( OPERATION == 'MINN')THEN
        !Local maximum tacking care of undef values.
        tmp_mean=UNDEF
        DO iin=1,NITEMS
          IF( tmp_mean == UNDEF .AND. tmp_field(iin) /= UNDEF )tmp_mean=tmp_field(iin)
          IF( tmp_field(iin)  < tmp_mean .AND. tmp_mean /= UNDEF .AND. tmp_field(iin) /= UNDEF )THEN
            tmp_mean = tmp_field(iin)
          ENDIF
        ENDDO
        dataout(ii,jj,kk)=tmp_mean
      ENDIF
  ENDDO
 ENDDO
ENDDO 
!$OMP END PARALLEL DO

RETURN
END SUBROUTINE BOX_FUNCTIONS_2D

SUBROUTINE  ECHO_TOP(reflectivity,heigth,rrange,na,nr,ne,nx,ny,nz,output_data_3d,output_data_2d)
!Curretnly this routine:
!Compute 3D echo top, echo base , echo depth , max dbz and max dbz z
!Performs interpolation from original radar grid (r,elevation) to an uniform (r,z) grid, where
!the parameters are computed.
!Then the result is interpolated back to the original radar grid.

use qc_const


IMPLICIT NONE
INTEGER     ,INTENT(IN)  :: na,nr,ne
INTEGER     ,INTENT(IN)  :: nx,ny,nz
REAL(r_size),INTENT(IN)  :: reflectivity(na,nr,ne) , heigth(nr,ne) , rrange(nr,ne)  
REAL(r_size),INTENT(OUT) :: output_data_3d(na,nr,ne,NPAR_ECHO_TOP_3D)  !Echo top , echo base , echo depth , max_dbz , maz_dbz_z , vertical_z_gradient
REAL(r_size),INTENT(OUT) :: output_data_2d(na,nr,NPAR_ECHO_TOP_2D)  !Max echo top, max_echo_base, max_echo_depth, col_max, height weighted col_max
!REAL(r_size)             :: tmp_output_data_3d(na,nr,ne,NPAR_ECHO_TOP_3D)
!-----> The following variables will be saved within calls to speed up the computation.
REAL(r_size), ALLOCATABLE,SAVE :: Z(:,:) , R(:,:) 
INTEGER, ALLOCATABLE,SAVE      :: REGJ(:,:,:) , REGI(:,:,:) , INVI(:,:,:) , INVJ(:,:,:) 
INTEGER, ALLOCATABLE,SAVE      :: NEARESTN(:,:) , INVNEARESTN(:,:)
REAL(r_size),ALLOCATABLE,SAVE  :: W(:,:,:),INVW(:,:,:)
INTEGER,SAVE                   :: REGNZ , REGNR
LOGICAL,SAVE                   :: INITIALIZED=.FALSE.
!----->
REAL(r_size), ALLOCATABLE      :: REGREF(:,:)
CHARACTER(4)                   :: OPERATION='MEAN'
INTEGER                        :: i, ii , jj , ia , ip 
REAL(r_size),ALLOCATABLE       :: tmp_data3d(:,:,:,:) , tmp_data3d_2(:,:,:,:) , tmp_data2d(:,:,:), tmp_data2d_2(:,:,:)


!WRITE(6,*)'HELLO FROM COMPUTE_ECHO_TOP'
IF( .NOT. INITIALIZED) THEN
!Perform this part only in the first call.

REGNZ=INT(MAX_Z_ECHO_TOP / DZ_ECHO_TOP)+1
REGNR=INT(MAX_R_ECHO_TOP / DX_ECHO_TOP)+1
!WRITE(6,*)'REGNZ = ',REGNZ,' REGNR = ',REGNR


  ALLOCATE( Z(REGNR,REGNZ) , R(REGNR,REGNZ) )
  ALLOCATE( REGI(REGNR,REGNZ,4) , REGJ(REGNR,REGNZ,4), NEARESTN(REGNR,REGNZ) )
  ALLOCATE( INVI(nr,ne,4) , INVJ(nr,ne,4), INVNEARESTN(nr,ne) )
  ALLOCATE( W(REGNR,REGNZ,4),INVW(nr,ne,4) )
  !Set interpolation from range-elevation to range-z grid

  !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ii,jj)
  DO ii=1,REGNR
   DO jj=1,REGNZ
     Z(ii,jj)=REAL(jj-1,r_size)*DZ_ECHO_TOP
     R(ii,jj)=REAL(ii-1,r_size)*DX_ECHO_TOP
     CALL com_xy2ij(nr,ne,rrange,heigth,R(ii,jj),Z(ii,jj),REGI(ii,jj,:),REGJ(ii,jj,:),W(ii,jj,:),NEARESTN(ii,jj))
   ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  !Set interpolation from range-z to range-elevation grid.

  !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ii,jj)
  DO ii=1,nr
   DO jj=1,ne
     CALL com_xy2ij(REGNR,REGNZ,R,Z,rrange(ii,jj),heigth(ii,jj),INVI(ii,jj,:),INVJ(ii,jj,:),INVW(ii,jj,:),INVNEARESTN(ii,jj)) 
   ENDDO
  ENDDO
  !$OMP END PARALLEL DO

ENDIF !End of first call only section.

ALLOCATE( tmp_data3d(na,REGNR,REGNZ,NPAR_ECHO_TOP_3D))
ALLOCATE( tmp_data3d_2(na,REGNR,REGNZ,NPAR_ECHO_TOP_3D))
ALLOCATE( tmp_data2d(na,REGNR,NPAR_ECHO_TOP_2D))
ALLOCATE( tmp_data2d_2(na,REGNR,NPAR_ECHO_TOP_2D))

ALLOCATE( REGREF(REGNR,REGNZ) )


tmp_data3d=UNDEF
tmp_data3d_2=UNDEF
!tmp_output_data_3d=UNDEF

output_data_3d=UNDEF

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ia,ii,jj,REGREF)
DO ia=1,na
 REGREF=UNDEF
 !Interp reflectivity from elevation-range grid to z-range grid. (nearest neighbor)
 DO ii=1,REGNR
  DO jj=1,REGNZ 
     IF( NEARESTN(ii,jj) > 0)THEN
        REGREF(ii,jj)=reflectivity(ia,REGI(ii,jj,NEARESTN(ii,jj)),REGJ(ii,jj,NEARESTN(ii,jj)))
     ENDIF
  ENDDO

  CALL ECHO_TOP_SUB(REGREF(ii,:),Z(ii,:),REGNZ,tmp_data3d(ia,ii,:,:),tmp_data2d(ia,ii,:),MAX_ECHO_TOP_LEVS,DBZ_THRESHOLD_ECHO_TOP)

 ENDDO

 !----> DEBUG
 !IF(ia == 67)THEN
 !WRITE(*,*)REGNR,REGNZ
 !OPEN(34,FILE='slice.grd',FORM='unformatted',access='sequential')
 !WRITE(34)REAL(tmp_data3d(ia,:,:,6),r_sngl)
 !WRITE(34)REAL(REGREF(:,:),r_sngl)
 !ENDIF
 !----> DEBUG


ENDDO
!$OMP END PARALLEL DO


!DO ip=1,NPAR_ECHO_TOP_3D
CALL BOX_FUNCTIONS_2D(tmp_data3d(:,:,:,1),na,REGNR,REGNZ,nx,ny,nz,'MEAN',0.0d0,tmp_data3d_2(:,:,:,1))
CALL BOX_FUNCTIONS_2D(tmp_data3d(:,:,:,2),na,REGNR,REGNZ,nx,ny,nz,'MEAN',0.0d0,tmp_data3d_2(:,:,:,2))
CALL BOX_FUNCTIONS_2D(tmp_data3d(:,:,:,3),na,REGNR,REGNZ,nx,ny,nz,'MEAN',0.0d0,tmp_data3d_2(:,:,:,3))
CALL BOX_FUNCTIONS_2D(tmp_data3d(:,:,:,4),na,REGNR,REGNZ,nx,ny,nz,'MEAN',0.0d0,tmp_data3d_2(:,:,:,4))
CALL BOX_FUNCTIONS_2D(tmp_data3d(:,:,:,5),na,REGNR,REGNZ,0,1,0,'MINN',0.0d0,tmp_data3d_2(:,:,:,5))
CALL BOX_FUNCTIONS_2D(tmp_data3d(:,:,:,6),na,REGNR,REGNZ,0,1,0,'MINN',0.0d0,tmp_data3d_2(:,:,:,6))
!END DO

DO ip=1,NPAR_ECHO_TOP_2D
CALL BOX_FUNCTIONS_2D(tmp_data2d(:,:,ip),na,REGNR,1,nx,ny,0,'MEAN',0.0d0,tmp_data2d_2(:,:,ip))
END DO


!Interpolate back to the elevation-range grid. (Using nearest neighbor)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ia,ii,jj)
DO ia=1,na
 DO ii=1,nr
  DO jj=1,ne
    IF( INVNEARESTN(ii,jj) .GT. 0 )THEN
      output_data_3d(ia,ii,jj,:)=tmp_data3d_2(ia,INVI(ii,jj,INVNEARESTN(ii,jj)),INVJ(ii,jj,INVNEARESTN(ii,jj)),:)
    ENDIF
    IF( output_data_3d(ia,ii,jj,5) /= UNDEF)THEN
      output_data_3d(ia,ii,jj,5)=output_data_3d(ia,ii,jj,5) !-topography(ia,ii,jj)  (dejo para mas adelante)
    ENDIF
  ENDDO
    IF( INVNEARESTN(ii,1) .GT. 0 )THEN !We interpolate the data to the lowest level.
      output_data_2d(ia,ii,:)=tmp_data2d_2(ia,INVI(ii,1,INVNEARESTN(ii,1)),:)
    ENDIF
 ENDDO

 !We will use the vertical reflectivity gradient and maximum ref height over the terrain only for heights between 0 and 2 km
 !WHERE( heigth - topography(ia,:,:) < 0 .OR. heigth - topography(ia,:,:) > 2000 )
 !     output_data_3d(ia,:,:,6)=UNDEF
 !     output_data_3d(ia,:,:,5)=UNDEF
 !ENDWHERE


ENDDO   
!$OMP END PARALLEL DO

DEALLOCATE( tmp_data3d_2 , tmp_data3d )
DEALLOCATE( REGREF )


INITIALIZED=.TRUE.
RETURN
END SUBROUTINE ECHO_TOP

SUBROUTINE ECHO_TOP_SUB(reflectivity,z,nz,output_3d,output_2d,max_levs,threshold)
!Vertical columns calculations
!Compute the possition of multiple echo tops in a single reflectivity column.
!Compute echo depth of each echo layer
!compute echo base
!compute max dbz
!This routine returns a vertical profile of echo base, echo top , echo depth , max dbz and max dbz a fore 
!each cloud layer.
!It also returns the vertical profile of the vertical gradient of the reflectivity field.

use qc_const

IMPLICIT NONE
INTEGER, INTENT(IN)  :: nz , max_levs
REAL(r_size),INTENT(IN) :: reflectivity(nz) , z(nz)
REAL(r_size),INTENT(OUT):: output_3d(nz,NPAR_ECHO_TOP_3D) !echo_top, echo_base, echo_depth , max_dbz , max_dbz_z , reflectivity gradient
REAL(r_size),INTENT(OUT):: output_2d(NPAR_ECHO_TOP_2D) !Max echo top, max_echo_base, max_echo_depth, col_max, height weighted col_max , first ref maximum height , intensity.

REAL(r_size),INTENT(IN) :: threshold    !Reflectivity threshold to detect echo top.
INTEGER, PARAMETER      :: Nlevelstop=2
REAL(r_size)            :: tmp(max_levs,5) !echo_top, echo_base, echo_depth , max_dbz , max_dbz_z
REAL(r_size)            :: ref(nz) , ave_ref , sum_z
INTEGER                 :: jj, iz , base_count , top_count , tmp_count , itop , imax
LOGICAL                 :: base_detected , top_detected
LOGICAL                 :: found_first_maximum
REAL(r_size), PARAMETER :: first_maximum_threshold = 10.0d0 
INTEGER     , PARAMETER :: NDELTAZ=5      ! NDELTAZ * dz is the distance used to estimate vertical reflectivity gradient
REAL(r_size), PARAMETER :: refmin =0.0d0  ! Reflectivity value that will be assumed for UNDEF values in gradient computation.
                                          

output_3d=UNDEF
output_2d=UNDEF
tmp=UNDEF

base_count=0
top_count=0

ref=reflectivity   !reflectivity is intent in.


base_detected=.false.
top_detected=.false.

!Before computation extend data one or to levels below the first echo. This is done to prevent the first level to fall outside the computation
!of these scores.
DO iz=1,nz
   IF( z(iz) > 3000 )EXIT

   IF( ref(iz) /= UNDEF )THEN
      IF( iz>= 2)THEN
       ref(iz-1)=ref(iz)
      ENDIF

     EXIT
   ENDIF
ENDDO


DO iz=1,nz
   !Look for an echo base
   IF( ref(iz) > threshold .AND.  .NOT. base_detected .AND. ref(iz) /= UNDEF )THEN
       !An echo base has been detected.
       IF( base_count < max_levs)THEN
       base_detected=.true.
       top_detected=.false.
       base_count=base_count+1
       tmp(base_count,2)=z(iz)   !Echo base
       tmp(base_count,4)=ref(iz) !Max dbz
       tmp(base_count,5)=z(iz)   !Max dbz_z
       ENDIF
   ENDIF
   !Look for an echo top.
   IF( iz > Nlevelstop )THEN
     tmp_count=0
     DO jj=iz-Nlevelstop+1,iz
        IF( ref(jj) < threshold .OR. ref(jj) == UNDEF )tmp_count=tmp_count+1
     ENDDO
     IF( tmp_count == Nlevelstop .AND. .NOT. top_detected .AND. base_detected )THEN
     !An echo top has been detected
        top_detected=.true.
        base_detected=.false.
        IF( base_count <= max_levs )THEN
           tmp(base_count,1)=z(iz-Nlevelstop)  !Echo top
        ENDIF
     ENDIF
   ENDIF
   !Echo top associated with top of the radar domain.
   IF( iz == nz .AND. base_detected .AND. .NOT. top_detected )THEN
   !Domain is over but echo top has not been found! :( 
   !Force echo top
       IF( base_count <= max_levs )THEN
           tmp(base_count,1)=z(iz)  !Echo top
       ENDIF
   ENDIF
   !Compute max dbz
   IF( base_detected .AND. .NOT. top_detected .AND. ref(iz) /=UNDEF )THEN
       !We are within a cloud or an echo region. Compute max dbz.
       IF( ref(iz) > tmp(base_count,4) )THEN  !Max dbz
           tmp(base_count,4)=ref(iz)  !Max dbz
           tmp(base_count,5)=z(iz)    !Max dbz z
       ENDIF
   ENDIF
   !Compute vertical gradient of reflectivity.
   IF( iz <= nz-NDELTAZ)THEN
   IF( ref(iz) /= UNDEF )THEN
    IF(  ref( iz + NDELTAZ ) /= UNDEF )THEN
        output_3d(iz,6)= ( ref(iz+NDELTAZ) - ref(iz) ) /( z(iz+NDELTAZ) - z(iz) )
    ELSE
        output_3d(iz,6)= ( refmin          - ref(iz) ) /( z(iz+NDELTAZ) - z(iz) )
    ENDIF
   ENDIF
   ENDIF


ENDDO !End for loop over levels

DO itop=1,max_levs
   IF( tmp(itop,1) .NE. UNDEF  .AND.  tmp(itop,2) .NE. UNDEF )THEN  !Echo top and echo base
       DO iz=1,nz-1
          IF( z(iz) >= tmp(itop,2) .AND. z(iz) <= tmp(itop,1))THEN
               output_3d(iz,1:2)=tmp(itop,1:2)
               output_3d(iz,3)  =tmp(itop,1)-tmp(itop,2)
               output_3d(iz,4:5)=tmp(itop,4:5)
          ENDIF
          IF( z(iz) > tmp(itop,1) )EXIT
       ENDDO
   ENDIF
   !Find maximum echo top
   IF ( tmp(itop,1) .NE. UNDEF )THEN
      IF( output_2d(1) .EQ. UNDEF )THEN
        output_2d(1)=tmp(itop,1)
        ELSE
         IF( tmp(itop,1) >= output_2d(1) )THEN
           output_2d(1) = tmp(itop,1)
         ENDIF
      ENDIF
   ENDIF

   !Find maximum echo base
   IF ( tmp(itop,2) .NE. UNDEF )THEN
      IF( output_2d(2) .EQ. UNDEF )THEN
        output_2d(2)=tmp(itop,2)
        ELSE
         IF( tmp(itop,2) >= output_2d(2) )THEN
           output_2d(2) = tmp(itop,2)
         ENDIF
      ENDIF
   ENDIF

   !Find maximum echo depth
   IF ( tmp(itop,3) .NE. UNDEF )THEN
      IF( output_2d(3) .EQ. UNDEF )THEN
        output_2d(3)=tmp(itop,3)
        ELSE
         IF( tmp(itop,3) >= output_2d(3) )THEN
           output_2d(3) = tmp(itop,3)
         ENDIF
      ENDIF
   ENDIF

   !Find maximum reflectivity (colmax)
   IF ( tmp(itop,4) .NE. UNDEF )THEN
      IF( output_2d(4) .EQ. UNDEF )THEN
        output_2d(4)=tmp(itop,4)
        ELSE
         IF( tmp(itop,4) >= output_2d(4) )THEN
           output_2d(4) = tmp(itop,4)
         ENDIF
      ENDIF
   ENDIF

   IF ( tmp(itop,4) .NE. UNDEF )THEN
      IF( output_2d(4) .EQ. UNDEF )THEN
        output_2d(4)=tmp(itop,4)
        ELSE
         IF( tmp(itop,4) >= output_2d(4) )THEN
           output_2d(4) = tmp(itop,4)
         ENDIF
      ENDIF
   ENDIF

ENDDO


!Compute heigh weigthed averaged reflectivity, the height of the first reflectivity maximum and its intensity.


ave_ref=0
sum_z=0
found_first_maximum=.FALSE.
imax=0
DO iz=1,nz

 IF( reflectivity(iz) .NE. UNDEF )THEN
   ave_ref = z(iz) * reflectivity(iz)
   sum_z   = z(iz)
 ENDIF

 IF( reflectivity(iz) /= UNDEF .AND. output_2d(7) == UNDEF )THEN
   output_2d(7) = reflectivity(iz)    !Intensity of first reflectivity maximun
   output_2d(6) = z(iz)               !Height of first reflectivity maximum
 ENDIF
 IF( reflectivity(iz) /= UNDEF .AND. (reflectivity(iz) - output_2d(7)) <  & 
     first_maximum_threshold .AND. .NOT. found_first_maximum )THEN
   found_first_maximum=.TRUE. 
 ELSE
   IF( reflectivity(iz) > output_2d(6) )THEN
     output_2d(7) = reflectivity(iz) !Keep updating the maximum until we reach the first maximum.
     output_2d(6) = z(iz)            !Keep updating the height of the maximum 
   ENDIF
 ENDIF

ENDDO
IF( sum_z .GT. 0 )THEN
  output_2d(5) = ave_ref / sum_z
ELSE
  output_2d(5) = 0
ENDIF

!Max echo top, max_echo_base, max_echo_depth, col_max, height weighted col_max
RETURN
END SUBROUTINE ECHO_TOP_SUB


!-----------------------------------------------------------------------
! (X,Y) --> (i,j) conversion (General pourpuse interpolation)
!   [ORIGINAL AUTHOR:] Masaru Kunii
!-----------------------------------------------------------------------
SUBROUTINE com_xy2ij(nx,ny,fx,fy,datax,datay,dist_min_x,dist_min_y,ratio,nearestn)

  use qc_const

  IMPLICIT NONE
  ! --- inout variables
  INTEGER,INTENT(IN) :: nx,ny !number of grid points
  REAL(r_size),INTENT(IN) :: fx(nx,ny),fy(nx,ny) !(x,y) at (i,j)
  REAL(r_size),INTENT(IN) :: datax,datay !target (lon,lat)
  ! --- local work variables
  LOGICAL,PARAMETER :: detailout = .FALSE.
  INTEGER,PARAMETER :: num_grid_ave = 4  ! fix
  INTEGER :: ix,jy,ip,wk_maxp
  INTEGER :: iorder_we,iorder_sn
  INTEGER :: nxp,nyp
  INTEGER,PARAMETER :: order = 2
  REAL(r_size),PARAMETER :: max_dist = 2.0e+6
  REAL(r_size) :: rxmax, rxmin, rymax, rymin   
  REAL(r_size) :: dist(num_grid_ave)  , tmp_dist(num_grid_ave)
  INTEGER,INTENT(OUT) :: dist_min_x( num_grid_ave)
  INTEGER,INTENT(OUT) :: dist_min_y( num_grid_ave) 
  INTEGER,INTENT(OUT) :: nearestn(1)
  REAL(r_size) :: wk_dist, sum_dist
  REAL(r_size),INTENT(OUT) :: ratio(num_grid_ave)
  
  IF(detailout) THEN
    WRITE(6,'(A)') '====================================================='
    WRITE(6,'(A)') '      Detailed output of SUBROUTINE com_pos2ij       '
    WRITE(6,'(A)') '====================================================='    
  END IF
  ! ================================================================
  !   Check the Order of fx,fy 
  ! ================================================================   
  iorder_we = 1
  iorder_sn = 1
  IF(fx(1,1) > fx(2,1)) THEN
    iorder_we = -1
  END IF
  IF(fy(1,1) > fy(1,2)) THEN
    iorder_sn = -1
  END IF
  IF(detailout) THEN  
    WRITE(6,'(3X,A,I5)') 'X Order (WE) :',iorder_we 
    WRITE(6,'(3X,A,I5)') 'Y Order (SN) :',iorder_sn 

  END IF
   
  ratio=UNDEF
  dist_min_x=0
  dist_min_y=0
  nearestn=0
    ! ================================================================
    !   Nearest 4 Grid Points Interpolation
    ! ================================================================   
      ! ------------------------------------------------------------
      !    Search 4-Grid Points
      ! ------------------------------------------------------------      
      dist(1:num_grid_ave) = 1.D+10
      DO jy=1,ny-1
        DO ix=1,nx-1
          rxmax = MAXVAL(fx(ix:ix+1, jy:jy+1))
          rxmin = MINVAL(fx(ix:ix+1, jy:jy+1))
          rymax = MAXVAL(fy(ix:ix+1, jy:jy+1))
          rymin = MINVAL(fy(ix:ix+1, jy:jy+1))
         IF(rxmin <= datax .AND. rxmax >= datax .AND. &
           & rymin <= datay .AND. rymax >= datay ) THEN
          tmp_dist(1)=( fx(ix,jy) - datax )** order + ( fy(ix,jy) - datay )** order
          tmp_dist(2)=( fx(ix+1,jy) - datax )** order + ( fy(ix+1,jy) - datay )** order
          tmp_dist(3)=( fx(ix+1,jy+1) - datax )** order + ( fy(ix+1,jy+1) - datay )** order
          tmp_dist(4)=( fx(ix,jy+1) - datax )** order + ( fy(ix,jy+1) - datay )** order

   
          IF( maxval(tmp_dist) <= maxval(dist) )THEN
            nearestn=minloc(tmp_dist)
            dist=tmp_dist
            dist_min_x(1)=ix
            dist_min_x(2)=ix+1
            dist_min_x(3)=ix+1
            dist_min_x(4)=ix
            dist_min_y(1)=jy
            dist_min_y(2)=jy
            dist_min_y(3)=jy+1
            dist_min_y(4)=jy+1
          ENDIF 
         ENDIF

        END DO
      END DO

      IF( dist_min_x(1) > 0)THEN
      sum_dist = dist(1) + dist(2) + dist(3) + dist(4)
      ratio(1) = dist(1)/sum_dist
      ratio(2) = dist(2)/sum_dist
      ratio(3) = dist(3)/sum_dist
      ratio(4) = dist(4)/sum_dist
      ENDIF
      !IF(detailout) WRITE(6,'(2X,A,5F15.5)') 'ratio      :',ratio(1:4),SUM(ratio(1:4))

        

  RETURN
END SUBROUTINE com_xy2ij



END MODULE QC
