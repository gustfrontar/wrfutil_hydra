MODULE common_met_em
!=======================================================================
!
! [PURPOSE:] Common Information for WRF
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   07/14/2010 Takemasa Miyoshi  modified for WRF
!   10/16/2023 Juan Ruiz         modified for met_em_perturbation
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE netcdf

  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------

  !MODULE VARIABLES (NOT READ FROM NAMELIST) 
  INTEGER, PARAMETER :: nvar = 28
  INTEGER            :: nv3d , nv2d , ns3d

  REAL(r_sngl) :: tmpatt
  INTEGER,SAVE :: dx,dy ! grid spacing [m]
  REAL(r_size),SAVE :: pdx,pdy
  INTEGER,SAVE :: nlon
  INTEGER,SAVE :: nlat
  INTEGER,SAVE :: nlev
  INTEGER,SAVE :: nlev_soil
  INTEGER,SAVE :: ntime  !Number of times in file.

  INTEGER,SAVE      :: nij0
  INTEGER,SAVE      :: nlevall
  INTEGER,SAVE      :: ngpv
  REAL(r_size),ALLOCATABLE,SAVE :: ri(:,:) , rj(:,:) 

  CHARACTER(19) :: hdate !Date in input files.

  !Variables to be read and perturbed.
  LOGICAL, PARAMETER :: OPTIONAL_VARIABLE(nvar)= &  !Some input variables are optional (eg. condensates)
  &(/ .false. , .false. , .true. , .false. , .true. , .false. , .false. , .true. , .true. , .true. , &
  !   U          V         W         T         P1         P2   ,   GHT       RH        QV        QC        
  & .true. , .true. , .true. , .true. , .true. , .false. , .false. ,  &
  !   QR       QCI      QS       QG       QH       PS       SLP 
  &  .true. , .true. , .true. , .true. , .true. , .true. , .true. , .true. , .true. , &
  ! SKINTEMP ST100200 ST040100 ST010040 ST000010 SM100200 SM040100 SM010040 SM000010 
  & .true. , .true. /)
  !   SM       ST 

 !Variables to be read and perturbed.
  LOGICAL :: PERTURB_VARIABLE(nvar)= &  !Some input variables are optional (eg. condensates)
  &(/ .true. , .true. , .false. , .true. , .true. , .true. , .true. , .true.  ,.true. , .false. , &
  !   U          V         W         T         P1     P2   ,   GHT      RH        QV        QC        
  & .false. , .false. , .false. , .false. , .false. , .true. , .true. , &
  !   QR       QCI          QS       QG       QH       PS       SLP 
  &  .true. , .true. , .true. , .true. , .true. , .true. , .true. , .true. , .true. , &
  ! SKINTEMP ST100200 ST040100 ST010040 ST000010 SM100200 SM040100 SM010040 SM000010
  & .true. , .true. /)
  !   SM       ST 

  CHARACTER(3), PARAMETER :: VAR_TYPE(nvar)= &  !Which type of variable (A3D - atmo 3d, A2D - atmo 2d , S3D - soil 3d)
  &(/ 'A3D' , 'A3D' , 'A3D' , 'A3D' , 'A3D' , 'A3D' , 'A3D' , 'A3D' , 'A3D' , 'A3D' ,   &
  !     U       V       W       T       P1      P2     GHT      RH      QV      QC        
  &  'A3D'  , 'A3D' , 'A3D' , 'A3D' , 'A3D' , 'A2D' , 'A2D' , &
  !   QR       QCI      QS      QG      QH     PS      SLP  
  &  'A2D' ,  'A2D' ,   'A2D' ,  'A2D' ,  'A2D' ,  'A2D' ,  'A2D' ,  'A2D' ,  'A2D' , &
  ! SKINTEMP ST100200 ST040100 ST010040 ST000010 SM100200 SM040100 SM010040 SM000010 
  &  'S3D' , 'S3D' /)
  !  SM    ST 

  CHARACTER(1), PARAMETER :: STAGGER_TYPE(nvar)= &       !Which time of staggering do we have (U,V,W,M,S)
  &(/ 'U' , 'V' , 'W' , 'M' , 'M' , 'M' , 'M' ,  'M'  , 'M' , 'M' ,   &
  !   U      V     W     T     P1   P2    GHT    RH     QV    QC        
  &  'M'  , 'M' , 'M' , 'M' , 'M' , 'M' , 'M' , &
  !   QR    QCI   QS    QG    QH   PS    SLP  
  &  'M'    ,  'M'   ,   'M'  ,  'M'   ,  'M'   ,  'M'   ,  'M'   ,  'M'   ,  'M'   , &
  ! SKINTEMP ST100200 ST040100 ST010040 ST000010 SM100200 SM040100 SM010040 SM000010 
  &  'S' , 'S' /)
  !  SM    ST 
  CHARACTER(20), PARAMETER :: VAR_NAME(nvar) = &
  &(/ 'UU' , 'VV' , 'WW' , 'TT' , 'PRESSURE' , 'PRES' , 'GHT' , 'RH' , 'SPECHUMD' , 'QC' ,  &
  !   U      V        W     T        P1          P2      GHT    RH        QV        QC        
  &  'QR'  , 'QI' , 'QS' , 'QG' , 'QH' ,  'PSFC' , 'PMSL' , &
  !   QR      QCI   QS      QG     QH      PS       SLP  
  &'SKINTEMP','ST100200','ST040100','ST010040','ST000010','SM100200','SM040100','SM010040','SM000010', &
  ! SKINTEMP   ST100200   ST040100   ST010040   ST000010   SM100200   SM040100   SM010040   SM000010
  &  'SM' , 'ST' /)
  !  SM      ST  

  !Variables to store the actual variable names that are included in the data structures.
  CHARACTER(len=20) , ALLOCATABLE :: VAR_NAME_3D(:) , VAR_NAME_2D(:) , VAR_NAME_S3D(:)  

  LOGICAL :: PRESENT_VARIABLE(nvar) = .false.  !For checking input variables.

  CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------

SUBROUTINE get_weigths( member_ori , member_tar , method , wmatrix , sigma_nml ) 
IMPLICIT NONE
REAL(r_sngl) , INTENT(OUT) :: wmatrix(member_ori,member_tar)
REAL(r_size)               :: bufr( member_ori )
REAL(r_sngl) , INTENT(IN)  :: sigma_nml
INTEGER      , INTENT(IN ) :: member_ori , member_tar , method
INTEGER                    :: ii , jj , im

!This routine generates a transformation matrix that produce the desired perturbation.


IF ( method == 1 .or. method == 2 ) THEN  !Particle rejuvenation.
   DO jj = 1 , member_tar
     IF ( method == 1 ) THEN 
        CALL com_randn( member_ori , bufr )  !With Gaussian distribution
     ELSE
        CALL com_rand(  member_ori , bufr )  !With Uniform distribution
     ENDIF
     wmatrix(:,jj) = REAL( bufr , r_sngl )
   ENDDO 
   wmatrix = wmatrix * ( sigma_nml / REAL( member_ori ) ) ** 0.5 
   im = 1
   DO jj = 1 , member_tar 
      wmatrix(im,jj) = wmatrix(im,jj) + 1.0e0
      im = im + 1
      IF ( im > member_ori ) THEN
         im = 1
      ENDIF
      WRITE(6,*)'Weigths for member ',jj 
      WRITE(6,*)wmatrix(:,jj) 
   ENDDO
ENDIF

IF ( method == 3 ) THEN   !Specular 

   IF ( member_tar /= 2 * member_ori ) THEN
      WRITE(*,*)'Error: In specular method the target ensemble size should be '
      WRITE(*,*)'       twice the original ensmble size'
      WRITE(*,*)'       Original ensemble size ',member_ori
      WRITE(*,*)'       Target ensemble size   ',member_tar
      STOP 1
   ENDIF
   wmatrix = 0.0e0 
   im = 0
   DO ii = 1 , member_ori
      wmatrix(ii,ii)=1.0e0
      wmatrix(ii,ii+member_ori)=-1.0e0
   ENDDO
ENDIF

END SUBROUTINE get_weigths

SUBROUTINE set_common_met_em(inputfile)
  IMPLICIT NONE
  INTEGER :: i , j 
  INTEGER(4) :: ncid,varid,dimids(4),shape(4)
  CHARACTER(*)           :: inputfile
  INTEGER :: ierr , iv3d , iv2d , is3d 

  WRITE(6,'(A)') 'Hello from set_common_wrf'
  !
  WRITE(6,'(A)') 'READING: ',inputfile
  CALL check_io(NF90_OPEN(inputfile,NF90_NOWRITE,ncid))
  CALL check_io(NF90_INQ_VARID(ncid,TRIM(VAR_NAME(1)),varid))
  CALL check_io(NF90_INQUIRE_VARIABLE(ncid,varid,dimids=dimids) )
  DO i=1,4
    CALL check_io(NF90_INQUIRE_DIMENSION(ncid,dimids(i),len=shape(i)))
  END DO
  nlon = shape(1)
  nlat = shape(2) + 1
  nlev = shape(3) + 1
  ntime= shape(4)

  CALL check_io(NF90_INQ_VARID(ncid,'SM',varid))
  CALL check_io(NF90_INQUIRE_VARIABLE(ncid,varid,dimids=dimids) )
  DO i=1,4
    CALL check_io(NF90_INQUIRE_DIMENSION(ncid,dimids(i),len=shape(i)))
  END DO
  nlev_soil = shape(3) 


  WRITE(6,'(A)') '*** grid information ***'
  WRITE(6,'(3(2X,A,I5))') 'nlon =',nlon,'nlat =',nlat,'nlev =',nlev,'nlev_soil=',nlev_soil
  nij0=nlon*nlat
  ngpv=nij0*nlevall

  WRITE(6,'(a)') '***  END  : READ NETCDF (CNST) ***'

  nv3d = 0
  nv2d = 0
  ns3d = 0

  !Check input variables
  DO i = 1,nvar
      !Check if the variable is present.
      ierr=NF90_INQ_VARID(ncid,TRIM(VAR_NAME(i)),varid)
      IF( ierr == NF90_NOERR )THEN
        PRESENT_VARIABLE(i) = .true.
        WRITE(6,*)"Variable ", VAR_NAME(i)," is present in input file."
      ELSE
        IF( OPTIONAL_VARIABLE(i) )THEN
          WRITE(6,*)"Warning: Variable ", VAR_NAME(i)," is not present in input file."
        ELSE 
          WRITE(6,*)"Error: Non optional variable ",VAR_NAME(i)," is not present in input file -> Aborting execution"
          STOP 1
        ENDIF
      ENDIF
      !Get the number of variables in the different variable types that will
      !be allocated.
      IF ( PRESENT_VARIABLE(i) .and. PERTURB_VARIABLE(i) ) THEN
        SELECT CASE ( VAR_TYPE(i) )
          CASE( 'A3D' )
             nv3d = nv3d + 1
          CASE( 'A2D' )
             nv2d = nv2d + 1
          CASE( 'S3D' )
             ns3d = ns3d + 1
          CASE DEFAULT
             WRITE(6,*)'Warning: not recognized variable type: ',VAR_TYPE(i)
        END SELECT
      ENDIF
  ENDDO
  CALL check_io(NF90_CLOSE(ncid))
  nlevall=nlev*nv3d+nv2d+ns3d*nlev_soil

  WRITE(6,*)'The number of 3D atmospheric variables to be perturbed is ',nv3d
  WRITE(6,*)'The number of 2D variables to be perturbed is             ',nv2d
  WRITE(6,*)'The number of 3D soil variables to be perturbed is        ',ns3d
  WRITE(6,*)'Nleval is ',nlevall

  ALLOCATE( VAR_NAME_3D(nv3d) , VAR_NAME_2d(nv2d) , VAR_NAME_S3D(ns3d) )

  iv3d = 0
  iv2d = 0
  is3d = 0
  !Check input variables
  DO i = 1,nvar
      IF ( PRESENT_VARIABLE(i) .and. PERTURB_VARIABLE(i) ) THEN
        SELECT CASE ( VAR_TYPE(i) )
          CASE( 'A3D' )
             iv3d = iv3d + 1
             VAR_NAME_3D(iv3d)=VAR_NAME(i)
          CASE( 'A2D' )
             iv2d = iv2d + 1
             VAR_NAME_2d(iv2d)=VAR_NAME(i)
          CASE( 'S3D' )
             is3d = is3d + 1
             VAR_NAME_s3d(is3d)=VAR_NAME(i)
        END SELECT
      ENDIF
  ENDDO


  !Generate the indices for the matrix
  ALLOCATE(ri(nlon,nlat))
  ALLOCATE(rj(nlon,nlat))
  DO j=1,nlat
    DO i=1,nlon
      ri(i,j) = REAL(i,r_sngl)
      rj(i,j) = REAL(j,r_sngl)
    END DO
  END DO

  RETURN
END SUBROUTINE set_common_met_em


SUBROUTINE read_grd(filename,v3d,v2d,vs3d)
IMPLICIT NONE
CHARACTER(len=*) :: filename
REAL(r_sngl)  :: v3d(nlon,nlat,nlev,nv3d) , v2d(nlon,nlat,nv2d) , vs3d(nlon,nlat,nlev_soil,ns3d)
INTEGER       :: ncid
INTEGER       :: iv , iv2d , iv3d , ivs3d
INTEGER       :: varid
INTEGER, ALLOCATABLE         :: start(:),count(:)
INTEGER                      :: nxvar,nyvar,nzvar
LOGICAL                      :: readvar

iv2d = 1
iv3d = 1
ivs3d =1
!Open the file.
CALL open_wrf_file(filename,'ro',ncid)

DO iv = 1 , nvar !Loop over the variables.

   readvar = present_variable(iv) .and. perturb_variable(iv)  
   nxvar=nlon-1
   nyvar=nlat-1
   nzvar=nlev-1
   SELECT CASE ( STAGGER_TYPE(iv) )
      CASE('U')
        nxvar=nlon
      CASE('V')
        nyvar=nlat
      CASE('W')
        nzvar=nlev
      CASE('S')
        nzvar=nlev_soil  
      CASE('M')  
      CASE DEFAULT
        WRITE(6,*)"ERROR: Stagger option not recognized: ",STAGGER_TYPE(iv)
        STOP 1
   END SELECT

   IF ( readvar ) THEN

      SELECT CASE ( VAR_TYPE(iv) )
         CASE('A3D')
           ALLOCATE(start(4),count(4))
           start = (/ 1,1,1,1 /)
           count = (/ nxvar,nyvar,nzvar,1 /)
           !READ VARIABLE
           CALL check_io(NF90_INQ_VARID(ncid,TRIM(VAR_NAME(iv)),varid))
           CALL check_io(NF90_GET_VAR(ncid,varid,v3d(1:nxvar,1:nyvar,1:nzvar,iv3d),start,count))
           iv3d=iv3d+1
           DEALLOCATE(start,count)

         CASE('S3D')
           ALLOCATE(start(4),count(4))
           start = (/ 1,1,1,1 /)
           count = (/ nxvar,nyvar,nzvar,1 /)      
           !READ VARIABLE
           CALL check_io(NF90_INQ_VARID(ncid,TRIM(VAR_NAME(iv)),varid))
           CALL check_io(NF90_GET_VAR(ncid,varid,vs3d(1:nxvar,1:nyvar,1:nzvar,ivs3d),start,count))
           ivs3d=ivs3d+1
           DEALLOCATE(start,count)

         CASE('A2D')
           ALLOCATE(start(3),count(3))
           start = (/ 1,1,1 /)
           count = (/ nxvar,nyvar,1,1 /) !READ ONE SINGLE LEVEL
           !READ VARIABLE
           CALL check_io(NF90_INQ_VARID(ncid,TRIM(VAR_NAME(iv)),varid))
           CALL check_io(NF90_GET_VAR(ncid,varid,v2d(1:nxvar,1:nyvar,iv2d),start,count))
           iv2d=iv2d+1         
           DEALLOCATE(start,count)  

         CASE DEFAULT
           WRITE(6,*)"Error: Not recognized variable type: ",VAR_TYPE(iv)
           STOP 1

      END SELECT

   END IF   ! End of wheter we are going to read the variable or not.

END DO !End do over variables   
CALL close_wrf_file(ncid)

RETURN

END SUBROUTINE read_grd

SUBROUTINE write_grd(filename,v3d,v2d,vs3d)
IMPLICIT NONE
CHARACTER(*) :: filename
REAL(r_sngl) :: v3d(nlon,nlat,nlev,nv3d) , v2d(nlon,nlat,nv2d) , vs3d(nlon,nlat,nlev_soil,ns3d)
INTEGER(4)   :: ncid
INTEGER      :: iv , iv2d , iv3d , ivs3d
INTEGER(4)   :: varid
INTEGER   , ALLOCATABLE :: start(:),count(:)
INTEGER                 :: nxvar,nyvar,nzvar
LOGICAL                 :: writevar

!Open the file
CALL open_wrf_file(filename,'rw',ncid)

!Read the data
iv2d = 1
iv3d = 1
ivs3d =1

DO iv = 1 , nvar !Main loop over variables

   writevar = present_variable(iv) .and. perturb_variable(iv)
   nxvar=nlon-1
   nyvar=nlat-1
   nzvar=nlev-1
   SELECT CASE ( STAGGER_TYPE(iv) )
      CASE('U')
        nxvar=nlon
      CASE('V')
        nyvar=nlat
      CASE('W')
        nzvar=nlev
      CASE('S')
        nzvar=nlev_soil
      CASE('M')
      CASE DEFAULT
        WRITE(6,*)"Error: Staggering type not recognized ",STAGGER_TYPE(iv)
        STOP 1
   END SELECT

   IF ( writevar ) THEN

      SELECT CASE ( VAR_TYPE(iv) )
         CASE('A3D')
           ALLOCATE(start(4),count(4))
           start = (/ 1,1,1,1 /)
           count = (/ nxvar,nyvar,nzvar,1 /)
           !WRITE VARIABLE
           CALL check_io(NF90_INQ_VARID(ncid,TRIM(VAR_NAME(iv)),varid))
           CALL check_io(NF90_PUT_VAR(ncid,varid,v3d(1:nxvar,1:nyvar,1:nzvar,iv3d),start,count))
           iv3d=iv3d+1
           DEALLOCATE( start , count )

         CASE('S3D') 
           ALLOCATE(start(4),count(4))
           start = (/ 1,1,1,1 /)
           count = (/ nxvar,nyvar,nzvar,1 /)
           !WRITE VARIABLE
           CALL check_io(NF90_INQ_VARID(ncid,TRIM(VAR_NAME(iv)),varid))
           CALL check_io(NF90_PUT_VAR(ncid,varid,vs3d(1:nxvar,1:nyvar,1:nzvar,ivs3d),start,count))
           ivs3d=ivs3d+1
           DEALLOCATE( start , count )

         CASE('A2D')
           ALLOCATE(start(3),count(3))
           start = (/ 1,1,1 /)
           count = (/ nxvar,nyvar,1,1 /) !READ ONE SINGLE LEVEL
           !WRITE VARIABLE
           CALL check_io(NF90_INQ_VARID(ncid,TRIM(VAR_NAME(iv)),varid))
           CALL check_io(NF90_PUT_VAR(ncid,varid,v2d(1:nxvar,1:nyvar,iv2d),start,count))
           iv2d=iv2d+1    
           DEALLOCATE( start , count )       

         CASE DEFAULT
           WRITE(6,*)"ERROR: Not recognized variable type: ",VAR_TYPE(iv)
           STOP 1

       END SELECT

   END IF   ! End of wheter we are going to write the variable or not.

END DO !End do over the model variables.

CALL close_wrf_file(ncid)


RETURN

END SUBROUTINE write_grd

!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
SUBROUTINE monit_grd(v3d,v2d,s3d)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: s3d(nlon,nlat,nlev_soil,ns3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: k,n

  DO k=1,nlev
    WRITE(6,'(I2,A)') k,'th level'
    DO n=1,nv3d
      WRITE(6,'(A,2ES10.2)') VAR_NAME(n),MAXVAL(v3d(:,:,k,n)),MINVAL(v3d(:,:,k,n))
    END DO
  END DO

  DO k=1,nlev_soil
    WRITE(6,'(I2,A)') k,'th level'
    DO n=1,ns3d
      WRITE(6,'(A,2ES10.2)') VAR_NAME(n),MAXVAL(s3d(:,:,k,n)),MINVAL(s3d(:,:,k,n))
    END DO
  END DO  

  DO n=1,nv2d
    WRITE(6,'(A,2ES10.2)') VAR_NAME(nv3d+n),MAXVAL(v2d(:,:,n)),MINVAL(v2d(:,:,n))
  END DO

  RETURN
END SUBROUTINE monit_grd

!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
SUBROUTINE ensmean_grd(member,nij,v3d,v2d,s3d,v3dm,v2dm,s3dm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member
  INTEGER,INTENT(IN) :: nij
  REAL(r_sngl),INTENT(IN)  :: v3d(nij,nlev,member,nv3d)
  REAL(r_sngl),INTENT(IN)  :: s3d(nij,nlev_soil,member,ns3d)
  REAL(r_sngl),INTENT(IN)  :: v2d(nij,member,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: s3dm(nij,nlev_soil,ns3d)
  REAL(r_sngl),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: m

  v3dm = v3d(:,:,1,:)
  DO m=2,member
     v3dm = v3dm + v3d(:,:,m,:)
  END DO
  v3dm = v3dm / REAL(member,r_sngl)

  s3dm = s3d(:,:,1,:)
  DO m=2,member
     s3dm = s3dm + s3d(:,:,m,:)
  END DO
  s3dm = s3dm / REAL(member,r_sngl)

  v2dm = v2d(:,1,:)
  DO m=2,member
     v2dm = v2dm + v2d(:,m,:)
  END DO
  v2dm = v2dm / REAL(member,r_sngl)

  RETURN
END SUBROUTINE ensmean_grd

SUBROUTINE write_date_met_em(inputfile,new_date)
IMPLICIT NONE
character(*)  , intent(in)  :: inputfile
character(19) , intent(in)  :: new_date
integer :: start(2) , count(2) , ncid
character(19) :: tmpdate
integer :: varid
  WRITE(6,*)'Updating the date' 
  CALL check_io(NF90_OPEN(inputfile,NF90_WRITE,ncid))
  CALL check_io(NF90_INQ_VARID(ncid,'Times',varid))
  start = (/ 1 ,1 /)
  count = (/ 19,1 /)
  !Read the current date in file to get the date format.
  CALL check_io(NF90_GET_VAR(ncid,varid,tmpdate,start,count))
  WRITE(6,*)'Previous date was ',tmpdate,' new date is ',new_date

  !Write the new date in the file.
  CALL check_io(NF90_PUT_VAR(ncid,varid,new_date,start,count))
  !CALL check_io(NF90_PUT_ATT(ncid,NF90_GLOBAL,'START_DATE',new_date))
  CALL check_io(NF90_PUT_ATT(ncid,NF90_GLOBAL,'SIMULATION_START_DATE',new_date))
  CALL check_io(NF90_CLOSE(ncid))

END SUBROUTINE write_date_met_em



SUBROUTINE check_io(status)
  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: status

  IF(status /= NF90_NOERR ) THEN
    WRITE(6,*) TRIM(NF90_STRERROR(status))
    STOP 10
  ENDIF

  RETURN
END SUBROUTINE check_io

SUBROUTINE open_wrf_file(filename,mode,ncid)
  IMPLICIT NONE
  INTEGER(4), INTENT(OUT)   :: ncid
  CHARACTER(*) , INTENT(IN) :: filename,mode

  IF(mode .eq. 'ro' )THEN
  CALL check_io(NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid))
  ELSEIF(mode .eq. 'rw' )THEN
  CALL check_io(NF90_OPEN(TRIM(filename),NF90_WRITE,ncid))
  END IF

  RETURN
END SUBROUTINE open_wrf_file

SUBROUTINE close_wrf_file(ncid)
  IMPLICIT NONE
  INTEGER(4), INTENT(IN)   :: ncid

  CALL  check_io(NF90_CLOSE(ncid))

  RETURN
END SUBROUTINE close_wrf_file


END MODULE common_met_em
