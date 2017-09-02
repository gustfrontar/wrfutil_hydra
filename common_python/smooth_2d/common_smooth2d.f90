MODULE common_smooth2d
  !!======================================================================
  !!                     ***  MODULE smooth2d  ***
  !!=====================================================================
  !!  ** Purpose :  perform a spatial filtering on input field.
  !!               - various filters are available :
  !!               1: Lanczos (default)
  !!               2: hanning
  !!               3: shapiro
  !!
  !!  ** Method  : input filed, output filtered field
  !!
  !! History : --   : 1995     : J.M. Molines : Original code for spem
  !!         : 2.1  : 07/2007  : J.M. Molines : port in cdftools
  !!         : 2.1  : 05/2010  : R. Dussin    : Add shapiro filter
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!                  05/2014  : J. Ruiz      : Implemented in LETKF-WRF system
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!  filterinit   : initialise weight
  !!  filter       : main routine for filter computation
  !!  initlanc     : initialise lanczos weights
  !!  inithann     : initialise hanning weights
  !!  initshap     : initialise shapiro routine
  !!  initbox      : initialize weight for box car average
  !!  lislanczos2d : Lanczos filter
  !!  lishan2d     : hanning 2d filter
  !!  lisshapiro1d : shapiro filter
  !!  lisbox       : box car filter
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  PUBLIC filter_2d

!-----------------------------------------------------------------------
! Variable size definitions
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------
  !
  INTEGER, PARAMETER                    :: jp_lanc=1         ! lancszos id
  INTEGER, PARAMETER                    :: jp_hann=2         ! hanning id
  INTEGER, PARAMETER                    :: jp_shap=3         ! shapiro id
  INTEGER, PARAMETER                    :: jp_boxc=4         ! box car id
  INTEGER                               :: jk, jt, jvar      ! dummy loop index
  INTEGER                               :: npiglo, npjglo    ! size of the domain
  INTEGER                               :: ncut, nband       ! cut period/ length, bandwidth
  INTEGER                               :: nfilter = jp_lanc ! default value
  INTEGER, DIMENSION(:,:),  ALLOCATABLE :: iw                ! flag for bad values (or land masked )
  REAL(8), DIMENSION(:,:),     ALLOCATABLE :: v2d, w2d          ! raw data,  filtered result
  REAL(8)                                  :: fn, rspval        ! cutoff freq/wavelength, spval

  REAL(8), DIMENSION(:,:),     ALLOCATABLE :: dec2d             ! working array
  REAL(8), DIMENSION(:),       ALLOCATABLE :: dec, de           ! weight in r8, starting index 0:nband

  !!----------------------------------------------------------------------
 
CONTAINS
  SUBROUTINE filter_2d(inputvar2d,outputvar2d,mask,lambdaf,ctyp,nx,ny)
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: inputvar2d(nx,ny)
  REAL(8), INTENT(OUT):: outputvar2d(nx,ny)
  INTEGER, INTENT(IN)      :: mask(nx,ny)
  INTEGER, INTENT(IN)      :: nx,ny
  INTEGER, INTENT(IN)      :: lambdaf
  CHARACTER(LEN=256), INTENT(IN) :: ctyp              ! filter type
  ncut=lambdaf
  !usage : cdfsmooth IN-file ncut [filter_type]'
  !   PRINT *,'      '
  !  PRINT *,'     PURPOSE :'
  !   PRINT *,'       Perform a spatial smoothing on the file using a particular'
  !   PRINT *,'       filter as specified in the option. Available filters' 
  !   PRINT *,'       are : Lanczos, Hanning, Shapiro, Box car average. Default'
  !   PRINT *,'       is Lanczos filter.'
  !   PRINT *,'      '
  !   PRINT *,'     ARGUMENTS :'
  !   PRINT *,'        ncut     : number of grid step to be filtered'
  !   PRINT *,'      '
  !   PRINT *,'     OPTIONS :'
  !   PRINT *,'        [filter_type] : Lanczos, L, l  (default)'
  !   PRINT *,'                        Hanning, H, h'
  !   PRINT *,'                        Shapiro, S, s'
  !   PRINT *,'                        Box    , B, b'
  !   PRINT *,'      
  !  remark: for a spatial filter, fn=dx/lambda where dx is spatial step, lamda is cutting wavelength
  fn    = 1./ncut
  nband = 2*ncut    ! Bandwidth of filter is twice the filter span

  ALLOCATE ( dec(0:nband) , de(0:nband) )

     SELECT CASE ( ctyp)
     CASE ( 'Lanczos','L','l') 
        nfilter=jp_lanc
     CASE ( 'Hanning','H','h')
        nfilter=jp_hann
        ALLOCATE ( dec2d(0:2,0:2) )
     CASE ( 'Shapiro','S','s')
        nfilter=jp_shap
     CASE ( 'Box','B','b')
        nfilter=jp_boxc
     CASE DEFAULT
        PRINT *, TRIM(ctyp),' : undefined filter ' ; STOP
     END SELECT

  CALL filterinit (nfilter, fn, nband)
  ! Look for input file and create outputfile
  npiglo = nx
  npjglo = ny
  

  ALLOCATE ( v2d(npiglo,npjglo),iw(npiglo,npjglo), w2d(npiglo,npjglo) )

   iw = mask !IF mask == 1 filter will be performed , if mask == 0 point will be skiped.
   v2d=inputvar2d
   IF ( ncut /= 0 ) CALL filter( nfilter, v2d, iw, w2d)
   IF ( ncut == 0 ) w2d = v2d
   WHERE( iw == 0 )w2d = v2d

   outputvar2d=w2d

   DEALLOCATE( v2d,iw,w2d,dec,de )



  END SUBROUTINE filter_2d



  SUBROUTINE filterinit(kfilter, pfn, kband)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE filterinit  ***
    !!
    !! ** Purpose :   initialise weight according to filter type
    !!
    !!----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(in) :: kfilter  ! filter number
    REAL(8),    INTENT(in) :: pfn      ! filter cutoff frequency/wavelength
    INTEGER, INTENT(in) :: kband    ! filter bandwidth
    !!----------------------------------------------------------------------
    SELECT CASE ( kfilter)
    CASE ( jp_lanc )
       CALL initlanc (pfn, kband)
    CASE ( jp_hann )
       CALL inithann (pfn, kband)
    CASE ( jp_shap )
       CALL initshap (pfn, kband)
    CASE ( jp_boxc )
       CALL initbox  (pfn, kband)
    END SELECT

  END SUBROUTINE filterinit


  SUBROUTINE filter (kfilter, px, kpx, py)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE filter  ***
    !!
    !! ** Purpose :  Call the proper filter routine according to filter type 
    !!
    !!----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER,                 INTENT(in) :: kfilter ! filter number
    REAL(8), DIMENSION(:,:),    INTENT(in) :: px      ! input data
    INTEGER, DIMENSION(:,:), INTENT(in) :: kpx     ! validity flag
    REAL(8), DIMENSION(:,:),   INTENT(out) :: py      ! output data
    !!----------------------------------------------------------------------
    SELECT CASE ( kfilter)
    CASE ( jp_lanc )
       CALL lislanczos2d (px, kpx, py, npiglo, npjglo, fn, nband)
    CASE ( jp_hann )
       CALL lishan2d     (px, kpx, py, ncut, npiglo, npjglo)
    CASE ( jp_shap )
       CALL lisshapiro1d (px, kpx, py, ncut, npiglo, npjglo)
    CASE ( jp_boxc )
       CALL lisbox       (px, kpx, py, npiglo, npjglo, fn, nband)
    END SELECT

  END SUBROUTINE filter


  SUBROUTINE initlanc(pfn, knj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE initlanc  ***
    !!
    !! ** Purpose : initialize lanczos weights
    !!
    !!----------------------------------------------------------------------
    IMPLICIT NONE



    REAL(8),    INTENT(in) :: pfn  ! cutoff freq/wavelength
    INTEGER, INTENT(in) :: knj  ! bandwidth
    INTEGER             :: ji   ! dummy loop index
    REAL(8)                :: dl_pi, dl_ey, dl_coef
    !!----------------------------------------------------------------------
    dl_pi   = ACOS(-1.d0)
    dl_coef = 2.0d0*dl_pi*pfn

    de(0) = 2.d0*pfn
    DO  ji=1,knj
       de(ji) = SIN(dl_coef*ji)/(dl_pi*ji)
    END DO
    !
    dec(0) = 2.d0*pfn
    DO ji=1,knj
       dl_ey   = dl_pi*ji/knj
       dec(ji) = de(ji)*SIN(dl_ey)/dl_ey
    END DO
  
  END SUBROUTINE initlanc


  SUBROUTINE inithann(pfn, knj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE inithann  ***
    !!
    !! ** Purpose : Initialize hanning weight 
    !!
    !!----------------------------------------------------------------------
    IMPLICIT NONE



    REAL(8),    INTENT(in) :: pfn  ! cutoff freq/wavelength
    INTEGER, INTENT(in) :: knj  ! bandwidth

    REAL(8)                :: dl_sum
    !!----------------------------------------------------------------------
    dec2d(:,:) = 0.d0 
    ! central point
    dec2d(1,1) = 4.d0
    ! along one direction
    dec2d(1,0) = 1.d0 ;  dec2d(1,2) = 1.d0
    ! and the other 
    dec2d(0,1) = 1.d0 ;  dec2d(2,1) = 1.d0

    ! normalize
    dl_sum     = SUM(dec2d)
    dec2d(:,:) = dec2d(:,:) / dl_sum

  END SUBROUTINE inithann


  SUBROUTINE initshap(pfn, knj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE initshap  ***
    !!
    !! ** Purpose :  Dummy routine to respect program structure 
    !!
    !!----------------------------------------------------------------------
    IMPLICIT NONE

    

    REAL(8),    INTENT(in) :: pfn  ! cutoff freq/wavelength
    INTEGER, INTENT(in) :: knj  ! bandwidth
    !!----------------------------------------------------------------------
    !   nothing to do 

  END SUBROUTINE initshap


  SUBROUTINE initbox(pfn, knj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE initbox  ***
    !!
    !! ** Purpose :  Init weights for box car 
    !!
    !!----------------------------------------------------------------------
    IMPLICIT NONE
 
    

 
    REAL(8),    INTENT(in) :: pfn  ! cutoff freq/wavelength
    INTEGER, INTENT(in) :: knj  ! bandwidth
    !!----------------------------------------------------------------------
    dec(:) = 1.d0

  END SUBROUTINE initbox


  SUBROUTINE lislanczos2d(px, kiw, py, kpi, kpj, pfn, knj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE lislanczos2d  ***
    !!
    !! ** Purpose : Perform lanczos filter
    !!
    !! ** Method  :   px      = input data
    !!                kiw     = validity of input data
    !!                py      = output filter
    !!                kpi,kpj = number of input/output data
    !!                pfn     = cutoff frequency
    !!                knj     = bandwith of the filter
    !!
    !! References : E. Blayo (1992) from CLS source and huge optimization 
    !!----------------------------------------------------------------------
    IMPLICIT NONE



    REAL(8),    DIMENSION(:,:), INTENT(in ) :: px               ! input array
    INTEGER, DIMENSION(:,:), INTENT(in ) :: kiw              ! flag input array
    REAL(8),    DIMENSION(:,:), INTENT(out) :: py               ! output array
    INTEGER,                 INTENT(in ) :: kpi, kpj         ! size of input/output
    REAL(8),                    INTENT(in ) :: pfn              ! cutoff frequency/wavelength
    INTEGER,                 INTENT(in ) :: knj              ! filter bandwidth

    INTEGER                              :: ji, jj, jmx, jkx !  dummy loop index
    INTEGER                              :: ik1x, ik2x, ikkx
    INTEGER                              :: ifrst=0
    INTEGER                              :: inxmin, inxmaxi
    INTEGER                              :: inymin, inymaxi
    REAL(8), DIMENSION(kpi,kpj)             :: dl_tmpx, dl_tmpy
    REAL(8)                                 :: dl_yy, dl_den
    !!----------------------------------------------------------------------
    inxmin   =  knj
    inxmaxi  =  kpi-knj+1
    inymin   =  knj
    inymaxi  =  kpj-knj+1

    !PRINT *,' filtering parameters'
    !PRINT *,'    nx    = ', kpi
    !PRINT *,'    nband = ', knj
    !PRINT *,'    fn    = ', pfn

    DO jj=1,kpj
       DO  jmx=1,kpi
          ik1x = -knj
          ik2x =  knj
          !
          IF (jmx <= inxmin ) ik1x = 1-jmx
          IF (jmx >= inxmaxi) ik2x = kpi-jmx
          !
          dl_yy  = 0.d0
          dl_den = 0.d0
          !
          DO jkx=ik1x,ik2x
             ikkx=ABS(jkx)
             IF (kiw(jkx+jmx,jj)  ==  1) THEN
                dl_den = dl_den + dec(ikkx)
                dl_yy  = dl_yy  + dec(ikkx)*px(jkx+jmx,jj)
             END IF
          END DO
          !
          dl_tmpx(jmx,jj)=dl_yy/dl_den
       END DO
    END DO

    DO ji=1,kpi
       DO  jmx=1,kpj
          ik1x = -knj
          ik2x =  knj
          !
          IF (jmx <= inymin ) ik1x = 1-jmx
          IF (jmx >= inymaxi) ik2x = kpj-jmx
          !
          dl_yy  = 0.d0
          dl_den = 0.d0
          !
          DO jkx=ik1x,ik2x
             ikkx=ABS(jkx)
             IF (kiw(ji,jkx+jmx)  ==  1) THEN
                dl_den = dl_den + dec(ikkx)
                dl_yy  = dl_yy  + dec(ikkx)*dl_tmpx(ji,jkx+jmx)
             END IF
          END DO
          py(ji,jmx)=0.
          IF (dl_den /=  0.) py(ji,jmx) = dl_yy/dl_den
       END DO
    END DO
    !
  END SUBROUTINE lislanczos2d

  SUBROUTINE lishan2d(px, kiw, py, korder, kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE lishan2d  ***
    !!
    !! ** Purpose : compute hanning filter at order korder
    !!
    !!----------------------------------------------------------------------
    IMPLICIT NONE



    REAL(8),    DIMENSION(:,:), INTENT(in ) :: px      ! input data
    INTEGER, DIMENSION(:,:), INTENT(in ) :: kiw     ! validity flags
    REAL(8),    DIMENSION(:,:), INTENT(out) :: py      ! output data
    INTEGER,                 INTENT(in ) :: korder  ! order of the filter
    INTEGER,                 INTENT(in ) :: kpi, kpj ! size of the data
    INTEGER                              :: jj, ji, jorder  ! loop indexes
    INTEGER                              :: iiplus1, iiminus1
    INTEGER                              :: ijplus1, ijminus1
    REAL(8), DIMENSION(:,:), ALLOCATABLE    :: ztmp
    !!----------------------------------------------------------------------
    ALLOCATE( ztmp(kpi,kpj) )

    py(:,:)   = 0.
    ztmp(:,:) = px(:,:)

    DO jorder = 1, korder
       DO jj   = 2, kpj-1
          DO ji = 2, kpi-1
             !treatment of the domain frontiers
             iiplus1 = MIN(ji+1,kpi) ; iiminus1 = MAX(ji-1,1) 
             ijplus1 = MIN(jj+1,kpj) ; ijminus1 = MAX(jj-1,1) 

             ! we don't compute in land
             IF ( kiw(ji,jj) == 1 ) THEN
                py(ji,jj) = SUM( dec2d(:,:) * ztmp(iiminus1:iiplus1,ijminus1:ijplus1) )
             ENDIF
          ENDDO
       ENDDO
       ! update the ztmp array
       ztmp(:,:) = py(:,:)
    ENDDO

  DEALLOCATE(ztmp)
  END SUBROUTINE lishan2d

  SUBROUTINE lisshapiro1d(px, kiw, py, korder, kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE lisshapiro1d  ***
    !!
    !! ** Purpose :  compute shapiro filter 
    !!
    !! References :  adapted from Mercator code
    !!----------------------------------------------------------------------
    IMPLICIT NONE


    REAL(8),    DIMENSION(:,:), INTENT(in ) :: px      ! input data
    INTEGER, DIMENSION(:,:), INTENT(in ) :: kiw     ! validity flags
    REAL(8),    DIMENSION(:,:), INTENT(out) :: py      ! output data
    INTEGER,                 INTENT(in ) :: korder  ! order of the filter
    INTEGER,                 INTENT(in ) :: kpi, kpj ! size of the data
    INTEGER                              :: jj, ji, jorder  ! loop indexes
    INTEGER                              :: imin, imax, ihalo=0
    REAL(8), PARAMETER                      :: rp_aniso_diff_XY = 2.25 !  anisotrope case
    REAL(8)                                 :: zalphax, zalphay, znum
    REAL(8), DIMENSION(:,:), ALLOCATABLE    :: ztmp , zpx , zpy, zkiw
    LOGICAL                                      :: ll_cycl = .TRUE.
    !!----------------------------------------------------------------------

    IF(ll_cycl) ihalo=1
    ! we allocate with an ihalo
    ALLOCATE( ztmp(0:kpi+ihalo,kpj) , zkiw(0:kpi+ihalo,kpj) )
    ALLOCATE( zpx (0:kpi+ihalo,kpj) , zpy (0:kpi+ihalo,kpj) )

    IF(ll_cycl) THEN
       zpx(1:kpi,:) = px(:  ,:) ;  zkiw(1:kpi,:) = kiw(:  ,:)
       zpx(0    ,:) = px(kpi,:) ;  zkiw(0    ,:) = kiw(kpi,:)
       zpx(kpi+1,:) = px(1  ,:) ;  zkiw(kpi+1,:) = kiw(1  ,:)
    ELSE
       zpx(:    ,:) = px(:  ,:)
    ENDIF

    zpy (:,:) = zpx(:,:)  ! init?
    ztmp(:,:) = zpx(:,:)  ! init

    zalphax=1./2.
    zalphay=1./2.

    !  Dx/Dy=rp_aniso_diff_XY  , D_ = vitesse de diffusion
    !  140 passes du fitre, Lx/Ly=1.5, le rp_aniso_diff_XY correspondant est:

    IF ( rp_aniso_diff_XY >=  1. ) zalphay=zalphay/rp_aniso_diff_XY
    IF ( rp_aniso_diff_XY <   1. ) zalphax=zalphax*rp_aniso_diff_XY

    DO jorder=1,korder
       imin = 2     - ihalo
       imax = kpi-1 + ihalo
       DO ji = imin,imax
          DO jj = 2,kpj-1
             ! We crop on the coast
             znum =    ztmp(ji,jj)                                                  &
                  &    + 0.25*zalphax*(ztmp(ji-1,jj  )-ztmp(ji,jj))*zkiw(ji-1,jj  ) &
                  &    + 0.25*zalphax*(ztmp(ji+1,jj  )-ztmp(ji,jj))*zkiw(ji+1,jj  ) &
                  &    + 0.25*zalphay*(ztmp(ji  ,jj-1)-ztmp(ji,jj))*zkiw(ji  ,jj-1) &
                  &    + 0.25*zalphay*(ztmp(ji  ,jj+1)-ztmp(ji,jj))*zkiw(ji  ,jj+1)
             zpy(ji,jj) = znum*zkiw(ji,jj)+zpx(ji,jj)*(1.-zkiw(ji,jj))
          ENDDO  ! end loop ji
       ENDDO  ! end loop jj

       IF ( ll_cycl ) THEN
          zpy(0    ,:) = zpy(kpi,:) 
          zpy(kpi+1,:) = zpy(1  ,:) 
       ENDIF

       ! update the tmp array
       ztmp(:,:) = zpy(:,:)

    ENDDO

    ! return this array
    IF( ll_cycl ) THEN
       py(:,:) = zpy(1:kpi,:)
    ELSE
       py(:,:) = zpy(:    ,:)
    ENDIF

  DEALLOCATE(ztmp,zpx,zkiw,zpy)

  END SUBROUTINE lisshapiro1d


  SUBROUTINE lisbox(px, kiw, py, kpi, kpj, pfn, knj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE lisbox  ***
    !!
    !! ** Purpose :  Perform box car filtering 
    !!
    !!----------------------------------------------------------------------
    IMPLICIT NONE




    REAL(8),    DIMENSION(:,:), INTENT(in ) :: px               ! input array
    INTEGER, DIMENSION(:,:), INTENT(in ) :: kiw              ! flag input array
    REAL(8),    DIMENSION(:,:), INTENT(out) :: py               ! output array
    INTEGER,                 INTENT(in ) :: kpi, kpj         ! size of input/output
    REAL(8),                    INTENT(in ) :: pfn              ! cutoff frequency/wavelength
    INTEGER,                 INTENT(in ) :: knj              ! filter bandwidth

    INTEGER                              :: ji, jj
    INTEGER                              :: ik1x, ik2x, ik1y, ik2y
    REAL(8)                                 :: dl_den
    LOGICAL, DIMENSION(kpi,kpj)                  :: ll_mask
    !!----------------------------------------------------------------------
    ll_mask=.TRUE.
    WHERE (kiw == 0 ) ll_mask=.FALSE.
    DO ji=1,kpi
       ik1x = ji-knj       ; ik2x = ji+knj
       ik1x = MAX(1,ik1x)  ; ik2x = MIN(kpi,ik2x)
       DO jj=1,kpj
          ik1y = jj-knj       ; ik2y = jj+knj
          ik1y = MAX(1,ik1y)  ; ik2y = MIN(kpj,ik2y)
          dl_den = SUM(kiw(ik1x:ik2x,ik1y:ik2y) )
          IF ( dl_den /= 0 ) THEN
             py(ji,jj) = SUM(px(ik1x:ik2x,ik1y:ik2y), mask=ll_mask(ik1x:ik2x,ik1y:ik2y) )/dl_den
          ELSE
             py(ji,jj) = rspval
          ENDIF
       END DO
    END DO

  END SUBROUTINE lisbox



END MODULE common_smooth2d
