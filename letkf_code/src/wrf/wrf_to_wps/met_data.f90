module met_data_module
   use common

   ! Derived types
   type met_data
      integer                       :: version, nx, ny, iproj, nlev
      real(r_sngl)                  :: xfcst, startlat, startlon, starti, startj, &
                                       deltalat, deltalon, dx, dy, xlonc, stdlon, &
                                       truelat1 , truelat2 , earth_radius , level
      real(r_sngl),allocatable      :: slab(:,:,:),xlvl(:)
      character (len=9)             :: field
      character (len=24)            :: hdate
      character (len=25)            :: units
      character (len=32)            :: map_source
      character (len=46)            :: desc
      character (len=8)             :: startloc
      logical                       :: is_wind_grid_rel
   end type met_data

   type(met_data)                   :: output_data

CONTAINS

subroutine WRITE_SLAB(imt_grd)
  IMPLICIT NONE
  INTEGER            :: i,j,k
  INTEGER,INTENT(IN) :: imt_grd

  DO k= 1,output_data % nlev 
  WRITE(imt_grd)output_data % version
  WRITE(imt_grd)output_data % hdate,     &
                output_data % xfcst,     &
                output_data % map_source,&
                output_data % field,     &
                output_data % units,     &
                output_data % desc,      &
                output_data % xlvl(k),   &
                output_data % nx,        &
                output_data % ny,        &
                output_data % iproj
  IF(output_data % iproj == 1)THEN !MERCATOR
                             WRITE(imt_grd) output_data % startloc, &
                                            output_data % startlat, &
                                            output_data % startlon, &
                                            output_data % dx,       &
                                            output_data % dy,       &
                                            output_data % truelat1, &
                                            output_data % earth_radius
  ENDIF

  IF( output_data % iproj == 4)THEN !LATLON
                             WRITE(imt_grd) output_data % startloc, &
                                            output_data % startlat, &
                                            output_data % startlon, &
                                            output_data % deltalat, &
                                            output_data % deltalon, &
                                            output_data % earth_radius
  ENDIF

  IF( output_data % iproj == 3)THEN !LAMBERT
                             WRITE(imt_grd) output_data % startloc, &
                                            output_data % startlat, &
                                            output_data % startlon, &
                                            output_data % dx,       &
                                            output_data % dy,       &
                                            output_data % stdlon,   &
                                            output_data % truelat1, &
                                            output_data % truelat2, &
                                            output_data % earth_radius
  ENDIF

  WRITE(imt_grd)output_data % is_wind_grid_rel
  WRITE(imt_grd)output_data % slab(:,:,k)
  END DO
  
  
                                                                                                                                                                
end subroutine WRITE_SLAB


subroutine READ_SLAB(imt_grd)
  IMPLICIT NONE
  INTEGER            :: i,j,k
  INTEGER,INTENT(IN) :: imt_grd
  REAL(r_sngl)       :: level
  INTEGER            :: reading_ierror

  !WRITE IN WPS FORMAT

  k=0
  DO 
  k=k+1
  WRITE(*,*)' Going to read record ',k
  READ(imt_grd, iostat = reading_ierror )output_data % version
  IF ( reading_ierror /= 0 )THEN
     WRITE(*,*)'It seems that we reached the end of the file ...'
     RETURN
  ENDIF
  READ(imt_grd)output_data % hdate,     &
               output_data % xfcst,     &
               output_data % map_source,&
               output_data % field,     &
               output_data % units,     &
               output_data % desc,      &
               output_data % level,     &
               output_data % nx,        &
               output_data % ny,        &
               output_data % iproj

     WRITE(*,*)'Header for record =',k
     WRITE(*,*)'WPSVersion=',output_data % version
     WRITE(*,*)'Date=',output_data % hdate
     WRITE(*,*)'FCST=',output_data % xfcst
     WRITE(*,*)'Map source=',output_data % map_source
     WRITE(*,*)'Field=',output_data % field
     WRITE(*,*)'Units=',output_data % units
     WRITE(*,*)'Desc=',output_data % desc
     WRITE(*,*)'Level=',output_data % level
     WRITE(*,*)'Nx=',output_data % nx
     WRITE(*,*)'Ny=',output_data % ny
     WRITE(*,*)'Proj type =',output_data % iproj

  IF(output_data % iproj == 1)THEN !MERCATOR
                             READ(imt_grd) output_data % startloc, &
                                            output_data % startlat, &
                                            output_data % startlon, &
                                            output_data % dx,       &
                                            output_data % dy,       &
                                            output_data % truelat1, &
                                            output_data % earth_radius

     WRITE(*,*)'Mercator projection detected for record =',k
     WRITE(*,*)output_data % startloc
     WRITE(*,*)output_data % startlat
     WRITE(*,*)output_data % startlon
     WRITE(*,*)output_data % dx
     WRITE(*,*)output_data % dy
     WRITE(*,*)output_data % truelat1
     WRITE(*,*)output_data % earth_radius


  ENDIF

  IF( output_data % iproj == 0 )THEN !LATLON
                             READ(imt_grd) output_data % startloc, &
                                            output_data % startlat, &
                                            output_data % startlon, &
                                            output_data % deltalat, &
                                            output_data % deltalon, &
                                            output_data % earth_radius

     WRITE(*,*)'Latlon projection detected for record =',k
     WRITE(*,*)output_data % startloc
     WRITE(*,*)output_data % startlat
     WRITE(*,*)output_data % startlon
     WRITE(*,*)output_data % deltalat
     WRITE(*,*)output_data % deltalon
     WRITE(*,*)output_data % earth_radius

  ENDIF

  IF( output_data % iproj == 3)THEN !LAMBERT
                             READ(imt_grd) output_data % startloc, &
                                            output_data % startlat, &
                                            output_data % startlon, &
                                            output_data % dx,       &
                                            output_data % dy,       &
                                            output_data % stdlon,   &
                                            output_data % truelat1, &
                                            output_data % truelat2, &
                                            output_data % earth_radius

     WRITE(*,*)'Lambert projection detected for record =',k
     WRITE(*,*)output_data % startloc
     WRITE(*,*)output_data % startlat
     WRITE(*,*)output_data % startlon
     WRITE(*,*)output_data % dx
     WRITE(*,*)output_data % dy
     WRITE(*,*)output_data % stdlon
     WRITE(*,*)output_data % truelat1
     WRITE(*,*)output_data % truelat2
     WRITE(*,*)output_data % earth_radius

  ENDIF

  READ(imt_grd)output_data % is_wind_grid_rel
  ALLOCATE( output_data % slab( output_data % nx , output_data % ny , 1 ) )
  WRITE(*,*)output_data % is_wind_grid_rel 
  READ(imt_grd) output_data % slab

  WRITE(*,*)'Data sample for record =',k

  WRITE(*,*)output_data % is_wind_grid_rel
  WRITE(*,*)'Min data ',MINVAL( output_data % slab ) , 'Max data ',MAXVAL( output_data % slab )

  DEALLOCATE( output_data % slab )

  ENDDO


end subroutine READ_SLAB

subroutine GET_SLAB_RECORD( imt_grd , slab_data , reading_ierr )
  IMPLICIT NONE
  INTEGER            :: i,j,k
  INTEGER,INTENT(IN) :: imt_grd     !An already opened file in WPS intermediate format
  REAL(r_sngl)       :: level
  TYPE(met_data), INTENT(OUT)      :: slab_data
  INTEGER       , INTENT(OUT)      :: reading_ierr

  !WRITE IN WPS FORMAT

  READ(imt_grd, iostat = reading_ierr )slab_data % version
  READ(imt_grd, iostat = reading_ierr )slab_data % hdate,     &
               slab_data % xfcst,     &
               slab_data % map_source,&
               slab_data % field,     &
               slab_data % units,     &
               slab_data % desc,      &
               slab_data % level,     &
               slab_data % nx,        &
               slab_data % ny,        &
               slab_data % iproj

  IF(slab_data % iproj == 1)THEN !MERCATOR
                             READ(imt_grd , iostat = reading_ierr )  slab_data % startloc, &
                                            slab_data % startlat, &
                                            slab_data % startlon, &
                                            slab_data % dx,       &
                                            slab_data % dy,       &
                                            slab_data % truelat1, &
                                            slab_data % earth_radius
  ENDIF

    IF( slab_data % iproj == 0 )THEN !LATLON
                             READ(imt_grd , iostat = reading_ierr ) slab_data % startloc, &
                                            slab_data % startlat, &
                                            slab_data % startlon, &
                                            slab_data % deltalat, &
                                            slab_data % deltalon, &
                                            slab_data % earth_radius
  ENDIF

  IF( slab_data % iproj == 3)THEN !LAMBERT
                             READ(imt_grd , iostat = reading_ierr ) slab_data % startloc, &
                                            slab_data % startlat, &
                                            slab_data % startlon, &
                                            slab_data % dx,       &
                                            slab_data % dy,       &
                                            slab_data % stdlon,   &
                                            slab_data % truelat1, &
                                            slab_data % truelat2, &
                                            slab_data % earth_radius
  ENDIF

  READ(imt_grd , iostat = reading_ierr )slab_data % is_wind_grid_rel
  ALLOCATE( slab_data % slab( slab_data % nx , slab_data % ny , 1 ) )
  READ(imt_grd , iostat = reading_ierr ) slab_data % slab

  ALLOCATE( slab_data % xlvl( 1 ) )
  slab_data % xlvl = slab_data % level


end subroutine GET_SLAB_RECORD


subroutine PUT_SLAB_RECORD(imt_grd , slab_data )
  IMPLICIT NONE
  INTEGER            :: i,j,k
  INTEGER,INTENT(IN) :: imt_grd
  TYPE(met_data), INTENT(IN)      :: slab_data

  DO k= 1,slab_data % nlev
  WRITE(imt_grd)slab_data % version
  WRITE(imt_grd)slab_data % hdate,     &
                slab_data % xfcst,     &
                slab_data % map_source,&
                slab_data % field,     &
                slab_data % units,     &
                slab_data % desc,      &
                slab_data % xlvl(k),   &
                slab_data % nx,        &
                slab_data % ny,        &
                slab_data % iproj
  IF(slab_data % iproj == 1)THEN !MERCATOR
                             WRITE(imt_grd) slab_data % startloc, &
                                            slab_data % startlat, &
                                            slab_data % startlon, &
                                            slab_data % dx,       &
                                            slab_data % dy,       &
                                            slab_data % truelat1, &
                                            slab_data % earth_radius
  ENDIF

  IF( slab_data % iproj == 4)THEN !LATLON
                             WRITE(imt_grd) slab_data % startloc, &
                                            slab_data % startlat, &
                                            slab_data % startlon, &
                                            slab_data % deltalat, &
                                            slab_data % deltalon, &
                                            slab_data % earth_radius
  ENDIF

  IF( slab_data % iproj == 3)THEN !LAMBERT
                             WRITE(imt_grd) slab_data % startloc, &
                                            slab_data % startlat, &
                                            slab_data % startlon, &
                                            slab_data % dx,       &
                                            slab_data % dy,       &
                                            slab_data % stdlon,   &
                                            slab_data % truelat1, &
                                            slab_data % truelat2, &
                                            slab_data % earth_radius
  ENDIF

  WRITE(imt_grd)slab_data % is_wind_grid_rel
  WRITE(imt_grd)slab_data % slab(:,:,k)
  END DO



end subroutine PUT_SLAB_RECORD

  
  
SUBROUTINE GET_SLP(NX,NY,NZ,HGT,P,T,Q,TOPO,SLP)

!SLP computation routine from UPP


!                .      .    .     
! SUBPROGRAM:    NGMSLP      NMC SEA LEVEL PRESSURE REDUCTION
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-02-02       
!     
! ABSTRACT:
!     
!     THIS ROUTINE COMPUTES SEA LEVEL PRESSURE USING THE
!     HYDROSTATIC EQUATION WITH THE SHUELL CORRECTION.  THE
!     FOLLOWING IS BASED ON DOCUMENTATION IN SUBROUTINE 
!     OUTHYDRO OF THE NGM:
!     
!     THE FUNDAMENTAL HYDROSTATIC EQUATION IS
!        D(HEIGHT)
!        ---------  =  TAU = VIRTUAL TEMPERATURE * (RGAS/GRAVITY)
!        D (Z)
!      WHERE
!        Z = MINUS LOG OF PRESSURE (-LN(P)).
!
!     SEA-LEVEL PRESSURE IS COMPUTED FROM THE FORMULA
!        PRESS(MSL) = PRESS(GROUND) * EXP( F)
!     WHERE
!        F        = HEIGHT OF GROUND / MEAN TAU
!        MEAN TAU = ( TAU(GRND) + TAU(SL) ) / 2
!     
!     IN THE NGM TAU(GRND) AND TAU(SL) ARE FIRST SET USING A 
!     6.5DEG/KM LAPSE RATE FROM LOWEST MDL LEVEL.  THIS IS MODIFIED
!     BY A CORRECTION BASED ON THE CRITICAL TAU OF THE SHUELL
!     CORRECTION:
!                  TAUCR=(RGASD/GRAVITY) * 290.66
!   
!     1) WHERE ONLY TAU(SL) EXCEEDS TAUCR, CHANGE TAU(SL) TO TAUCR.
!
!     2) WHERE BOTH TAU(SL) AND TAU(GRND) EXCEED TAUCR,
!        CHANGE TAU(SL) TO TAUCR-CONST*(TAU(GRND)-TAUCR  )**2
!        WHERE CONST = .005 (GRAVITY/RGASD)
!   
!     THE AVERAGE OF TAU(SL) AND TAU(GRND) IS THEN USED TOGETHER
!     WITH THE GROUND HEIGHT AND PRESSURE TO DERIVE THE PRESSURE
!     AT SEA LEVEL. 
!     
!     HEIGHT OF THE 1000MB SURFACE IS COMPUTED FROM THE MSL PRESSURE
!     FIELD USING THE FORMULA:
!     
!       P(MSL) - P(1000MB) = MEAN DENSITY * GRAVITY * HGT(1000MBS)
!     
!     WHERE P(MSL) IS THE SEA LEVEL PRESSURE FIELD WE HAVE JUST
!     COMPUTED.
!     
!
!     MEB 6/13/02: THIS CODE HAS BEEN SIMPLIFIED CONSIDERABLY FROM
!     THE ONE USED IN ETAPOST.  HORIZONTAL SMOOTHING HAS BEEN
!     REMOVED AND THE FIRST MODEL LEVEL IS USED RATHER
!     THAN THE MEAN OF THE VIRTUAL TEMPERATURES IN
!     THE LOWEST 30MB ABOVE GROUND TO COMPUTE TAU(GRND).
!     
!   .     
!     
! PROGRAM HISTORY LOG:
!   93-02-02  RUSS TREADON
!   98-06-08  T BLACK - CONVERSION FROM 1-D TO 2-D
!   00-01-04  JIM TUCCILLO - MPI VERSION
!   01-10-25  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   01-11-02  H CHUANG - MODIFIED LINE 234 FOR COMPUTATION OF 
!                         SIGMA/HYBRID SLP
!   01-12-18  H CHUANG - INCLUDED SMOOTHING ALONG BOUNDARIES TO BE
!                         CONSISTENT WITH MESINGER SLP
!   02-06-13  MIKE BALDWIN - WRF VERSION
!   06-12-18  H CHUANG - BUG FIX TO CORRECT TAU AT SFC
!   19-09-23  J RUIZ - COUPLED WITH WRF_TO_WPS TOOL
! USAGE:    CALL NGMSLP
!   INPUT ARGUMENT LIST:
!     NONE     
!
!   OUTPUT ARGUMENT LIST: 
!     NONE
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       NONE
!     LIBRARY:
!       COMMON   - CTLBLK
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!     
      !use vrbls3d,    only: zint, pint, t, q, zmid
      !use vrbls2d,    only: slp, fis, z1000 (sea level pressure , surfac height, 1000 hPa geopotential heigth)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!     
      real(r_sngl),PARAMETER :: ZSL=0.0
      real(r_sngl),PARAMETER :: TAUCR=rd*290.66/gg,CONST=0.005*gg/rd
      real(r_sngl),PARAMETER :: GORD=gg/rd,DP=60.E2
      real(r_sngl),PARAMETER :: D608=0.608 , D50=0.5 , P1000=1000.E2
      real(r_sngl),PARAMETER :: GAMMAL=6.5E-3 , H1=1.0 
     
      integer, INTENT(IN)     :: NX,NY,NZ 
      real(r_sngl),INTENT(IN) :: HGT(nx,ny,nz),TOPO(nx,ny)
      real(r_sngl),INTENT(IN) :: P(nx,ny,nz),T(nx,ny,nz),Q(nx,ny,nz)
      real(r_sngl),INTENT(OUT):: SLP(nx,ny)
      
      integer I,J
      real ZSFC,PSFC,TVRT,TAU,TVRSFC,TAUSFC,TVRSL,TAUSL,TAUAVG
           
!     
!**********************************************************************
!     START NGMSLP HERE.
!     
!     LOOP OVER HORIZONTAL GRID.
!
!!$omp  parallel do
!!$omp& private(psfc,tau,tauavg,tausfc,tausl,tvrsfc,tvrsl,
!!$omp&         tvrt,zsfc)
       DO J=1,NY
       DO I=1,NX
         PSFC = P(I,J,1)
         SLP(I,J) = PSFC 
         IF ( ABS(TOPO(I,J)) .GT. 1.0 ) THEN
            ZSFC = HGT(I,J,1)
!           COMPUTE LAYER TAU (VIRTUAL TEMP*RD/G).
            TVRT = T(I,J,1)*(H1+D608*Q(I,J,1))
            TAU  = TVRT*rd/gg
!     
!           COMPUTE TAU AT THE GROUND (Z=ZSFC) AND SEA LEVEL (Z=0)
!           ASSUMING A CONSTANT LAPSE RATE OF GAMMAL=6.5DEG/KM.
            TVRSFC = TVRT + (HGT(I,J,1) - ZSFC)*GAMMAL ! Chuang
            TAUSFC = TVRSFC*rd/gg
            TVRSL  = TVRT + (HGT(I,J,1) - ZSL)*GAMMAL
            TAUSL  = TVRSL*rd/gg
!     
!           IF NEED BE APPLY SHEULL CORRECTION.
            IF ((TAUSL.GT.TAUCR).AND.(TAUSFC.LE.TAUCR)) THEN
               TAUSL=TAUCR
            ELSEIF ((TAUSL.GT.TAUCR).AND.(TAUSFC.GT.TAUCR)) THEN
               TAUSL = TAUCR-CONST*(TAUSFC-TAUCR)**2
            ENDIF
!     
!           COMPUTE MEAN TAU.
            TAUAVG = D50*(TAUSL+TAUSFC)
!     
!           COMPUTE SEA LEVEL PRESSURE.
            IF (ABS(TOPO(I,J)).GT.1.0)SLP(I,J) = PSFC*EXP(ZSFC/TAUAVG)
            !WRITE(*,*)SLP(I,J),PSFC,EXP(ZSFC/TAUAVG),ZSFC
!     
        ENDIF    
!     MOVE TO NEXT HORIZONTAL GRIDPOINT.
      ENDDO
      ENDDO
!!$omp end parallel do      
!     
!     
!     END OF ROUTINE.
!     

      RETURN
END SUBROUTINE GET_SLP
                                                                                                                     

end module met_data_module

