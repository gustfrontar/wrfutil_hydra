PROGRAM dec_prepbufr
!
! NOTE: output goes to fort.90
!
!  USE common
!  USE common_obs_wrf

  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  REAL(r_size),PARAMETER :: pi=3.1415926535d0
  REAL(r_size),PARAMETER :: gg=9.81d0
  REAL(r_size),PARAMETER :: rd=287.0d0
  REAL(r_size),PARAMETER :: cp=1005.7d0
  REAL(r_size),PARAMETER :: re=6371.3d3
  REAL(r_size),PARAMETER :: r_omega=7.292d-5
  REAL(r_size),PARAMETER :: t0c=273.15d0
  REAL(r_size),PARAMETER :: undef=9.99d33

  CHARACTER(20)  :: minlonchar , maxlonchar , minlatchar , maxlatchar


  REAL :: minlon , maxlon , minlat , maxlat
  INTEGER,PARAMETER :: maxlev = 255     !Maximum number of BUFR levels
  INTEGER,PARAMETER :: maxevn = 10      !Maximum number of BUFR event sequences
!  CHARACTER(MXSTRL) :: head = 'SID XOB YOB DHR ELV TYP T29 ITP'
!  CHARACTER(MXSTRL) :: ostr(MXR8VT) = (/'POB PQM PPC PRC PFC PAN CAT',&
!                                    &   'QOB QQM QPC QRC QFC QAN CAT',&
!                                    &   'TOB TQM TPC TRC TFC TAN CAT',&
!                                    &   'ZOB ZQM ZPC ZRC ZFC ZAN CAT',&
!                                    &   'UOB WQM WPC WRC UFC UAN CAT',&
!                                    &   'VOB WQM WPC WRC VFC VAN CAT'/)
  CHARACTER(8) :: obtype

  INTEGER,PARAMETER :: nid_obs=7
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
  INTEGER,PARAMETER :: id_tv_obs=3079
!
! surface observations codes > 9999
!
  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_rain_obs=19999
  INTEGER,PARAMETER :: id_tclon_obs=99991
  INTEGER,PARAMETER :: id_tclat_obs=99992
  INTEGER,PARAMETER :: id_tcmip_obs=99993

!  CHARACTER :: var(8) = (/'P','Q','T','Z','U','V'/)
  INTEGER,PARAMETER :: nobtype = 20
  INTEGER :: IREADNS
  INTEGER :: idate,idummy,nilev,n
  INTEGER :: ilev,ievn
  INTEGER :: iobs(nobtype+1),iobs_out(nobtype+1)
  CHARACTER(6) :: obtypelist(nobtype)
  REAL(r_dble) :: station(5)
  CHARACTER(8) :: cs
  REAL(r_dble) :: prs(4,maxlev,maxevn)
  REAL(r_dble) :: obs(4,maxlev,maxevn)
  REAL(r_sngl) :: wk(7)
  INTEGER :: iunit
!  INTEGER,PARAMETER :: id_tv_obs=3079
  real(r_sngl) :: vtcd
  integer :: vtcdi  

  character(32) :: primera
  character(62) :: segunda
  character(6) :: ot

  ! Get domain limits

  CALL getarg(1,minlonchar)
  CALL getarg(2,maxlonchar)
  CALL getarg(3,minlatchar)
  CALL getarg(4,maxlatchar)

  WRITE(*,*)"Domain limits are:"
  WRITE(*,*)"MINLON=",minlonchar
  WRITE(*,*)"MAXLON=",maxlonchar
  WRITE(*,*)"MINLAT=",minlatchar
  WRITE(*,*)"MAXLAT=",maxlatchar


  READ(minlonchar,*)minlon
  READ(maxlonchar,*)maxlon
  READ(minlatchar,*)minlat
  READ(maxlatchar,*)maxlat

  if(  maxlat <= minlat .OR. minlon < 0 .OR. maxlon < 0 .OR. minlon > 360 .OR. maxlon > 360  & 
       .OR. maxlat > 90 .OR. minlat < -90 )then
    WRITE(*,*)"Error!: Input domain is not correct."
    WRITE(*,*)"Please insert input domain as dec_prepbufr minlon maxlon minlat maxlat" 
    WRITE(*,*)"minlon and maxlon should be from 0 to 360 with the westernmost longitude first"
    STOP
  endif


  !
  ! Open the input file
  !
  OPEN(11,FILE='prepbufr.in',FORM='unformatted',CONVERT='little_endian')
!  print *, 'abro el prepbufr.in, llamo a OPENBF'
  CALL OPENBF(11,'IN',11)
!  print *, 'llamo a DATELEN'
  CALL DATELEN(10)    ! YYYYMMDDHH
!  print *, 'llamo a UFBQCD'
  call UFBQCD(11,'VIRTMP',vtcd)
  vtcdi = nint(vtcd)  
  
  !
  ! Main loop
  !
!  print *, 'empiezo main loop'
  obtypelist = (/'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR',&
               & 'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG',&
               & 'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND',&
               & 'MSONET', 'GPSIPW', 'RASSDA', 'WDSATR', 'ASCATW'/)
  iobs = 0
  iobs_out = 0
  DO
    !
    ! next record
    !
    IF(IREADNS(11,obtype,idate) /= 0) EXIT
    print *, 'obtype= , idate= ',obtype,idate ! for debugging
!    READ(*,*)            ! for debugging
    DO n=1,nobtype+1
!    print *, 'n=', n
      IF(obtype == obtypelist(n) .OR. n > nobtype) THEN
        iobs(n) = iobs(n)+1
        EXIT
      END IF
    END DO
    wk(7) = REAL(n)
!    IF(n /= 15) CYCLE      ! for debugging
    !
    ! station location (lon, lat)
    !
    CALL UFBINT(11,station,5,1,idummy,'SID XOB YOB ELV DHR')
    WRITE(cs(1:8),'(A8)') station(1)
!    IF(n < 2) CYCLE                                       ! for debugging
!    print '(A,4F10.3,I)',cs,station(2:5),NINT(station(5)) ! for debugging
!    READ(*,*)                                             ! for debugging
    wk(2:3) = station(2:3)
 IF( maxlon > minlon )THEN 
    IF(wk(2) < minlon .OR. maxlon < wk(2) .OR. &
     & wk(3) < minlat .OR. maxlat < wk(3)) CYCLE ! domain check
 ELSE
     IF( ( wk(2) < minlon .AND.  wk(2) > maxlon ) .OR. &
     & wk(3) < minlat .OR. maxlat < wk(3)) CYCLE ! domain check   
 ENDIF
    wk(4) = station(4)
    IF(NINT(station(5)) < -3 .OR. 3 < NINT(station(5))) CYCLE
    iunit = 90+NINT(station(5))
    !
    ! obs
    !
    CALL UFBEVN(11,prs,4,maxlev,maxevn,nilev,'POB POE PQM PPC')
    IF(obtype == 'ADPSFC' .OR.&
     & obtype == 'SFCSHP' .OR.&
     & obtype == 'SFCBOG') CALL output_ps ! surface pressure report
    CALL UFBEVN(11,obs,4,maxlev,maxevn,ilev,'QOB QOE QQM QPC')
    CALL output(id_q_obs)
    CALL UFBEVN(11,obs,4,maxlev,maxevn,ilev,'TOB TOE TQM TPC')
    CALL output(id_t_obs)
    CALL UFBEVN(11,obs,4,maxlev,maxevn,ilev,'UOB WOE WQM WPC')
    CALL output(id_u_obs)
    CALL UFBEVN(11,obs,4,maxlev,maxevn,ilev,'VOB WOE WQM WPC')
    CALL output(id_v_obs)
  END DO

  PRINT '(A)','================================================================================'
  PRINT '(A)','                  SYNOPSIS OF PREPBUFR DECODER using BUFRLIB'
  PRINT '(A)','                              by TAKEMASA MIYOSHI'
  PRINT '(A)','--------------------------------------------------------------------------------'
  PRINT '(A,I10)',' TOTAL NUMBER OF READ-IN RECORDS:',SUM(iobs)
  PRINT '(A,I10,A)',' TOTAL NUMBER OF WRITTEN RECORDS:',SUM(iobs_out),' (levels/variables separately)'
  PRINT '(A)','--------------------------------------------------------------------------------'
  PRINT '(10(2X,A))',obtypelist(1:10)
  PRINT '(10I8)',iobs(1:10)
  PRINT '(10I8)',iobs_out(1:10)
  PRINT '(A)','--------------------------------------------------------------------------------'
  PRINT '(10(2X,A))',obtypelist(11:20)
  PRINT '(10I8)',iobs(11:20)
  PRINT '(10I8)',iobs_out(11:20)
  PRINT '(A)','--------------------------------------------------------------------------------'
  PRINT '(2X,A)','OTHERS'
  PRINT '(I8)',iobs(21)
  PRINT '(I8)',iobs_out(21)
  PRINT '(A)','================================================================================'

!----------------------------------------------------------------------
! to save this information in a txt

! (al principio defino estas variables)

 primera='TOTAL NUMBER OF READ-IN RECORDS:'
 segunda=' TOTAL NUMBER OF WRITTEN RECORDS(levels/variables separately):'
 ot='OTHERS'

 open(10,status='UNKNOWN',form='formatted',file='info')

 WRITE(10,100) primera,SUM(iobs),segunda,SUM(iobs_out)
 WRITE(10,200) obtypelist(1:20),ot
 WRITE(10,300) iobs(1:21)
 WRITE(10,300) iobs_out(1:21)

100     format (A32,I10,A61,I10)
200     format (21(2X,A))
300     format (21(I8))
!----------------------------------------------------------------------

  STOP
CONTAINS
SUBROUTINE output(id)
  INTEGER,INTENT(IN) :: id
  INTEGER :: iqm, iseq

  IF(ilev /= nilev) THEN
    PRINT *,'FATAL ERROR, nilev /= ilev',nilev,ilev
    STOP
  END IF
!  wk(1) = id

  DO ilev=1,nilev
    wk(1) = id    
    wk(4) = prs(1,ilev,1) ! hPa
    iqm = NINT(prs(3,ilev,1))
    IF(iqm < 0 .OR. 2 < iqm) CYCLE
    wk(5:6) = obs(1:2,ilev,1)
    IF(id == id_q_obs) THEN
      wk(5) = wk(5) * 1.E-6 ! mg/kg -> kg/kg
      wk(6) = MAX(wk(5)*wk(6)*0.15,1.0E-7)
    END IF

! ---> VT & T
!    IF(id == id_t_obs) wk(5) = wk(5) + t0c
    IF(id == id_t_obs) then
      wk(5) = wk(5) + t0c
      do iseq = 1,10
        if (nint(obs(4,ilev,iseq)) == vtcdi) then
!  write (*, '(A,2I6,2F12.3)') obtype, ilev, iseq, obs(4,ilev,iseq), wk(5)	
          wk(1) = id_tv_obs
          exit
        end if
      end do
    END IF
! <--- VT & T    

    iqm = NINT(obs(3,ilev,1))
    IF(iqm < 0 .OR. 2 < iqm) CYCLE
    IF(wk(6) > 1.E10) CYCLE
    WRITE(iunit) wk
    iobs_out(n) = iobs_out(n) + 1
  END DO

  RETURN
END SUBROUTINE output

SUBROUTINE output_ps
  INTEGER :: iqm

  wk(1) = id_ps_obs
  iqm = NINT(prs(3,1,1))
  IF(iqm < 0 .OR. 2 < iqm) RETURN
  wk(5:6) = prs(1:2,1,1)
  IF(wk(6) > 1.E10) RETURN
  WRITE(iunit) wk
  iobs_out(n) = iobs_out(n) + 1

  RETURN
END SUBROUTINE output_ps

END PROGRAM dec_prepbufr
