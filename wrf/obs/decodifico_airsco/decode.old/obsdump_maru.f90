PROGRAM obsdump
  IMPLICIT NONE
  REAL(4) :: wk(7)
  INTEGER :: ios
  CHARACTER(1) :: S
  INTEGER :: ii
  OPEN(UNIT=23,FILE='./2012122205.dat',FORM='unformatted')

  DO
 1111   READ(23,end=1000) ( wk(ii),ii=1,7 )
    PRINT *,NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7)
      goto 1111
!    PRINT '(I6,2F7.2,F10.2,3ES12.2)',NINT(wk(1)),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7)
    PRINT '(A)','PRESS "S" TO STOP'
    READ(5,'(A1)') S
    IF(S == 'S') EXIT
  END DO

 1000  CLOSE(23)
  STOP
END PROGRAM obsdump
