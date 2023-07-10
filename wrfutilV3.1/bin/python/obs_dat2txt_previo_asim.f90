PROGRAM obsdump
  IMPLICIT NONE
  REAL(4) :: wk(7)      ! los archivos de obs.dat previos a asimilacion solo
                        ! tienen 7 columnas
  INTEGER :: ios
  CHARACTER(1) :: S
  INTEGER :: ii

  character(100) :: dir
! ******** Leer argumentos de entrada
  integer :: num_args, ix
  character(len=200), dimension(:), allocatable :: args ! Tamano maximo de cada argumento

  num_args = command_argument_count()
  call get_command_argument(1,dir)





  OPEN(UNIT=23,FILE=dir,FORM='unformatted')

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
