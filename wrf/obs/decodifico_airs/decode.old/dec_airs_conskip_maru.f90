  PROGRAM dec_airsret
  USE common
  USE common_obs_wrf

  IMPLICIT NONE

!  REAL :: pepe
  LOGICAL,PARAMETER :: verbose = .TRUE.
  INTEGER,PARAMETER :: airstype = 21 ! observation type record for AIRS
! dominio que usamos para decodificar los prepbufr:
! poner valores de longitud de 0 a 360, despues hace una conversion para que
! sean compatibles con los airs que vienen de -180 a 180
  REAL(r_size),PARAMETER :: lon1 =  250.d0 ! minimum longitude
  REAL(r_size),PARAMETER :: lon2 =  340.d0 ! maximum longitude
  REAL(r_size),PARAMETER :: lat1 =  -70.d0 ! minimum latitude
  REAL(r_size),PARAMETER :: lat2 =   20.d0 ! maximum latitude

! usamos el mismo dominio
  ! use full data in the following domain.
  ! (use thinned data in other domain (iskip)).
  REAL(r_size),PARAMETER :: loni1 =  250.0d0 ! minimum longitude (inner)
  REAL(r_size),PARAMETER :: loni2 =  340.0d0 ! maximum longitude (inner)
  REAL(r_size),PARAMETER :: lati1 =  -70.d0 ! minimum latitude (inner)
  REAL(r_size),PARAMETER :: lati2 =   20.d0 ! maximum latitude (inner)

  REAL(r_size),PARAMETER :: coef_error = 1.0d0  
  INTEGER,PARAMETER :: iskip = 2  ! terminas tomando un dato cada 4 (ver
!                                   cuaderno 2015)
  INTEGER,PARAMETER :: nunit=70
  INTEGER,PARAMETER :: nx=30    ! dimension de los airs (footprints cross track)
  INTEGER,PARAMETER :: ny=45    ! dimension de los airs (lines along track)
  INTEGER,PARAMETER :: nzt=28   ! niveles de los airs
  INTEGER,PARAMETER :: nzq=14   ! niveles de los airs
  INTEGER :: statn
  INTEGER :: fid
  INTEGER :: nswath
  INTEGER :: nchar
  INTEGER :: swid
  INTEGER :: swopen,swinqswath,swattach,swrdattr,swrdfld,swdetach,swclose
  INTEGER :: pepei,pepej
  CHARACTER(256) :: swathname
  CHARACTER(100) :: filename
  REAL(r_dble) :: rlon(nx,ny)
  REAL(r_dble) :: rlat(nx,ny)
  REAL(r_dble) :: time(nx,ny),time1,time2
  REAL(r_sngl) :: levt(nzt)
  REAL(r_sngl) :: levq0(nzq+1)
  REAL(r_sngl) :: levq(nzq)
  REAL(r_sngl) :: t(nzt,nx,ny)
  REAL(r_sngl) :: te(nzt,nx,ny)
  REAL(r_sngl) :: q(nzq,nx,ny)
  REAL(r_sngl) :: qe(nzq,nx,ny)
  REAL(r_sngl) :: psurf(nx,ny)
  REAL(r_sngl) :: pbest(nx,ny)
  REAL(r_sngl) :: pgood(nx,ny)
  INTEGER(2) :: qq(nx,ny)
  INTEGER :: i,j,k
  REAL(r_sngl) :: wk(7)
  INTEGER :: iy,im,id,ih,imin
  INTEGER :: iy1,im1,id1,ih1,imin1
  INTEGER :: iy2,im2,id2,ih2,imin2
  REAL(r_sngl) :: sec4
  REAL(r_size) :: sec
  LOGICAL :: ex
  CHARACTER(14) :: outfile='yyyymmddhh.dat'
!                           123456789012345
!  character(10) :: fecha1
!  character(10) :: fecha2
  INTEGER :: io
  INTEGER :: ii
  INTEGER :: n1,n2,n,nt,nq
  INTEGER :: ntot,npar

!=======================================================================
! READING HDF
!=======================================================================
  !
  ! open
  !
  READ(5,'(A)') filename
  PRINT *,filename
  fid = swopen(filename,1)
   print*,'fid',fid
  IF(fid == -1) THEN
    PRINT *,'IO ERROR: opening file ',filename
    STOP
  END IF
  nswath = swinqswath(filename,swathname,nchar)
  print*,'nswath',nswath
  IF(nswath /= 1) THEN
    PRINT *,'FILE ERROR: bad nswath ',nswath
    STOP
  END IF
  IF(swathname /= 'L2_Standard_atmospheric&surface_product') THEN
    PRINT *,'FILE ERROR: bad swath name ',swathname
    STOP
  END IF
  swid = swattach(fid, swathname)
  print*,'swid',swid
  IF(swid == -1) THEN
    PRINT *,'FILE ERROR: failed to attach to swath ',swathname
    STOP
  END IF
  !
  ! read
  !
  statn = swrdattr(swid,'start_year',iy)
  statn = swrdattr(swid,'start_month',im)
  statn = swrdattr(swid,'start_day',id)
  statn = swrdattr(swid,'start_hour',ih)
  print *, 'esta leyendo ih=', ih
  statn = swrdattr(swid,'start_minute',imin)
  statn = swrdattr(swid,'start_sec',sec4)
  sec = sec4

! cuando los statn dan -1 es porque no esta leyendo bien
! cuando dan 0 es porque esta todo ok
  statn = swrdattr(swid,'start_Time',time1)
  statn = swrdattr(swid,'end_Time',time2)
  statn = swrdfld(swid,'Longitude',(/0,0/),(/1,1/),(/nx,ny/),rlon)
  statn = swrdfld(swid,'Latitude',(/0,0/),(/1,1/),(/nx,ny/),rlat)
  statn = swrdfld(swid,'Time',(/0,0/),(/1,1/),(/nx,ny/),time)
  statn = swrdfld(swid,'pressStd',0,1,nzt,levt)
  statn = swrdfld(swid,'pressH2O',0,1,nzq+1,levq0)

  DO k=1,nzq
    levq(k) = exp(0.5d0 * (log(levq0(k)) + log(levq0(k+1))))
!    levq(k) = 0.5d0 * ((levq0(k)) + (levq0(k+1)))
  END DO

  statn = swrdfld(swid,'TAirStd',(/0,0,0/),(/1,1,1/),(/nzt,nx,ny/),t)
  statn = swrdfld(swid,'TAirStdErr',(/0,0,0/),(/1,1,1/),(/nzt,nx,ny/),te)
  statn = swrdfld(swid,'PSurfStd',(/0,0/),(/1,1/),(/nx,ny/),psurf)
  statn = swrdfld(swid,'PBest',(/0,0/),(/1,1/),(/nx,ny/),pbest)
  statn = swrdfld(swid,'PGood',(/0,0/),(/1,1/),(/nx,ny/),pgood)
  statn = swrdfld(swid,'H2OMMRStd',(/0,0,0/),(/1,1,1/),(/nzq,nx,ny/),q)
  q = q*0.001d0 ! g/kg -> kg/kg
  statn = swrdfld(swid,'H2OMMRStdErr',(/0,0,0/),(/1,1,1/),(/nzq,nx,ny/),qe)  
  qe = qe*0.001d0 ! g/kg -> kg/kg
  statn = swrdfld(swid,'Qual_H2O',(/0,0/),(/1,1/),(/nx,ny/),qq)

  !
  ! close
  !
  statn = swdetach(swid)
  statn = swclose(fid)
!=======================================================================
! QC -> OUTPUT LETKF OBS FORMAT
!=======================================================================
  wk(7) = REAL(airstype) ! TYPE CODE for AIRS retrieval

  CALL com_tai2utc(time1,iy1,im1,id1,ih1,imin1,sec)
  IF(imin1 > 30) CALL com_tai2utc(time1+1800.d0,iy1,im1,id1,ih1,imin1,sec)
  WRITE(outfile(1:10),'(I4.4,3I2.2)') iy1,im1,id1,ih1
  OPEN(nunit+ih1,FILE=outfile,FORM='unformatted',ACCESS='sequential')
  
  n1=0
  DO
    READ(nunit+ih1,IOSTAT=io) ( wk(ii),ii=1,7 )
    if (io>0) then
!     WRITE(*,*) 'Something wrong with output.dat'
     exit
   else if (io==0) then
     n1=n1+1
    end if
   ENDDO
print *, 'n1=', n1

  CALL com_tai2utc(time2,iy2,im2,id2,ih2,imin2,sec)
  IF(imin2 > 30) CALL com_tai2utc(time2+1800.d0,iy2,im2,id2,ih2,imin2,sec)
  IF(ih2 /= ih1) THEN
    WRITE(outfile(1:10),'(I4.4,3I2.2)') iy2,im2,id2,ih2
    OPEN(nunit+ih2,FILE=outfile,FORM='unformatted',ACCESS='sequential')
  END IF

   n2=0
   DO
    READ(nunit+ih2,IOSTAT=io) ( wk(ii),ii=1,7 )
    if (io>0) then
!     WRITE(*,*) 'Something wrong with output.dat'
     exit
    else if (io==0) then
     n2=n2+1
    end if
   ENDDO
print *, 'n2=', n2

!  open(2,status='scratch')
  n=max(n1,n2)
print*, 'n=', n

  nt=0
  nq=0
  npar=0

  DO j=1,ny
    DO i=1,nx
      IF(rlon(i,j) < 0) rlon(i,j) = rlon(i,j) + 360.0d0

      IF(rlon(i,j) < lon1) CYCLE
      IF(rlon(i,j) > lon2) CYCLE
      IF(rlat(i,j) < lat1) CYCLE
      IF(rlat(i,j) > lat2) CYCLE

! hago el thinning:
       ! IF(mod(i,iskip) /= 1) then
       ! print*,'salteo pepei=', pepei,'i=', i,'j=', j
       ! else
       ! print*,'NO salteo pepei=', pepei,'i=', i,'j=', j
       ! endif
       ! IF(mod(j,iskip) /= 1) then
       ! print*,'salteo pepej=', pepej,'i=', i,'j=', j
       ! else
       ! print*,'NO salteo pepej=', pepej,'i=', i,'j=', j
       ! endif
        IF(mod(i,iskip) /= 1) CYCLE
        IF(mod(j,iskip) /= 1) CYCLE

                 
      CALL com_tai2utc(time(i,j),iy,im,id,ih,imin,sec)
      IF(imin > 30) CALL com_tai2utc(time(i,j)+1800.d0,iy,im,id,ih,imin,sec)
 !veo que contador tengo que usar:     
!     if (ih==ih1) then
!     n=n1
!     else if (ih==ih2) then
!     n=n2
!     end if 
!    print*, 'n=', n

! primero leo el archivo y lo reescribo, si es que ya existe     
!  n=0 !ojo xq aca estoy adentro del do i,j, entonces va punto por punto. ah, tal
!vez por eso da los io>0, sera que busca el i,j pero no lo encuentra en el .dat,
!seguramente porque no tenia datos ahi
!DO
!    READ(nunit+ih,IOSTAT=io) ( wk(ii),ii=1,7 )
!     READ(nunit+ih) ( wk(ii),ii=1,7 ) !si pongo solo esto me da error xq no
!     sabe el largo
! si lo vuelvo a escribir en realidad despues del primer punto ya estaria
! sobreescrito el archivo, entonces pierdo info igual que antes
! tendria que buscar la forma de contar las lineas, y escribir lo nuevo a partir
! de la ultima linea
!   if (io>0) then
!     WRITE(*,*) 'Something wrong with output.dat'
!     exit
!    else if (io<0) then
!     GOTO 100
!    else if (io==0) then
!     n=n+1
!         WRITE(nunit+ih) ,wk(1),wk(2),wk(3),wk(4),wk(5),wk(6),wk(7)
!    end if

!  END DO
!print *, 'n=', n

! ahora me ocupo de traducir el nuevo hdf:
         wk(2) = rlon(i,j)
         wk(3) = rlat(i,j)

      ntot=n+npar
    print*, 'ntot antes de loop T:', ntot
      !
      ! T
      !
      wk(1) = id_t_obs
      DO k=1,nzt
        IF(levt(k) > psurf(i,j)) CYCLE !below surface
        IF(levt(k) > pbest(i,j)) CYCLE !QC quality control
        wk(4) = levt(k)
        wk(5) = t(k,i,j)
        wk(6) = te(k,i,j) * sqrt(coef_error)
     !   if (nt==0) then
     !   WRITE(nunit+ih,line(n+1)) wk  !aca poner formato para q escriba desde n
     !   else
        WRITE(nunit+ih) wk
 !       endif
!      print *, 'escribiendo T, nunit=', nunit, 'ih=', ih
        IF(verbose) WRITE(6,'(F5.0,7F8.2)') wk(1:6),pbest(i,j),pgood(i,j)
       nt=nt+1
      END DO
      !
      ! Q
      !
      wk(1) = id_q_obs
      DO k=1,nzq
        IF(qq(i,j) == 2) CYCLE !QC
        IF(levq0(k) > psurf(i,j)) CYCLE !below surface
        IF(qq(i,j) == 1 .AND. levq0(k) > pbest(i,j)) CYCLE !QC
        wk(4) = levq(k)
        wk(5) = q(k,i,j)
        wk(6) = qe(k,i,j) * sqrt(coef_error)
        WRITE(nunit+ih) wk
!      print *, 'escribiendo Q, nunit=', nunit, 'ih=', ih
        IF(verbose) WRITE(6,'(F5.0,3F8.2,2ES11.2,2F8.2,I2)') wk(1:6),levq0(k),pbest(i,j),qq(i,j)
      nq=nq+1
      END DO
 npar=nt+nq
 print *, 'nt=', nt
 print *, 'nq=', nq    
 print *, 'nparcial=', npar
    END DO
  END DO

  CLOSE(nunit+ih1)
  IF(ih2 /= ih1) CLOSE(nunit+ih2)

!  STOP
END PROGRAM dec_airsret
