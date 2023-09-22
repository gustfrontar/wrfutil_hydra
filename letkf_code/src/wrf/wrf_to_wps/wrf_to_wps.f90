PROGRAM wrf_to_wps
!Write WRF output in WPS format.
!This version reads the input from the LETKF sigma binary files and generates a file in the WPS binary intermediate format.
!The data can only be input in mercator or lat-lon proyections. (metgrid does not accept other proyections as input)
!WARNING! Surface data is initialized with dummy values 
!QCLOUD, QRAIN, etc are included in the file. Namelist variables has to be set in order to have it vertically interpolated by real routine.

  USE met_data_module
  USE common_wrf
  USE common

  IMPLICIT NONE
  INTEGER(r_sngl),PARAMETER :: grdfid = 50

  INTEGER :: nlev_ns,nlon_ns,nlat_ns,nlev_soil_ns
  INTEGER :: i,j,k
  INTEGER :: ncid
  !REAL(4),ALLOCATABLE :: dnw(:) 
  
  REAL(r_sngl),ALLOCATABLE :: ps(:,:)
  REAL(r_sngl),ALLOCATABLE :: u(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: v(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: w(:,:,:)  
  REAL(r_sngl),ALLOCATABLE :: p(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: t(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: ght(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: qv(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: qc(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: qr(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: qci(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: qs(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: qg(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: rh(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: smois(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: tslb(:,:,:)
  
  REAL(r_sngl),ALLOCATABLE :: t2(:,:)
  REAL(r_sngl),ALLOCATABLE :: q2(:,:)  
  REAL(r_sngl),ALLOCATABLE :: rh2(:,:)
  REAL(r_sngl),ALLOCATABLE :: u10(:,:)
  REAL(r_sngl),ALLOCATABLE :: v10(:,:)
  REAL(r_sngl),ALLOCATABLE :: landsea(:,:)
  REAL(r_sngl),ALLOCATABLE :: skintemp(:,:)
  REAL(r_sngl),ALLOCATABLE :: snow(:,:)
  REAL(r_sngl),ALLOCATABLE :: snowh(:,:)
  REAL(r_sngl),ALLOCATABLE :: seaice(:,:)
  REAL(r_sngl),ALLOCATABLE :: canwat(:,:)

  REAL(r_sngl),ALLOCATABLE :: topo(:,:)
  REAL(r_sngl),ALLOCATABLE :: slp(:,:)

  REAL(r_sngl),ALLOCATABLE :: AUX(:,:,:)
  CHARACTER(9) :: grdin = 'input.nc'       !Meteorological fields input.
  CHARACTER(10) :: grdout= 'output.grd'

  include 'netcdf.inc'

!
! On input nlon,nlat,nlev are the same as nx,ny,nz in namelist.wps and namelist.input. Internally the size
! of the arrays is referenced to nlon-1, nlat-1 and nlev-1 (to be consistent with ncio routines)
!
  call set_common_wrf(grdin)

! Get non-staggered dimensions

  nlon_ns = nlon-1
  nlat_ns = nlat-1
  nlev_ns = nlev-1
  nlev_soil_ns = nlev_soil

 ALLOCATE(output_data % xlvl(nlev))

  output_data % nx = nlon_ns
  output_data % ny = nlat_ns
  output_data % nlev = nlev_ns
  if     ( pcode == 3 )then !Mercator
    output_data % iproj = 1
  elseif ( pcode == 1 )then !Lambert
    output_data % iproj = 3
  elseif ( pcode == 0 )then !Latlon
    output_data % iproj = 4
  endif
  output_data % deltalat = 0e0
  output_data % deltalon = 0e0
  output_data % dx = pdx !20
  output_data % dy = pdy !20
  output_data % truelat1 = ptruelat1 !22.5
  output_data % truelat2 = ptruelat2 
  output_data % stdlon   = pstdlon
  output_data % hdate =  hdate ! WRF INPUT/OUTPUT FORMAT
  output_data % xfcst = 0.0
  DO k=1,nlev_ns
  output_data % xlvl(k) = k !znu(k) !1e5-REAL(k)   !Tengo que ver que hay que poner aca.
  ENDDO
  output_data % version = 5  !WPS
  output_data % map_source = 'WRF MODEL OUTPUT'
  output_data % earth_radius =  6371.229004
  output_data % startloc = 'SWCORNER'
  output_data % is_wind_grid_rel = .TRUE.
  output_data % startlat = lat(1,1)
  output_data % startlon = lon(1,1)

  ALLOCATE(output_data % slab (nlon_ns,nlat_ns,nlev_ns))

!  WRITE(6,'(A)') '*** grid information ***'
!  WRITE(6,'(3(2X,A,I5))') 'nlon =',nlon,'nlat =',nlat,'nlev =',nlev 
!  WRITE(6,'(3(2X,A,I5))') 'nlon_ns =',nlon_ns,'nlat_ns =',nlat_ns,'nlev_ns =',nlev_ns
!
! ALLOCATE valiables
!
  ALLOCATE(ps(nlon_ns,nlat_ns))
  ALLOCATE(slp(nlon_ns,nlat_ns))
  ALLOCATE(u(nlon_ns,nlat_ns,nlev_ns))
  ALLOCATE(v(nlon_ns,nlat_ns,nlev_ns))
  ALLOCATE(w(nlon_ns,nlat_ns,nlev_ns))  
  ALLOCATE(t(nlon_ns,nlat_ns,nlev_ns))
  ALLOCATE(p(nlon_ns,nlat_ns,nlev_ns))
  ALLOCATE(ght(nlon_ns,nlat_ns,nlev_ns))
  ALLOCATE(qv(nlon_ns,nlat_ns,nlev_ns))
  ALLOCATE(qc(nlon_ns,nlat_ns,nlev_ns))
  ALLOCATE(qr(nlon_ns,nlat_ns,nlev_ns))
  ALLOCATE(qci(nlon_ns,nlat_ns,nlev_ns))
  ALLOCATE(qs(nlon_ns,nlat_ns,nlev_ns))
  ALLOCATE(qg(nlon_ns,nlat_ns,nlev_ns))
  ALLOCATE(rh(nlon_ns,nlat_ns,nlev_ns))
  ALLOCATE(smois(nlon_ns,nlat_ns,nlev_soil_ns))
  ALLOCATE(tslb(nlon_ns,nlat_ns,nlev_soil_ns))  
  ALLOCATE(t2(nlon_ns,nlat_ns))
  ALLOCATE(q2(nlon_ns,nlat_ns))  
  ALLOCATE(rh2(nlon_ns,nlat_ns))
  ALLOCATE(landsea(nlon_ns,nlat_ns))
  ALLOCATE(skintemp(nlon_ns,nlat_ns))
  ALLOCATE(snow(nlon_ns,nlat_ns))
  ALLOCATE(snowh(nlon_ns,nlat_ns))
  ALLOCATE(seaice(nlon_ns,nlat_ns))
  ALLOCATE(canwat(nlon_ns,nlat_ns))
  

  ALLOCATE(topo(nlon_ns,nlat_ns))
  ALLOCATE(AUX(nlon,nlat,nlev))
!
! READ WRF DATA AND DE-STAGGER DATA
!
  WRITE(6,'(A)') 'READING: ', grdin
  CALL check_io(NF_OPEN(grdin,NF_NOWRITE,ncid))

  !!! U
  write(*,*)'Reading U'
  call read_var_wrf(ncid,iv3d_u,1,nlev,1,AUX,'3d')
  u=(AUX(1:nlon_ns,1:nlat_ns,1:nlev_ns)+AUX(2:nlon_ns+1,1:nlat_ns,1:nlev_ns))*0.5

  !!! V
  write(*,*)'Reading V'
  call read_var_wrf(ncid,iv3d_v,1,nlev,1,AUX,'3d')
  v=(AUX(1:nlon_ns,1:nlat_ns,1:nlev_ns)+AUX(1:nlon_ns,2:nlat_ns+1,1:nlev_ns))*0.5

  !!! W
  write(*,*)'Reading W'
  call read_var_wrf(ncid,iv3d_w,1,nlev,1,AUX,'3d')
  w=(AUX(1:nlon_ns,1:nlat_ns,1:nlev_ns)+AUX(1:nlon_ns,1:nlat_ns,2:nlev_ns+1))*0.5

  !!! T
  write(*,*)'Reading T'
  call read_var_wrf(ncid,iv3d_t,1,nlev,1,AUX,'3d')
  t=AUX(1:nlon_ns,1:nlat_ns,1:nlev_ns)

  !!! P
  write(*,*)'Reading P'
  call read_var_wrf(ncid,iv3d_p,1,nlev,1,AUX,'3d')
  p=AUX(1:nlon_ns,1:nlat_ns,1:nlev_ns)

  !!! GHT
  write(*,*)'Reading GHT'
  call read_var_wrf(ncid,iv3d_ph,1,nlev,1,AUX,'3d')
  ght=(AUX(1:nlon_ns,1:nlat_ns,1:nlev_ns)+AUX(1:nlon_ns,1:nlat_ns,2:nlev_ns+1))*0.5 / gg

  !!! QVAPOR
  write(*,*)'Reading QVAPOR'
  call read_var_wrf(ncid,iv3d_qv,1,nlev,1,AUX,'3d')
  qv=AUX(1:nlon_ns,1:nlat_ns,1:nlev_ns)

  !!! QCLOUD
  write(*,*)'Reading QCLOUD'
  call read_var_wrf(ncid,iv3d_qc,1,nlev,1,AUX,'3d')
  qc=AUX(1:nlon_ns,1:nlat_ns,1:nlev_ns) 

  !!! QRAIN
  write(*,*)'Reading QRAIN'
  call read_var_wrf(ncid,iv3d_qr,1,nlev,1,AUX,'3d')
  qr=AUX(1:nlon_ns,1:nlat_ns,1:nlev_ns)

  !!! QICE
  write(*,*)'Reading QICE'
  call read_var_wrf(ncid,iv3d_qci,1,nlev,1,AUX,'3d')
  qci=AUX(1:nlon_ns,1:nlat_ns,1:nlev_ns)

  !!! QSNOW
  write(*,*)'Reading QSNOW'
  call read_var_wrf(ncid,iv3d_qs,1,nlev,1,AUX,'3d')
  qs=AUX(1:nlon_ns,1:nlat_ns,1:nlev_ns)

  !!! QGRAUPEL
  write(*,*)'Reading QGRAUPEL'
  call read_var_wrf(ncid,iv3d_qg,1,nlev,1,AUX,'3d')
  qg=AUX(1:nlon_ns,1:nlat_ns,1:nlev_ns)

  !!! SMOIS
  write(*,*)'Reading SMOIS'
  call read_var_wrf(ncid,is3d_sm,1,nlev_soil_ns,1,AUX(:,:,1:nlev_soil_ns),'so')
  smois=AUX(1:nlon_ns,1:nlat_ns,1:nlev_soil_ns)

  !!! TSLB
  write(*,*)'Reading TSLB'
  call read_var_wrf(ncid,is3d_st,1,nlev_soil_ns,1,AUX(:,:,1:nlev_soil_ns),'so')
  tslb=AUX(1:nlon_ns,1:nlat_ns,1:nlev_soil_ns)

  !!! PS
  write(*,*)'Reading PS'
  call read_var_wrf(ncid,iv2d_ps,1,1,1,AUX(:,:,1),'2d')
  ps=AUX(1:nlon_ns,1:nlat_ns,1)

  !!! T2
  write(*,*)'Reading T2'
  call read_var_wrf(ncid,iv2d_t2,1,1,1,AUX(:,:,1),'2d')
  t2=AUX(1:nlon_ns,1:nlat_ns,1)

  !!! Q2
  write(*,*)'Reading Q2'
  call read_var_wrf(ncid,iv2d_q2,1,1,1,AUX(:,:,1),'2d')
  q2=AUX(1:nlon_ns,1:nlat_ns,1)

  !!! U10
  write(*,*)'Reading U10'
  call read_var_wrf(ncid,iv2d_u10,1,1,1,AUX(:,:,1),'2d')
  u10=(AUX(1:nlon_ns,1:nlat_ns,1)+AUX(2:nlon_ns+1,1:nlat_ns,1))*0.5

  !!! V10
  write(*,*)'Reading V10'
  call read_var_wrf(ncid,iv2d_v10,1,1,1,AUX(:,:,1),'2d')
  v10=(AUX(1:nlon_ns,1:nlat_ns,1)+AUX(1:nlon_ns,2:nlat_ns+1,1))*0.5

  !!! LANDSEA
  write(*,*)'Reading LANDSEA'
  call read_var_wrf(ncid,iv2d_ls,1,1,1,AUX(:,:,1),'2d')
  landsea=AUX(1:nlon_ns,1:nlat_ns,1)

  !!! SKINTEMP
  write(*,*)'Reading SKINTEMP'
  call read_var_wrf(ncid,iv2d_st,1,1,1,AUX(:,:,1),'2d')
  skintemp=AUX(1:nlon_ns,1:nlat_ns,1)

  !!! SNOW
  write(*,*)'Reading SNOW'
  call read_var_wrf(ncid,iv2d_sn,1,1,1,AUX(:,:,1),'2d')
  snow=AUX(1:nlon_ns,1:nlat_ns,1)

  !!! SNOWH
  write(*,*)'Reading SNOWH'
  call read_var_wrf(ncid,iv2d_sh,1,1,1,AUX(:,:,1),'2d')
  snowh=AUX(1:nlon_ns,1:nlat_ns,1)

  !!! SEAICE
  write(*,*)'Reading SEAICE'
  call read_var_wrf(ncid,iv2d_si,1,1,1,AUX(:,:,1),'2d')
  seaice=AUX(1:nlon_ns,1:nlat_ns,1)

  !!! CANWAT
  write(*,*)'Reading CANWAT'
  call read_var_wrf(ncid,iv2d_cw,1,1,1,AUX(:,:,1),'2d')
  canwat=AUX(1:nlon_ns,1:nlat_ns,1)

  !!! TOPOGRAPHY
  write(*,*)'Reading TOPOGRAPHY'
  topo=phi0(1:nlon_ns,1:nlat_ns)

  CALL check_io(NF_CLOSE(ncid))
  
  !!! COMPUTE SEA LEVEL PRESSURE
 
  CALL GET_SLP(nlon_ns,nlat_ns,nlev_ns,ght,p,t,qv,topo,slp)
 
  !WRITE IN WPS FORMAT
  open(unit=grdfid, file=trim(grdout), status='unknown', form='unformatted')
  !OPEN(grdfid,FILE=grdout,FORM='unformatted',ACCESS='sequential')


  !UU
  output_data % field ='UU'
  output_data % units ='m s-1'
  output_data % desc ='U'
  output_data % slab =  u
  output_data % nlev = nlev_ns
  CALL WRITE_SLAB(grdfid)

  !VV
  output_data % field ='VV'
  output_data % units ='m s-1'
  output_data % desc ='V'
  output_data % slab =  v 
  output_data % nlev = nlev_ns
  CALL WRITE_SLAB(grdfid)

  !WW
  output_data % field ='WW'
  output_data % units ='m s-1'
  output_data % desc ='W'
  output_data % slab =  w
  output_data % nlev = nlev_ns
  CALL WRITE_SLAB(grdfid)

  !TT
  output_data % field ='TT'
  output_data % units ='K'
  output_data % desc ='Temperature'
  output_data % slab =  t
  output_data % nlev = nlev_ns
  CALL WRITE_SLAB(grdfid)

  !P
  output_data % field ='PRESSURE'
  output_data % units ='Pa'
  output_data % desc ='none'
  output_data % slab =  p
  output_data % nlev = nlev_ns
  CALL WRITE_SLAB(grdfid)

  !HGT
  output_data % field ='HGT'
  output_data % units ='m'
  output_data % desc ='none'
  output_data % slab =  ght / gg
  output_data % nlev = nlev_ns
  CALL WRITE_SLAB(grdfid)

  !QV
  output_data % field ='SPECHUMD'
  output_data % units ='Kg Kg-1'
  output_data % desc ='Specific Humidity'
  output_data % slab =  qv
  output_data % nlev = nlev_ns
  CALL WRITE_SLAB(grdfid)
  
  !QC
  output_data % field ='QC'
  output_data % units ='Kg Kg-1'
  output_data % desc ='none'
  output_data % slab =  qc
  output_data % nlev = nlev_ns
  CALL WRITE_SLAB(grdfid)

  !QR
  output_data % field ='QR'
  output_data % units ='Kg Kg-1'
  output_data % desc ='none'
  output_data % slab =  qr
  output_data % nlev = nlev_ns
  CALL WRITE_SLAB(grdfid)

  !QI
  output_data % field ='QI'
  output_data % units ='Kg Kg-1'
  output_data % desc ='none'
  output_data % slab =  qci
  output_data % nlev = nlev_ns
  CALL WRITE_SLAB(grdfid)

  !QS
  output_data % field ='QS'
  output_data % units ='Kg Kg-1'
  output_data % desc ='none'
  output_data % slab =  qs
  output_data % nlev = nlev_ns
  CALL WRITE_SLAB(grdfid)

  !QG
  output_data % field ='QG'
  output_data % units ='Kg Kg-1'
  output_data % desc ='none'
  output_data % slab =  qg
  output_data % nlev = nlev_ns
  CALL WRITE_SLAB(grdfid)

  !PS
  output_data % field ='PSFC'
  output_data % units ='Pa'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  ps
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !Q2
  output_data % field ='SPECHUMD'
  output_data % units ='Kg Kg-1'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  q2
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !T2
  output_data % field ='TT'
  output_data % units ='K'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  t2
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !SOILHGT
  output_data % field ='SOILHGT'
  output_data % units ='m'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  topo
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !U10
  output_data % field ='UU'
  output_data % units ='m s-1'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  u10(:,:)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !V10
  output_data % field ='VV'
  output_data % units ='m s-1'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  v10(:,:)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !WW
  output_data % field ='WW'
  output_data % units ='m s-1'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  w(:,:,1)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)


  !LANDSEA
  output_data % field ='LANDSEA'
  output_data % units ='fraction'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  landsea(:,:)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !SKINTEMP
  output_data % field ='SKINTEMP'
  output_data % units ='K'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  skintemp(:,:)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !SNOW
  output_data % field ='SNOW'
  output_data % units ='Kg m-2'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  snow(:,:)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !SNOWH
  output_data % field ='SNOWH'
  output_data % units ='m'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  snowh(:,:)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !SEAICE
  output_data % field ='SEAICE'
  output_data % units ='proprtn'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  seaice(:,:)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !SLP
  output_data % field ='PMSL'
  output_data % units ='Pa'
  output_data % desc ='Sea Level Pressure'
  output_data % slab(:,:,1) =  slp(:,:)
  output_data % nlev =  1
  output_data % xlvl(1) = 201300
  CALL WRITE_SLAB(grdfid)

  !CANWAT
  output_data % field ='CANWAT'
  output_data % units ='kg m-2'
  output_data % desc ='Plant Canopy Surface Water'
  output_data % slab(:,:,1) =  canwat(:,:)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !Los niveles de suelo son diferentes en la salida del UPP a 
  !los niveles nativos del WRF. Ensayo una interpolacion muy rudimentaria.
  !quiza esto se pueda mejorar.
   
  !SM000010
  output_data % field ='SM000010'
  output_data % units ='m3 m-3'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  smois(:,:,1)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !SOILM005
  output_data % field ='SM010040'
  output_data % units ='m3 m-3'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  smois(:,:,2)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !SOILM020
  output_data % field ='SM040100'
  output_data % units ='m3 m-3'
  output_data % desc ='none'
  output_data % slab(:,:,1) = smois(:,:,3) 
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !SOILM040
  output_data % field ='SM100200'
  output_data % units ='kg m-3'
  output_data % desc ='none'
  output_data % slab(:,:,1) = smois(:,:,4) 
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !ST000010
  output_data % field ='ST000010'
  output_data % units ='K'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  tslb(:,:,1)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !ST010040
  output_data % field ='ST010040'
  output_data % units ='K'
  output_data % desc ='none'
  output_data % slab(:,:,1) =  tslb(:,:,2)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !ST040100
  output_data % field ='ST040100'
  output_data % units ='K'
  output_data % desc ='none'
  output_data % slab(:,:,1) = tslb(:,:,3) 
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  !ST100200
  output_data % field ='ST100200'
  output_data % units ='K'
  output_data % desc ='none'
  output_data % slab(:,:,1) = tslb(:,:,4)
  output_data % nlev =  1
  output_data % xlvl(1) = 200100
  CALL WRITE_SLAB(grdfid)

  CLOSE(grdfid)

  STOP
END PROGRAM wrf_to_wps



