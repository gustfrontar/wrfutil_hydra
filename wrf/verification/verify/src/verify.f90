!===============================================================================
! verify: main program
! This program compares two grid datasets.
! Grid information is read from a ctl file 
! An alternative (usually lower resolution) regular grid can also be defined
! to perform verification at larger scales (eg. to compare against low
! resolution analysis)
!-------------------------------------------------------------------------------
!===============================================================================

program verify
!-------------------------------------------------------------------------------

  use common
  use verify_tools
  implicit none

  !Original grid files and variables
  character(20) :: fcstfile = 'fcst00000.grd'
  character(20) :: analfile = 'anal.grd'
  character(20) :: meanfile = 'mean.grd' !Ensemble mean
  character(20) :: sprdfile = 'sprd.grd' !Ensemble spread
  character(20) :: merrfile = 'merr.grd' !Ensemble mean error
  character(20) :: ctl_file_for = 'inputfor.ctl'
  character(20) :: ctl_file_anl = 'inputanl.ctl'

  character(20) :: areafile = 'area.txt'

  real(r_size) ,ALLOCATABLE :: vf(:,:,:) , vfm(:,:,:) , vfs(:,:,:) , ve(:,:,:)
  real(r_size) ,ALLOCATABLE :: va(:,:,:) 
  real(r_size) ,ALLOCATABLE :: tmpgrid(:,:,:)

  !Regrid grid files and variables
  character(20) :: meanfilereg = 'meanreg.grd' !Ensemble mean
  character(20) :: sprdfilereg = 'sprdreg.grd' !Ensemble spread
  character(20) :: merrfilereg = 'merrreg.grd' !Ensemble mean error

  character(20) :: areafilereg = 'areareg.txt'

  real(r_size)    ,ALLOCATABLE :: vfreg(:,:,:) , vfmreg(:,:,:) , vfsreg(:,:,:) , vereg(:,:,:) , vareg(:,:,:)
  integer         ,ALLOCATABLE :: numreg(:,:,:)

  integer,ALLOCATABLE  ::     num_ana(:,:)  !narea,nv 
  real(r_size),ALLOCATABLE :: wei_ana(:,:)  !narea,nv
  real(r_size),ALLOCATABLE :: bias_ana(:,:) !narea,nv
  real(r_size),ALLOCATABLE :: abse_ana(:,:) !narea,nv
  real(r_size),ALLOCATABLE :: rmse_ana(:,:) !narea,nv
  real(r_size),ALLOCATABLE :: sprd_ana(:,:) !narea,nv

  integer,ALLOCATABLE  ::     num_anareg(:,:)  !narea,nv
  real(r_size),ALLOCATABLE :: wei_anareg(:,:)  !narea,nv
  real(r_size),ALLOCATABLE :: bias_anareg(:,:) !narea,nv
  real(r_size),ALLOCATABLE :: abse_anareg(:,:) !narea,nv
  real(r_size),ALLOCATABLE :: rmse_anareg(:,:) !narea,nv
  real(r_size),ALLOCATABLE :: sprd_anareg(:,:) !narea,nv

  real(r_size), ALLOCATABLE :: rmse(:,:,:),bias(:,:,:),abse(:,:,:),coun(:,:,:)
  real(r_size), ALLOCATABLE :: rmsereg(:,:,:),biasreg(:,:,:),absereg(:,:,:),counreg(:,:,:)

  LOGICAL , ALLOCATABLE :: undefmaskf(:,:,:) , totalundefmaskf(:,:,:) , tmpmask(:,:,:)
  LOGICAL , ALLOCATABLE :: undefmaska(:,:,:) , totalundefmaska(:,:,:) 
  LOGICAL , ALLOCATABLE :: undefmaskreg(:,:,:) , totalundefmaskreg(:,:,:)

  real(r_size) :: dz, hx, dep, wei, latm1, latm2
  integer :: n, nobs,  k, vk, vid, iarea
  integer :: i, i1, i2, j, j1, j2 , ii , jj , kk 

  type(ctl_info) :: ctlanl , ctlfor
  integer , allocatable :: mapfor(:) , mapanl(:)
  
!  integer :: nobsused

  integer :: iunit, iolen, irec

!-------------------------------------------------------------------------------

!===============================================================================

  call read_namelist !ensemble size and general options.

  call parse_ctl(ctl_file_anl,ctlanl) !Get ctl dimmensions and variables for the verification dataset.
  call parse_ctl(ctl_file_for,ctlfor) !Get ctl dimmensions and variables for the forecast dataset.

  !nlon=ctlanl%nlon
  !nlat=ctlanl%nlat

  !Check if the grids are equivalent (if they are not we can only compare these
  !two grids in the common regrid grid)
  if( ctlanl%nlon /= ctlfor%nlon .or. ctlanl%nlat /= ctlfor%nlat )then
    samegrid_output=.false. !We cannot compare this two grids directly.
  else 
   do ii=1,ctlanl%nlon
     do jj=1,ctlanl%nlat
        if( ctlanl%lon(ii,jj) /= ctlfor%lon(ii,jj) )samegrid_output=.false.
        if( ctlanl%lat(ii,jj) /= ctlfor%lat(ii,jj) )samegrid_output=.false.
     enddo
   enddo
  endif

  if( .not. samegrid_output )then
    WRITE(*,*)"Grids are not equivalent we can only compare regrided grids"
  else
    WRITE(*,*)"Grids are equivalent we will produce errors at the original model grid"
  endif

  ALLOCATE( lon(ctlanl%nlon,ctlanl%nlat) , lat(ctlanl%nlon,ctlanl%nlat) )
 
  lon=ctlanl%lon
  lat=ctlanl%lat

  ALLOCATE(mapanl(ctlanl%nfields),mapfor(ctlfor%nfields))
  call match_ctl_variables(ctlanl,ctlfor,mapanl,mapfor)  !See which are the common variables and levels in both ctls.
                                           !verification will be performed only on the common variables and levels.

  !Allocate arrays to store the verification for the common fields.
  ALLOCATE( num_ana(narea,commonnfields) )
  ALLOCATE( wei_ana(narea,commonnfields) )
  ALLOCATE( bias_ana(narea,commonnfields) )
  ALLOCATE( abse_ana(narea,commonnfields) )
  ALLOCATE( rmse_ana(narea,commonnfields) )
  ALLOCATE( sprd_ana(narea,commonnfields) )
  ALLOCATE( vf(ctlfor%nlon,ctlfor%nlat,commonnfields) , va(ctlanl%nlon,ctlanl%nlat,commonnfields) )
  ALLOCATE( vfm(ctlfor%nlon,ctlfor%nlat,commonnfields) , vfs(ctlfor%nlon,ctlfor%nlat,commonnfields) )

  if( samegrid_output )then
     ALLOCATE( ve(ctlfor%nlon,ctlfor%nlat,commonnfields) )
  endif
  ALLOCATE( undefmaskf(ctlfor%nlon,ctlfor%nlat,commonnfields) )
  ALLOCATE( totalundefmaskf(ctlfor%nlon,ctlfor%nlat,commonnfields) )
  ALLOCATE( undefmaska(ctlanl%nlon,ctlanl%nlat,commonnfields) )
  ALLOCATE( totalundefmaska(ctlanl%nlon,ctlanl%nlat,commonnfields) )


  if( regrid_output )then !We will generate also a regrid output in a regular grid.
    call get_regrid_grid(ctlanl%minlon,ctlanl%maxlon,ctlanl%minlat,ctlanl%maxlat,ctlanl%minlev, &
                         ctlanl%maxlev,commonvarnameall,commonlevall,commonnfields)
    !call get_regrid_grid(ctlanl)   !Get the dimensions of the regrid grid.

    !Allocate arrays for the regrided fields.
    ALLOCATE( vfreg(nlonreg,nlatreg,nfieldsreg) )
    ALLOCATE( vareg(nlonreg,nlatreg,nfieldsreg) )
    ALLOCATE( vfmreg(nlonreg,nlatreg,nfieldsreg) )
    ALLOCATE( vfsreg(nlonreg,nlatreg,nfieldsreg) )
    ALLOCATE( vereg(nlonreg,nlatreg,nfieldsreg) )
    ALLOCATE( numreg(nlonreg,nlatreg,nfieldsreg) )
    ALLOCATE( undefmaskreg(nlonreg,nlatreg,nfieldsreg) )
    ALLOCATE( totalundefmaskreg(nlonreg,nlatreg,nfieldsreg) )

    ALLOCATE( num_anareg(narea,nfieldsreg) )
    ALLOCATE( wei_anareg(narea,nfieldsreg) )
    ALLOCATE( bias_anareg(narea,nfieldsreg) )
    ALLOCATE( abse_anareg(narea,nfieldsreg) )
    ALLOCATE( rmse_anareg(narea,nfieldsreg) )
    ALLOCATE( sprd_anareg(narea,nfieldsreg) )
  endif

!-------------------------------------------------------------------------------

! Comput the ensemble mean 

totalundefmaskf=.true.
totalundefmaska=.true.

vfm=0.0d0
vfs=0.0d0

DO i=1,nbv

  WRITE( fcstfile(5:9), '(I5.5)' )i

  ALLOCATE(tmpgrid(ctlfor%nlon,ctlfor%nlat,ctlfor%nfields),tmpmask(ctlfor%nlon,ctlfor%nlat,ctlfor%nfields))
 
  call read_grd(fcstfile,ctlfor%nlon,ctlfor%nlat,ctlfor%nfields,tmpgrid,tmpmask,ctlfor%undefbin)

    DO ii=1,ctlfor%nfields
       if( ctlfor%readvar(ii) )then
        vf(:,:,mapfor(ii))=tmpgrid(:,:,ii)
        undefmaskf(:,:,mapfor(ii))=tmpmask(:,:,ii)
       endif
    ENDDO
  DEALLOCATE( tmpgrid , tmpmask )

  WHERE( .NOT. undefmaskf )
    totalundefmaskf=.false.
  END WHERE
    

  vfm=vfm+vf
  vfs=vfs+vf**2

ENDDO


  ALLOCATE( tmpgrid(ctlanl%nlon,ctlanl%nlat,ctlanl%nfields),tmpmask(ctlanl%nlon,ctlanl%nlat,ctlanl%nfields) )

  call read_grd(analfile,ctlanl%nlon,ctlanl%nlat,ctlanl%nfields,tmpgrid,tmpmask,ctlanl%undefbin)

      DO ii=1,ctlanl%nfields
       if( ctlanl%readvar(ii) )then
        va(:,:,mapanl(ii))=tmpgrid(:,:,ii)
        undefmaska(:,:,mapanl(ii))=tmpmask(:,:,ii)
       endif
      ENDDO
  DEALLOCATE( tmpgrid , tmpmask )
  

  WHERE( .NOT. undefmaska )
   totalundefmaska=.false.
  END WHERE

  !At this point total undef mask is false
  !if any of the ensemble members is undef 
  !or if the analysis is undef.

  !Finish computing the analysis ensemble mean and spread.

  vfm=vfm/REAL(nbv,r_size)
  vfs=vfs/REAL(nbv,r_size)
  vfs=vfs-vfm**2
  WHERE(vfs < 0)
    vfs=0.0d0
  ENDWHERE
  vfs=SQRT(vfs)

  !Writing the ensemble mean and spread in the postproc grid.

  CALL write_grd(meanfile,ctlfor%nlon,ctlfor%nlat,commonnfields,vfm,totalundefmaskf,ctlanl%undefbin)
  CALL write_grd(sprdfile,ctlfor%nlon,ctlfor%nlat,commonnfields,vfs,totalundefmaskf,ctlanl%undefbin)

  if ( samegrid_output )then
  !If analysis and forecast share the same original grid, then coordinate the
  !masks so the error is computed where both analysis and forecast are
  !available.
    WHERE( .NOT. totalundefmaskf )
     totalundefmaska=.false.
    ENDWHERE
    WHERE( .NOT. totalundefmaska )
     totalundefmaskf=.false.
    ENDWHERE
  !Compute the error
    ve=vfm-va
  !Write it out
    CALL write_grd(merrfile,ctlfor%nlon,ctlfor%nlat,commonnfields,ve ,totalundefmaskf,ctlanl%undefbin)
  endif


  !Computeing and writing the regrid output
 if ( regrid_output ) then
  totalundefmaskreg=.true.
  CALL regrid_var(ctlfor%nlon,ctlfor%nlat,commonnfields,nlonreg,nlatreg,nfieldsreg,levregmap,ctlfor%lon,ctlfor%lat, &
                  lonreg,latreg,vfm,totalundefmaskf,vfmreg,totalundefmaskreg,numreg)
  CALL regrid_var(ctlfor%nlon,ctlfor%nlat,commonnfields,nlonreg,nlatreg,nfieldsreg,levregmap,ctlfor%lon,ctlfor%lat, &
                  lonreg,latreg,vfs,totalundefmaskf,vfsreg,totalundefmaskreg,numreg)
  CALL regrid_var(ctlanl%nlon,ctlanl%nlat,commonnfields,nlonreg,nlatreg,nfieldsreg,levregmap,ctlanl%lon,ctlanl%lat, &
                  lonreg,latreg,va ,totalundefmaska,vareg ,totalundefmaskreg,numreg)
  !Compute error in the regrid grid
  vereg=vfmreg-vareg  

  CALL write_grd(meanfilereg,nlonreg,nlatreg,nfieldsreg,vfmreg,totalundefmaskreg,ctlanl%undefbin)
  CALL write_grd(sprdfilereg,nlonreg,nlatreg,nfieldsreg,vfsreg,totalundefmaskreg,ctlanl%undefbin)
  CALL write_grd(merrfilereg,nlonreg,nlatreg,nfieldsreg,vereg ,totalundefmaskreg,ctlanl%undefbin)

 endif

  !Compute rmse over selected regions and generate text ouput. 

 if( samegrid_output )then
    num_ana = 0
    wei_ana = 0.0d0
    bias_ana = 0.0d0
    abse_ana = 0.0d0
    rmse_ana = 0.0d0
    sprd_ana = 0.0d0

    do iarea = 1, narea
     do i=1,ctlfor%nlon
      do j=1,ctlfor%nlat
       if( lon(i,j) < vlon1(iarea) )cycle
       if( lon(i,j) > vlon2(iarea) )cycle
       if( lat(i,j) < vlat1(iarea) )cycle
       if( lat(i,j) > vlat2(iarea) )cycle

        do vid = 1, commonnfields 
         if ( totalundefmaska(i,j,vid) ) then
                dep = vfm(i,j,vid) - va(i,j,vid)
                num_ana(iarea,vid) = num_ana(iarea,vid) + 1
                bias_ana(iarea,vid) = bias_ana(iarea,vid) + dep
                abse_ana(iarea,vid) = abse_ana(iarea,vid) + abs(dep)
                rmse_ana(iarea,vid) = rmse_ana(iarea,vid) + dep*dep
                sprd_ana(iarea,vid) = sprd_ana(iarea,vid) + vfs(i,j,vid)**2
          end if
         end do
      end do
     end do
    end do ! [ iarea = 1, narea ]

    do vid = 1, commonnfields
        do iarea = 1, narea
          if (num_ana(iarea,vid) == 0) then
            bias_ana(iarea,vid)=ctlfor%undefbin
            rmse_ana(iarea,vid)=ctlfor%undefbin
            abse_ana(iarea,vid)=ctlfor%undefbin
          else
            bias_ana(iarea,vid)=bias_ana(iarea,vid)/real(num_ana(iarea,vid),r_size)
            rmse_ana(iarea,vid)=SQRT( rmse_ana(iarea,vid)/real(num_ana(iarea,vid),r_size) )
            abse_ana(iarea,vid)=abse_ana(iarea,vid)/real(num_ana(iarea,vid),r_size)
            sprd_ana(iarea,vid)=SQRT( sprd_ana(iarea,vid)/real(num_ana(iarea,vid),r_size) )
          end if
        end do
    end do

    iunit=101
    open(iunit,file=areafile,form='formatted',recl=1000)
      WRITE(iunit,*)'AREAL ERROR STATISTICS COMPUTED FROM THE ORIGINAL GRID'
      WRITE(iunit,*)'UNDEF ',ctlfor%undefbin,' NVARS ',commonnfields,' NAREA ',narea
      WRITE(iunit,*)'VLON1',vlon1(1:narea)
      WRITE(iunit,*)'VLON2',vlon2(1:narea)
      WRITE(iunit,*)'VLAT1',vlat1(1:narea)
      WRITE(iunit,*)'VLAT2',vlat2(1:narea)
      WRITE(iunit,*)'RMSE SECTION'
      WRITE(iunit,*)'VAR','  LEV (Pa) ',(' AREA',i,i=1,narea)
      DO i=1,commonnfields
        WRITE(iunit,*)commonvarnameall(i),commonlevall(i)*100.0d0,(REAL(rmse_ana(iarea,i),r_sngl),iarea=1,narea)
      ENDDO
      WRITE(iunit,*)'SPRD SECTION'
      WRITE(iunit,*)'VAR','  LEV (Pa) ',(' AREA',i,i=1,narea)
      DO i=1,commonnfields
        WRITE(iunit,*)commonvarnameall(i),commonlevall(i)*100.0d0,(REAL(sprd_ana(iarea,i),r_sngl),iarea=1,narea)
      ENDDO
      WRITE(iunit,*)'BIAS SECTION'
      WRITE(iunit,*)'VAR','  LEV (Pa) ',(' AREA',i,i=1,narea)
      DO i=1,commonnfields
        WRITE(iunit,*)commonvarnameall(i),commonlevall(i)*100.0d0,(REAL(bias_ana(iarea,i),r_sngl),iarea=1,narea)
      ENDDO
      WRITE(iunit,*)'ABSERROR SECTION'
      WRITE(iunit,*)'VAR','  LEV (Pa) ',(' AREA',i,i=1,narea)
      DO i=1,commonnfields
        WRITE(iunit,*)commonvarnameall(i),commonlevall(i)*100.0d0,(REAL(abse_ana(iarea,i),r_sngl),iarea=1,narea)
      ENDDO
      WRITE(iunit,*)'COUNT SECTION'
      WRITE(iunit,*)'VAR','  LEV (Pa) ',(' AREA',i,i=1,narea)
      DO i=1,commonnfields
        WRITE(iunit,*)commonvarnameall(i),commonlevall(i)*100.0d0,(REAL(num_ana(iarea,i),r_sngl),iarea=1,narea)
      ENDDO
    close(iunit)

  endif

  if( regrid_output )then

    num_anareg = 0
    wei_anareg = 0.0d0
    bias_anareg = 0.0d0
    abse_anareg = 0.0d0
    rmse_anareg = 0.0d0
    sprd_anareg = 0.0d0

    do iarea = 1, narea
     do i=1,nlonreg
      do j=1,nlatreg
       if( lonreg(i,j) < vlon1(iarea) )cycle
       if( lonreg(i,j) > vlon2(iarea) )cycle
       if( latreg(i,j) < vlat1(iarea) )cycle
       if( latreg(i,j) > vlat2(iarea) )cycle

        do vid = 1, nfieldsreg
         if ( totalundefmaskreg(i,j,vid) ) then
                dep = vereg(i,j,vid) 
                num_anareg(iarea,vid) = num_anareg(iarea,vid) + 1
                bias_anareg(iarea,vid) = bias_anareg(iarea,vid) + dep
                abse_anareg(iarea,vid) = abse_anareg(iarea,vid) + abs(dep)
                rmse_anareg(iarea,vid) = rmse_anareg(iarea,vid) + dep*dep
                sprd_anareg(iarea,vid) = sprd_anareg(iarea,vid) + vfsreg(i,j,vid) ** 2
          end if
         end do
      end do
     end do
    end do ! [ iarea = 1, narea ]



    do vid = 1, nfieldsreg
        do iarea = 1, narea
          if (num_anareg(iarea,vid) == 0) then
            bias_anareg(iarea,vid)=ctlfor%undefbin
            rmse_anareg(iarea,vid)=ctlfor%undefbin
            abse_anareg(iarea,vid)=ctlfor%undefbin
            sprd_anareg(iarea,vid)=ctlfor%undefbin
          else
            bias_anareg(iarea,vid)=bias_anareg(iarea,vid)/real(num_anareg(iarea,vid),r_size)
            rmse_anareg(iarea,vid)=SQRT( rmse_anareg(iarea,vid)/real(num_anareg(iarea,vid),r_size) )
            abse_anareg(iarea,vid)=abse_anareg(iarea,vid)/real(num_anareg(iarea,vid),r_size)
            sprd_anareg(iarea,vid)=SQRT( sprd_anareg(iarea,vid)/real(num_anareg(iarea,vid),r_size) )
          end if
        end do
    end do

    iunit=101
    open(iunit,file=areafilereg,form='formatted',recl=1000)
      WRITE(iunit,*)'AREAL ERROR STATISTICS COMPUTED FROM THE REGRID GRID'
      WRITE(iunit,*)'UNDEF ',ctlfor%undefbin,' NVARS ',nfieldsreg,' NAREA ',narea
      WRITE(iunit,*)'VLON1',vlon1(1:narea)
      WRITE(iunit,*)'VLON2',vlon2(1:narea)
      WRITE(iunit,*)'VLAT1',vlat1(1:narea)
      WRITE(iunit,*)'VLAT2',vlat2(1:narea)
      WRITE(iunit,*)'RMSE SECTION'
      WRITE(iunit,*)'VAR','  LEV (Pa) ',(' AREA',i,i=1,narea)
      DO i=1,nfieldsreg
        WRITE(iunit,*)varnameallreg(i),levallreg(i)*100.0d0,(REAL(rmse_anareg(iarea,i),r_sngl),iarea=1,narea)
      ENDDO
      WRITE(iunit,*)'SPRD SECTION'
      WRITE(iunit,*)'VAR','  LEV (Pa)  ',(' AREA',i,i=1,narea)
      DO i=1,nfieldsreg
        WRITE(iunit,*)varnameallreg(i),levallreg(i)*100.0d0,(REAL(sprd_anareg(iarea,i),r_sngl),iarea=1,narea)
      ENDDO
      WRITE(iunit,*)'BIAS SECTION'
      WRITE(iunit,*)'VAR','  LEV (Pa) ',(' AREA',i,i=1,narea)
      DO i=1,nfieldsreg
        WRITE(iunit,*)varnameallreg(i),levallreg(i)*100.0d0,(REAL(bias_anareg(iarea,i),r_sngl),iarea=1,narea)
      ENDDO
      WRITE(iunit,*)'ABSERROR SECTION'
      WRITE(iunit,*)'VAR','  LEV (Pa) ',(' AREA',i,i=1,narea)
      DO i=1,nfieldsreg
        WRITE(iunit,*)varnameallreg(i),levallreg(i)*100.0d0,(REAL(abse_anareg(iarea,i),r_sngl),iarea=1,narea)
      ENDDO
      WRITE(iunit,*)'COUNT SECTION'
      WRITE(iunit,*)'VAR','  LEV (Pa) ',(' AREA',i,i=1,narea)
      DO i=1,nfieldsreg
        WRITE(iunit,*)varnameallreg(i),levallreg(i)*100.0d0,(REAL(num_anareg(iarea,i),r_sngl),iarea=1,narea)
      ENDDO
    close(iunit)

  endif

STOP

!-------------------------------------------------------------------------------
end program
!===============================================================================
