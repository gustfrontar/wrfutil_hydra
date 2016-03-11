module met_data_module
   use common

   ! Derived types
   type met_data
      integer                       :: version, nx, ny, iproj, nlev
      real(r_sngl)                  :: xfcst, startlat, startlon, starti, startj, &
                                       deltalat, deltalon, dx, dy, xlonc, stdlon, &
                                       truelat1 , truelat2 , earth_radius 
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

  DO k=1,output_data % nlev
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



end module met_data_module

