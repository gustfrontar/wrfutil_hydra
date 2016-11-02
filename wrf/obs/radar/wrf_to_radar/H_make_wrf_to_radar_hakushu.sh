#!/bin/sh
set -ex
PGM=./wrf_to_radar.exe
F90=ifort  

LIB_NETCDF="-L/opt/netcdf-fortran/4.2/lib/ -lnetcdff"
INC_NETCDF="-I/opt/netcdf-fortran/4.2/include/"


OMP=
F90OPT=' -convert big_endian -O3 -openmp'

cd ./src

#PRE CLEAN
rm -f *.mod
rm -f *.o

ln -sf ../../common/common_radar_tools.f90     .
ln -sf ../../../../common/common_wrf.f90       .
ln -sf ../../../../common/common_obs_wrf.f90   .
ln -sf ../../../../common/module_map_utils.f90 .
#ln -sf ../../../../common/common_namelist.f90  .
ln -sf ../../../../../common/common.f90        .
ln -sf ../../../../../common/SFMT.f90          .

#COMPILIN
$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c module_map_utils.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_wrf.f90
$F90 $OMP $F90OPT -c common_namelist_wrf_to_radar.f90
$F90 $OMP $F90OPT -c common_obs_wrf.f90
$F90 $OMP $F90OPT -c common_radar_tools.f90
$F90 $OMP $F90OPT -c common_wrf_to_radar.f90
$F90 $OMP $F90OPT -c main_wrf_to_radar.f90
$F90 $OMP $F90OPT -o ${PGM} *.o  ${LIB_NETCDF}

mv *.exe ../

#CLEAN UP
rm -f *.mod
rm -f *.o
rm common_wrf.f90
#rm common_namelist.f90
rm module_map_utils.f90
rm common_radar_tools.f90
rm common_obs_wrf.f90
rm common.f90
rm SFMT.f90


echo "NORMAL END"
