#!/bin/sh
set -ex
PGM=./wrf_to_radar.exe
F90=ifort  

LIB_NETCDF="-L/home/jruiz/netcdf4.intel/lib/ -lnetcdff"
INC_NETCDF="-I/home/jruiz/netcdf4.intel/include/"


OMP=
F90OPT='-g -traceback ' # -convert big_endian -O3 -openmp'

cd ./src

#PRE CLEAN
rm -f *.mod
rm -f *.o

ln -sf ../../common/common_radar_tools_cfradial.f90     .
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
$F90 $OMP $F90OPT $INC_NETCDF -c common_namelist_wrf_to_radar.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_obs_wrf.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_radar_tools_cfradial.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_wrf_to_radar.f90
$F90 $OMP $F90OPT $INC_NETCDF -c main_wrf_to_radar.f90
$F90 $OMP $F90OPT $INC_NETCDF -o ${PGM} *.o  ${LIB_NETCDF}

mv *.exe ../

#CLEAN UP
rm -f *.mod
rm -f *.o
#rm common_wrf.f90
#rm common_namelist.f90
rm module_map_utils.f90
rm common_radar_tools_cfradial.f90
rm common_obs_wrf.f90
rm common.f90
rm SFMT.f90


echo "NORMAL END"
