#!/bin/sh
set -ex
PGM=./scale_to_radar.exe
F90=ifort  

LIB_NETCDF="-L/usr/local/netcdf4.intel/lib/ -lnetcdff"
INC_NETCDF="-I/usr/local/netcdf4.intel/include/"


OMP=
F90OPT='-g -traceback'

cd ./src

#PRE CLEAN
rm -f *.mod
rm -f *.o

ln -sf ../../common/common_radar_tools_cfradial.f90     .
ln -sf ../../../../../common/common.f90        .
ln -sf ../../../../../common/SFMT.f90          .

#COMPILIN
$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c common_namelist_scale_to_radar.f90
$F90 $OMP $F90OPT -c scale_const.f90
$F90 $OMP $F90OPT -c scale_mapproj.F90
$F90 $OMP $F90OPT $INC_NETCDF -c common_ncio.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_scale_hist.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_obs_scale.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_radar_tools_cfradial.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_scale_to_radar.f90
$F90 $OMP $F90OPT $INC_NETCDF -c main_scale_to_radar.f90
$F90 $OMP $F90OPT $INC_NETCDF -o ${PGM} *.o  ${LIB_NETCDF}

mv *.exe ../

#CLEAN UP
rm -f *.mod
rm -f *.o
rm common_radar_tools_cfradial.f90
rm common.f90
rm SFMT.f90


echo "NORMAL END"
