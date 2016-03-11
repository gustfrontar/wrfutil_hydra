#!/bin/sh
set -ex
PGM=./model_to_radar_invap_lessvertres.exe
F90=ifort  #mpif90

LIB_NETCDF="-L/usr/local/netcdf.intel/lib/ -lnetcdff"
INC_NETCDF="-I/usr/local/netcdf.intel/include/"


OMP=
F90OPT='-convert big_endian -O3 -openmp'

cd ./src

#PRE CLEAN
rm -f *.mod
rm -f *.o


#COMPILIN
$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c module_map_utils.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_wrf_memnc.f90
$F90 $OMP $F90OPT -c common_radar_tools.f90
$F90 $OMP $F90OPT -c common_wrf_to_radar.f90
$F90 $OMP $F90OPT -c main_invap_lessvertres.f90
$F90 $OMP $F90OPT -o ${PGM} *.o  ${LIB_NETCDF}

mv *.exe ../

#CLEAN UP
rm -f *.mod
rm -f *.o


echo "NORMAL END"
