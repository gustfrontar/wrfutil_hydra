#!/bin/sh
set -ex
PGM=./test.exe
F90=ifort  #mpif90

LIB_NETCDF="-L//data2/opt/netcdf/4.3.2/lib/ -lnetcdff"
INC_NETCDF="-I//data2/opt/netcdf/4.3.2/include/"


OMP=
F90OPT='-convert big_endian -O3 -openmp'

cd ./src

#PRE CLEAN
rm -f *.mod
rm -f *.o


#COMPILIN
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_wrf_memnc.f90
$F90 $OMP $F90OPT -c common_radar_tools.f90
$F90 $OMP $F90OPT -c common_wrf_to_radar.f90
$F90 $OMP $F90OPT -c test.f90
$F90 $OMP $F90OPT -o ${PGM} *.o  ${LIB_NETCDF}

mv *.exe ../

#CLEAN UP
rm -f *.mod
rm -f *.o


echo "NORMAL END"
