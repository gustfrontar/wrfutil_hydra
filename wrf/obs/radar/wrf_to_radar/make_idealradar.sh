#!/bin/sh
set -ex
PGM=./model_to_radar_ideal.exe
F90=ifort  #mpif90

LIB_NETCDF="-L//data1/opt/netcdf/4.3.2/lib/ -lnetcdff"
INC_NETCDF="-I//data1/opt/netcdf/4.3.2/include/"


OMP=
F90OPT='-convert big_endian -O3 '

cd ./src

#PRE CLEAN
rm -f *.mod
rm -f *.o


#COMPILIN
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_wrf_memnc.f90
$F90 $OMP $F90OPT -c main_idealradar.f90
$F90 $OMP $F90OPT -o ${PGM} *.o  ${LIB_NETCDF}

mv $PGM ../

#CLEAN UP
rm -f *.mod
rm -f *.o


echo "NORMAL END"
