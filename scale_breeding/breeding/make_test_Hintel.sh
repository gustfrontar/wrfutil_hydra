#!/bin/sh
set -ex
PGM=../test
F90=ifort
OMP=

F90OPT='-g -traceback' #-Kfast,parallel' # -Hs'

cd ./src

rm -f *.mod
rm -f *.o

LIB_NETCDF="-L/data/opt/netcdf-fortran/4.4.1_intel/lib/ -lnetcdff"
INC_NETCDF="-I/data/opt/netcdf-fortran/4.4.1_intel/include/"

$F90 $OMP $F90OPT -c $INC_NETCDF test.f90
$F90 $OMP $F90OPT -o ${PGM} *.o $LIB_NETCDF 

rm -f *.mod
rm -f *.o

echo "NORMAL END"
