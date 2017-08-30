#!/bin/sh
set -ex
PGM=../breeding
F90=ifort
OMP=

F90OPT='-O3' #-Kfast,parallel' # -Hs'

cd ./src

rm -f *.mod
rm -f *.o

LIB_NETCDF="-L/data/opt/netcdf-fortran/4.4.1_intel/lib/ -lnetcdff"
INC_NETCDF="-I/data/opt/netcdf-fortran/4.4.1_intel/include/"


$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c $INC_NETCDF common_ncio.f90
$F90 $OMP $F90OPT -c scale_const.f90
$F90 $OMP $F90OPT -c common_smooth2d.f90
$F90 $OMP $F90OPT -c common_namelist.f90
$F90 $OMP $F90OPT -c $INC_NETCDF common_scale.f90
$F90 $OMP $F90OPT -c $INC_NETCDF common_breeding.f90
$F90 $OMP $F90OPT -c $INC_NETCDF main_breeding.f90
$F90 $OMP $F90OPT -o ${PGM} *.o $LIB_NETCDF 

rm -f *.mod
rm -f *.o

echo "NORMAL END"
