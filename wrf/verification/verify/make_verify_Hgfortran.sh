#!/bin/sh
set -ex
PGM=verify.exe
F90=gfortran
OMP=

F90OPT='-O3' # -g -traceback' #

BLAS=1 #0: no blas 1: using blas

rm -f *.mod
rm -f *.o

cd ./src/

COMMONDIR=../../../../common/
ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/common.f90 ./
#ln -fs $COMMONDIR/common_smooth2d.f90 ./

ln -fs ../../../common/module_map_utils.f90 ./
ln -fs ../../common/common_verification.f90 ./


$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c module_map_utils.f90
$F90 $OMP $F90OPT -c common_verification.f90
#$F90 $OMP $F90OPT -c common_smooth2d.f90
$F90 $OMP $F90OPT -c verify_tools.f90
$F90 $OMP $F90OPT -c verify.f90
$F90 $OMP $F90OPT -o ${PGM} *.o  $LIB_NETCDF 

rm -f *.mod
rm -f *.o
rm -f netlib2.f

rm SFMT.f90 
rm common.f90 
rm module_map_utils.f90 
#rm common_smooth2d.f90

echo "NORMAL END"
