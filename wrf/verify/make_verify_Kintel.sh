#!/bin/sh
set -ex
PGM=verify.exe
F90=mpif90
OMP=

F90OPT='-convert big_endian -O3' # -g -traceback' # -Hs' -Kfast,parallel

BLAS=1 #0: no blas 1: using blas
BASEDIR=${HOME}/libintel/netcdf-3.6.3/

rm -f *.mod
rm -f *.o

COMMONDIR=../../common/
#LIB_NETCDF="-L${BASEDIR}/netcdf/4.1.1/lib -lnetcdf"
#INC_NETCDF="-I${BASEDIR}/netcdf/4.1.1/include/"


COMMONDIR=../../common/
ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/common.f90 ./
#ln -fs $COMMONDIR/common_smooth2d.f90 ./

ln -fs ../common/module_map_utils.f90 ./


$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c module_map_utils.f90
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
