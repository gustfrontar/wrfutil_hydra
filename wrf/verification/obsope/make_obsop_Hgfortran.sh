#!/bin/sh
set -ex
PGM=obsope.exe
F90=/home/jruiz/mpich.gfortran/bin/mpif90
OMP=

F90OPT='-O3 -fconvert=big-endian' #-convert big_endian -O3 ' # -Hs' -Kfast,parallel

BLAS=1 #0: no blas 1: using blas
BASEDIR=${HOME}/libintel/netcdf-3.6.3/

rm -f *.mod
rm -f *.o

cd ./src/

COMMONDIR=../../../../common/
COMMONDIRWRF=../../../common/
LIB_NETCDF="-L/home/jruiz/netcdf.gfortran/lib -lnetcdff"
INC_NETCDF="-I/home/jruiz/netcdf.gfortran/include/"


ln -fs $COMMONDIR/SFMT.f90            ./
ln -fs $COMMONDIR/common.f90          ./
ln -fs $COMMONDIR/common_mpi.f90      ./
ln -fs $COMMONDIR/common_smooth2d.f90 ./

ln -fs $COMMONDIRWRF/common_wrf.f90              ./
ln -fs $COMMONDIRWRF/common_mpi_wrf.f90          ./
ln -fs $COMMONDIRWRF/common_obs_wrf.f90          ./
ln -fs $COMMONDIRWRF/module_map_utils.f90        ./


$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c common_mpi.f90
$F90 $OMP $F90OPT -c module_map_utils.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_wrf.f90
$F90 $OMP $F90OPT -c common_namelist.f90
$F90 $OMP $F90OPT -c common_smooth2d.f90
$F90 $OMP $F90OPT -c common_obs_wrf.f90
$F90 $OMP $F90OPT -c common_mpi_wrf.f90
$F90 $OMP $F90OPT -c obsope_tools.f90
$F90 $OMP $F90OPT -c obsope.f90
$F90 $OMP $F90OPT -o ${PGM} *.o  $LIB_NETCDF 

mv *.exe ../

rm -f *.mod
rm -f *.o
rm -f netlib2.f

rm SFMT.f90 
rm common.f90 
rm common_mpi.f90 
rm common_wrf.f90 
rm common_mpi_wrf.f90 
rm common_obs_wrf.f90
rm module_map_utils.f90 
rm common_smooth2d.f90

echo "NORMAL END"
