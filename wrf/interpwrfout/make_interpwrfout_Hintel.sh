#!/bin/sh
set -ex

#. /usr/share/modules/init/sh
#module unload pgi-12.10
#module load intel-2013.1.117

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH

LIB_NETCDF="-L/opt/netcdf-fortran/4.2/lib/ -lnetcdff"
INC_NETCDF="-I/opt/netcdf-fortran/4.2/include/ "


PGM=./interpwrfout.exe
F90=ifort  #mpif90



OMP=
F90OPT='-O3 '

cd ./src

#PRE CLEAN
rm -f *.mod
rm -f *.o

ln -sf ../../../common/common.f90        ./
ln -sf ../../../common/SFMT.f90          ./
ln -sf ../../common/module_map_utils.f90 ./

#COMPILIN
$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c module_map_utils.f90
$F90 $OMP $F90OPT -c ${INC_NETCDF} ${LIB_NETCDF} common_wrf.f90
$F90 $OMP $F90OPT -c main_interp_wrfout.f90
$F90 $OMP $F90OPT -o ${PGM} *.o  ${LIB_NETCDF}

mv *.exe ../

#CLEAN UP
rm -f *.mod
rm -f *.o
rm module_map_utils.f90
rm common.f90
rm SFMT.f90


echo "NORMAL END"
