#!/bin/sh
set -ex

#. /usr/share/modules/init/sh
#module unload pgi-12.10
#module load intel-2013.1.117

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH

#LIB_NETCDF="-L/usr/local/lib/ -lnetcdff"
#INC_NETCDF="-I/usr/local/include/ "

LIB_NETCDF="-L/opt/netcdf-fortran/4.2/lib/ -lnetcdff"
INC_NETCDF="-I/opt/netcdf-fortran/4.2/include/ "


PGM=../compute_pert_metem.exe
F90=ifort  #mpif90



OMP=
F90OPT='-convert big_endian -O3 '

cd ./src

#PRE CLEAN
rm -f *.mod
rm -f *.o


#COMPILIN
$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c common_namelist.f90
$F90 $OMP $F90OPT -c common_smooth2d.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_metem_memnc.f90
$F90 $OMP $F90OPT -c common_perturb_ensemble_metem.f90
$F90 $OMP $F90OPT -c main_compute_pert_metem.f90
$F90 $OMP $F90OPT -o ${PGM} *.o  ${LIB_NETCDF}

#mv *.exe ../

#CLEAN UP
rm -f *.mod
rm -f *.o


echo "NORMAL END"
