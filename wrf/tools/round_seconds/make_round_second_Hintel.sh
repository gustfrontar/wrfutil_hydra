#!/bin/sh
set -ex

#. /usr/share/modules/init/sh
#module unload pgi-12.10
#module load intel-2013.1.117

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH

LIB_NETCDF="-L/usr/local/lib/ -lnetcdff"
INC_NETCDF="-I/usr/local/include/ "

PGM=./round_second.exe
F90=ifort  #mpif90

OMP=
F90OPT='-O3 '

#PRE CLEAN
rm -f *.mod
rm -f *.o

cd ./src

#COMPILIN
$F90 $OMP $F90OPT $INC_NETCDF -c common_round_seconds.f90
$F90 $OMP $F90OPT -c main_round_seconds.f90
$F90 $OMP $F90OPT -o ${PGM} *.o  ${LIB_NETCDF}

#CLEAN UP
rm -f *.mod
rm -f *.o

mv $PGM ../


echo "NORMAL END"
