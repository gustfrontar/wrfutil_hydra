#!/bin/sh
set -ex

#Este script sirve para debugear los modulos usando el compilador.
#

#PGM=./scale_to_radar.exe
F90=ifort  

LIB_NETCDF="-L/usr/local/netcdf4.intel/lib/ -lnetcdff"
INC_NETCDF="-I/usr/local/netcdf4.intel/include/"


OMP=
F90OPT='-g -traceback'

#PRE CLEAN
rm -f *.mod
rm -f *.o


#COMPILIN
$F90 $OMP $F90OPT $INC_NETCDF -c common_qc_tools.f90

#CLEAN UP
rm -f *.mod
rm -f *.o

echo "NORMAL END"
