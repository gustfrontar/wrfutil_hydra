#!/bin/bash
PGM=merge_wrfinput.exe
F90=mpifrtpx

FLAGS="-O2"

BASEDIR=${HOME}/share/Libs/

LIB_NETCDF="-L${BASEDIR}/netcdf/4.1.1/lib -lnetcdf"
INC_NETCDF="-I${BASEDIR}/netcdf/4.1.1/include/"

#Build wrf_to_wps

$F90 $FLAGS $INC_NETCDF -c merge_wrfinput.f90
$F90 $FLAGS -o ${PGM} *.o $LIB_NETCDF 

rm *.o

echo "NORMAL END"
