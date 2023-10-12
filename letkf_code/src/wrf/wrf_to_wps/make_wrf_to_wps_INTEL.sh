#!/bin/bash
source /opt/load-libs.sh 1              #HYDRA  Load intel libraries and paths
. /opt/intel/oneapi/setvars.sh intel64  #Fugaku Load intel libraries and paths


PGM=wrf_to_wps.exe
F90=ifort

FLAGS="-convert big_endian -FR"

#NETCDF=/opt/netcdf/netcdf-4/intel/2021.4.0/            #Netcdf HYDRA
NETCDF=/home/ra000007/a04037/data/comp_libs/netcdf4/    #Fugaku - intel

LIB_NETCDF="-L${NETCDF}/lib -lnetcdff -lnetcdf"
INC_NETCDF="-I${NETCDF}/include/"


#Build wrf_to_wps
ln -sf ../../common/SFMT.f90       ./SFMT.f90
ln -sf ../common/module_map_utils.f90  ./module_map_utils.f90
ln -sf ../common/common_wrf.f90 ./common_wrf.f90
ln -sf ../../common/common.f90  ./common.f90

$F90 $FLAGS             -c SFMT.f90
$F90 $FLAGS             -c common.f90
$F90 $FLAGS             -c module_map_utils.f90
$F90 $FLAGS $INC_NETCDF -c common_wrf.f90
$F90 $FLAGS             -c met_data.f90
$F90 $FLAGS $INC_NETCDF -c wrf_to_wps.f90
$F90 $FLAGS -o ${PGM} *.o $LIB_NETCDF 

rm ./SFMT.f90
rm ./module_map_utils.f90
rm ./common_wrf.f90
rm ./common.f90
rm *.o *.mod


#Build read_intermediate
ln -sf ../../common/SFMT.f90    ./SFMT.f90
ln -sf ../../common/common.f90  ./common.f90

$F90 $FLAGS             -c SFMT.f90
$F90 $FLAGS             -c common.f90
$F90 $FLAGS             -c met_data.f90
$F90 $FLAGS             -c read_intermediate.f90
$F90 $FLAGS -o read_intermediate.exe *.o

rm ./SFMT.f90
rm ./common.f90 
rm *.o *.mod

tar -cvf ../../../wrf_to_wps_INTEL.tar ./*.exe

echo "NORMAL END"
