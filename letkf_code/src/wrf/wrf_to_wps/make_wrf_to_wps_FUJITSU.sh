#!/bin/bash
PGM=wrf_to_wps.exe
F90=mpifrtpx
FLAGS='-O3 -Nalloc_assign' 

#NETCDF=/opt/netcdf/netcdf-4/intel/2021.4.0/            #Netcdf HYDRA
NETCDF=/home/ra000007/a04037/data/comp_libs_fujitsu/netcdf/    #Fugaku - FUJITSU
LIB_NETCDF="-L${NETCDF}/lib -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5_hl -lhdf5 "
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

tar -cvf ../../../wrf_to_wps_FUJITSU.tar ./*.exe 

echo "NORMAL END"
