#!/bin/bash
# This script uses the environmented defined in compile_all.sh
#######
#set -ex

rm -f *.mod
rm -f *.o
rm -f *.exe

#Build wrf_to_wps
ln -sf ../../common/SFMT.f90       ./SFMT.f90
ln -sf ../common/module_map_utils.f90  ./module_map_utils.f90
ln -sf ../common/common_wrf.f90 ./common_wrf.f90
ln -sf ../../common/common.f90  ./common.f90

$PARF90 $FFLAGS             -c SFMT.f90
$PARF90 $FFLAGS             -c common.f90
$PARF90 $FFLAGS             -c module_map_utils.f90
$PARF90 $FFLAGS $INC_NETCDF -c common_wrf.f90
$PARF90 $FFLAGS             -c met_data.f90
$PARF90 $FFLAGS $INC_NETCDF -c wrf_to_wps.f90
$PARF90 $FFLAGS -o wrf_to_wps.exe *.o $LIB_NETCDF 

rm ./SFMT.f90
rm ./module_map_utils.f90
rm ./common_wrf.f90
rm ./common.f90
rm *.o *.mod

#Build interpolate_intermediate
ln -sf ../../common/SFMT.f90       ./SFMT.f90
ln -sf ../common/module_map_utils.f90  ./module_map_utils.f90
ln -sf ../common/common_wrf.f90 ./common_wrf.f90
ln -sf ../../common/common.f90  ./common.f90

$F90 $FFLAGS             -c SFMT.f90
$F90 $FFLAGS             -c common.f90
$F90 $FFLAGS             -c met_data.f90
$F90 $FFLAGS             -c interpolate_intermediate.f90
$F90 $FFLAGS -o interpolate_intermediate.exe *.o

rm ./SFMT.f90
rm ./module_map_utils.f90
rm ./common_wrf.f90
rm ./common.f90
rm *.o *.mod

#Build read_intermediate
ln -sf ../../common/SFMT.f90    ./SFMT.f90
ln -sf ../../common/common.f90  ./common.f90

$F90 $FFLAGS             -c SFMT.f90
$F90 $FFLAGS             -c common.f90
$F90 $FFLAGS             -c met_data.f90
$F90 $FFLAGS             -c read_intermediate.f90
$F90 $FFLAGS -o read_intermediate.exe *.o

rm ./SFMT.f90
rm ./common.f90 
rm *.o *.mod

tar -cvf ${current_dir}/bin/wrf_to_wps_${COMPILATION_NAME}.tar ./*.exe 

GREEN=$'\e[0;32m'
RED=$'\e[0;31m'
NC=$'\e[0m'
if [ -e ./wrf_to_wps.exe ] ; then
   echo "${GREEN}#######################################################"	
   echo "${GREEN} SUCCESSFULLY COMPILED WRF TO WPS!!!!"
   echo "${GREEN}#######################################################"
   echo "${NC}"
else
   echo "${RED}#########################################################"       
   echo "${RED} ERROR: COMPILATION OF WRF TO WPS FAILED!!!!"
   echo "${RED}#########################################################"	
   echo "${NC}"
   exit 1
fi



