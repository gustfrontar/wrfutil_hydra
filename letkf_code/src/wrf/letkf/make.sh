#!/bin/bash
# This script uses the environmented defined in compile_all.sh
#######
#set -ex

rm -f *.mod
rm -f *.o
rm -f *.exe

COMMONDIR=../../common/
ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/netlib.f ./
ln -fs $COMMONDIR/common.f90 ./
ln -fs $COMMONDIR/common_mpi.f90 ./
ln -fs $COMMONDIR/common_mtx.f90 ./
ln -fs $COMMONDIR/common_letkf.f90 ./
ln -fs $COMMONDIR/common_smooth2d.f90 ./

ln -fs ../common/common_wrf.f90 ./
ln -fs ../common/common_mpi_wrf.f90 ./
ln -fs ../common/common_obs_wrf.f90 ./
ln -fs ../common/module_map_utils.f90 ./
ln -fs ../common/common_namelist.f90 ./


$PARF90 $OMP $FFLAGS -c SFMT.f90
$PARF90 $OMP $FFLAGS -c common.f90
$PARF90 $OMP $FFLAGS -c common_mpi.f90
$PARF90 $OMP $FFLAGS -c common_mtx.f90
$PARF90 $OMP -O3     -c netlib.f
$PARF90 $OMP $FFLAGS -c common_letkf.f90
$PARF90 $OMP $FFLAGS -c module_map_utils.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c common_wrf.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c common_namelist.f90
$PARF90 $OMP $FFLAGS -c common_smooth2d.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c common_obs_wrf.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c common_mpi_wrf.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c letkf_obs.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c letkf_tools.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c letkf.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -o letkf.exe *.o $LBLAS $LIB_NETCDF

rm -f *.mod
rm -f *.o

#Build update_wrf_time.f90
$PARF90 $OMP $FFLAGS -c SFMT.f90
$PARF90 $OMP $FFLAGS -c common.f90
$PARF90 $OMP $FFLAGS -c module_map_utils.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c common_wrf.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c update_wrf_time.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -o update_wrf_time.exe *.o $LBLAS $LIB_NETCDF

rm -f *.mod
rm -f *.o

rm SFMT.f90
rm common.f90
rm common_mpi.f90
rm common_mtx.f90
rm common_letkf.f90
rm common_wrf.f90
rm common_mpi_wrf.f90
rm common_obs_wrf.f90
rm module_map_utils.f90
rm common_smooth2d.f90
rm common_namelist.f90
rm netlib.f

tar -cvf ${current_dir}/bin/letkf_${COMPILATION_NAME}.tar ./*.exe ./letkf.namelist

GREEN=$'\e[0;32m'
RED=$'\e[0;31m'
NC=$'\e[0m'
if [ -e ./letkf.exe ] ; then
   echo "${GREEN}#######################################################"       
   echo "${GREEN} SUCCESSFULLY COMPILED LETKF!!!!"
   echo "${GREEN}#######################################################"
   echo "${NC}"
else
   echo "${RED}#########################################################"       
   echo "${RED} ERROR: COMPILATION OF WRF LETKF!!!!"
   echo "${RED}#########################################################"       
   echo "${NC}"
   exit 1
fi


