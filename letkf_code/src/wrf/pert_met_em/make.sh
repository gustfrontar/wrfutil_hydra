#!/bin/bash
# This script uses the environmented defined in compile_all.sh
#######
rm -f *.mod
rm -f *.o
rm -f *.exe

COMMONDIR=../../common/

ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/common.f90 ./
ln -fs $COMMONDIR/common_mpi.f90 ./
ln -fs $COMMONDIR/netlib.f ./

ln -fs ../common/common_met_em.f90 ./
ln -fs ../common/common_mpi_met_em.f90 ./
ln -fs ../common/common_namelist_met_em.f90 ./

$PARF90 $OMP $FFLAGS -c SFMT.f90
$PARF90 $OMP $FFLAGS -c common.f90
$PARF90 $OMP $FFLAGS -c common_mpi.f90
$PARF90 $OMP         -c netlib.f
$PARF90 $OMP $FFLAGS $INC_NETCDF -c common_met_em.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c common_namelist_met_em.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c common_mpi_met_em.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c met_em_tools.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -c pert_met_em.f90
$PARF90 $OMP $FFLAGS $INC_NETCDF -o pert_met_em.exe *.o $LBLAS $LIB_NETCDF

rm -f *.mod
rm -f *.o

$PARF90 $FFLAGS -c dummy_mpi.f90 
$PARF90         -o dummy_mpi.exe *.o

tar -cvf ${current_dir}/bin/pert_met_em_${COMPILATION_NAME}.tar ./*.exe 

GREEN=$'\e[0;32m'
RED=$'\e[0;31m'
NC=$'\e[0m'
if [ -e ./pert_met_em.exe ] ; then
   echo "${GREEN}#######################################################"
   echo "${GREEN} SUCCESSFULLY COMPILED PERT MET_EM!!!!"
   echo "${GREEN}#######################################################"
   echo "${NC}"
else
   echo "${RED}#########################################################"
   echo "${RED} ERROR: COMPILATION OF PERT MET_EM FAILED!!!!"
   echo "${RED}#########################################################"
   echo "${NC}"
   exit 1
fi

