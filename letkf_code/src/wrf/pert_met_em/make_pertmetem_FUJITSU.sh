#!/bin/bash
#######
set -ex
F90=mpifrtpx
OMP='' #'-Kopenmp'

F90OPT='-O3 -Nalloc_assign' 
BLAS=1 #0: no blas 1: using blas

rm -f *.mod
rm -f *.o

COMMONDIR=../../common/

cat $COMMONDIR/netlib.f > netlib2.f
if [ $BLAS -eq 1 ] ;then
   LBLAS="-SSL2BLAMP"
   LSCALAPACK="-SCALAPACK"
else
   cat $COMMONDIR/netlibblas.f >> netlib2.f
   LBLAS=""
fi

NETCDF=/home/ra000007/a04037/data/comp_libs_fujitsu/netcdf/        #Fugaku

LIB_NETCDF="-L$NETCDF/lib/ -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5_hl -lhdf5  -lfjprofmpi -lmpi_cxx"
INC_NETCDF="-I$NETCDF/include/"

COMMONDIR=../../common/
ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/common.f90 ./
ln -fs $COMMONDIR/common_mpi.f90 ./

ln -fs ../common/common_met_em.f90 ./
ln -fs ../common/common_mpi_met_em.f90 ./
ln -fs ../common/common_namelist_met_em.f90 ./

$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c common_mpi.f90
$F90 $OMP $F90OPT -c netlib2.f
$F90 $OMP $F90OPT $INC_NETCDF -c common_met_em.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_namelist_met_em.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_mpi_met_em.f90
$F90 $OMP $F90OPT $INC_NETCDF -c met_em_tools.f90
$F90 $OMP $F90OPT $INC_NETCDF -c pert_met_em.f90
$F90 $OMP $F90OPT $INC_NETCDF -o pert_met_em.exe *.o $LBLAS $LIB_NETCDF

rm -f *.mod
rm -f *.o
rm -f netlib.f

tar -cvf ../../../pertmetem_INTEL.tar ./*.exe ./pertmetem.namelist

echo "NORMAL END"
