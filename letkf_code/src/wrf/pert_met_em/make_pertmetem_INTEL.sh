#!/bin/bash
#######
source /opt/load-libs.sh 1              #HYDRA  Load intel libraries
. /opt/intel/oneapi/setvars.sh intel64  #FUGAKU Load intel libraries

set -ex
F90=mpiifort
OMP=''

F90OPT='-O3 -xHost ' #-convert big_endian' #-O3 -convert big_endian' #-Kfast,parallel' # -Hs'


BLAS=1 #0: no blas 1: using blas

rm -f *.mod
rm -f *.o

COMMONDIR=../../common/

cat $COMMONDIR/netlib.f > netlib2.f
if test $BLAS -eq 1
then
LBLAS="-mkl"

else
cat $COMMONDIR/netlibblas.f >> netlib2.f
LBLAS=""
fi

NETCDF=/home/ra000007/a04037/data/comp_libs/netcdf4/                #Fugaku - intel


LIB_NETCDF="-L$NETCDF/lib/ -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5_hl -lhdf5 "
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

$F90 $OMP $F90OPT -c dummy_mpi.f90 
$F90              -o dummy_mpi.exe *.o

tar -cvf ../../../pertmetem_INTEL.tar ./*.exe ./pertmetem.namelist


echo "NORMAL END"
