#!/bin/sh
set -ex
PGM=letkf_inteleigenexa.exe
F90=mpif90
OMP= #'-openomp'

F90OPT='-convert big_endian' # -Hs'

BLAS=0 #0: no blas 1: using blas
BASEDIR=${HOME}/libintel/netcdf-3.6.3/

INC_EIGEN="-I${HOME}/share/EigenExa-2.3cintel/"
LIB_EIGEN="-L${HOME}/share/EigenExa-2.3cintel/ -lEigenExa"

echo $LIB_EIGEN

rm -f *.mod
rm -f *.o

COMMONDIR=../../common/

cat $COMMONDIR/netlib.f > netlib2.f
if test $BLAS -eq 1
then
LSCALAPACK="-SCALAPACK"

else
cat $COMMONDIR/netlibblas.f >> netlib2.f
LBLAS=""
fi

LIB_NETCDF="-L${BASEDIR}/lib -lnetcdf"
INC_NETCDF="-I${BASEDIR}/include/"


COMMONDIR=../../common/
ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/common.f90 ./
ln -fs $COMMONDIR/common_mpi.f90 ./
ln -fs $COMMONDIR/common_mtx_eigenexa.f90 ./
ln -fs $COMMONDIR/common_letkf.f90 ./
ln -fs $COMMONDIR/common_smooth2d.f90  ./

ln -fs ../common/common_wrf.f90 ./
ln -fs ../common/common_mpi_wrf.f90 ./
ln -fs ../common/common_obs_wrf.f90 ./
ln -fs ../common/module_map_utils.f90 ./
ln -fs ../common/common_namelist.f90 ./


$F90 $F90OPT -c SFMT.f90
$F90 $F90OPT -c common.f90
$F90 $F90OPT -c common_mpi.f90
$F90 $F90OPT $INC_EIGEN $LIB_EIGEN -c common_mtx_eigenexa.f90
$F90 $F90OPT -c netlib2.f
$F90 $F90OPT -c common_letkf.f90
$F90 $F90OPT -c module_map_utils.f90
$F90 $F90OPT -c common_wrf.f90
$F90 $F90OPT -c common_namelist.f90
$F90 $F90OPT -c common_smooth2d.f90
$F90 $F90OPT -c common_obs_wrf.f90
$F90 $F90OPT -c common_mpi_wrf.f90
$F90 $F90OPT -c letkf_obs.f90
$F90 $F90OPT -c letkf_tools.f90
$F90 $F90OPT $INC_EIGEN $LIB_EIGEN -c letkf_eigenexa.f90
$F90 $OMP $F90OPT -o ${PGM} *.o $LBLAS $LSCALAPACK $LIB_NETCDF $LIB_EIGEN

rm -f *.mod
rm -f *.o
rm -f netlib2.f

rm SFMT.f90 
rm common.f90 
rm common_mpi.f90 
rm common_mtx_eigenexa.f90 
rm common_letkf.f90 
rm common_wrf.f90 
rm common_mpi_wrf.f90 
rm common_obs_wrf.f90
rm module_map_utils.f90 
rm common_smooth2d.f90
rm common_namelist.f90

echo "NORMAL END"
