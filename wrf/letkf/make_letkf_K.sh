#!/bin/sh
set -ex
PGM=letkf_noeigenexa.exe
F90=mpifrtpx
OMP=

F90OPT='-Kfast' #-Kfast,parallel' # -Hs'

BLAS=0 #0: no blas 1: using blas
BASEDIR=${HOME}/share/Libs/

rm -f *.mod
rm -f *.o

COMMONDIR=../../common/

cat $COMMONDIR/netlib.f > netlib2.f
if test $BLAS -eq 1
then
LBLAS="-SSL2BLAMP"
LSCALAPACK="-SCALAPACK"
INCLUDE= #"-I${BASEDIR}/EigenExa-2.3c/"
Leigen= #"-L${BASEDIR}/EigenExa-2.3c/ -lEigenExa"

else
cat $COMMONDIR/netlibblas.f >> netlib2.f
LBLAS=""
LSCALAPACK=""
fi

LIB_NETCDF="-L${BASEDIR}/netcdf/4.1.1/lib -lnetcdff"
INC_NETCDF="-I${BASEDIR}/netcdf/4.1.1/include/"


COMMONDIR=../../common/
ln -fs $COMMONDIR/SFMT.f90 ./
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


$F90 $OMP $F90OPT $INLINE -c SFMT.f90
$F90 $OMP $F90OPT $INLINE -c common.f90
$F90 $OMP $F90OPT -c common_mpi.f90
$F90 $OMP $F90OPT $INLINE $INCLUDE $Leigen -c common_mtx.f90
$F90 $OMP $F90OPT $INLINE -c netlib2.f
$F90 $OMP $F90OPT $INCLUDE $Leigen -c common_letkf.f90
$F90 $OMP $F90OPT $INLINE -c module_map_utils.f90
$F90 $OMP $F90OPT $INLINE $INC_NETCDF -c common_wrf.f90
$F90 $OMP $F90OPT $INLINE -c common_namelist.f90
$F90 $OMP $F90OPT $INLINE -c common_smooth2d.f90
$F90 $OMP $F90OPT -c common_obs_wrf.f90
$F90 $OMP $F90OPT -c common_mpi_wrf.f90
$F90 $OMP $F90OPT $INCLUDE $Leigen -c letkf_obs.f90
$F90 $OMP $F90OPT $INCLUDE $Leigen -c letkf_tools.f90
$F90 $OMP $F90OPT $INCLUDE $Leigen -c letkf.f90
$F90 $OMP $F90OPT -o ${PGM} *.o $LBLAS $LSCALAPACK $LIB_NETCDF $Leigen

rm -f *.mod
rm -f *.o
rm -f netlib2.f

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

echo "NORMAL END"
