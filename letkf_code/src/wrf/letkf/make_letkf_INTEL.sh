#!/bin/bash

#####
# Nota: Correr el source envars.sh , ppor ejemplo del /share/apps/Build_WRF_INTEL -- primncipalemte para setear el path de los compiladore y porque NECESITA !!!!!! Setear la variable NETCDF con el path al netcdf.
# Note2: Como debe tomar del entorno la variable NETCDF, hay que correr este script con source 
#######
set -ex
PGM=letkf.exe
F90=mpiifort
#OMP=''
OMP='-qopenmp'

F90OPT='-O3 -xHost ' #-convert big_endian' #-O3 -convert big_endian' #-Kfast,parallel' # -Hs'


BLAS=1 #0: no blas 1: using blas

rm -f *.mod
rm -f *.o

COMMONDIR=../../common/

cat $COMMONDIR/netlib.f > netlib2.f
if test $BLAS -eq 1
then
LBLAS="-mkl"
#INCLUDE= #"-I${BASEDIR}/EigenExa-2.3c/"
#Leigen= #"-L${BASEDIR}/EigenExa-2.3c/ -lEigenExa"

else
cat $COMMONDIR/netlibblas.f >> netlib2.f
LBLAS=""
fi

#NETCDF=/opt/netcdf/netcdf_c-4.8.1_fortran-4.5.3/intel/2021.4.0/   #Hydra
NETCDF=/home/ra000007/a04037/data/comp_libs/netcdf/                #Fugaku - intel


LIB_NETCDF="-L$NETCDF/lib/ -lnetcdff -lnetcdf"
INC_NETCDF="-I$NETCDF/include/"


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


$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c common_mpi.f90
$F90 $OMP $F90OPT -c common_mtx.f90
$F90 $OMP $F90OPT -c netlib2.f
$F90 $OMP $F90OPT -c common_letkf.f90
$F90 $OMP $F90OPT -c module_map_utils.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_wrf.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_namelist.f90
$F90 $OMP $F90OPT -c common_smooth2d.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_obs_wrf.f90
$F90 $OMP $F90OPT $INC_NETCDF -c common_mpi_wrf.f90
$F90 $OMP $F90OPT $INC_NETCDF -c letkf_obs.f90
$F90 $OMP $F90OPT $INC_NETCDF -c letkf_tools.f90
$F90 $OMP $F90OPT $INC_NETCDF -c letkf.f90
$F90 $OMP $F90OPT $INC_NETCDF -o ${PGM} *.o $LBLAS $LIB_NETCDF

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


#cd ../wrf_to_wps

#PGM=merge_wrfinput.exe

#FLAGS="-O2"

##Build wrf_to_wps

#$F90 $FLAGS $INC_NETCDF -c merge_wrfinput.f90
#$F90 $FLAGS -o ${PGM} *.o $LIB_NETCDF

#rm *.o

tar -cvf ../../../letkf.tar ./*.exe ./letkf.namelist

echo "NORMAL END"
