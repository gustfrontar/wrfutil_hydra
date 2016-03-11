#!/bin/sh
set -ex
PGM=./radar_qc.exe
F90=ifort  #mpif90


OMP=
F90OPT='-O3 -openmp '  #-O3 -openmp'

cd ./src

#cd ./src

#PRE CLEAN
rm -f *.mod
rm -f *.o

ln -sf ../../common/common_radar_tools.f90        .
ln -sf ../../../../../common/SFMT.f90             .

#COMPILIN
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common_radar_tools.f90
$F90 $OMP $F90OPT -c common_smooth2d.f90
$F90 $OMP $F90OPT -c common_qc_tools.f90
$F90 $OMP $F90OPT -c main_qc.f90
$F90 $OMP $F90OPT -o ${PGM} *.o  

#CLEAN UP
rm -f *.mod
rm -f *.o

rm -f common_radar_tools.f90
rm -f SFMT.f90

mv *.exe ../

echo "NORMAL END"
