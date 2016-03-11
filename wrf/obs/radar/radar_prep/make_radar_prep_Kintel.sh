#!/bin/sh
set -ex
PGM=./radar_prep.exe

F90=ifort  #mpif90

OMP=
F90OPT='-O3 -convert big_endian'

cd ./src

#PRE CLEAN
rm -f *.mod
rm -f *.o


ln -sf ../../../../../common/common_smooth2d.f90 .
ln -sf ../../../../../common/common.f90       .
ln -sf ../../../../../common/SFMT.f90         .
ln -sf ../../common/common_radar_tools.f90    .

#COMPILING
$F90 $OMP $F90OPT -c common.f90
$F90 $OMP $F90OPT -c common_smooth2d.f90
$F90 $OMP $F90OPT -c SFMT.f90
$F90 $OMP $F90OPT -c common_radar_tools.f90
$F90 $OMP $F90OPT -c common_namelist.f90
$F90 $OMP $F90OPT -c common_superobbing.f90
$F90 $OMP $F90OPT -c main_radar_prep.f90
$F90 $OMP $F90OPT -o ${PGM} *.o 

mv *.exe ../



#CLEAN UP
rm -f *.mod
rm -f *.o
rm -f common_smooth2d.f90
rm -f common.f90
rm -f SFMT.f90
rm -f common_radar_tools.f90


echo "NORMAL END"
