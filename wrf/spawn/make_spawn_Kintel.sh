#!/bin/sh
set -ex
PGM=spawn.exe
F90=/usr/local/MPICH2/gnu/bin/mpif90

F90OPT='-g  ' # -Hs'


$F90 $F90OPT -c spawn.f90
$F90 $F90OPT -o ${PGM} *.o 

rm *.o

$F90 $F90OPT -c child.f90
$F90 $F90OPT -o child.exe *.o

rm *.o

