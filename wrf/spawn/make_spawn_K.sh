#!/bin/sh
set -ex
PGM=spawn.exe
F90=mpifrtpx


F90OPT=' ' # -Hs'


$F90 $F90OPT -c spawn.f90
$F90 $F90OPT -o ${PGM} *.o 

rm *.o

$F90 $F90OPT -c child.f90
$F90 $F90OPT -o child.exe *.o

rm *.o

