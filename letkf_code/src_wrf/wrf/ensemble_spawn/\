#!/bin/bash
PGM=spawn.exe
F90=mpif90

FLAGS=""

#Build  the program
$F90 $FLAGS            -c spawn.f90
$F90 $FLAGS -o ${PGM} *.o 

rm *.o 

tar -cvf ../../spawn.tar ./*.exe

echo "NORMAL END"
