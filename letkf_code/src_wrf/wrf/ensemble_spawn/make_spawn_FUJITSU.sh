#!/bin/bash
#Env set
########

F90="mpifrtpx"

FLAGS='' #"-traceback"


#Build  the program
echo "Building spawn.exe"
$F90 $FLAGS -c spawn.f90
$F90 $FLAGS -o spawn.exe *.o 
rm *.o 

#Build  the program
echo "Building spawn_serial.exe"
$F90 $FLAGS -c spawn_serial.f90
$F90 $FLAGS -o spawn_serial.exe *.o
rm *.o


#Build the program
echo "Building mpi_serial_wrapper" 
$F90 $FLAGS -c mpi_serial_wrapper.f90
$F90 $FLAGS -o mpi_serial_wrapper.exe *.o
rm *.o

#Build the program
echo "Building ungrib_wrapper"
$F90 $FLAGS -c mpi_ungrib_wrapper.f90
$F90 $FLAGS -o mpi_ungrib_wrapper.exe *.o
rm *.o


tar -cvf ../../../spawn_FUJITSU.tar ./*.exe

echo "NORMAL END"
