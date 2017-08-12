
cd ./bufrlib
rm *.o *.a
echo "Building C-routines"
gcc -DUNDERSCORE -O3 -c *.c
echo "Bouilding F-routines"
gfortran -O3 -c *.f
echo "Creating the library"
ar crv libbufr.a *.o
rm *.o

