
cd ./bufrlib
rm *.o *.a
echo "Building C-routines"
icc -DUNDERSCORE -O3 -c *.c
echo "Bouilding F-routines"
ifort -O3 -c *.f
echo "Creating the library"
ar crv libbufr.a *.o
rm *.o

