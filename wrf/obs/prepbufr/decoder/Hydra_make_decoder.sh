set -x

F90=gfortran
FFLAGS=

cd ./src

#
# Build grabbufr
#
$F90 $FFLAGS -O3 -o ../grabbufr grabbufr.f spbufr.f

#
# Build decoder
#
$F90 $FFLAGS -O3 -fconvert='big-endian' dec_prepbufr.f90 -L../bufrlib -lbufr -o ../decoder




