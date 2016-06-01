set -x

F90=ifort
FFLAGS=

cd ./src

#
# Build grabbufr
#
$F90 $FFLAGS -O3 -o ../grabbufr grabbufr.f spbufr.f

#
# Build decoder
#
$F90 $FFLAGS -O3 -DUNDERSCORE -convert big_endian -fPIC -ffree-form dec_prepbufr.f90 -L../bufrlib -lbufr -o ../decoder




