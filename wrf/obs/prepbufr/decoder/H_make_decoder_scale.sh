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
$F90 $FFLAGS -O3 -fconvert='big-endian' dec_prepbufr_scale_1hr.f90 -L../bufrlib -lbufr -o ../decoder_scale_1hr

$F90 $FFLAGS -O3 -fconvert='big-endian' dec_prepbufr_scale_3hr.f90 -L../bufrlib -lbufr -o ../decoder_scale_3hr

$F90 $FFLAGS -O3 -fconvert='big-endian' dec_prepbufr_scale_6hr.f90 -L../bufrlib -lbufr -o ../decoder_scale_6hr




