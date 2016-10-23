export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
F90=gfortran
F90OPT="-fconvert=big-endian -fno-underscoring"
LIB="-lblas -lhdfeos -lGctp -lmfhdf -ldf -lhdf5  -ljpeg -lsz  -lz"
PGM=decoder


$F90 $F90OPT -c common.f90
$F90 $F90OPT -c common_mtx.f90
$F90 $F90OPT -c common_airs.f90
$F90 $F90OPT -c netlib.f
$F90 $F90OPT -c dec_airsco.f90
$F90 $F90OPT -o ${PGM} *.o $LIB
rm *.o

