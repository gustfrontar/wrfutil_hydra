#!/bin/bash

F90='ifort'

FFLAGS='-static-intel -assume underscore'

DIRLIB1=" -L$HDFEOS/lib \
          -L$HDF5/lib \
          -L$HDF4/lib \
          -L$ZLIB/lib \
	  -L$JASPER/lib \
          -L$SZLIB/lib"

#$F90  $FFLAGS -o binarios_airs dec_airsret_smn_rra.f90   $DIRLIB1 -lhdfeos -lGctp  -lhdf5  -lmfhdf -ldf  -lz -lm -ljpeg  -lsz
#echo "$F90 $FFLAGS -o decoder.airs dec_airsret_smn_rra_select_level.f90   $DIRLIB1 -lhdfeos -lGctp  -lhdf5  -lmfhdf -ldf  -lz -lm -ljpeg  -lsz"
$F90 $FFLAGS -o decoder.airs dec_airsret_smn_rra_select_level.f90   $DIRLIB1 -lm -Bstatic -lhdfeos -lGctp  -lhdf5  -lmfhdf -ldf  -lz -ljpeg  -lsz


####3
## Lo siguiente es un un test de que todo este funcionando correctamente
export OBSWIN=60
export DESPLAZAMIENTO=60
export OBSFREC=10
./decoder.airs AIRS.2019.09.05.053.L2.RetStd_IR.v6.0.31.0.R19248025727.hdf 200 2019 09 05 05 17

if [ -s AIRSRT_20190905052000.dat ]
then
     echo "Compilacion EXITOSA !!!!!"
else
     echo "ERROR !!! Verificar las librerias EOS2 libjpeg jasper hdf4"
fi
