#!/bin/bash
export PATH=$PATH:/home/devel/QC_radar/qc-radar
set -ex

export SO_DATAROOT='/data/OBS/qc-radar/asimilacion/'
export SO_INSTRUMENT_LIST='ANG,RMA1,RMA2,RMA3,RMA4,RMA5,RMA6,RMA7,RMA8,RMA10,RMA11,PAR'
export SO_REFERENCE_DATE='20200918_120000'
export SO_OUTPATH='./'
export DIROBSBIN=$SO_OUTPATH
export SO_INSTRUMENT_LIST="ANG,RMA1,RMA2,RMA3,RMA4,RMA5,RMA6,RMA7,RMA8,RMA10,RMA11,PAR"
export SO_GRID='10000,1000,15e3,240e3'     #resolucion horizontal, resolucion vertical, altura maxima, distancia horizontal maxima (en metros)
export SO_VARS='cref'                      #Variables corregidas a las que se les aplica el superobbing
export SO_CREF_OPT='4001,5,-0.1'           #Error, valor minimo de la variable
export SO_WINDOW_LENGTH='10'               #Longitud de la ventana del slot (en minutos)
export SO_MIN_NOBS='10'                    #Cantidad minima de observaciones en una caja para calcular el promedio por cajas

source activate wrfutil
export OMP_NUM_THREADS=12
export GOMP_STACKSIZE=512m
export OMP_STACKSIZE=512m
ulimit -s unlimited
python  so-radar.py
source deactivate 
