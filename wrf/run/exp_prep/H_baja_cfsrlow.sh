#!/bin/bash

#BAJA LOS CFSR DE BAJA RESOLUCION.

ITIME=$1   #Start Time
ETIME=$2   #End Time
INT=$3     #Frequency
DESTDIR=$4

source ../util.sh


PROXY="proxy.fcen.uba.ar:8080"
CURL="/usr/bin/curl"

CTIME=$ITIME
mkdir -p $DESTDIR


while [ $CTIME -le $ETIME ]
do
echo "Downloading the following file: $CTIME"
FECHA=`echo $CTIME | cut -c1-8`
ANIO=`echo $CTIME | cut -c1-4`
MES=`echo $CTIME | cut -c5-6`
DIA=`echo $CTIME | cut -c7-8`
HORA=`echo $CTIME | cut -c9-10`

curl --proxy ${PROXY}  http://nomads.ncdc.noaa.gov/modeldata/cmd_grblow/$ANIO/${ANIO}${MES}/${ANIO}${MES}${DIA}/pgbl00.gdas.${ANIO}${MES}${DIA}${HORA}.grb2 -o $DESTDIR/pgbl00.gdas.${ANIO}${MES}${DIA}${HORA}.grb2

CTIME=`date_edit2  $CTIME $INT `
echo $CTIME
done
