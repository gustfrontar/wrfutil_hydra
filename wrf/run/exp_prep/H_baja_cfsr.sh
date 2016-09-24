#!/bin/bash

#BAJA LOS CFSR DE BAJA RESOLUCION.

ITIME=$1   #Start Time (yyyymmddhhMMss)
ETIME=$2   #End Time
INT=$3     #Frequency
DESTDIR=$4

source ../util.sh


CTIME=$ITIME
while [ $CTIME -le $ETIME ]
do
echo "Voy a bajar el GDAS correspondiente a la fecha: $CTIME"
FECHA=`echo $CTIME | cut -c1-8`
ANIO=`echo $CTIME | cut -c1-4`
MES=`echo $CTIME | cut -c5-6`
DIA=`echo $CTIME | cut -c7-8`
HORA=`echo $CTIME | cut -c9-10`

mkdir -p  $DESTDIR

curl http://nomads.ncdc.noaa.gov/modeldata/cmd_pgbh/$ANIO/${ANIO}${MES}/${ANIO}${MES}${DIA}/pgbh00.gdas.${ANIO}${MES}${DIA}${HORA}.grb2 -o $DESTDIR/pgbh00.gdas.${ANIO}${MES}${DIA}${HORA}.grb2


CTIME=`date_edit2  $CTIME $INT `

echo $CTIME
done

