#!/bin/bash

#CONVIERTE A UNA REGION MAS PEQUENIA

ITIME=20060101000000
ETIME=20091231180000
INT=21600
PARALLEL=8

REGNAME="ARGENTINA"
LATRANGE="-50:-20"
LONRANGE="270:310"

#REGNAME="JAPAN"
#LATRANGE="25:50"
#LONRANGE="120:150"


source ../util.sh

DESTDIR=$HOME/share/DATA/GRIB/CFSR/HIRES/$REGNAME



CTIME=$ITIME
while [ $CTIME -le $ETIME ]
do

  echo "Voy a convertir el CFSR correspondiente a la fecha: $CTIME"
  FECHA=`echo $CTIME | cut -c1-8`
  ANIO=`echo $CTIME | cut -c1-4`
  MES=`echo $CTIME | cut -c5-6`
  DIA=`echo $CTIME | cut -c7-8`
  HORA=`echo $CTIME | cut -c9-10`

  mv $DESTDIR/argentina_pgbh00.gdas.${ANIO}${MES}${DIA}${HORA}.grb2 $DESTDIR/${ANIO}${MES}${DIA}${HORA}0000.grb  


  CTIME=`date_edit2 $CTIME $INT`
echo $CTIME
done

