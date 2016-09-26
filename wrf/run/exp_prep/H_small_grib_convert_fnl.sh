#!/bin/bash

#CONVIERTE A UNA REGION MAS PEQUENIA

ITIME=20140120000000
ETIME=20140129180000
INT=21600
PARALLEL=8
export OMP_NUM_THREADS=1

REGNAME="ARGENTINA"
LATRANGE="-50:-20"
LONRANGE="270:310"

#REGNAME="JAPAN"
#LATRANGE="25:50"
#LONRANGE="120:150"


source ../util.sh

ORIGINALDIR=$HOME/datos/DATA/GRIB/FNL/HIRES/GLOBAL/

DESTDIR=$HOME/datos/DATA/GRIB/FNL/HIRES/$REGNAME

mkdir -p $DESTDIR

CTIME=$ITIME
while [ $CTIME -le $ETIME ]
do

 cparallel=1
 while [ $cparallel -le $PARALLEL -a $CTIME -le $ETIME  ]
 do 
  echo "Voy a convertir el CFSR correspondiente a la fecha: $CTIME"
  FECHA=`echo $CTIME | cut -c1-8`
  ANIO=`echo $CTIME | cut -c1-4`
  MES=`echo $CTIME | cut -c5-6`
  DIA=`echo $CTIME | cut -c7-8`
  HORA=`echo $CTIME | cut -c9-10`

  wgrib2 $ORIGINALDIR/fnl_${ANIO}${MES}${DIA}_${HORA}_00.grib2 -set_grib_type jpeg -match ":(UGRD|VGRD|TMP|HGT|RH|PRES|SOILW|PRMSL|LAND|ICEC|WEASD):" -small_grib $LONRANGE $LATRANGE $DESTDIR/${ANIO}${MES}${DIA}${HORA}0000.grb  > ./tmp.log  &


  CTIME=`date_edit2 $CTIME $INT`
  cparallel=`expr $cparallel + 1 `
#echo $CTIME
  done
  time wait
done

