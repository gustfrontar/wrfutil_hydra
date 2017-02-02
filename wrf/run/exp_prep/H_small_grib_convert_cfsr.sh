#!/bin/bash

#CONVIERTE A UNA REGION MAS PEQUENIA

ITIME=20060101000000
ETIME=20091231180000
INT=21600
PARALLEL=8

#Select a region or create a new region.
REGNAME="EUROPE"
LATRANGE="30:80"
LONRANGE="-30:40"

#REGNAME="ARGENTINA"
#LATRANGE="-50:-20"
#LONRANGE="270:310"

#REGNAME="JAPAN"
#LATRANGE="25:50"
#LONRANGE="120:150"

source ../util.sh

ORIGINALDIR=$HOME/share/DATA/GRIB/CFSR/HIRES/GLOBAL/

DESTDIR=$HOME/share/DATA/GRIB/CFSR/HIRES/$REGNAME


CTIME=$ITIME
while [ $CTIME -le $ETIME ]
do

 cparallel=1
 while [ $cparallel -le $PARALLEL -a $CTIME -le $ETIME  ]
 do 
  echo "Reducing the file corresponding to the following time: $CTIME"
  YEAR=`echo $CTIME  | cut -c1-4`
  MONTH=`echo $CTIME | cut -c5-6`
  DAY=`echo $CTIME   | cut -c7-8`
  HOUR=`echo $CTIME  | cut -c9-10`

  wgrib2 $ORIGINALDIR/pgbh00.gdas.${YEAR}${MONTH}${DAY}${HOUR}.grb2 -set_grib_type jpeg -match ":(UGRD|VGRD|TMP|HGT|RH|PRES|SOILW|TSOIL|PRMSL|LAND|ICEC|WEASD):" -small_grib $LONRANGE $LATRANGE $DESTDIR/${YEAR}${MONTH}${DAY}${HOUR}0000.grb  > ./tmp.log  &


  CTIME=`date_edit2 $CTIME $INT`
  cparallel=`expr $cparallel + 1 `

  done
  time wait
done

