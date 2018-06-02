#!/bin/bash

#KEEP ONLY A SMALL PART OF A GRIB.

ITIME=19980401000000
ETIME=19980430180000

INT=21600
PARALLEL=1
export OMP_NUM_THREADS=1

#REGNAME="EUROPE"
#LATRANGE="30:70"
#LONRANGE="-20:60"

REGNAME="SA"
LATRANGE="-90:20"
LONRANGE="200:360"

#REGNAME="JAPAN"
#LATRANGE="25:50"
#LONRANGE="120:150"


source ../util.sh

ORIGINALDIR=$HOME/share/DATA/GRIB/CFSR/HIRES/GLOBAL/

DESTDIR=$HOME/share/DATA/GRIB/CFSR/HIRES/$REGNAME
#DESTDIRSFC=$HOME/datos/DATA/GRIB/CFSR/HIRES/${REGNAME}SFC

mkdir -p $DESTDIR
#mkdir -p $DESTDIRSFC

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
  wgrib2 $ORIGINALDIR/pgbh06.gdas.${ANIO}${MES}${DIA}${HORA}.grb2 -set_ftime "anl" -set_grib_type jpeg -match ":(UGRD|VGRD|TMP|HGT|RH|PRES|PRMSL|SOILW|TSOIL|LAND|ICEC|WEASD):" -small_grib $LONRANGE $LATRANGE ./grib1  > ./tmp.log  
 
  wgrib2 $ORIGINALDIR/flxf06.gdas.${ANIO}${MES}${DIA}${HORA}.grb2 -set_ftime "anl" -set_grib_type jpeg -match ":(UGRD|VGRD|TMP|HGT|RH|PRES|PRMSL|SOILW|TSOIL|LAND|ICEC|WEASD):" -small_grib $LONRANGE $LATRANGE ./grib2 > ./tmp.log


  wgrib2 ./grib2  -append -set_grib_type jpeg -grib  ./grib1 > ./tmp.log

  mv ./grib1 $DESTDIR/${ANIO}${MES}${DIA}${HORA}0000.grb

  CTIME=`date_edit2 $CTIME $INT`
  cparallel=`expr $cparallel + 1 `
#echo $CTIME
  done
  time wait
done

