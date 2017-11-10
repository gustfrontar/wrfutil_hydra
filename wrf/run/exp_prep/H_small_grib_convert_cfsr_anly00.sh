#!/bin/bash

#KEEP ONLY A SMALL PART OF A GRIB.

ITIME=19980417000000
ETIME=19980430000000

INT=21600
PARALLEL=1
export OMP_NUM_THREADS=1

#REGNAME="EUROPE"
#LATRANGE="30:70"
#LONRANGE="-20:60"

REGNAME="ARGENTINA"
LATRANGE="-50:-20"
LONRANGE="270:310"

#REGNAME="SA"
#LATRANGE="-90:20"
#LONRANGE="200:360"

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
  echo $ORIGINALDIR/pgbh00.gdas.${ANIO}${MES}${DIA}${HORA}.grb2
  wgrib2 $ORIGINALDIR/pgbh00.gdas.${ANIO}${MES}${DIA}${HORA}.grb2 -set_grib_type jpeg -match ":(UGRD:10 m above ground|VGRD:10 m above ground|TMP:2 m above ground|RH:2 m above ground|SOILW|TSOIL|LAND|ICEC|WEASD|TMP:0-0.1 m below ground|TMP:0.1-0.4 m below ground|TMP:0.4-1 m below ground|TMP:1-2 m below ground):" -small_grib $LONRANGE $LATRANGE ./grib1 > log1.log

  echo $ORIGINALDIR/pgbhnl.gdas.${ANIO}${MES}${DIA}${HORA}.grb2 
  wgrib2 $ORIGINALDIR/pgbhnl.gdas.${ANIO}${MES}${DIA}${HORA}.grb2 -set_grib_type jpeg -match ":(UGRD|VGRD|TMP|HGT|HGT|RH|PRES|PRMSL):" -small_grib $LONRANGE $LATRANGE ./grib2  > log2.log 

  echo $DESTDIR/${ANIO}${MES}${DIA}${HORA}0000.grb
  wgrib2 ./grib2  -append -set_grib_type jpeg -grib ./grib1.grb  > ./log3.log
  mv grib1.grb  $DESTDIR/${ANIO}${MES}${DIA}${HORA}0000.grb 

  CTIME=`date_edit2 $CTIME $INT`
  cparallel=`expr $cparallel + 1 `
#echo $CTIME
  done
  time wait
done

