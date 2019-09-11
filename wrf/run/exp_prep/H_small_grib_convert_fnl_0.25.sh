#!/bin/bash
#CONVIERTE A UNA REGION MAS PEQUENIA

ITIME=20181211000000
ETIME=20181212120000
INT=21600
PARALLEL=1
export OMP_NUM_THREADS=1

REGNAME="SA"
LATRANGE="-90:20"
LONRANGE="200:360"

WGRIB2=wgrib2

GRIBDIR=$HOME/share/DATA/GRIB/

source ../util.sh

ORIGINALDIR=$GRIBDIR/FNL/HIRES/GLOBAL/

DESTDIR=$GRIBDIR/FNL/HIRES/$REGNAME

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

  #gdas1.fnl0p25.2018121100.f00.grib2
  $WGRIB2 $ORIGINALDIR/gdas1.fnl0p25.${ANIO}${MES}${DIA}${HORA}.f00.grib2 -set_grib_type jpeg -match ":(UGRD|VGRD|TMP|HGT|RH|PRES|SOILW|TSOIL|PRMSL|LAND|ICEC|WEASD):" -small_grib $LONRANGE $LATRANGE $DESTDIR/${ANIO}${MES}${DIA}${HORA}0000.grb  > ./tmp.log  &

  CTIME=`date_edit2 $CTIME $INT`
  cparallel=`expr $cparallel + 1 `

  done

  time wait

done

