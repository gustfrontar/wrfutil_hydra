#!/bin/bash
#Toma un experimento en donde habilitamos la salida del UPP del WRF y renombra los datos de forma tal que puedan ser utilizados como input para otro experimento.

ITIME=20091117010000
ETIME=20091118000000
INT=3600
ENSSIZE=60

EXPNAME="ANALYSIS_PARANA_10KM_control_paranafnl10k_60m_prepbufrandsurface_1hr_grib_Hydra"


GRIBDIR=$HOME/share/DATA/GRIB/
WRFOUTDIR=$HOME/salidas/EXPERIMENTS/

source ../util.sh

operation="ln -sf"  #La operacion puede ser ln -sf o cp (si es ln -sf no duplicamos espacio en disco)


DESTDIR=$GRIBDIR/WRF/$EXPNAME/

mkdir -p $DESTDIR

CTIME=$ITIME
while [ $CTIME -le $ETIME ] ; do

  echo "Renaming time $CTIME"

  member=1
  while [ $member -le $ENSSIZE ] ; do 

    members=`add_zeros $member 5 `
    mkdir  -p $DESTDIR/$members
    ln -sf $WRFOUTDIR/$EXPNAME/anal/$CTIME/slev${members}.grib $DESTDIR/$members/${CTIME}.grib

    member=`expr $member + 1`

  done

  CTIME=`date_edit2 $CTIME $INT  `
  

done

