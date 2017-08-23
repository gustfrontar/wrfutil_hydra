#!/bin/bash

OBS='ANGUIL'
CASE='20100111'
DATA='240'
INPUT='region_dealiased_data'

BASEDIR='/home/paula.maldonado/datosmate/RADAROBS'
FILEDIR="$BASEDIR/$OBS/$CASE/$DATA/$INPUT"

cd $FILEDIR

for i in "$FILEDIR"/*.nc; do
      nc3="${i}3"
      echo $nc3
      nccopy -k classic $i ${nc3}
done

if [ "$DATA" = "240" ];then
  rm *.nc
fi


