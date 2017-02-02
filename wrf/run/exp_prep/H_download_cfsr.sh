#!/bin/bash

#DOWNLOAD HI RESOLUTION CFSR ANALYSIS

ITIME=$1   #Start Time (yyyymmddhhMMss)
ETIME=$2   #End Time
INT=$3     #Frequency
DESTDIR=$4

source ../util.sh


CTIME=$ITIME
while [ $CTIME -le $ETIME ]
do
echo "Downloading file corresponding to the following date: $CTIME"
YEAR=`echo $CTIME | cut -c1-4`
MONTH=`echo $CTIME | cut -c5-6`
DAY=`echo $CTIME | cut -c7-8`
TIME=`echo $CTIME | cut -c9-10`

mkdir -p  $DESTDIR

curl http://nomads.ncdc.noaa.gov/modeldata/cmd_pgbh/$YEAR/${YEAR}${MONTH}/${YEAR}${MONTH}${DAY}/pgbh00.gdas.${YEAR}${MONTH}${DAY}${TIME}.grb2 -o $DESTDIR/pgbh00.gdas.${YEAR}${MONTH}${DAY}${TIME}.grb2


CTIME=`date_edit2  $CTIME $INT `

echo $CTIME
done

