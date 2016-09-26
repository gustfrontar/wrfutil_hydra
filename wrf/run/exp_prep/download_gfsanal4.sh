#!/bin/bash

#DESCARGA DE LOS ANALISIS OPERATIVOS DEL GFS DE 0.5°

IDATE=20160215060000
EDATE=20160216180000
FREQ=21600
DESTDIR=/home/paula.maldonado/share/DATA/GFSANL_4/

source ./dateedit.sh

CDATE=$IDATE
while [ $CDATE -le $EDATE ]
do

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Downloading GFS operative analysis: $CDATE"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo " "

DATE=`echo $CDATE | cut -c1-8`
YEAR=`echo $CDATE | cut -c1-4`
MONTH=`echo $CDATE | cut -c5-6`
DAY=`echo $CDATE | cut -c7-8`
HOUR=`echo $CDATE | cut -c9-10`
wget -c http://nomads.ncdc.noaa.gov/data/gfsanl/${YEAR}${MONTH}/${YEAR}${MONTH}${DAY}/gfsanl_4_${YEAR}${MONTH}${DAY}_${HOUR}00_000.grb2 -a $DESTDIR/wget-report_${IDATE}

CDATE=`date_edit2  $CDATE $FREQ `
done
