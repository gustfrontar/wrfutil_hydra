#!/bin/bash

source util.sh 

EXPNAME=OsakaPAR_1km_control1000m_smallrandompert_new
INIDATE=20130713051030  #Assimilation start date
ENDDATE=20130713053930  #Assimilation end date.
FREQ=30                 #Assimilation frequency (seconds)
MAX_DOWNLOAD=1          #Maximum number of simultaneous downloads.


USER=a03094

REMOTEPATH=/home/ra000015/a03094/data/exp/
LOCALPATH=/home/jruiz/share/scale_input_data/

MEMBER=0001

GET_ANALYSIS=1
GET_GUES=1

mkdir -p $LOCALPATH


CDATE=$INIDATE

while [ $CDATE -le $ENDDATE  ] ; do

 N_DOWNLOADS=1
 while [ $N_DOWNLOADS -le $MAX_DOWNLOAD -a $CDATE -le $ENDDATE ] ; do

 if [ $GET_ANALYSIS -eq 1 ] ; then
   mkdir -p $LOCALPATH/$EXPNAME/$CDATE/anal/$MEMBER
 fi
 if [ $GET_GUES -eq 1 ] ; then
   mkdir -p $LOCALPATH/$EXPNAME/$CDATE/gues/$MEMBER
 fi

 if [ $GET_ANALYSIS -eq 1 ] ; then
      scp -r ${USER}@k.aics.riken.jp:$REMOTEPATH/$EXPNAME/$CDATE/anal/$MEMBER/ $LOCALPATH/$EXPNAME/$CDATE/anal/  &
 fi
 if [ $GET_GUES -eq 1 ] ; then
      scp -r ${USER}@k.aics.riken.jp:$REMOTEPATH/$EXPNAME/$CDATE/gues/$MEMBER/ $LOCALPATH/$EXPNAME/$CDATE/gues/  &
 fi

 CDATE=`date_edit2 $CDATE $FREQ `
 N_DOWNLOADS=`expr $N_DOWNLOADS + 1 `
 done
 time wait

done





