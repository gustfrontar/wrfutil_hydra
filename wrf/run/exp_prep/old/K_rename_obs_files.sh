#!/bin/bash
#=======================================================================
#   This script renames obs files for conventional observations and airs
#   data. Original files are stores as t-3, t-2, t-1 whit respectto the 
#   correspinding center of the assimilation window. 
#   This scripts generates one file per hour with the format
#   yyyymmddhhmmss.dat 
#   Observations at -3 and +3 from consecutive windows are combined in a
#   single file.
#=======================================================================
OBSDIR=/home/hp150019/k02016/share/OBS/ucar_airs_th3x3_corrected/
OBSFREC=3600        #Observation files frequency (secs)
IE=20080806000000   #Experiment initial date.
EE=20080931180000   #Experiment end date.

CWD=`pwd`

source ../util.sh

CDATE=$IE

while [  $CDATE -le ${EE}  ] ; do
# --- set file date ---
# make namelist (model level) ---
echo "PROCESING ${CDATE} "

        cy=`echo $CDATE | cut -c1-4`
        cm=`echo $CDATE | cut -c5-6`
        cd=`echo $CDATE | cut -c7-8`
        ch=`echo $CDATE | cut -c9-10`
        cn=`echo $CDATE | cut -c11-12`
        cs=`echo $CDATE | cut -c13-14`
echo $ch
if [ $ch == 00 -o $ch == 06 -o $ch == 12 -o $ch == 18 ] ; then
DIRDATE=$CDATE
OBSFILE=t.dat
fi

if [ $ch == 23 -o $ch == 05 -o $ch == 11 -o $ch == 17 ] ; then
DIRDATE=`date_edit2 $CDATE 3600 `
OBSFILE=t-1.dat
fi

if [ $ch == 22 -o $ch == 04 -o $ch == 10 -o $ch == 16 ] ; then
DIRDATE=`date_edit2 $CDATE 7200 `
OBSFILE=t-2.dat
fi

if [ $ch == 01 -o $ch == 07 -o $ch == 13 -o $ch == 19 ] ; then
DIRDATE=`date_edit2 $CDATE -3600 `
OBSFILE=t+1.dat
fi

if [ $ch == 02 -o $ch == 08 -o $ch == 14 -o $ch == 20 ] ; then
DIRDATE=`date_edit2 $CDATE -7200 `
OBSFILE=t+2.dat
fi
DIRDATE=`echo $DIRDATE | cut -c1-10`
cp $OBSDIR/obs${DIRDATE}/$OBSFILE $OBSDIR/${CDATE}.dat


if [ $ch == 21 -o $ch == 03 -o $ch == 09 -o $ch == 15 ] ; then
DIRDATE=`date_edit2 $CDATE -10800 `
DIRDATE=`echo $DIRDATE | cut -c1-10`
cp $OBSDIR/obs${DIRDATE}/t+3.dat $OBSDIR/tmp1.dat
DIRDATE=`date_edit2 $CDATE 10800 `
DIRDATE=`echo $DIRDATE | cut -c1-10`
cp $OBSDIR/obs${DIRDATE}/t+3.dat $OBSDIR/tmp2.dat
 
cat $OBSDIR/tmp1.dat $OBSDIR/tmp2.dat > $OBSDIR/${CDATE}.dat

rm $OBSDIR/tmp1.dat $OBSDIR/tmp2.dat

fi

#Update the date
CDATE=`date_edit2 ${CDATE} ${OBSFREC}  `

done



























