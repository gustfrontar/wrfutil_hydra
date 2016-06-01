#!/bin/bash
# Runing on LINUX
SYSTEM=1 #We need this for some bash functions.
# Initial time
IDATE=20100801000000
# Final time
EDATE=20101101000000

#OBSNAME
OBSNAME=PREPBUFRSA

OBSDIR=$HOME/share/OBS/$OBSNAME
TMPDIR=$HOME/data/TMP/PREP_$OBSNAME

source $HOME/share/LETKF_WRF/wrf/run/util.sh #Load time functions.

DECODER=$HOME/share/LETKF_WRF/wrf/obs/prepbufr/decoder/

#
# WKDIR
#

safe_rm_tmpdir $TMPDIR
mkdir -p $TMPDIR
cd $TMPDIR

ln -sf $DECODER/decoder  .
ln -sf $DECODER/grabbufr .
#
# Settings for DATASRC
#

#
# Cycle run # MAIN LOOP #
#

CDATE=$IDATE

while [ $CDATE -le $EDATE ]
do
echo " >>>"
echo " >>> BEGIN DECODIFICATION OF $CDATE "
echo " >>>"
#
# GET NCEP PREPBUFR
#

#
# DECODE
#
#test -f prepbufr.tmp && rm -f prepbufr.tmp
AUXDATE=`echo $CDATE | cut -c1-10`
ln -sf $OBSDIR/prepbufr/${AUXDATE}.prepbufr.nr prepbufr.tmp
wc -c prepbufr.tmp | ./grabbufr prepbufr.tmp prepbufr.in  > info

time ./decoder >> info

#Correspondance between files and relative time.
#  fort.87 t-3.dat
#  fort.88 t-2.dat
#  fort.89 t-1.dat
#  fort.90 t.dat
#  fort.91 t+1.dat
#  fort.92 t+2.dat
#  fort.93 t+3.dat

touch fort.87
touch fort.88
touch fort.89
touch fort.90
touch fort.91
touch fort.92
touch fort.93
touch t-3bis.dat

#T-3 relative time.
AUXDATE=`date_edit2 $CDATE -10800 `
cat t-3bis.dat fort.87 > $OBSDIR/OBS${AUXDATE}.dat

#T-2 relative time
AUXDATE=`date_edit2 $CDATE -7200 `
mv fort.88 $OBSDIR/OBS${AUXDATE}.dat

#T-1 relative time
AUXDATE=`date_edit2 $CDATE -3600 `
mv fort.89 $OBSDIR/OBS${AUXDATE}.dat

#T0 relative time
AUXDATE=`date_edit2 $CDATE 0 `
mv fort.90 $OBSDIR/OBS${AUXDATE}.dat

#T+1 relative time
AUXDATE=`date_edit2 $CDATE 3600 `
mv fort.91 $OBSDIR/OBS${AUXDATE}.dat

#T+2 relative time
AUXDATE=`date_edit2 $CDATE 7200 `
mv fort.92 $OBSDIR/OBS${AUXDATE}.dat

#T+3 relative time (store it to merge it with t-3 of the next bufr file)
mv fort.93 t-3bis.dat

mv info $OBSDIR/info$CDATE

#
# Date change ### MAIN LOOP END ###
#

CDATE=`date_edit2 $CDATE 21600 `

done
echo "NORMAL END"
