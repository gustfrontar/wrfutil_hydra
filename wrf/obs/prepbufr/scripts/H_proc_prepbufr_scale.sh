#!/bin/bash
#Hourly prepbufr output 
# Runing on LINUX

SYSTEM=1 #We need this for some bash functions.


# Initial time
IDATE=20160518000000
# Final time
EDATE=20160610000000

MINLON=300
MAXLON=80
MINLAT=10
MAXLAT=85

OBSNAME=prepbufr_europe_scale

BUFRDIR=$HOME/share/DATA/OBS/prepbufr/
OBSDIR=$HOME/share/DATA/OBS/${OBSNAME}/

TMPDIR=$HOME/data/TMP/PREP_$OBSNAME

source $HOME/share/LETKF_WRF/wrf/run/util.sh #Load time functions.

DECODER=$HOME/share/LETKF_WRF/wrf/obs/prepbufr/decoder/

mkdir -p $OBSDIR/scale_1hr

mkdir -p $OBSDIR/scale_3hr

mkdir -p $OBSDIR/scale_6hr

#
# WKDIR
#

safe_rm_tmpdir $TMPDIR
mkdir -p $TMPDIR
cd $TMPDIR

ln -sf $DECODER/decoder_scale_*  .
ln -sf $DECODER/grabbufr .

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
# DECODE
#
#test -f prepbufr.tmp && rm -f prepbufr.tmp
AUXDATE=`echo $CDATE | cut -c1-10`
ln -sf $BUFRDIR/${AUXDATE}.prepbufr.nr prepbufr.tmp
wc -c prepbufr.tmp | ./grabbufr prepbufr.tmp prepbufr.in  > outinfo

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Run decoder for 1 hr analysis
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

echo "Running decoder for 1 hr analysis "

rm -f fort.*

time ./decoder_scale_1hr $MINLON $MAXLON $MINLAT $MAXLAT  >> outinfo

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

#T-3 relative time (merge with t+3 of previous time if exists)

AUXDATE=`date_edit2 $CDATE -10800 `

echo "Generating $OBSDIR/scale_1hr/obs_${AUXDATE}.dat"

touch $OBSDIR/scale_1hr/obs_${AUXDATE}.dat
mv $OBSDIR/scale_1hr/obs_${AUXDATE}.dat ./tmp.dat
cat tmp.dat fort.87 > $OBSDIR/scale_1hr/obs_${AUXDATE}.dat

#T-2 relative time
AUXDATE=`date_edit2 $CDATE -7200 `

echo "Generating $OBSDIR/scale_1hr/obs_${AUXDATE}.dat"

mv fort.88 $OBSDIR/scale_1hr/obs_${AUXDATE}.dat

#T-1 relative time
AUXDATE=`date_edit2 $CDATE -3600 `

echo "Generating $OBSDIR/scale_1hr/obs_${AUXDATE}.dat"

mv fort.89 $OBSDIR/scale_1hr/obs_${AUXDATE}.dat

#T0 relative time
AUXDATE=`date_edit2 $CDATE 0 `

echo "Generating $OBSDIR/scale_1hr/obs_${AUXDATE}.dat"

mv fort.90 $OBSDIR/scale_1hr/obs_${AUXDATE}.dat

#T+1 relative time
AUXDATE=`date_edit2 $CDATE 3600 `

echo "Generating $OBSDIR/scale_1hr/obs_${AUXDATE}.dat"

mv fort.91 $OBSDIR/scale_1hr/obs_${AUXDATE}.dat

#T+2 relative time
AUXDATE=`date_edit2 $CDATE 7200 `

echo "Generating $OBSDIR/scale_1hr/obs_${AUXDATE}.dat"

mv fort.92 $OBSDIR/scale_1hr/obs_${AUXDATE}.dat

#T+3 relative time
AUXDATE=`date_edit2 $CDATE 10800 `

echo "Generating $OBSDIR/scale_1hr/obs_${AUXDATE}.dat"

mv fort.93 $OBSDIR/scale_1hr/obs_${AUXDATE}.dat

cat info outinfo > $OBSDIR/scale_1hr/info$CDATE

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Run decoder for 3 hr analysis
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

echo "Running decoder for 3 hr analysis "

rm -f fort.*

time ./decoder_scale_3hr $MINLON $MAXLON $MINLAT $MAXLAT  >> outinfo

touch fort.87
touch fort.88
touch fort.89
touch fort.90
touch fort.91
touch fort.92
touch fort.93

AUXDATE=`date_edit2 $CDATE -10800 `

echo "Generating $OBSDIR/scale_1hr/obs_${AUXDATE}.dat"

touch $OBSDIR/scale_3hr/obs_${AUXDATE}.dat
mv $OBSDIR/scale_3hr/obs_${AUXDATE}.dat tmp.dat
cat tmp.dat fort.87 fort.88 > $OBSDIR/scale_3hr/obs_${AUXDATE}.dat

AUXDATE=`date_edit2 $CDATE 0 `

echo "Generating $OBSDIR/scale_1hr/obs_${AUXDATE}.dat"

cat fort.89 fort.90 fort.91 > $OBSDIR/scale_3hr/obs_${AUXDATE}.dat

AUXDATE=`date_edit2 $CDATE 10800 `

echo "Generating $OBSDIR/scale_1hr/obs_${AUXDATE}.dat"

cat fort.92 fort.93         > $OBSDIR/scale_3hr/obs_${AUXDATE}.dat

cat info outinfo > $OBSDIR/scale_3hr/info$CDATE


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Run decoder for 6 hr analysis
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

echo "Running decoder for 6 hr analysis "

rm -f fort.*

time ./decoder_scale_6hr $MINLON $MAXLON $MINLAT $MAXLAT  >> outinfo

touch fort.87
touch fort.88
touch fort.89
touch fort.90
touch fort.91
touch fort.92
touch fort.93

AUXDATE=`date_edit2 $CDATE 0 `

echo "Generating $OBSDIR/scale_1hr/obs_${AUXDATE}.dat"

cat fort.87 fort.88 fort.89 fort.90 fort.91 fort.92  > $OBSDIR/scale_6hr/obs_${AUXDATE}.dat


cat info outinfo > $OBSDIR/scale_6hr/info$CDATE

#
# Date change ### MAIN LOOP END ###
#

CDATE=`date_edit2 $CDATE 21600 `

done
echo "NORMAL END"
