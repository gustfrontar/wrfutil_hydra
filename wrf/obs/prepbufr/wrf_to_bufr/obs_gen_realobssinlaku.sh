#!/bin/bash
#=======================================================================
# letkf_cycle.sh
# To generate input and boundary files as a prior step to run LETKF cycle
# in the K computer. 
#=======================================================================
#set -x
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
### parallel computation setting
### LETKF setting
WINDOW=360            # assimilation window length (min)
GUESFT=540            # first guess forecast length (min)
INITIME=180           # first observation time.
TRUE_RUN=TRUE_RUN_SINLAKU_CONSTANT_QFX0.8_60KM
TRUE_PATH=/home/jruiz/salidas/EXPERIMENTS/${TRUE_RUN}/  # name of experiment
OBSTYPE=3             # Types of observing network.(1 random, 2 regular, 3 input from file)
OBSDENSITY=0.005      # OBSDENSITY (for random network only)
SKIP=1                # SKIP in X and Y (for regular networks)
SKIPZ=3               # SKIP in Z for all types of networks except 3
OBSINPUT=/home/jruiz/DATA/OBS/ucar_airs_th3x3_corrected/  #If OBSTYPE == 3 provide the path from where to read the file.

EXPNAME=OBS_$TRUE_RUN
### initial date setting
IY=2008
IM=08
ID=07
IH=00
IMN=00
### final date setting
EY=2008
EM=09
ED=29
EH=00
EMN=00

### directory settings
CDIR=`pwd`
cd ../..
WRF=`pwd`
NCIO=$WRF/ncio/
TMPDIR=/home/jruiz/TMP/TMP_$EXPNAME # work directory
OUTPUT=/home/jruiz/DATA/OBS/

ulimit -s unlimited
#-----------------------------------------------------------------------
# Usually do not modify below
#-----------------------------------------------------------------------
echo '>>>'
echo ">>> STARTING OBS  GENERATION CYCLE FROM $IY/$IM/$ID/$IH TO $EY/$EM/$ED/$EH"
echo '>>>'

cd $CDIR
source util.sh
#
# Work directory
#
echo "Setting up work directories.."
echo " >> removing old work directories.."
rm -rf $TMPDIR
mkdir -p $TMPDIR

cp $CDIR/obs_gen.exe  $TMPDIR/obs_gen.exe

# 
# Cycle run ### MAIN LOOP ###
#
while test $IY$IM$ID$IH -le $EY$EM$ED$EH
do

echo '>>>'
echo ">>> GENERATING OBS AROUND $IY/$IM/$ID/$IH"
echo '>>>'

date_edit $IY $IM $ID $IH 00 $WINDOW > anal_time.txt
read AY AM AD AH AMN < anal_time.txt
rm -f anal_time.txt

mkdir -p $OUTPUT/$EXPNAME/obs${AY}${AM}${AD}${AH}

cd $TMPDIR

#CREATE NAMELIST
echo "$OBSDENSITY (OBSDENSITY) #1 would be observations at each grid point."                                 > parameters.in
echo "$OBSTYPE    (OBSTYPE)    #1 random, 2 regular, 3 input file"                                          >> parameters.in
echo "$SKIP       (SKIP)       #Regular observation distribution only. How many points to skip in X and Y." >> parameters.in
echo "$SKIPZ      (SKIPZ)      #SKIP for vertical distribution of observations. (OBSTYPE 1 and 2 only)"     >> parameters.in

DT=$INITIME
date_edit $IY $IM $ID $IH 00 $DT > tmpfile
read OY OM OD OH OMN < tmpfile
rm -f tmpfile

ln -sf $TRUE_PATH/wrfout_d01_${OY}-${OM}-${OD}_${OH}:${OMN}:00 gs01001
ln -sf $OUTPUT/$EXPNAME/obs${AY}${AM}${AD}${AH}/t-3.dat out.dat
if  [ $OBSTYPE -eq 3  ]
then
ln -sf $OBSINPUT/obs${AY}${AM}${AD}${AH}/t-3.dat $TMPDIR/obs.dat
fi
#Run the program.
./obs_gen.exe


DT=`expr $INITIME + 60 `
date_edit $IY $IM $ID $IH 00 $DT > tmpfile
read OY OM OD OH OMN < tmpfile
rm -f tmpfile

ln -sf $TRUE_PATH/wrfout_d01_${OY}-${OM}-${OD}_${OH}:${OMN}:00 gs01001
ln -sf $OUTPUT/$EXPNAME/obs${AY}${AM}${AD}${AH}/t-2.dat out.dat
if  [ $OBSTYPE -eq 3  ]
then
ln -sf $OBSINPUT/obs${AY}${AM}${AD}${AH}/t-2.dat $TMPDIR/obs.dat
fi
#Run the program.
./obs_gen.exe

DT=`expr $INITIME + 120 `
date_edit $IY $IM $ID $IH 00 $DT > tmpfile
read OY OM OD OH OMN < tmpfile
rm -f tmpfile

ln -sf $TRUE_PATH/wrfout_d01_${OY}-${OM}-${OD}_${OH}:${OMN}:00 gs01001
ln -sf $OUTPUT/$EXPNAME/obs${AY}${AM}${AD}${AH}/t-1.dat out.dat
if  [ $OBSTYPE -eq 3  ]
then
ln -sf $OBSINPUT/obs${AY}${AM}${AD}${AH}/t-1.dat $TMPDIR/obs.dat
fi
#Run the program.
./obs_gen.exe

DT=`expr $INITIME + 180 `
date_edit $IY $IM $ID $IH 00 $DT > tmpfile
read OY OM OD OH OMN < tmpfile
rm -f tmpfile

ln -sf $TRUE_PATH/wrfout_d01_${OY}-${OM}-${OD}_${OH}:${OMN}:00 gs01001
ln -sf $OUTPUT/$EXPNAME/obs${AY}${AM}${AD}${AH}/t.dat out.dat
if  [ $OBSTYPE -eq 3  ]
then
ln -sf $OBSINPUT/obs${AY}${AM}${AD}${AH}/t.dat $TMPDIR/obs.dat
fi
#Run the program.
./obs_gen.exe

DT=`expr $INITIME + 240 `
date_edit $IY $IM $ID $IH 00 $DT > tmpfile
read OY OM OD OH OMN < tmpfile
rm -f tmpfile

ln -sf $TRUE_PATH/wrfout_d01_${OY}-${OM}-${OD}_${OH}:${OMN}:00 gs01001
ln -sf $OUTPUT/$EXPNAME/obs${AY}${AM}${AD}${AH}/t+1.dat out.dat
if  [ $OBSTYPE -eq 3  ]
then
ln -sf $OBSINPUT/obs${AY}${AM}${AD}${AH}/t+1.dat $TMPDIR/obs.dat
fi
#Run the program.
./obs_gen.exe

DT=`expr $INITIME + 300 `
date_edit $IY $IM $ID $IH 00 $DT > tmpfile
read OY OM OD OH OMN < tmpfile
rm -f tmpfile

ln -sf $TRUE_PATH/wrfout_d01_${OY}-${OM}-${OD}_${OH}:${OMN}:00 gs01001
ln -sf $OUTPUT/$EXPNAME/obs${AY}${AM}${AD}${AH}/t+2.dat out.dat
if  [ $OBSTYPE -eq 3  ]
then
ln -sf $OBSINPUT/obs${AY}${AM}${AD}${AH}/t+2.dat $TMPDIR/obs.dat
fi
#Run the program.
./obs_gen.exe

DT=`expr $INITIME + 360 `
date_edit $IY $IM $ID $IH 00 $DT > tmpfile
read OY OM OD OH OMN < tmpfile
rm -f tmpfile

ln -sf $TRUE_PATH/wrfout_d01_${OY}-${OM}-${OD}_${OH}:${OMN}:00 gs01001
ln -sf $OUTPUT/$EXPNAME/obs${AY}${AM}${AD}${AH}/t+3.dat out.dat
if  [ $OBSTYPE -eq 3  ]
then
ln -sf $OBSINPUT/obs${AY}${AM}${AD}${AH}/t+3.dat $TMPDIR/obs.dat
fi
#Run the program.
./obs_gen.exe


#
IY=$AY
IM=$AM
ID=$AD
IH=$AH
IMN=$AMN
ITER=`expr $ITER + 1`

done
echo "NORMAL END"
