#!/bin/bash

#Pourpose: Run the scale model in Hakushu

CONFIGURATION="breeding_osaka_pawr_1km_bip_global_1000m_UVT"

RESTART=1
RESTARTDATE=20130713053130
RESTARTITER=10

SRCDIR=`pwd`

source $SRCDIR/../configuration/${CONFIGURATION}/config.sh

source $SRCDIR/util.sh 

#Create tmpdir            #######################################

safe_init_tmpdir $TMPDIR

set_tmpdir_breeding $TMPDIR

#Create outputdir         #######################################

if [ $RESTART -eq 0 ] ; then

safe_init_outputdir $OUTPUTDIR

else

mkdir $TMPDIR/out/$RESTARTDATE
cp -r $OUTPUTDIR/$RESTARTDATE/* $TMPDIR/out/$RESTARTDATE

fi

##################################################
# START CYCLE IN TIME
##################################################

echo '>>>'
echo ">>> STARTING SCALE BREEDING CYCLE FROM $IDATE TO $EDATE"
echo '>>>'

if [ $RESTART -eq 0 ] ; then
CDATE=$INIDATE
ITER=1
else
CDATE=$RESTARTDATE
ITER=$RESTARTITER
fi


while [ $CDATE -le $ENDDATE ] ; do

#DEFINE IMORTANT DATES FOR THE CURRENT CYCLE (YYYYMMDDHHNNSS)M
FDATE=`date_edit2 $CDATE $BVFREQ `           #Forecast end date

if [ $ITER -eq 1 ]; then

 echo "Initialize bred vector "
 init_bv $CDATE

fi

evolve_rescale_bv $CDATE


CDATE=`date_edit2 $CDATE $BVFREQ `
ITER=`expr $ITER + 1`
RESTART=0    #TURN RESTART FLAG TO 0 IN CASE IT WAS GREATHER THAN 0.

done    ### MAIN LOOP ###

mv $TMPDIR/out/* $OUTPUTDIR/
mv $TMPDIR/out/logs/* $OUTPUTDIR/logs/*
mv $TMPDIR/out/$RESTARTDATE/* $OUTPUTDIR/$RESTARTDATE/


echo "NORMAL END"
