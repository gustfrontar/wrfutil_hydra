#!/bin/bash
#=======================================================================
# K_driver.sh
#   To run the WRF-LETKF in the K computer.
#   This scripts prepares all the data needed to run the letkf cycle 
#   then submmits the job to the K computer.
#   Based on the script developed by Shigenori Otsuka
#   Diferences from V3:
#   -Boundary perturbations, met em are perturbed and real is run in 
#    K computer nodes.
#   -Letkf incorporates relaxation to prior inflation and eigen exa
#    matrix computations.
#   -LETKF uses NETCDF IO.
#=======================================================================
#set -x
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
#Get root dir
CDIR=`pwd`

#TODO(some of these variables should be included in a configuration file)
ITER_BEGIN_STATS=40  #How many cycles before start the computation of accumulated rmse and bias.
NAREA=3              #Number of verification regions where statistics will be averaged
VLON1="0,110,110,"   #Limits for the verifcation regions.
VLON2="360,160,170"
VLAT1="-90,5,25,"
VLAT2="90,25,50,"
REGRID_OUTPUT=".true."        #Wheter regrid ouput will be generated or not
REGRID_RES="1.0d0"            #Resolution of regrided output.
REGRID_VERT_RES="10000.0d0"     #Vertical resolution of regrided output (Pa)
REGRID_VERT_ZRES="1000.0d0"   #Vertical resolution of regrided output (m)
FILTER_INPUT=".false."        #Not coded yet.
TMPDIRVERIF=$HOME/data/TMP/VERIFICATION_512M60K
DO_OBSGRID=".true."
SLOTSTEP="0"                  #Slot step for obsgrid verifcation
SLOTOFFSET="1"                #Slot offset for obsgrid verification
OBSOPENAMELIST=control
LETKFNAMELIST=control #[Not used]
GROSSERROROBSOP="10.0d0"
GLOBALANALYSIS_DATA_VERIFICATION_FREQ=21600

#CONFIGURATION
DOMAINCONF=SINLAKU_60K               #Define a domain
CONFIGURATION=control512m            #Define a experiment configuration
MCONFIGURATION=machine_60k_K         #Define a machine configuration [System type]

VINIDATE=
VENDDATE=


MYHOST=`hostname`
PID=$$
MYSCRIPT=${0##*/}
echo ">>>> I'm RUNNING IN $MYHOST and my PID is $PID"
echo ">>>> My config file is $CONFIGURATION         "
echo ">>>> My domain is $DOMAINCONF                 "
echo ">>>> I' am $CDIR/$MYSCRIPT                    "

if [ -e $CDIR/../configuration/analysis_conf/${CONFIGURATION}.sh ];then
source $CDIR/../configuration/analysis_conf/${CONFIGURATION}.sh
else
echo "CAN'T FIND CONFIGURATION FILE $CDIR/../configuration/analysis_conf/${CONFIGURATION}.sh"
exit
fi

if [ -e $CDIR/../configuration/machine_conf/${MCONFIGURATION}.sh ];then
source $CDIR/../configuration/machine_conf/${MCONFIGURATION}.sh
else
echo "CAN'T FIND MACHINE CONFIGURATION FILE $CDIR/../configuration/machine_conf/${MCONFIGURATION}.sh"
exit
fi

NSLOTS=1 #Override the number of time slots for the observation operatior.
NBSLOT=1
GROSS_ERROR=$GROSSERROROBSOP #Override obserror setting for the verification.

echo ">>>> My LETKFNAMELIST is $NAMELISTLETKF       "

source $UTIL #Load all the functions that will be used in this script.

echo '>>>'
echo ">>> INITIALIZING WORK DIRECTORY"
echo '>>>'

TMPDIR=$TMPDIRVERIF

echo $TMPDIR

safe_init_tmpdir $TMPDIR

echo '>>>'
echo ">>> COPYING DATA TO WORK DIRECTORY "
echo '>>>'

copy_data

##################################################
# START CYCLE IN TIME
##################################################
#If VINIDATE and VENDATE are set then use them to control the 
#start and end of the verification period.
if [ -n "$VINIDATE" ] ; then
  IDATE=$VINIDATE
  echo "[Warning]: Original configuration initial date overrid by VINIDATE value: $VINIDATE"
fi
if [ -n "$VENDDATE" ] ; then
  EDATE=$VENDDATE
  echo "[Warning]: Original configuration initial date overrid by VENDDATE value: $VENDDATE"
fi

echo '>>>'
echo ">>> STARTING ANALYSIS VERIFICATION FROM $IDATE TO $EDATE"
echo '>>>'

CDATE=$IDATE
ITER=1

while test $CDATE -le $EDATE
do

echo '>>>'
echo ">>> BEGIN COMPUTATION OF $CDATE  ITERATION: $ITER"
echo '>>>'

set_cycle_dates

#CREATE OUTPUT DIRECTORIES FOR THE CURRENT CYCLE (YYYYMMDDHHNNSS)
RESULTDIRA=$OUTPUTDIR/anal/$CDATE/
RESULTDIRG=$OUTPUTDIR/gues/$CDATE/

echo '>>>'
echo ">>> PROCESSING GLOBAL ANALYSIS DATA "
echo '>>>'

time postproc_global_analysis $CDATE $CDATE

echo '>>>'
echo ">>> ANALYSIS VERIFICATION USING GLOBAL ANALYSIS "
echo '>>>'

time verification_against_global_analysis $CDATE $CDATE

echo " >>"
echo " >> ANALYSIS VERIFICATION USING OBSERVATIONS "
echo " >>"

time verification_against_obs $CDATE $CDATE


CDATE=$ADATE
ITER=`expr $ITER + 1`

done	### MAIN LOOP ###

echo "NORMAL END"
