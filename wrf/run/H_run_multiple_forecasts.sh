#PBS -l nodes=@@NODES@@:ppn=@@PPN@@
#PBS -S /bin/bash

#=======================================================================
# This script runs multiple data assimilation cycles.
# The driver script sets the TMP directory with all the required data
# and this scripts executes multiple data assimilation cycle and the output
# goes to the TMP directory.
#=======================================================================
#set -x
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
#Get root dir
CDIR=`pwd`

CONFIGURATION=@@CONFIGURATION@@

source @@CONFFILE@@
source @@MCONFFILE@@

source $TMPDIR/SCRIPTS/util.sh

RESTART=@@RESTART@@
RESTARTDATE=@@RESTARTDATE@@
RESTARTITER=10

PPSSERVER=`hostname`

#In this case we will use the same number of nodes for the forecast and for LETKF.
TOTAL_NODES_LETKF=$TOTAL_NODES_FORECAST

#Override value of GRIBDIR and PERTGRIBDIR
GRIBDIR=$TMPDIR/GRIB/
PERTGRIBDIR=$TMPDIR/PERTGRIB/
MYHOST=@@MYHOST@@
PID=@@PID@@
MYSCRIPT=@@MYSCRIPT@@

#Overrride the value of OBSDIR and RADAROBSDIR.
OBSDIR=$TMPDIR/OBS/
RADAROBSDIR=$TMPDIR/RADAROBS/

#Overrride the value of output dir.
OUTPUTDIR=$TMPDIR/output/

#Overrride analyisis source
ANALYSIS_SOURCE=$TMPDIR/input/


#Create a copy of the OUTPUTDIR structure in the TMPDIR.
safe_init_outputdir $OUTPUTDIR

#Set output log
set_my_log

{

echo ">>>> I'm RUNNING IN $MYHOST and my PID is $PID" 
echo ">>>> My config file is $CONFIGURATION         " 
echo ">>>> My domain is $DOMAINCONF                 " 
echo ">>>> My machine is $MCONFIGURATION            " 
echo ">>>> I' am $CDIR/$MYSCRIPT                    "


echo '>>>'                                           
echo ">>> SET METEM DATA FREQ "          
echo '>>>' 

set_pre_processing_intervals

echo '>>>'                                           
echo ">>> GET INITIAL RANDOM DATES "          
echo '>>>' 

get_random_dates 

##################################################
# START CYCLE IN TIME
##################################################

echo ">>>"                                                         
echo ">>> STARTING WRF-ENSEMBLE-FORECAST FROM $IDATE TO $EDATE"         
echo ">>>"                                                          
 
if [ $RESTART -eq 0 ] ; then
CDATE=$IDATE
ITER=1
else
CDATE=$RESTARTDATE
ITER=$RESTARTITER
fi

while test $CDATE -le $EDATE
do

echo '>>>'                                                           
echo ">>> BEGIN COMPUTATION OF $CDATE  ITERATION: $ITER"             
echo '>>>'                                                           

set_cycle_dates   

echo " >>"                                                           
echo " >> GET THE UNPERTURBED INITIAL AND BOUNDARY CONDITIONS"                                  
echo " >>" 

#RUN MODEL PRE-PROCESSING FROM GLOBAL ANALYSIS OR FORECASTS (run in PPS)

get_met_em_from_grib_noqueue

echo " >>"                                                           
echo " >> GENERATING PERTURBATIONS"                                  
echo " >>"                                                           

#PERTURB MET_EM FILES USING RANDOM BALANCED OR RANDOM SMOOTHED PERTURBATIONS (run in PPS)

perturb_met_em_from_grib_noqueue


echo " >>"                                                           
echo " >> ENSEMBLE FORECASTS"
echo " >>"

#CREATE OUTPUT DIRECTORIES FOR THE CURRENT FORECAST (YYYYMMDDHHNNSS)
RESULTDIRG=$OUTPUTDIR/forecast/$CDATE/
RESULTDIRA=$OUTPUTDIR/anal/$ADATE/


mkdir -p $RESULTDIRG
mkdir -p $RESULTDIRA

echo " >>"
echo " >> RUNNING THE ENSEMBLE"
echo " >>"

run_ensemble_forecast_noqueue

echo " >>"
echo " >> DOING POST PROCESSING"
echo " >>"

arw_postproc_forecast

echo " >>"
echo " >> CHECKING CYCLE"
echo " >>"

check_postproc


CDATE=$ADATE
ITER=`expr $ITER + 1`
RESTART=0    #TURN RESTART FLAG TO 0 IN CASE IT WAS GREATHER THAN 0.

done	### MAIN LOOP ###

echo "NORMAL END"

} > $my_log 2>&1
