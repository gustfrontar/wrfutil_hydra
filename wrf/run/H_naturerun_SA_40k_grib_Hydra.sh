#!/bin/bash
#=======================================================================
# Main driver for nature run (i.e. free WRF run) experiments.
#   To run the WRF-LETKF in a cluster with qsub.
#=======================================================================
#set -x
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
#Get root dir
CDIR=`pwd`


#CONFIGURATION
CONFIGURATION=control_run_SA40k        #Define a experiment configuration
MCONFIGURATION=machine_nature_run      #Define a machine configuration (number of nodes, etc)

RESTART=0
RESTARTDATE=20080823060000
RESTARTITER=10

MYHOST=`hostname`
PID=$$
MYSCRIPT=${0##*/}

if [ -e $CDIR/configuration/forecast_conf/${CONFIGURATION}.sh ];then
source $CDIR/configuration/forecast_conf/${CONFIGURATION}.sh
else
echo "CAN'T FIND CONFIGURATION FILE $CDIR/configuration/forecast_conf/${CONFIGURATION}.sh"
exit
fi


if [ -e $CDIR/configuration/machine_conf/${MCONFIGURATION}.sh ];then
source $CDIR/configuration/machine_conf/${MCONFIGURATION}.sh
else
echo "CAN'T FIND MACHINE CONFIGURATION FILE $CDIR/configuration/machine_conf/${MCONFIGURATION}.sh"
exit
fi

source $UTIL

echo ">>>"
echo ">>> INITIALIZING WORK DIRECTORY AND OUTPUT DIRECTORY"
echo ">>>"

safe_init_outputdir $OUTPUTDIR

set_my_log

{

echo ">>>> I'm RUNNING IN $MYHOST and my PID is $PID"
echo ">>>> My config file is $CONFIGURATION         "
echo ">>>> My domain is $DOMAINCONF                 "
echo ">>>> My machine is $MCONFIGURATION            "
echo ">>>> I' am $CDIR/$MYSCRIPT                    "
echo ">>>> My LETKFNAMELIST is $LETKFNAMELIST       "

echo '>>>'                                           
echo ">>> INITIALIZING TMPDIR "          
echo '>>>'

safe_init_tmpdir $TMPDIR

echo '>>>'                                           
echo ">>> SAVING CONFIGURATION "          
echo '>>>'

save_configuration $CDIR/$MYSCRIPT

echo '>>>'                                           
echo ">>> COPYING DATA TO WORK DIRECTORY "          
echo '>>>'                                         

copy_data

echo '>>>'                                           
echo ">>> GENERATING DOMAIN "          
echo '>>>' 

get_domain

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

echo '>>>'
echo ">>> STARTING WRF SIMULATIONS FROM $IDATE TO $EDATE"
echo '>>>'

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
echo ">>> BEGIN SIMULATION STARTING ON $CDATE  ( ITERATION: $ITER )"             
echo '>>>'                                                           

set_cycle_dates

echo " >>"                                                           
echo " >> GET THE UNPERTURBED INITIAL AND BOUNDARY CONDITIONS"                                  
echo " >>" 

#RUN MODEL PRE-PROCESSING FROM GLOBAL ANALYSIS OR FORECASTS (run in PPS)

get_met_em_from_grib

echo " >>"                                                           
echo " >> GENERATING PERTURBATIONS"                                  
echo " >>"                                                           

#PERTURB MET_EM FILES USING RANDOM BALANCED OR RANDOM SMOOTHED PERTURBATIONS (run in PPS)

perturb_met_em_from_grib

echo " >>"
echo " >> ENSEMBLE FORECASTS"
echo " >>"

#CREATE OUTPUT DIRECTORIES FOR THE CURRENT CYCLE (YYYYMMDDHHNNSS)
RESULTDIRG=$OUTPUTDIR/forecast/$CDATE/
mkdir -p $RESULTDIRG

run_ensemble_forecast 

echo " >>"
echo " >> DOING POST PROCESSING"
echo " >>"

arw_postproc_forecast 

check_postproc


CDATE=$ADATE
ITER=`expr $ITER + 1`
RESTART=0    #TURN RESTART FLAG TO 0 IN CASE IT WAS GREATHER THAN 0.


done	### MAIN LOOP ###

echo "NORMAL END"  
}  > $my_log 2>&1
