#!/bin/bash
#=======================================================================
# Main driver for nature run (i.e. free WRF run) experiments.
#   To run a WRF ensemble.
#=======================================================================
#set -x
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
#Get root dir
CDIR=`pwd`

#CONFIGURATION
CONFIGURATION=gsi_boundary_prep             #Define a experiment configuration
MCONFIGURATION=machine_scale_boundary_prep    #Define a machine configuration (number of nodes, etc)

RESTART=0
RESTARTDATE=20080823060000
RESTARTITER=10

MYHOST=`hostname`
PID=$$
MYSCRIPT=${0##*/}

if [ -e $CDIR/configuration/analysis_conf/${CONFIGURATION}.sh ];then
source $CDIR/configuration/analysis_conf/${CONFIGURATION}.sh
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

copy_data_multiplecycles 

PERTGRIBDIR=$TMPDIR/PERTGRIB/

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

echo '>>>'                                                           
echo ">>> SET SOME MACHINE DEPENDENT PARAMETERS"             
echo '>>>' 

get_node_list

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

#echo " >>"
#echo " >> GENERATE WRINPUT FROM PERTURBED MET_EM"
#echo " >>"

#CREATE OUTPUT DIRECTORIES FOR THE CURRENT CYCLE (YYYYMMDDHHNNSS)

#get_wrfinput_from_met_em


CDATE=$ADATE
ITER=`expr $ITER + 1`
RESTART=0    #TURN RESTART FLAG TO 0 IN CASE IT WAS GREATHER THAN 0.


done	### MAIN LOOP ###

echo "NORMAL END"  
}  > $my_log 2>&1
