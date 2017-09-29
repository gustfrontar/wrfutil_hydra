#!/bin/bash
#=======================================================================
# Driver Script
#   To run the WRF-LETKF in a qsub cluster.
#   This scripts prepares all the data needed to run the letkf cycle 
#   then submmits the job to the torque queu
#   New developments:
#   -WPS is run in pre/post nodes and real and wrf are run in computation nodes.
#   -Letkf incorporates relaxation to prior inflation and eigen exa
#    matrix computations. 
#=======================================================================
#set -x
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
#Get root dir
CDIR=`pwd`

#CONFIGURATION
CONFIGURATION=ensemble_forecast_exprtps09_60m_radar_grib_Hakushu #Define a experiment configuration
MCONFIGURATION=machine_radar60m_Hakushu_multiple                 #Define a machine configuration (number of nodes, etc)

RESTART=0
RESTARTDATE=20140122175000
RESTARTITER=10

MYHOST=`hostname`
PID=$$
MYSCRIPT=${0##*/}

if [ -e $CDIR/configuration/forecast_conf/${CONFIGURATION}.sh ];then
source $CDIR/configuration/forecast_conf/${CONFIGURATION}.sh
else
echo "CAN'T FIND CONFIGURATION FILE $CDIR/configuration/analysis_conf/${CONFIGURATION}.sh"
exit
fi

if [ -e $CDIR/configuration/machine_conf/${MCONFIGURATION}.sh ];then
source $CDIR/configuration/machine_conf/${MCONFIGURATION}.sh
else
echo "CAN'T FIND MACHINE CONFIGURATION FILE $CDIR/configuration/machine_conf/${MCONFIGURATION}.sh"
exit
fi


source $UTIL

echo '>>>'
echo ">>> INITIALIZING WORK DIRECTORY AND OUTPUT DIRECTORY"
echo '>>>'

safe_init_outputdir $OUTPUTDIR

set_my_log

echo '>>>'
echo ">>> OUTPUT WILL BE REDIRECTED TO THE FOLLOWING FILE"
echo '>>>'

echo $my_log

#Start of the section that will be output to my log.
{
safe_init_tmpdir $TMPDIR

save_configuration $CDIR/$MYSCRIPT

echo ">>>> I'm RUNNING IN $MYHOST and my PID is $PID" 
echo ">>>> My config file is $CONFIGURATION         " 
echo ">>>> My domain is $DOMAINCONF                 " 
echo ">>>> My machine is $MCONFIGURATION            " 
echo ">>>> I' am $CDIR/$MYSCRIPT                    "

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

echo ">>>"                                                         
echo ">>> STARTING WRF-LETKF CYCLE FROM $IDATE TO $EDATE"         
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

get_met_em_from_grib

echo " >>"                                                           
echo " >> GENERATING PERTURBATIONS"                                  
echo " >>"                                                           

#PERTURB MET_EM FILES USING RANDOM BALANCED OR RANDOM SMOOTHED PERTURBATIONS (run in PPS)

perturb_met_em_from_grib

echo " >>"                                                           
echo " >> ENSEMBLE FORECASTS AND LETKF"
echo " >>"

#CREATE OUTPUT DIRECTORIES FOR THE CURRENT CYCLE (YYYYMMDDHHNNSS)
RESULTDIRG=$OUTPUTDIR/forecast/$CDATE/
RESULTDIRA=$OUTPUTDIR/anal/$ADATE/

mkdir -p $RESULTDIRG
mkdir -p $RESULTDIRA

echo " >>"
echo " >> RUNNING THE ENSEMBLE"
echo " >>"

run_ensemble_forecast

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
