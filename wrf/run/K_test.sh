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


#CONFIGURATION
DOMAINCONF=OSAKA_1KM_NESTED            #Define a domain
CONFIGURATION=osaka_nested             #Define a experiment configuration
MCONFIGURATION=machine_osaka_nested    #Define a machine configuration (number of nodes, etc)
LETKFNAMELIST=control                  #Define a letkf namelist template

RESTART=0
RESTARTDATE=20080904120000
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

#safe_init_outputdir $OUTPUTDIR

set_my_log

{

echo ">>>> I'm RUNNING IN $MYHOST and my PID is $PID"
echo ">>>> My config file is $CONFIGURATION         "
echo ">>>> My domain is $DOMAINCONF                 "
echo ">>>> My machine is $MCONFIGURATION            "
echo ">>>> I' am $CDIR/$MYSCRIPT                    "
echo ">>>> My LETKFNAMELIST is $LETKFNAMELIST       "

safe_init_tmpdir $TMPDIR

echo '>>>'
echo ">>> COPYING DATA TO WORK DIRECTORY "
echo '>>>'

copy_data

##################################################
# START CYCLE IN TIME
##################################################

echo '>>>'
echo ">>> STARTING WRF-LETKF CYCLE FROM $IDATE TO $EDATE"
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
echo ">>> BEGIN COMPUTATION OF $CDATE  ITERATION: $ITER"
echo '>>>'

set_cycle_dates

echo " >>"
echo " >> GENERATING PERTURBATIONS"
echo " >>"

#PERTURB MET_EM FILES USING RANDOM BALANCED OR RANDOM SMOOTHED PERTURBATIONS (run in PPS)

METEMDIR=$METEMDIRFOR/$CDATE/
run_script=$TMPDIR/SCRIPTS/perturb_met_em.sh
#perturb_met_em $run_script
run_script=$TMPDIR/SCRIPTS/wrf_to_met_em.sh
#wrf_to_met_em $run_script   #Only works if interp ana == 1 (i.e. letkf analysis grid is different from forecast grid)


echo " >>"
echo " >> ENSEMBLE FORECASTS AND LETKF"
echo " >>"

#CREATE OUTPUT DIRECTORIES FOR THE CURRENT CYCLE (YYYYMMDDHHNNSS)
RESULTDIRG=$OUTPUTDIR/forecast/$CDATE/

mkdir -p $RESULTDIRG

#PREPARE THE SCRIPT
cd $TMPDIR/SCRIPTS/

run_forecast_script=$TMPDIR/SCRIPTS/run_forecast_script.sh
generate_run_forecast_script_k $run_forecast_script

#sub_and_wait $run_forecast_script $ITER

#check_forecast

echo " >>"
echo " >> DOING POST PROCESSING"
echo " >>"

#arw_postproc_forecast 

echo " >>"
echo " >> CHECKING CYCLE"
echo " >>"

#check_postproc


#echo " >>" (moving it to an independent verifcation script)
#echo " >> DOING VERIFICATION WITH OBS"
#echo " >>"
#
#
#verfication_against_obs


CDATE=$ADATE
ITER=`expr $ITER + 1`
RESTART=0    #TURN RESTART FLAG TO 0 IN CASE IT WAS GREATHER THAN 0.

done	### MAIN LOOP ###

echo "NORMAL END"  
}  > $my_log 2>&1
