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
CONFIGURATION=relampago_nature_run   #Define a experiment configuration
MCONFIGURATION=machine_relampago_nature_run    #Define a machine configuration (number of nodes, etc)

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

#Load functions for LETKF-WRF cycle.
source $UTIL

echo '>>>'
echo ">>> INITIALIZING WORK DIRECTORY AND OUTPUT DIRECTORY"
echo '>>>'

safe_init_outputdir $OUTPUTDIR

#Start of the section that will be output to my log.
#{

safe_init_tmpdir $TMPDIR

save_configuration $CDIR/$MYSCRIPT

echo '>>>'                                           
echo ">>> COPYING DATA TO WORK DIRECTORY "          
echo '>>>'  

copy_data

copy_data_multiplecycles

#Generating the domain requires acces to GEOG database.
echo '>>>'                                           
echo ">>> GENERATING DOMAIN "          
echo '>>>' 

get_domain


edit_multiplecycle $TMPDIR/SCRIPTS/H_run_multiple_forecasts.sh

#Run multiple cycles with only one QSUB
sub_and_wait $TMPDIR/SCRIPTS/H_run_multiple_forecasts.sh

#Move experiment data to OUTPUTDIR
mv $TMPDIR/output/* $OUTPUTDIR

echo "NORMAL END"


