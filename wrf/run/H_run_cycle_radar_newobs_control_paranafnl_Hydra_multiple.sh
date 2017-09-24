#!/bin/bash
#=======================================================================
# Driver Script for multiple cycle qsub script.
#   To run the WRF-LETKF in a qsub cluster.
#   This scripts prepares all the data needed to run the letkf cycle 
#   then submmits the job to the torque queu
#   New developments:
#   -WPS is run in pre/post nodes and real and wrf are run in computation nodes.
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
CONFIGURATION=control_paranafnl_newobs_60m_radar_grib_Hydra    #Define a experiment configuration
MCONFIGURATION=machine_radar60m_multiple_Hydra          #Define a machine configuration (number of nodes, etc)

RESTART=0
RESTARTDATE=20091117220000
RESTARTITER=10

MYHOST=`hostname`
PID=$$
MYSCRIPT=${0##*/}

CONFFILE=$CDIR/configuration/analysis_conf/${CONFIGURATION}.sh   
MCONFFILE=$CDIR/configuration/machine_conf/${MCONFIGURATION}.sh

if [ -e $CONFFILE ];then
source $CONFFILE
else
echo "CAN'T FIND CONFIGURATION FILE $CONFFILE "
exit
fi

if [ -e $MCONFFILE ];then
source $MCONFFILE
else
echo "CAN'T FIND MACHINE CONFIGURATION FILE $MCONFFILE "
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

edit_multiplecycle $TMPDIR/SCRIPTS/H_run_multiple_cycles.sh

#Run multiple cycles with only one QSUB
sub_and_wait $TMPDIR/SCRIPTS/H_run_multiple_cycles.sh  

#Move experiment data to OUTPUTDIR
mv $TMPDIR/output/* $OUTPUTDIR

echo "NORMAL END"

#} > $my_log 2>&1
