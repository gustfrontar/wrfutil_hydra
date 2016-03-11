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
#CONFIGURATION
DOMAINCONF=CORDOBA_2K               #Define a domain
CONFIGURATION=control60m_radar      #Define a experiment configuration
MCONFIGURATION=machine_radar60m_K   #Define a machine configuration (number of nodes, etc)
LETKFNAMELIST=control               #Define a letkf namelist template

RESTART=0
RESTARTDATE=20080810000000
RESTARTITER=10


MYHOST=`hostname`
PID=$$
MYSCRIPT=${0##*/}

if [ -e $CDIR/configuration/analysis_conf/${CONFIGURATION}.sh ];then
source $CDIR/configuration/analysis_conf/${CONFIGURATION}.sh
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
echo ">>>> My LETKFNAMELIST is $LETKFNAMELIST       " 

echo '>>>'                                           
echo ">>> COPYING DATA TO WORK DIRECTORY "          
echo '>>>'                                         

copy_data   


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
echo " >> GENERATING PERTURBATIONS"                                  
echo " >>"                                                           

#PERTURB MET_EM FILES USING RANDOM BALANCED OR RANDOM SMOOTHED PERTURBATIONS (run in PPS)
run_script=$TMPDIR/SCRIPTS/perturb_met_em.sh                         
perturb_met_em $run_script                                           


echo " >>"                                                           
echo " >> ENSEMBLE FORECASTS AND LETKF"
echo " >>"

#CREATE OUTPUT DIRECTORIES FOR THE CURRENT CYCLE (YYYYMMDDHHNNSS)
RESULTDIRG=$OUTPUTDIR/gues/$ADATE/
RESULTDIRA=$OUTPUTDIR/anal/$ADATE/

mkdir -p $RESULTDIRG
mkdir -p $RESULTDIRA

#GET THE OBSERVATIONS
get_conventional_observations 

get_radar_observations_oldformat

cd $TMPDIR
#SUBMITT AND WAIT FOR THE JOB
echo " >>"
echo " >> WAITING FOR ENSEMBLE RUN"
echo " >>"

run_ensemble_forecast

echo " >>"
echo " >> WAITING FOR LETKF"
echo " >>"

run_letkf_script=$TMPDIR/SCRIPTS/run_letkf_script.sh
generate_run_letkf_script_k $run_letkf_script

MEM=`ens_member $MM `
MEMEAN=`ens_member $MEANMEMBER `
cp $TMPDIR/CURRENT_LETKF/gs${NBSLOT}${MEM} $TMPDIR/CURRENT_LETKF/gs${NBSLOT}${MEMEAN}
cp $TMPDIR/CURRENT_LETKF/gs${NBSLOT}${MEM} ${RESULTDIRG}/gues${MEMEAN}
ln -sf ${RESULTDIRG}/gues${MEMEAN} $TMPDIR/CURRENT_LETKF/gues${MEMEAN}

sub_and_wait $run_letkf_script $ITER

check_analysis

echo " >>"
echo " >> COPYING THE ENSEMBLE GUES "
echo " >>"

M=1
while [ $M -le $MEANMEMBER ] ; do

  RUNNING=0
  while [ $RUNNING -le $MAX_RUNNING -a $M -le $MEANMEMBER ] ; do
    MEM=`ens_member $M `
    cp ${TMPDIR}/CURRENT_LETKF/gs${NBSLOT}${MEM}   ${RESULTDIRG}/gues${MEM} &
    RUNNING=`expr $RUNNING + 1 `
    M=`expr $M + 1 `
  done
  time wait
done

echo " >>"
echo " >> DOING POST PROCESSING"
echo " >>"

arw_postproc 

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
