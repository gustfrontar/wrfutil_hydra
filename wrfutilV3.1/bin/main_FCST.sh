#!/bin/bash 
#The main pourpose of this script is to run forecasts initialized from the LETKF analysis.
#######################################################
if [ ! -z ${PBS_O_WORKDIR}    ]; then cd ${PBS_O_WORKDIR}   ;fi
if [ ! -z ${PJM_O_WORKDIR}    ]; then cd ${PJM_O_WORKDIR}   ;fi
if [ ! -z ${SLURM_SUBMIT_DIR} ]; then cd ${SLURM_SUBMIT_DIR};fi

#############
# Servicio Meteorologico Nacional
# Autor: Maximiliano A. Sacco y tantos otros! (Yani, Maru, Cyn, Juan)
# Date: 01/2018
# Ported to hydra
# Date: 07/2023
# Ported to FUGAKU
# Date: 11/2023
#############

#Load experiment configuration
BASEDIR=$(pwd)/../
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/$EXPTYPE.conf
source $BASEDIR/conf/machine.conf

rm $BASEDIR/PROCS/*_ERROR
rm $BASEDIR/PROCS/*_ENDOK

###################################
#We compute the number of forecasts (steps)
####################################
STEP=$INI_STEP

echo "Running DA-FCST cycle starting at STEP: $STEP"

REMAINING_STEPS=$((( ($(date -d "$FCST_END_DATE" +%s) - $(date -d "$FCST_INI_DATE" +%s)) )/$FCST_INI_FREQ + 1 ))
REMAINING_STEPS=$((10#$REMAINING_STEPS-10#$STEP))

echo "The first forecast starts at $FCST_INI_DATE "
echo "The last forecast starts at $FCST_END_DATE"
echo "We did $STEP forecast and we need to perform $REMAINING_STEPS forecast"

rm -f $PROCSDIR/*_ENDOK #Remove control files from previous runs.

####################################
#Main loop over steps (forecast initialization)
####################################

while [ $REMAINING_STEPS -gt 0 ] ; do

   ###### 1st forecast only
   if [[ $STEP == 0 ]]; then
      echo " Step | TimeStamp" > $LOGDIR/da_forecasts.log
   fi
   #####  all forecasts cycles
   echo "Running forecast for initialization: $STEP"
   echo "$(printf "%02d" $STEP)  | $(date +'%T')" >>  $LOGDIR/cycles.log

   echo "Seting important dates: $STEP"
   write_step_conf "FORECAST" #Generate step.conf
   if [[ $RUN_WPS -eq 1 && -z ${EXTWPSPATH} ]] ; then 
      write_step_conf "WPS"
      echo "Running WPS" 
      time $BASEDIR/bin/run_WPS.sh >> $LOGDIR/wps_${STEP}.log   2>&1
      if [ $? -ne 0 ] ; then
         echo "Error: run_WPS finished with errors!"
         echo "Aborting this step"
         exit 1 
      fi
   else 
      if [ ! -z ${EXTWPSPATH} ] ; then
         echo "We will use met_em files from the following path: ${EXTWPSPATH} " 
         echo "WPS is set to 0: I'm not going to run WPS" 
      fi
   fi
   if [[ $BDY_PERT -eq 1 && $RUN_BDY_PERT -eq 1 && -z ${EXTWPSPATH} ]] ; then
      echo "Running Pert met em" 
      time $BASEDIR/bin/run_Pert.sh >> $LOGDIR/pert_met_em_${STEP}.log   2>&1
      if [ $? -ne 0 ] ; then
         echo "Error: run_Pert finished with errors!"
         echo "Aborting this step"
         exit 1 
      fi
   elif [[ $BDY_PERT -eq 0 && -z ${EXTWPSPATH} && ${STEP} -eq 0 ]] ; then
      echo "Linking met_em directory" 
      ln -sf $HISTDIR/WPS/met_em_ori $HISTDIR/WPS/met_em >> $LOGDIR/pert_met_em_${STEP}.log  2>&1
   elif [[ ! -z ${EXTWPSPATH} && ${STEP} -eq 0 ]] ; then 
      mkdir $HISTDIR/WPS/                      >> $LOGDIR/pert_met_em_${STEP}.log  2>&1
      rm -fr $HISTDIR/WPS/met_em               >> $LOGDIR/pert_met_em_${STEP}.log  2>&1  
      ln -sf ${EXTWPSPATH} $HISTDIR/WPS/met_em >> $LOGDIR/pert_met_em_${STEP}.log  2>&1  #We will use existing met_ems from a previous experiment      
   fi

   if [ $RUN_FCST -eq 1 ] ; then
      echo "Running forecast for initialization: $STEP"
      echo "$(printf "%02d" $STEP)  | $(date +'%T')" >>  $LOGDIR/da_forecasts.log
      echo "Running the model" 
      time $BASEDIR/bin/run_ensfcst.sh >> $LOGDIR/dafcst_${STEP}.log  2>&1
      if [ $? -ne 0 ] ; then
         echo "Error: run_DAFcst finished with errors!"
         echo "Aborting STEP "$STEP
         exit 1 
      fi
   fi
   REMAINING_STEPS=$((10#$REMAINING_STEPS-1))
   STEP=$((10#$STEP+1))

done

echo "Exiting... hasta la proxima!"
echo "We finish running da forecasts @ "$(date)
exit 0


