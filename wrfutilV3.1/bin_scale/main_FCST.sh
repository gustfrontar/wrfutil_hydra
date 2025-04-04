#!/bin/bash 
#The main purpose of this script is to run forecasts nested in the GFS or a global model.

#######################################################
# ESTE SCRIPT PUEDE SER ENVIADO DIRECTAMENTE A LA COLA
# O EJECUTADO DESDE EL HEAD NODE PARA QUE ENCOLE LOS 
# DIFERENTES PASOS
# En el primer caso la funcion de encolar debe ser SSH
# en el segundo caso PBH_block
#######################################################
if [ ! -z ${PBS_O_WORKDIR}    ]; then cd ${PBS_O_WORKDIR}   ;fi
if [ ! -z ${PJM_O_WORKDIR}    ]; then cd ${PJM_O_WORKDIR}   ;fi
if [ ! -z ${SLURM_SUBMIT_DIR} ]; then cd ${SLURM_SUBMIT_DIR};fi


#############
# Servicio Meteorologico Nacional
# Autor: Maximiliano A. Sacco y tantos otros! (Yani, Maru, Cyn, Juan)
# Fecha: 01/2018
# Readaptado a hydra
# Fecha: 07/2023
# Readaptado a SCALE por Arata Amemiya
# Fecha: 03/2025
#############

#Load experiment configuration
BASEDIR=$(pwd)/../
source $BASEDIR/conf/config.env
source $BASEDIR/conf/exp.conf
source $BASEDIR/conf/model.conf
source $BASEDIR/conf/machine.conf
if [ $MACHINE == "FUGAKU" ] ;then
  source $TOPDIR/setup_spack.sh
fi

rm -f $BASEDIR/PROCS/*_ERROR
rm -f $BASEDIR/PROCS/*_ENDOK

###################################
#We compute the number of forecasts (steps)
####################################
STEP=$INI_STEP

echo "Running DA-FCST cycle starting at STEP: $STEP"

REMAINING_STEPS=$((( ($(date -d "$FCST_END_DATE" +%s) - $(date -d "$FCST_INI_DATE" +%s)) )/$FCST_INI_FREQ + 1 ))
REMAINING_STEPS=$(( 10#$REMAINING_STEPS - 10#$STEP ))

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

   echo "Setting important dates: $STEP"
   write_step_conf "FORECAST" #Generate step.conf

   if [[ $RUN_WPS -eq 1 ]] ; then 
      time $BASEDIR/bin/run_prep.sh >> $LOGDIR/wps_${STEP}.log   2>&1
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

