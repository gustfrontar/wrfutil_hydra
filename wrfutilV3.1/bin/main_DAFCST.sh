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

REMAINING_STEPS=$((( ($(date -d "$FORECAST_END_DATE" +%s) - $(date -d "$FORECAST_INI_DATE" +%s)) )/$FORECAST_INI_FREQ + 1 ))
REMAINING_STEPS=$((10#$REMAINING_STEPS-10#$FORECAST_STEP))

echo "The first forecast starts at $FORECAST_INI_DATE "
echo "The last forecast starts at $FORECAST_END_DATE"
echo "We did $FORECAST_STEP forecast and we need to perform $REMAINING_STEPS forecast"

rm -f $PROCSDIR/*_ENDOK #Remove control files from previous runs.

####################################
#Main loop over steps (forecast initialization)
####################################

while [ $REMAINING_STEPS -gt 0 ] ; do
   ###### 1st forecast only
   if [[ $FORECAST_STEP == 0 ]]; then
      echo " Step | TimeStamp" > $LOGDIR/da_forecasts.log
   fi
   #####  all forecasts cycles
   echo "Running forecast for initialization: $FORECAST_STEP"
   echo "$(printf "%02d" $FORECAST_STEP)  | $(date +'%T')" >>  $LOGDIR/da_forecasts.log

   echo "Running the model" > $LOGDIR/dafcst_${PASO}.log
   time $BASEDIR/bin/run_DAFcst.sh >> $LOGDIR/dafcst_${FORECAST_STEP}.log  2>&1
   REMAINING_STEPS=$((10#$REMAINING_STEPS-1))
   FORECAST_STEP=$((10#$FORECAST_STEP+1))
   #Update PASO in the configuration file.
   sed -i -e "/export FORECAST_STEP=/c\\export FORECAST_STEP=$FORECAST_STEP" $BASEDIR/conf/${EXPTYPE}.conf

done

echo "Exiting. . .  . hasta la proxima!"
echo "We finish running da forecasts @ "$(date)
exit 0


