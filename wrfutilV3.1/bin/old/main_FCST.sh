#!/bin/bash 
#The main pourpose of this script is to run forecasts nested in the GFS or a global model.

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
#Calculamos la cantidad de pasos
####################################

REMAINING_STEPS=$((( ($(date -d "$FCST_END_DATE" +%s) - $(date -d "$FCST_INI_DATE" +%s)) - $FORECAST_LEAD_TIME )/$FORECAST_INI_FREQ + 1 ))
REMAINING_STEPS=$((10#$REMAINING_STEPS-10#$STEP))

echo "El experimento abarca desde $FCST_INI_DATE hasta $FCST_END_DATE"
echo "Se hicieron $STEP pronosticos y resta hacer $REMAINING_STEPS"

rm -f $PROCSDIR/*_ENDOK #Remove control files from previous runs.

####################################
#Main loop over steps (forecast initialization)
####################################

while [ $REMAINING_STEPS -gt 0 ] ; do
   echo "Inicializando el primer pronostico"
   ###### 1st assimilation cycle only
   if [[ $STEP == 0 ]]; then
      echo " Step | TimeStamp" > $LOGDIR/cycles.log
   fi
   #####  all forecasts cycles
   echo "Running forecast for initialization: $STEP"
   echo "$(printf "%02d" $STEP)  | $(date +'%T')" >>  $LOGDIR/cycles.log

   if [[ $RUN_WPS -eq 1 && -z ${EXTWPSPATH} ]] ; then 
      echo "Running WPS" > $LOGDIR/wps_${STEP}.log
      time $BASEDIR/bin/run_WPS.sh >> $LOGDIR/wps_${STEP}.log   2>&1
      if [ $? -ne 0 ] ; then
         echo "Error: run_WPS finished with errors!"
         echo "Aborting this step"
         exit 1 
      fi
   else 
      [[ ! -z ${EXTWPSPATH} ]] && echo "We will used met_em files from the following path: ${EXTWPSPATH} " >> $LOGDIR/wps_${STEP}.log   2>&1
      [[ $RUN_WPS -eq 1  ]] && echo "RUN_WPS is set to 0" >> $LOGDIR/wps_${STEP}.log   2>&1
      echo "I'm not going to run WPS" >> $LOGDIR/wps_${STEP}.log   2>&1
   fi
   if [[ $BDY_PERT -eq 1 && $RUN_BDY_PERT -eq 1 && -z ${EXTWPSPATH} ]] ; then
      echo "Running Pert met em" > $LOGDIR/pert_met_em_${STEP}.log
      time $BASEDIR/bin/run_Pert.sh >> $LOGDIR/pert_met_em_${STEP}.log   2>&1
      if [ $? -ne 0 ] ; then
         echo "Error: run_Pert finished with errors!"
         echo "Aborting this step"
         exit 1 
      fi
   elif [[ $BDY_PERT -eq 0 && -z ${EXTWPSPATH} && ${STEP} -eq 0 ]] ; then
      echo "Linking met_em directory" >> $LOGDIR/pert_met_em_${STEP}.log
      ln -sf $HISTDIR/WPS/met_em_ori $HISTDIR/WPS/met_em >> $LOGDIR/pert_met_em_${STEP}.log  2>&1
   elif [[ ! -z ${EXTWPSPATH} && ${STEP} -eq 0 ]] ; then 
      mkdir $HISTDIR/WPS/                      >> $LOGDIR/pert_met_em_${STEP}.log  2>&1
      rm -fr $HISTDIR/WPS/met_em               >> $LOGDIR/pert_met_em_${STEP}.log  2>&1  
      ln -sf ${EXTWPSPATH} $HISTDIR/WPS/met_em >> $LOGDIR/pert_met_em_${STEP}.log  2>&1  #We will use existing met_ems from a previous experiment      
   fi

   echo "Running the model" > $LOGDIR/fcst_${STEP}.log
   time $BASEDIR/bin/run_ensfcst.sh >> $LOGDIR/fcst_${STEP}.log  2>&1
   if [ $? -ne 0 ] ; then
      echo "Error: run_FCST finished with errors!"
      echo "Aborting STEP "$STEP
      exit 1 
   fi

   REMAINING_STEPS=$((10#$REMAINING_STEPS-1))
   STEP=$((10#$STEP+1))
   #Update STEP in the configuration file.
   #sed -i -e "/export PASO=/c\\export PASO=$PASO" $BASEDIR/conf/${EXPTYPE}.conf

done

echo "Exiting... hasta la proxima!"
echo "We finish running fcst @ "$(date )
exit 0


