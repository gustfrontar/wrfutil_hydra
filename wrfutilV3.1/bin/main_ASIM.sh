#!/bin/bash 
#######################################################
# ESTE SCRIPT PUEDE SER ENVIADO DIRECTAMENTE A LA COLA
# O EJECUTADO DESDE EL HEAD NODE PARA QUE ENCOLE LOS 
# DIFERENTES PASOS
# En el primer caso la funcion de encolar debe ser SSH
# en el segundo caso PBS_block
#######################################################
# Servicio Meteorologico Nacional
# Autor: Maximiliano A. Sacco y tantos otros! (Yani, Maru, Cyn, Juan)
# Fecha: 01/2018
# Readaptado a hydra
# Fecha: 07/2023
# Portado a Fugaku
# Fecha: 10/2023
#######################################################

if [ ! -z ${PBS_O_WORKDIR}    ]; then cd ${PBS_O_WORKDIR}   ;fi
if [ ! -z ${PJM_O_WORKDIR}    ]; then cd ${PJM_O_WORKDIR}   ;fi
if [ ! -z ${SLURM_SUBMIT_DIR} ]; then cd ${SLURM_SUBMIT_DIR};fi

#Load experiment configuration
BASEDIR=$(pwd)/../
source $BASEDIR/conf/config.env
source $BASEDIR/conf/exp.conf
source $BASEDIR/conf/machine.conf

rm $BASEDIR/PROCS/*_ERROR
rm $BASEDIR/PROCS/*_ENDOK

####################################
#Calculamos la cantidad de pasos
####################################
STEP=$INI_STEP

REMAINING_STEPS=$((($(date -d "$DA_END_DATE" +%s) - $(date -d "$DA_INI_DATE" +%s) - $SPIN_UP_LENGTH )/$ANALYSIS_FREQ ))
REMAINING_STEPS=$((10#$REMAINING_STEPS-10#$STEP))

echo "The first assimilation cycle starts at $DA_INI_DATE "
echo "The last assimilation cycle starts at $DA_END_DATE"
echo "We did $STEP assimilation cycles and we need to perform $REMAINING_STEPS cycles"

rm -f $PROCSDIR/*_ENDOK #Remove control files from previous runs.

if [ ! -z ${PJM_SHAREDTMP} ] && [  ${USETMPDIR} -eq 1 ] ; then 
   echo "We will use Fugaku's temporary directory to speed up IO"
   echo "Copying the data to ${PJM_SHAREDTMP}"
   mkdir -p ${PJM_SHAREDTMP}/HIST   
   cp -r ${HISTDIR}/WPS        ${PJM_SHAREDTMP}/HIST/  &
   cp -r ${HISTDIR}/ANAL       ${PJM_SHAREDTMP}/HIST/  &
   cp -r ${HISTDIR}/GUES       ${PJM_SHAREDTMP}/HIST/  &
   mkdir -p ${PJM_SHAREDTMP}/WRF    ; cp -r $WRFDIR/* ${PJM_SHAREDTMP}/WRF      &
   mkdir -p ${PJM_SHAREDTMP}/LETKF  ; cp -r $LETKFDIR/* ${PJM_SHAREDTMP}/LETKF  &
   mkdir -p ${PJM_SHAREDTMP}/WPS    ; cp -r $WPSDIR/* ${PJM_SHAREDTMP}/WPS      &
   time wait
   echo "Finish copying the data"
   echo $( ls ${PJM_SHAREDTMP}/WRF/ )
fi


####################################
#Main loop over steps (data assimilation cycles)
####################################

while [ $REMAINING_STEPS -gt 0 ] ; do
   ###### 1st assimilation cycle only
   if [ $STEP == 0 ]; then
      echo " Step | TimeStamp" > $LOGDIR/cycles.log
      echo "$(printf "%02d" $STEP)  | $(date +'%T')"  >>  $LOGDIR/cycles.log
      if [ ! -z ${EXTWPSPATH} && $RUN_WPS -eq 0 ] ; then 
	 mkdir $HISTDIR/WPS/                           > $LOGDIR/pert_met_em_${STEP}.log  2>&1
	 rm -fr $HISTDIR/WPS/met_em                   >> $LOGDIR/pert_met_em_${STEP}.log  2>&1  
         ln -sf ${EXTWPSPATH} $HISTDIR/WPS/met_em_ori >> $LOGDIR/pert_met_em_${STEP}.log  2>&1  #We will use existing met_ems from a previous experiment      
      fi
   fi  #End of the special case in which WPS is run at the begining of the experiment

   #####  all forecasts cycles
   echo "Running data assimilation cycle: $STEP"
   echo "$(printf "%02d" $STEP)  | $(date +'%T')" >>  $LOGDIR/cycles.log

   if [ $RUN_WPS == 1 ] ; then 
      write_step_conf "WPS" #Generate step.conf
      echo "Running WPS" > $LOGDIR/wps_${STEP}.log
      time $BASEDIR/bin/run_WPS.sh                        >> $LOGDIR/wps_${STEP}.log
      if [ $? -ne 0 ] ; then
         echo "Error: run_WPS finished with errors!"
         echo "Aborting this step"
         exit 1 
      fi
      echo "Succesfully run WPS"
   fi

   if [[ $BDY_PERT -eq 1 && $RUN_BDY_PERT -eq 1 ]] ; then
      echo "Running Pert met em"                          >> $LOGDIR/pert_met_em_${STEP}.log
      time $BASEDIR/bin/run_Pert.sh                       >> $LOGDIR/pert_met_em_${STEP}.log   2>&1
      if [ $? -ne 0 ] ; then
         echo "Error: run_Pert finished with errors!"
         echo "Aborting this step"
         exit 1 
      fi
      echo "Succesfully run Pert met em"

   elif [ $BDY_PERT == 0 ] && [ $STEP == 0 ] ; then
      echo "Linking met_em directory"                     >> $LOGDIR/pert_met_em_${STEP}.log
      ln -sf $HISTDIR/WPS/met_em_ori $HISTDIR/WPS/met_em  >> $LOGDIR/pert_met_em_${STEP}.log  2>&1
   fi

   if [ $RUN_FCST -eq 1 ] ; then
      echo "Running forecast for STEP: $STEP " > $LOGDIR/guess_${STEP}.log
      write_step_conf "GUESS" #Generate step.conf
      time $BASEDIR/bin/run_ensfcst.sh         >> $LOGDIR/guess_${STEP}.log  2>&1
      if [ $? -ne 0 ] ; then
         echo "Error: run_Guess finished with errors!"
         echo "Aborting this step"
         exit 1 
      fi
   fi

   if [ $RUN_LETKF -eq 1 ] && [ $STEP -gt 0 ] ; then 
      echo "Running LETKF" > $LOGDIR/letkf_${STEP}.log
      time $BASEDIR/bin/run_LETKF.sh  >> $LOGDIR/letkf_${STEP}.log  2>&1
      if [ $? -ne 0 ] ; then
         echo "Error: run_LETKF finished with errors!"
         echo "Aborting this step"
         exit 1 
      fi
   fi

   REMAINING_STEPS=$((10#$REMAINING_STEPS-1))
   STEP=$((10#$STEP+1))
   #echo "Update STEP in the configuration file."
   #sed -i -e "/export STEP=/c\\export STEP=$STEP" $BASEDIR/conf/${EXPTYPE}.conf

done

echo "Exiting... hasta la proxima!"
echo "We finished running @:"$(date )
exit 0


