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
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/$EXPTYPE.conf
source $BASEDIR/conf/machine.conf

####################################
#Calculamos la cantidad de pasos
####################################

PASOS_RESTANTES=$((($(date -d "$FECHA_FIN" +%s) - $(date -d "$FECHA_INI" +%s) - $SPIN_UP_LENGTH )/$ANALISIS_FREC ))
PASOS_RESTANTES=$((10#$PASOS_RESTANTES-10#$PASO))

echo "El experimento abarca desde $FECHA_INI hasta $FECHA_FIN"
echo "Se hicieron $PASO pasos de asimilacion y resta hacer $PASOS_RESTANTES"

rm -f $PROCSDIR/*_ENDOK

if [ ! -z ${PJM_SHAREDTMP} -a  ${USETMPDIR} -eq 1 ] ; then 
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


while [ $PASOS_RESTANTES -gt 0 ] ; do
   ###### 1st assimilation cycle only
   if [ $PASO == 0 ]; then
      echo " Step | TimeStamp" > $LOGDIR/cycles.log
      echo "$(printf "%02d" $PASO)  | $(date +'%T')" >>  $LOGDIR/cycles.log
      if [[ $RUN_WPS -eq 1 && -z ${EXTWPSPATH} ]] ; then 
         echo "Running WPS" > $LOGDIR/wps_${PASO}.log
         time $BASEDIR/bin/run_WPS.sh >> $LOGDIR/wps_${PASO}.log
      else 
	 [[ ! -z ${EXTWPSPATH} ]] && echo "We will used met_em files from the following path: ${EXTWPSPATH} " >> $LOGDIR/wps_${PASO}.log
	 [[ $RUN_WPS -eq 1   ]] && echo "RUN_WPS is set to 0" >> $LOGDIR/wps_${PASO}.log
	 echo "I'm not going to run WPS" >> $LOGDIR/wps_${PASO}.log
      fi
      if [[ $BDY_PERT -eq 1 && $RUN_BDY_PERT -eq 1 && -z ${EXTWPSPATH} ]] ; then
	 echo "Running Pert met em" >> $LOGDIR/pert_met_em__${PASO}.log
         time $BASEDIR/bin/run_Pert.sh >> $LOGDIR/pert_met_em_${PASO}.log   2>&1
      elif [[ $BDY_PERT -eq 0 && -z ${EXTWPSPATH} ]] ; then
         echo "Linking met_em directory" >> $LOGDIR/pert_met_em__${PASO}.log
	 ln -sf $HISTDIR/WPS/met_em_ori $HISTDIR/WPS/met_em >> $LOGDIR/pert_met_em_${PASO}.log  2>&1
      elif [ ! -z ${EXTWPSPATH} ] ; then 
	 mkdir $HISTDIR/WPS/                      >> $LOGDIR/pert_met_em_${PASO}.log  2>&1
	 rm -fr $HISTDIR/WPS/met_em               >> $LOGDIR/pert_met_em_${PASO}.log  2>&1  
         ln -sf ${EXTWPSPATH} $HISTDIR/WPS/met_em >> $LOGDIR/pert_met_em_${PASO}.log  2>&1  #We will use existing met_ems from a previous experiment      
      fi
   fi

   #####  all assimilation cycles
   echo "Running cycle: $PASO"
   echo "$(printf "%02d" $PASO)  | $(date +'%T')" >>  $LOGDIR/cycles.log
   echo "Vamos a ejecutar el real, el da_upbdate_bc y el wrf" > $LOGDIR/guess_${PASO}.log
   time $BASEDIR/bin/run_Guess.sh    >> $LOGDIR/guess_${PASO}.log  2>&1
   if [ $PASO -gt 0 ] ; then 
      echo "Vamos a ejecutar el LETKF" > $LOGDIR/letkf_${PASO}.log
      time $BASEDIR/bin/run_LETKF.sh  >> $LOGDIR/letkf_${PASO}.log  2>&1
   fi
   PASOS_RESTANTES=$((10#$PASOS_RESTANTES-1))
   PASO=$((10#$PASO+1))
   echo "Update PASO in the configuration file."
   sed -i -e "/export PASO=/c\\export PASO=$PASO" $BASEDIR/conf/${EXPTYPE}.conf

done

echo "Exiting. . .  . hasta la proxima!"
echo "We finished running @:"$(date )
exit 0


