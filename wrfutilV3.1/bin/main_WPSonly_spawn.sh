#!/bin/bash 
#######################################################
# ESTE SCRIPT PUEDE SER ENVIADO DIRECTAMENTE A LA COLA
# O EJECUTADO DESDE EL HEAD NODE PARA QUE ENCOLE LOS 
# DIFERENTES PASOS
# En el primer caso la funcion de encolar debe ser SSH
# en el segundo caso PBH_block
#############
# Servicio Meteorologico Nacional
# Autor: Maximiliano A. Sacco y tantos otros! (Yani, Maru, Cyn, Juan)
# Fecha: 01/2018
# Readaptado a hydra
# Fecha: 07/2023

#Get the working directory
if [ ! -z ${PBS_O_WORKDIR}    ]; then cd ${PBS_O_WORKDIR}   ;fi
if [ ! -z ${PJM_O_WORKDIR}    ]; then cd ${PJM_O_WORKDIR}   ;fi
if [ ! -z ${SLURM_SUBMIT_DIR} ]; then cd ${SLURM_SUBMIT_DIR};fi

#Load experiment configuration
BASEDIR=$(pwd)/../
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/assimilation.conf
source $BASEDIR/conf/machine.conf 
#Set some environmental parameters
eval "$ENVSET"

####################################
#Calculamos la cantidad de pasos
####################################
PASOS_RESTANTES=$((($(date -d "$FECHA_FIN" +%s) - $(date -d "$FECHA_INI" +%s))/$ANALISIS_FREC ))
PASOS_RESTANTES=$((10#$PASOS_RESTANTES-10#$PASO))

echo "El experimento abarca desde $FECHA_INI hasta $FECHA_FIN"
echo "Se hicieron $PASO pasos de asimilacion y resta hacer $PASOS_RESTANTES"

rm -f $PROCSDIR/*_ENDOK

###### 1st assimilation cycle only
echo " Step | TimeStamp" > $LOGDIR/cycles.log
if [ $RUN_WPS -eq 1 ] ; then 
   echo "Corriendo el WPS"
   time $BASEDIR/bin/correr_WPS_spawn.sh > $LOGDIR/log_wps_${PASO}.log
fi
if [[ $BDY_PERT -eq 1 && $RUN_BDY_PERT -eq 1 ]] ; then
   echo "Vamos a perturbar los met_em"
   time $BASEDIR/bin/correr_Pert_spawn.sh > $LOGDIR/pert_met_em_${PASO}.log   2>&1
elif [ $BDY_PERT -eq 0 ] ; then
   echo "Linking met_em directory"
   time ln -sf $HISTDIR/WPS/met_em_ori $HISTDIR/WPS/met_em > $LOGDIR/pert_met_em_${PASO}.log  2>&1
fi

echo "Saliendo. . .  . hasta la proxima!"
echo "Termiamos de correr a:"$(date )
exit 0


