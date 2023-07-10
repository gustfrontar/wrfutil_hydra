#!/bin/bash 

#############
# Servicio Meteorologico Nacional
# Autor: Maximiliano A. Sacco y tantos otros! (Yani, Maru, Cyn, Juan)
# Fecha: 01/2018
# Readaptado a hydra
# Fecha: 07/2023
#############

### PARAMETROS
BASEDIR=$(pwd)/../
source $BASEDIR/lib/errores.env
CONFIG=$BASEDIR/conf/config.env
[ ! -e "$CONFIG" ] && dispararError 4 "Error: No encontre config.env"
source $CONFIG

### CONFIGURACION
[ ! -f "$BASEDIR/conf/$EXPMACH" ] && dispararError 4 "$BASEDIR/conf/$EXPMACH"
source $BASEDIR/conf/$EXPMACH   #Cargo las variables de la cola.
[ ! -f "$BASEDIR/conf/$EXPCONF" ] && dispararError 4 "$BASEDIR/conf/$EXPCONF"
source $BASEDIR/conf/$EXPCONF   #Cargo la configuracion del experimento.
[ ! -f "$BASEDIR/conf/$EXPDEP" ] && dispararError 4 "$BASEDIR/conf/$EXPDEP"
source $BASEDIR/conf/$EXPDEP    #Cargo la dependencia de los scripts.
[ ! -f "$BASEDIR/conf/$EXPASIM" ] && dispararError 4 "$BASEDIR/conf/$EXPASIM"
source $BASEDIR/conf/$EXPASIM   #Cargo la configuracion del exp de asimilacion

####################################
#Calculamos la cantidad de pasos
####################################

PASOS_RESTANTES=$((($(date -d "$FECHA_FIN" +%s) - $(date -d "$FECHA_INI" +%s))/$ANALISIS_FREC))
PASOS_RESTANTES=$((10#$PASOS_RESTANTES-10#$PASO))

echo "El experimento abarca desde $FECHA_INI hasta $FECHA_FIN"
echo "Se hicieron $PASO pasos de asimilacion y resta hacer $PASOS_RESTANTES"


while [ $PASOS_RESTANTES -gt 0 ] ; do

   ###### 1st assimilation cycle
   if [[ $PASO -lt 0 ]]; then
      echo "Reinicializando asimilacion"
      exit
   elif [[ $PASO == 0 ]]; then
      rm -f $PROCSDIR/*_ENDOK
      echo "Corriendo primer paso de asimilacion"
      echo " Stage |  Anl Date  | Step | TimeStamp" > $BASEDIR/LOGS/log_cycles.txt

      #DEBUG DEBUG $BASEDIR/bin/correr_WPS.sh

      if [ $PERTURBAR -eq 1 ] ; then
	 echo "Vamos a perturbar los met_em"
         #DEBUG DEBUG TODO $BASEDIR/bin/correr_Pert.sh
      else 

	 ln -sf $HISTDIR/WPS/met_em_ori $HISTDIR/WPS/met_em
      fi

      #TODO $BASEDIR/bin/correr_OBS.sh  #Armar el script que genera las observaciones para el letkf.

   fi

   ##### 2nd-inf assimilation cycle
   echo "Ejecutando paso: $PASO"
   echo "  INI  | $FECHA_CICLO |  $(printf "%02d" $PASO)  | $(date +'%s')" >>  $BASEDIR/LOGS/log_cycles.txt

   echo "Vamos a ejecutar el real, el da_upbdate_bc y el wrf"
   $BASEDIR/bin/correr_Guess.sh

   echo "Actualizando la fecha "
   FECHA_CICLO=$(date -u --date "$FECHA_INI UTC + $((10#$PASO*10#$ANALISIS_FREC)) seconds " +"%Y-%m-%d %T")

   #$BASEDIR/bin/correr_LETKF.sh  
   #echo "LETKF encolado a punto de lanzar el POST:" $(date )
   #$BASEDIR/bin/correr_LETKFpost.sh  

   PASOS_RESTANTES=$((10#$PASOS_RESTANTES-1))
   PASO=$((10#$PASO+1))
   #TODO Tengo que actualizar el paso en el archivo de configuracion

done


#echo 'LETKF MONIT'
#$WRFUTILDIR/bin/correr_LETKFmonit.sh 12

echo "Saliendo. . .  . hasta la proxima!"
echo "Termiamos de correr a:"$(date )
exit 0


