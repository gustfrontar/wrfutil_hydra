#!/bin/bash 
################################################################
# Author: DMSR, Servicio Meteorologico Nacional                #
# Create: 12/2021 - M. Sacco; P. Maldonado                     #
################################################################

###  Leyenda de USO
read -r -d '' USO << EOF
   Ud. deberia usar este escript de la siguiente manera:
      $0 <nombre entorno >
   Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${1:?"$USO"}  ${WRFUTILDIR:?"$USO"}

### PARAMETROS
export NOMBRE=$1

### CONFIGURACION
[ ! -f $WRFUTILDIR/lib/errores.env ] && exit 1
source $WRFUTILDIR/lib/errores.env
CONFIG=$WRFUTILDIR/config.env
[ -f $WRFUTILDIR/RUNs/$NOMBRE/config.env ] && CONFIG=$WRFUTILDIR/RUNs/$NOMBRE/config.env
CONFIG=$(readlink -e ${CONFIG} 2> /dev/null)
[ -z "$CONFIG" ] && dispararError 4 "config.env"
source $CONFIG  $NOMBRE
[ ! -d "$BASEDIR" ] && dispararError 7 "$BASEDIR"
[ ! -f "$BASEDIR/$EXPCONF" ] && dispararError 4 "$BASEDIR/$EXPCONF"
source $BASEDIR/$EXPCONF
[ ! -f "$BASEDIR/$EXPASIM" ] && dispararError 4 "$BASEDIR/$EXPASIM"
source $BASEDIR/$EXPASIM

### Update PASO
PASO=$(($PASO+1))
sed -i -e "/export PASO=/c\\export PASO=$PASO" $BASEDIR/experimento.asim

### Check reinitialization 
[[ ! -z $REINIT ]] && [[ $PASO -eq $REINIT ]] && PASO=-1 

### Select case base on PASO
##### Reinitialization
if [[ $PASO -lt 0 ]]; then

   echo "Reinicializando asimilacion"

   # Remove LETKF error control file
   rm -f $BASEDIR/FATAL.ERROR
  
   # Update PASO, FECHA_INI, CICLO 
   PASO=$(($PASO+1))
   sed -i -e "/export PASO=/c\\export PASO=$PASO" $BASEDIR/experimento.asim

   FECHA_INI=$(date -u --date "$FECHA_INI UTC + $(($PASO*$ANALISIS)) seconds + $CICLO_MIN minutes" +"%Y/%m/%d")
   CICLO=$(printf "%02d" $(((10#$CICLO + $ANALISIS/60)%24)))
   CICLO_MIN=$(printf "%02d" $(((10#$CICLO_MIN + $ANALISIS)%60)))
   sed -i -e "/export FECHA_INI=/c\\export FECHA_INI=\'$FECHA_INI\'" $BASEDIR/experimento.conf

###### 1st assimilation cycle
elif [[ $PASO == 0 ]]; then
   echo "Corriendo primer paso de asimilacion"
   echo " Stage |  Anl Date  | Cycle | Step | TimeStamp" > $BASEDIR/LOGS/log_cycles.txt

##### 2nd-inf assimilation cycle
else

   echo "Actualizando la fecha de inicio y el ciclo"
   FECHA_INI=$(date -u --date "$FECHA_INI UTC + $((10#$CICLO+$ANALISIS/60)) hours + $CICLO_MIN minutes" +"%Y/%m/%d")
   CICLO=$(printf "%02d" $(((10#$CICLO + $ANALISIS/60)%24)))
   CICLO_MIN=$(printf "%02d" $(((10#$CICLO_MIN + $ANALISIS)%60)))
   sed -i -e "/export FECHA_INI=/c\\export FECHA_INI=\'$FECHA_INI\'" $BASEDIR/experimento.conf
   sed -i -e "/export CICLO=/c\\export CICLO=$CICLO" $BASEDIR/experimento.conf
   sed -i -e "/export CICLO_MIN=/c\\export CICLO_MIN=$CICLO_MIN" $BASEDIR/experimento.conf

fi

echo "Ejecutando paso: $PASO"
echo "  INI  | $FECHA_INI |   $CICLO  |  $(printf "%02d" $PASO)  | $(date +'%s')" >>  $BASEDIR/LOGS/log_cycles.txt

