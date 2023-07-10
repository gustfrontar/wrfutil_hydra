#!/bin/bash 
#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF

	Ud. deberia usar este escript de la siguiente manera:
		$0 <nombre entonro > 
	Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${1:?"$USO"}  ${WRFUTILDIR:?"$USO"}


### PARAMETROS

NOMBRE=$1

### CONFIGURACION

[ ! -f $WRFUTILDIR/lib/errores.env ] && exit 1
source $WRFUTILDIR/lib/errores.env
CONFIG=$WRFUTILDIR/config.env
[ -f $WRFUTILDIR/RUNs/$NOMBRE/config.env ] && CONFIG=$WRFUTILDIR/RUNs/$NOMBRE/config.env
CONFIG=$(readlink -e ${CONFIG} 2> /dev/null)
[ -z "$CONFIG" ] && dispararError 4 "config.env"
source $CONFIG  $NOMBRE
[ ! -d "$BASEDIR" ] && dispararError 7 "$NOMBRE"
[ ! -f "$BASEDIR/$EXPCONF" ] && dispararError 4 "$BASEDIR/$EXPCONF"
source $BASEDIR/$EXPCONF
[ ! -f "$BASEDIR/$EXPDEP" ] && dispararError 4 "$BASEDIR/$EXPDEP"
source $BASEDIR/$EXPDEP
[ ! -f "$BASEDIR/$EXPASIM" ] && dispararError 4 "$BASEDIR/$EXPASIM"
source $BASEDIR/$EXPASIM


##### FIN INICIALIZACION ######


mkdir -p $DIROBSREPO $OBSprocesados 

# Fecha del analisis 
read -r anio mes dia hora min  <<< $(date -u -d "$FECHA_OBS +$CICLO_OBS hours +$CICLO_OBS_MIN minutes" +"%Y %m %d %H %M") 

# fecha fin de la ventana de observaciones
read -r Fanio Fmes Fdia Fhora Fmin  <<< $(date -u -d "$FECHA_OBS +$CICLO_OBS hours +$((10#$CICLO_OBS_MIN-10#$DESPLAZAMIENTO+10#$OBSWIN)) minutes" +"%Y %m %d %H %M")  

#fecha de inicio de la ventana de observaciones
read -r Ianio Imes Idia Ihora Imin  <<< $(date -u -d "$FECHA_OBS +$CICLO_OBS hours +$((10#$CICLO_OBS_MIN-10#$DESPLAZAMIENTO)) minutes" +"%Y %m %d %H %M")

fin="$(date -u -d "${Fanio}/${Fmes}/${Fdia} ${Fhora}:${Fmin} "  +"%Y-%m-%d %H:%M %z")"
ini="$(date -u -d "${Ianio}/${Imes}/${Idia} ${Ihora}:${Imin} "  +"%Y-%m-%d %H:%M %z")"
inc=$OBSFREC     ## Este incremento es en minutos                               

Ijuliano="$(date -d "${Ianio}/${Imes}/${Idia}" +"%j")"
Fjuliano="$(date -d "${Fanio}/${Fmes}/${Fdia}" +"%j")"
juliano="$(date -d "${anio}/${mes}/${dia}" +"%j")"

#####
## Def direcorios de obs
####
amdarDir=$DIROBSREPO
amvDir=$DIROBSREPO
emasDir=$DIROBSREPO/EMAs/
sfcDir=$DIROBSREPO
sondeosDir=$DIROBSREPO
shipDir=$DIROBSREPO
boyasDir=$DIROBSREPO

OBSBIN=${DIROBSBIN}/${anio}${mes}${dia}_${hora}${min}00/
#OBSASIM=${DIROBSASIM}/${anio}${mes}${dia}_${hora}${min}00/
OBSASIM=${DIROBSASIM}/

cd $OBSASIM

read -r -d '' QSCRIPTCMD << EOF

mkdir -p $OBSBIN $OBSASIM

######
## Procesamos los airs
######

#Fede: Agregue este source porque el decoder necesita algunas librerias
source $BASEDIR/LETKF/envvars.sh

cd $DIROBSREPO/AIRS/
for FILE in \$(ls *.hdf)
do
        $SCRIPTDIR/decoder.airs \$FILE $NIVELAIRS $anio $mes $dia $hora $min
done
mv AIRSRT_* $OBSBIN

#####
## Procesamos las Observaciones
#####
cd $SCRIPTDIR
export DATOS_DIR=$DIROBSREPO
export DATOS_OUT=$OBSBIN
iter="$ini"

#Fede: Le puse un -u a los date para que siga tomando la hora UTC
# le agregamos /home/fcutraro/.conda/envs/RRA/bin/  porque sino no encontraba el python desde el crontab
# OJO si cambiamos de usuario !!!
PATH=$CONDADIR:$PATH
source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil

while [[ "\$(date -u -d "\$iter" +"%Y%m%d%H%M")" -le "$(date -u -d "$fin"  +"%Y%m%d%H%M")"  ]]
do
	cd $SCRIPTDIR
        read -r anio mes dia hora min  <<< "\$(date -u -d "\$iter" +"%Y %m %d %H %M")"  #con el formato deseado
	echo "Procesando observaciones de :  \$anio/\$mes/\$dia \$hora:\$min"
	export DATOS_DIR=$amdarDir
 	python $SCRIPTDIR/amdar_RRA.py \$anio \$mes \$dia \$hora \$min
	export DATOS_DIR=$amvDir
	python $SCRIPTDIR/amv_RRA.py \$anio \$mes \$dia \$hora \$min
	export DATOS_DIR=$emasDir
        python $SCRIPTDIR/emas_HM.py \$anio \$mes \$dia \$hora \$min
	export DATOS_DIR=$sfcDir
	python $SCRIPTDIR/sfc_RRA.py \$anio \$mes \$dia \$hora \$min
	export DATOS_DIR=$sondeosDir
	python $SCRIPTDIR/sondeos_RRA.py \$anio \$mes \$dia \$hora \$min
	export DATOS_DIR=$shipDir
	python $SCRIPTDIR/ship_RRA.py \$anio \$mes \$dia \$hora \$min
	export DATOS_DIR=$boyasDir
	python $SCRIPTDIR/boyas_RRA.py \$anio \$mes \$dia \$hora \$min
	cd $OBSBIN
	cat *_\${anio}\${mes}\${dia}\${hora}\${min}00.dat > $OBSASIM/obs_\${anio}\${mes}\${dia}\${hora}\${min}.dat
        iter="\$(date  -u -d "\$iter + $inc minutes" +"%Y-%m-%d %H:%M %z")"
done
source $CONDADIR/deactivate || $CONDADIR/conda deactivate 

EOF
# Parametros de encolamiento
QDEPEND_NAME=""
QDEPENDTYPE=""
QPROC_NAME=PROCOBS_${1}_$(date +"%s")
QPROC=$OBSPROC
QTHREADS=$OBSTHREADS
QARRAY=""
QWALLTIME=$OBSWALLTIME
QEXCLU=0
queue

#OLD_IFS="$IFS"
#IFS=","
#export SO_INSTRUMENT_LIST=$(echo "'${RADAR_LIST[*]}'")
#IFS="$OLD_IFS"


########
## Procesando Radares
########

read -r -d '' QSCRIPTCMD << EOF
#####
## Procesamos las Observaciones
#####
cd $RADSCRIPTDIR
export SO_OUTPUTROOT=$OBSBIN
mkdir -p \$SO_OUTPUTROOT
export RADAR_TIME=10  # tiempo en minutos que un radar tarda en generar un volumen
export ADD_SLOTS=\$((((10#\$RADAR_TIME+10#\$OBSFREC)-1)/10#\$OBSFREC))

PATH=$CONDADIR:$RADSCRIPTDIR:$PATH
source $CONDADIR/activate wrfutil2 || $CONDADIR/conda activate wrfutil2
#source $RADSCRIPTDIR/experimento.qcr
export OMP_NUM_THREADS=$RADOBSTHREADS
export GOMP_STACKSIZE=512m
export OMP_STACKSIZE=512m
export SO_INSTRUMENT_LIST=$SO_INSTRUMENT_LIST

ulimit -s unlimited
delta=$(($OBSFREC/2))  # observaciones centradas el slot

slots=()
for slot in \${RADAR_SLOTSNUM[@]}
do 
        for i in \$(seq \$((\$slot-10#\$ADD_SLOTS)) \$slot)
        do 
                slots+=(\$((10#\$i*\$OBSFREC)) )
        done
done

# LA linea de abajo quita repetidos del array 
slots=($(echo "\${slots[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

#ini=\$(date -u -d "$ini +\$((\${slots[0]} - \$delta)) minutes" +"%Y%m%d_%H%M00")
echo "Se procesaran los siguientes SLOTS: \${slots[@]}"
for slot in \${slots[@]}
do
        export SO_REFERENCE_DATE=\$(date -u -d "$ini +\$slot minutes" +"%Y%m%d_%H%M00")
        echo "Procesando radar de :  \$SO_REFERENCE_DATE"
        python $SCRIPTDIR/radar_so/so-radar.py
done
source $CONDADIR/deactivate || $CONDADIR/conda deactivate 
EOF

# Parametros de encolamiento
QDEPEND_NAME=""
QDEPENDTYPE=""
QPROC_NAME=PROCOBSRAD_${1}_$(date +"%s")
QPROC=$RADOBSPROC
QTHREADS=$RADOBSTHREADS
QARRAY=""
QWALLTIME=$RADOBSWALLTIME
QEXCLU=0
queue
