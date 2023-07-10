#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF

        Ud. deberia usar este escript de la siguiente manera:
                $0 <nombre entorno > <cantidad de ciclos previos>
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${1:?"$USO"} ${2:?"$USO"} ${WRFUTILDIR:?"$USO"}


### PARAMETROS
export NOMBRE=$1
export NTIMES=$2
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
######### FIN VERIFICACION 

cd $BASEDIR/LETKF
#script de ejecucion
read -r -d '' QSCRIPTCMD << EOF
ulimit -s unlimited
ulimit -l unlimited
source envvars.sh
source $CONDADIR/deactivate || $CONDADIR/conda deactivate
source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil

export PATH_OBS=$BASEDIR/HIST/OBS/$NOMBRE
export PATH_PLOT=$BASEDIR/HIST/PLOT/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/
mkdir -p $PATH_PLOT

NSLOTS=$(($OBSWIN/$OBSFREC+1))
FECHA_ASSIM="$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y %m %d %H %M")"

#--------------------#
# Multiple cycles
#--------------------#
python $WRFUTILDIR/bin/python/plot_obsdat_obsmean.py \$FECHA_ASSIM \$NTIMES 

python $WRFUTILDIR/bin/python/plot_obsdat_dastats.py \$FECHA_ASSIM \$NTIMES

python $WRFUTILDIR/bin/python/plot_obsdat_obsstats.py \$FECHA_ASSIM \$NTIMES

#python $WRFUTILDIR/bin/python/plot_obsdat_obsstats_from_pkl.py \$FECHA_ASSIM \$NTIMES


EOF

# Parametros de encolamiento
QDEPEND_NAME=""
QDEPENDTYPE=""
QPROC_NAME=LETKFMONIT_$1_$PASO
QPROC=$LETKFPOSTPROC
QTHREADS=$LETKFPOSTTHREADS
QARRAY=""
QWALLTIME=$LETKFPOSTWALLTIME
QEXCLU=1

# Encolar
queue

