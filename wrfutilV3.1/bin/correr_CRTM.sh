#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

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
[ ! -d "$BASEDIR" ] && dispararError 7 "$NOMBRE"
[ ! -f "$BASEDIR/$EXPCONF" ] && dispararError 4 "$BASEDIR/$EXPCONF"
source $BASEDIR/$EXPCONF
[ ! -f "$BASEDIR/$EXPDEP" ] && dispararError 4 "$BASEDIR/$EXPDEP"
source $BASEDIR/$EXPDEP
[ ! -f "$BASEDIR/$EXPPLGIN" ] && dispararError 4 "$BASEDIR/$EXPPLGIN"
source $BASEDIR/$EXPPLGIN


##### FIN INICIALIZACION ######

ulimit -s unlimited

cd $BASEDIR/CRTM/
mkdir -p run
tar xf CRTM.tar -C "$BASEDIR/CRTM/run/"

read -r -d '' QSCRIPTCMD << "EOF"
plazoID=0
test $ARRAYID && plazoID=$ARRAYID
ulimit -s unlimited
ulimit -l unlimited
echo "Procesando hora $plazoID"
read -r IY IM ID IH Im  <<< $(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN)) minutes" +"%Y %m %d %H %M")
read -r WY WM WD WH Wm  <<< $(date -u -d "$FECHA_INI UTC +$((10#$CICLO+$plazoID)) hours+$((10#$CICLO_MIN)) minutes" +"%Y %m %d %H %M")
read -r FY FM FD FH Fm  <<< $(date -u -d "$FECHA_INI UTC +$((10#$CICLO+$PLAZO)) hours+$((10#$CICLO_MIN+10#$PLAZO_MIN)) minutes" +"%Y %m %d %H %M")
cd $WRFDIR

PATH=$CONDADIR:$PATH

RUNDIR="$BASEDIR/WRF/00/WRFOUT/"
FECHA=`echo $FECHA_INI$CICLO| cut -d'/' -s -f1,2,3 --output-delimiter=''`
OUT_PATH=$BASEDIR/HIST/CRTM/${IY}${IM}${ID}_${IH}${Im}00/$NOMBRE/DATA/

# Parametros de entrada del run_CRTM
export COEFF_PATH="$BASEDIR/CRTM/run/"
export FILEIN="$RUNDIR/wrfout_d01_${WY}-${WM}-${WD}_${WH}_${Wm}_00"		
export OUT_PATH	
export BASEPLOTDIR=$WRFUTILDIR/bin/python/ploteos/CRTM/

cd $COEFF_PATH
source envvars.sh

mkdir -p $OUT_PATH
echo "Procesando el archivo $FILEIN" 

$BASEDIR/CRTM/run/run_CRTM 


export CRTM_FIG=$BASEDIR/HIST/CRTM/${IY}${IM}${ID}_${IH}${Im}00/$NOMBRE/PLOT/  # Figuras del historico
export DIRSALIDA=/data/salida/$AMBIENTE/intra/   # Figuras para la intramet

mkdir -p $CRTM_FIG

cd $BASEPLOTDIR

source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil
param="$OUT_PATH/G16.${IY}${IM}${ID}_${IH}${Im}00.$(printf "%03g" $plazoID).nc"
python $BASEPLOTDIR/grafica_CRTM.py  $param

EOF

# Parametros de encolamiento
QDEPEND_NAME=$DEPCRTM
QDEPENDTYPE=$TYPECRTM
QPROC_NAME=CRTM_$1_$PASO
QPROC=$CRTMPROC
QTHREADS=$CRTMTHREADS
QARRAY="1-$((10#$PLAZO))"
QWALLTIME=$CRTMWALLTIME
QEXCLU=1

# Encolar
queue


