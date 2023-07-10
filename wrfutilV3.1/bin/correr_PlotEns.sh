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
source $CONFIG  $NOMBRE
[ -f $WRFUTILDIR/RUNs/$NOMBRE/config.env ] && CONFIG=$WRFUTILDIR/RUNs/$NOMBRE/config.env
CONFIG=$(readlink -e ${CONFIG} 2> /dev/null)
[ -z "$CONFIG" ] && dispararError 4 "config.env"
source $CONFIG  $NOMBRE
[ ! -d "$BASEDIR" ] && dispararError 7 "$NOMBRE"
[ ! -f "$BASEDIR/$EXPCONF" ] && dispararError 4 "$BASEDIR/$EXPCONF"
source $BASEDIR/$EXPCONF
[ ! -f "$BASEDIR/$EXPDEP" ] && dispararError 4 "$BASEDIR/$EXPDEP"
source $BASEDIR/$EXPDEP

##### FIN INICIALIZACION ######

cd $BASEDIR
#cp $WRFUTILDIR/templates/PInteres.dat /data/salida/$AMBIENTE/intra/
###
DIRPLOT=$WRFUTILDIR/RUNs/HIST/PLOT/
STARTDATE="${FECHA_INI//\/}_${CICLO}0000"
mkdir -p ${DIRPLOT}/${STARTDATE}/${NOMBRE}/
# Creacion y en colamiento del script
###
read -r -d '' QSCRIPTCMD << "EOF"
export DIRPOST=$WRFUTILDIR/RUNs/HIST/POST/
export DIRWRF=$WRFUTILDIR/RUNs/HIST/FCST/
export DIRPLOT=$WRFUTILDIR/RUNs/HIST/PLOT/
export DIRSALIDA=/data/salida/$AMBIENTE/
export WEBSALIDA=/data/salida/testweb/
export STARTDATE="${FECHA_INI//\/}"
export BASEPLOTDIR=$WRFUTILDIR/bin/python/ploteos/
export MIEMBRO_MIN=10
source $CONDADIR/activate wrfutil  || $CONDADIR/conda activate wrfutil 
cd $BASEPLOTDIR
mkdir -p $WEBSALIDA

	python -W ignore horariosPronoEns.py $STARTDATE $CICLO $ARRAYID
if [[ $ARRAYID -gt 0 ]]
then 
	python -W ignore acumPronoEns.py $STARTDATE $CICLO $ARRAYID
	python -W ignore acumPronoEnsWEB.py $STARTDATE $CICLO $ARRAYID
fi

EOF

# Parametros de encolamiento
QDEPEND_NAME=$DEPPLOT
QDEPENDTYPE=$TYPEPLOT
QPROC_NAME="PLOTEH_${NOMBRE}_${PASO}"
QPROC=32
QTHREADS=1
QARRAY="0-$PLAZO"
QWALLTIME="5:00"
QEXCLU=1

## Encolar
queue


read -r -d '' QSCRIPTCMD << "EOF"
export DIRPOST=$WRFUTILDIR/RUNs/HIST/POST/
export DIRWRF=$WRFUTILDIR/RUNs/HIST/FCST/
export DIRPLOT=$WRFUTILDIR/RUNs/HIST/PLOT/
export DIRSALIDA=/data/salida/$AMBIENTE/
export WEBSALIDA=/data/salida/testweb/
export STARTDATE="${FECHA_INI//\/}"
export OMP_NUM_THREAD=$PLOTTHREADS
PATH=$PATH:$CONDADIR
export MIEMBRO_MIN=10
export BASEPLOTDIR=$WRFUTILDIR/bin/python/ploteos/
export CONCURR=10

mkdir -p $WEBSALIDA
source $CONDADIR/activate wrfutil  || $CONDADIR/conda activate wrfutil 
cd $BASEPLOTDIR/
#param=$(sed -n ${ARRAYID}p $WRFUTILDIR/templates/PInteres.dat)
## La siguiente linea selecciona de todos los pinteres por CIPIs unicos (Si, hace eso !!!)
param=$(awk -F, '{ if (!a[$3]++ ) print ;}' $WRFUTILDIR/templates/PInteres.dat |sed -n $(($ARRAYID))p)
python -W ignore meteogramasEns.py $STARTDATE $CICLO $param
#param=$(cat $WRFUTILDIR/templates/PInteres.dat |grep -U 'WEB,' |sed -n $(($ARRAYID-1))p)
param=$(awk -F, ' /^WEB,/{ if (!a[$3]++ ) print ;}' $WRFUTILDIR/templates/PInteres.dat |sed -n $(($ARRAYID-1))p)
if [[ -n $param ]]
then
	python -W ignore meteogramasEnsWeb.py $STARTDATE $CICLO "$param"
fi
EOF
num=$(awk -F, '{ if (!a[$3]++ ) print ;}' $WRFUTILDIR/templates/PInteres.dat | wc -l)
#num=$((num+1))
# Parametros de encolamiento
QDEPEND_NAME=$DEPPLOT
QDEPENDTYPE=$TYPEPLOT
QPROC_NAME="PLOTM_${NOMBRE}_${PASO}"
QPROC=1
QTHREADS=1
QARRAY="2-$num"
QWALLTIME="20:00"
QEXCLU=0

# Encolar
queue
