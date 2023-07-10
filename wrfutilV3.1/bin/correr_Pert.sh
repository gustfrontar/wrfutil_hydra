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

##### FIN INICIALIZACION ######


###
# Creacion y en colamiento del script
###
cd $BASEDIR/WPS

### Verifico que esten todos los METGRIDs del DETERMINISTICO
### ES MUY DIFICL VERIFICAR, porque hay que hacerlo en el script que se encola
### y se encolan por hora y no por miembro. . . y los 'ls' tardan mucho como para repetirlos
FECHA_PERT=$(date -u -d "$FECHA_INI UTC +$CICLO hours" +"%Y%m%d_%H%M%S")
MET_ORI="$BASEDIR/HIST/WPS/$FECHA_PERT/$NOMBRE/originales/"
MET_PERT="$BASEDIR/HIST/WPS/$FECHA_PERT/$NOMBRE/perturbados/"
MET_CNT=$((10#$PLAZO * 10#$WPSPROC))
#cantidad= $(ls $MET_ORI/*|wc -w)
#[[ $MET_CNT -ge $cantidad ] && dispararError 8 "Se necesitan al menos $MET_CNT archivos en $MET_ORI y hay $cantidad"
#for MIEM in $(seg -f "%02g" $MIEMBRO_FIN $MIEMBRO_INI)
#do
#	cantidad= $(ls $MET_DET/$MIEM/*|wc -w)
#	[[ $MET_CNT -ge $cantidad ]] && dispararError 8 "Se necesitan al menos $MET_CNT archivos en $MET_DET/$MIEM/ y hay $cantidad"

#Creo los directorios para guardar los perturbados
eval mkdir -p $MET_PERT/{$MIEMBRO_INI..$MIEMBRO_FIN}

read -r -d '' QSCRIPTCMD << "EOF"
MIEM=$(printf "%02g" $((10#$MIEMBRO_FIN-10#$MIEMBRO_INI)))
test $ARRAYID && MIEM=$(printf "%02g" $ARRAYID)
ulimit -s unlimited
sleep  $((1*10#$MIEM))
export FECHA_PERT=$(date -u -d "$FECHA_INI UTC +$CICLO hours" +"%Y%m%d_%H%M%S")
export MIEMBROS=$((10#$MIEMBRO_FIN-10#$MIEMBRO_INI+1))
export MET_ORI="$BASEDIR/HIST/WPS/$FECHA_PERT/$NOMBRE/originales/"
export MET_PERT="$BASEDIR/HIST/WPS/$FECHA_PERT/$NOMBRE/perturbados/"
export ARRAYCNT=$SLURM_ARRAY_TASK_COUNT
export ARRAYID=$SLURM_ARRAY_TASK_ID
PATH=$CONDADIR:$PATH

source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil
python $WRFUTILDIR/bin/python/perturbar.py 
EC=$?
[[ $EC -ne 0 ]] && rm  $BASEDIR/HIST/WPS/last/${NOMBRE}_perturbados.txt
[[ $EC -ne 0 ]] && dispararError 9 "perturbar.py"

echo $FECHA_PERT > $BASEDIR/HIST/WPS/last/${NOMBRE}_perturbados.txt

EOF

# Parametros de encolamiento
QDEPEND_NAME=${DEPPERT}
QDEPENDTYPE=${TYPEPERT}
QPROC_NAME="PERT_${NOMBRE}_${PASO}"
QPROC=${PERTPROC}
QTHREADS=${PERTTHREADS}
QARRAY="0-$PLAZO"
QWALLTIME=${PERTWALLTIME}
QEXCLU=1

# Encolar
queue

