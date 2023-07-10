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


###
# Creacion y en colamiento del script
###
cd $WRFUTILDIR/RUNs/$NOMBRE/PostDF/
read -r -d '' QSCRIPTCMD << "EOF"
export MIEM=$(printf "%02g" $((10#$MIEMBRO_FIN-10#$MIEMBRO_INI)))
test $ARRAYID && MIEM=$(printf "%02g" $ARRAYID)
export DIRPOST=$WRFUTILDIR/RUNs/HIST/POST/
export DIRPLOT=$WRFUTILDIR/RUNs/HIST/PLOT/
export DIRPOSTDF=$BASEDIR/PostDF/
PATH=$PATH:$CONDADIR

[ ! -f "$BASEDIR/$EXPPLGIN" ] && dispararError 4 "$BASEDIR/$EXPPLGIN"
source $BASEDIR/$EXPPLGIN

for dataframe in ${POSTDFLIST[@]}
do 
	[[ ! -f "$DIRPOSTDF/$dataframe" ]] && dispararError 4 "$DIRPOSTDF/$dataframe"
	source $DIRPOSTDF/$dataframe
	source $CONDADIR/activate $entorno || $CONDADIR/conda activate $entorno
	[[ -z "$script" ]] && dispararError 4 "$script"
	python $script
	[[ $? -ne 0 ]] && dispararError 9 "$script"
	postscript    # ejecuta las rutinas definidas para finalizar la corrida
	source $CONDADIR/deactivate 
done


EOF

# Parametros de encolamiento
QDEPEND_NAME=${DEPPOSTDF}
QDEPENDTYPE=${TYPEPOSTDF}
QPROC_NAME="POSTDF_${NOMBRE}_${PASO}"
QPROC=1
QTHREADS=1
QARRAY=""
QWALLTIME=${POSTDFWALLTIME}
QEXCLU=1

# Encolar
queue


