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
cd $WRFUTILDIR/RUNs/$NOMBRE/CALIB/
read -r -d '' QSCRIPTCMD << "EOF"
export MIEM=$(printf "%02g" $((10#$MIEMBRO_FIN-10#$MIEMBRO_INI)))
test $ARRAYID && MIEM=$(printf "%02g" $ARRAYID)
export DIRPOST=$WRFUTILDIR/RUNs/HIST/POST/
export DIRPLOT=$WRFUTILDIR/RUNs/HIST/PLOT/
export DIRVERIF=$BASEDIR/CALIB/
PATH=$PATH:$CONDADIR
source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil


read -r IY IM ID IH Im  <<< "$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN)) minutes" +"%Y %m %d %H %M")"
if [[ $FILEIN_TYPE == "WRF" ]]
then
        export DIRIN=$BASEDIR/HIST/POST/${IY}${IM}${ID}_${IH}${Im}00/$NOMBRE/
elif [[ $FILEIN_TYPE == "GEFS" ]]
then
        export DIRIN=$BASEDIR/HIST/GFS_CAL/${IY}${IM}${ID}_${IH}/ENS/
else
        export DIRIN=$BASEDIR/HIST/GFS_CAL/${IY}${IM}${ID}_${IH}/DET/
fi


python $WRFUTILDIR/bin/python/estima_QR.py
[[ $? -ne 0 ]] && dispararError 9 "estima_QR.py"
python $WRFUTILDIR/bin/python/calib_temp.py
[[ $? -ne 0 ]] && dispararError 9 "calib_temp.py"

python $WRFUTILDIR/bin/python/estima_QR_minmax.py
[[ $? -ne 0 ]] && dispararError 9 "estima_QRminmax.py"
python $WRFUTILDIR/bin/python/calib_minmax.py
[[ $? -ne 0 ]] && dispararError 9 "calib_minmax.py"

rm ${OUTDIR_CALIB}/${IY}${IM}${ID}_${IH}${Im}00/$NOMBRE/*.npy

python $WRFUTILDIR/bin/python/estima_QR_viento.py
[[ $? -ne 0 ]] && dispararError 9 "estima_QR_viento.py"
python $WRFUTILDIR/bin/python/calib_viento.py
[[ $? -ne 0 ]] && dispararError 9 "calib_viento.py"

#·sed -i -e "/export FECHA_CALIB=/c\\export FECHA_CALIB=\'$FECHA_INI\'" $BASEDIR/experimento.plgin
#sed -i -e "/export CICLO_CALIB=/c\\export CICLO_CALIB=$CICLO" $BASEDIR/experimento.plgin


EOF

# Parametros de encolamiento
QDEPEND_NAME=${DEPCALIB}
QDEPENDTYPE=${TYPECALIB}
QPROC_NAME="CALIB_${NOMBRE}_${PASO}"
QPROC=1
QTHREADS=1
QARRAY=""
QWALLTIME=${PLOTWALLTIME}
QEXCLU=1

# Encolar
queue


