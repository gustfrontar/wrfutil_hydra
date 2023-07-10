#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF

        Ud. deberia usar este escript de la siguiente manera:
                $0 < nuevo nombre entonro > 
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${1?"$USO"}  ${WRFUTILDIR:?"$USO"}


### PARAMETROS

export NOMBRE=$1

### CONFIGURACION

source $WRFUTILDIR/lib/errores.env
CONFIG=$WRFUTILDIR/config.env
[ -f $WRFUTILDIR/RUNs/$NOMBRE/config.env ] && CONFIG=$WRFUTILDIR/RUNs/$NOMBRE/config.env
CONFIG=$(readlink -e ${CONFIG} 2> /dev/null)
[ -z "$CONFIG" ] && dispararError 4 "config.env"
source $CONFIG $NOMBRE


##### FIN INICIALIZACION ######

EXPDIR=$BASEDIR
[ ! -d "$EXPDIR" ]  && dispararError 7 "$EXPDIR"

###  Creando Entorno 

#mkdir -p $EXPDIR/CALIB/00
#cp  $TEMPLATE/calibracion/00/*.h5 $EXPDIR/CALIB/00/
#mkdir -p $EXPDIR/CALIB/06
#cp  $TEMPLATE/calibracion/06/*.h5 $EXPDIR/CALIB/06/
#mkdir -p $EXPDIR/CALIB/12
#cp  $TEMPLATE/calibracion/12/*.h5 $EXPDIR/CALIB/12/
#mkdir -p $EXPDIR/CALIB/18
#cp  $TEMPLATE/calibracion/18/*.h5 $EXPDIR/CALIB/18/
#cp $TEMPLATE/estaciones.csv $EXPDIR/CALIB/
#cp $TEMPLATE/ciudades.csv $EXPDIR/CALIB/


echo "### Calibracion y Verificacion" >> $EXPDIR/experimento.plgin 
echo "export FILEIN_TYPE='WRF'                                                 # tipo = "nc" -> WRF o "grib" ->GFS  " >> $EXPDIR/experimento.plgin
echo 'export OUTDIR_CALIB="$BASEDIR/HIST/POST/"' >> $EXPDIR/experimento.plgin
echo "### FIN Calibracion y Verificacion " >> $EXPDIR/experimento.plgin
