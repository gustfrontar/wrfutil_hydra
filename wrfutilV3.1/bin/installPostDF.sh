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

mkdir -p $EXPDIR/PostDF/
#cp $TEMPLATE/estaciones.csv $EXPDIR/PostDF/
#cp $TEMPLATE/ciudades.csv $EXPDIR/PostDF/

echo "DEPPOSTDF='' ">>$EXPDIR/experimento.dep
echo "TYPEPOSTDF'' ">>$EXPDIR/experimento.dep

echo "### INI PostDF  -- Armado de DataFrames" >> $EXPDIR/experimento.plgin 
echo "export POSTDFWALLTIME=" >> $EXPDIR/experimento.plgin
echo "export POSTDFLIST=( 'superficie' 'nubosidad' )" >> $EXPDIR/experimento.plgin
echo "### FIN PostDF  -- Armado de DataFrames " >> $EXPDIR/experimento.plgin
