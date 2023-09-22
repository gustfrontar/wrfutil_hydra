#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
# Adaptado al cluster hydra en CIMA
# Fecha: 06/2023
#############

###  Leyenda de USO

read -r -d '' USO << EOF

        Ud. deberia usar este escript de la siguiente manera:
                $0 < nuevo nombre entorno >
EOF
: ${1?"$USO"}  

WRFUTILDIR=$(pwd)/../

### PARAMETROS

export NOMBRE=$1

### CONFIGURACION

source $WRFUTILDIR/lib/errores.env
CONFIG=$WRFUTILDIR/conf/${NOMBRE}/config.env
[ ! -e "$CONFIG" ] && dispararError 4 "Error: No encontre config.env"
source $CONFIG

[ -f $WPSPATH ] || dispararError 7 "$WPSPATH ; Error: No se encontro el WPS en el path indicado" 
[ -f $WRFPATH ] || dispararError 7 "$WRFPATH ; Error: No se encontro el WRF en el path indicado"


##### FIN INICIALIZACION ######

EXPDIR=$BASEDIR
[ -d "$EXPDIR" ]  && dispararError 6 "$EXPDIR"
mkdir -p $EXPDIR || dispararError 5 "$EXPDIR"
cd $EXPDIR
[ -d "$TMPDIR" ]  && dispararError 6 "$TMPDIR"
mkdir -p $TMPDIR || dispararError 5 "$TMPDIR"
cd $TMPDIR


### Copio la configuracion y genero los archivos de configuracion para el experimento. 
mkdir $EXPDIR/conf 
cp $WRFUTILDIR/conf/${NOMBRE}/*  $EXPDIR/conf/ 


###  Creando Entorno 
ln -s $BDYPATH $BDYDIR
mkdir -p $LOGDIR
mkdir -p $WPSDIR
mkdir -p $PROCSDIR
cp    $WPSPATH $WPSDIR/wps.tar
cp    $SPAWNPATH $WPSDIR/spawn.tar
cp    $WRFUTILDIR/vtables/* $WPSDIR/
mkdir -p $WRFDIR
cp    $WRFPATH $WRFDIR/wrf.tar
cp    $SPAWNPATH $WRFDIR/spawn.tar
cp -r $WRFUTILDIR/bin $EXPDIR
cp -r $WRFUTILDIR/lib $EXPDIR
mkdir -p $NAMELISTDIR
cp    $WRFUTILDIR/namelists/$NOMBRE/* $NAMELISTDIR/
mkdir -p $HISTDIR

