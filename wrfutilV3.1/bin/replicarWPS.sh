#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF
        Ud. deberia usar este escript de la siguiente manera:
                $0 <nuevo nombre entonro WPS > </path/al/directorio/a/replicar/> <entorno.env>
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF

: ${1:?"$USO"} ${2:?"$USO"} ${3:?"$USO"} ${WRFUTILDIR:?"$USO"}

[ ! -f $3 ] && echo "archivo de configuracion de entorno inexistente" && exit 1 

### PARAMETROS

WPSENV=$1
WPSPATH=$2
ENTORNO=$3

### CONFIGURACION

source $WRFUTILDIR/lib/errores.env
CONFIG=$WRFUTILDIR/config.env
#[ -f $WRFUTILDIR/RUNs/$NOMBRE/config.env ] && CONFIG=$WRFUTILDIR/RUNs/$NOMBRE/config.env 
#CONFIG=$(readlink -e ${CONFIG} 2> /dev/null)
[ -z "$CONFIG" ] && dispararError 4 "config.env"
source $CONFIG

##### FIN INICIALIZACION ######


## Verifica que no exista la replica
[ -f "$WRFUTILDIR/WPSs/${WPSENV}.tar" ] && dispararError 3 $WPSENV
## Verifica que exista el original
[ ! -d "$WPSPATH" ] && dispararError 4 $WPSPATH

### WPS 

mkdir -p $TMPDIR/WPS
[ $? -ne 0 ] && dispararError 2 "$TMPDIR/WPS"
cd $TMPDIR/WPS
ln -s  $WPSPATH/geogrid.exe
ln -s  $WPSPATH/geogrid
ln -s  $WPSPATH/ungrib.exe
ln -s  $WPSPATH/ungrib
ln -s  $WPSPATH/metgrid.exe
ln -s  $WPSPATH/metgrid
ln -s  $WPSPATH/link_grib.csh
ln -sf ../Vtable Vtable
#cp $WPSPATH/../envvars.sh .
cp $ENTORNO .
#tar cf ../${WPSENV}.tar --exclude="$geos" .
tar cf ../${WPSENV}.tar  .
mv ../${WPSENV}.tar $WRFUTILDIR/WPSs/
cd ..
rm -r $TMPDIR/WPS

