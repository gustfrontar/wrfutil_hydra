#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF
        Ud. deberia usar este escript de la siguiente manera:
                $0 <nuevo nombre entonro LETKF > </path/al/letkkf.exe> </directorio/base/del/WRFDA/> <entorno.env>
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF

: ${1:?"$USO"} ${2:?"$USO"} ${3:?"$USO"} ${4:?"$USO"}  ${WRFUTILDIR:?"$USO"}

[ ! -f $4 ] && echo "archivo de configuracion de entorno inexistente" && exit 1 

### PARAMETROS

LETKFENV=$1
LETKFPATH=$2
WRFDA=$3
ENTORNO=$4

### CONFIGURACION

source $WRFUTILDIR/lib/errores.env
CONFIG=$WRFUTILDIR/config.env
[ -z "$CONFIG" ] && dispararError 4 "config.env"
source $CONFIG

##### FIN INICIALIZACION ######


## Verifica que no exista la replica
[ -f "$WRFUTILDIR/LETKFs/${LETKFENV}.tar" ] && dispararError 3 $LETKFENV
## Verifica que exista el original
[ ! -f "$LETKFPATH" ] && dispararError 4 $LETKFPATH
[ ! -d "$WRFDA" ] && dispararError 4 $WRFDA

### LETKF 

mkdir -p $TMPDIR/LETKF
[ $? -ne 0 ] && dispararError 2 "$TMPDIR/LETKF"
cd $TMPDIR/LETKF
cp  $LETKFPATH .
ln -s  $WRFDA/var/da/da_update_bc.exe .
cp $ENTORNO .

tar cf ../${LETKFENV}.tar  .
mv ../${LETKFENV}.tar $WRFUTILDIR/LETKFs/
cd ..
rm -r $TMPDIR/LETKF

