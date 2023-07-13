#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF

        Ud. deberia usar este escript de la siguiente manera:
                $0 < nuevo nombre entonro > < entorno CRTM > 
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${1?"$USO"}  ${2?"$USO"}   ${WRFUTILDIR:?"$USO"}


### PARAMETROS

export NOMBRE=$1

### CONFIGURACION

source $WRFUTILDIR/lib/errores.env
CONFIG=$WRFUTILDIR/config.env
[ -f $WRFUTILDIR/RUNs/$NOMBRE/config.env ] && CONFIG=$WRFUTILDIR/RUNs/$NOMBRE/config.env
CONFIG=$(readlink -e ${CONFIG} 2> /dev/null)
[ -z "$CONFIG" ] && dispararError 4 "config.env"
source $CONFIG $NOMBRE

[ -f $WRFUTILDIR/CRTMs/$2 ] || dispararError 7 "$2 ; Elija uno de los siguientes \n $(ls  $WRFUTILDIR/CRTMs/)" 


##### FIN INICIALIZACION ######

EXPDIR=$BASEDIR
[ ! -d "$EXPDIR" ]  && dispararError 7 "$EXPDIR"

###  Creando Entorno 

mkdir -p $EXPDIR/CRTM || dispararError 2 "$EXPDIR/CRTM"
cd $EXPDIR/CRTM
cp  $WRFUTILDIR/CRTMs/$2 $EXPDIR/CRTM/CRTM.tar

echo '### CRTM' >> $EXPDIR/experimento.plgin
echo 'export CRTMPROC=1                                                       # Cantidad de procesos total para correr el CRTM' >> $EXPDIR/experimento.plgin
echo 'export CRTMTHREADS=32                                                   # Cantidad de threads por porceso del CRTM' >> $EXPDIR/experimento.plgin
echo 'export CRTMPPN=1                                                        # Cantidad de procesos que se correran por nodo  del CRTM' >> $EXPDIR/experimento.plgin
echo 'export CRTMWALLTIME="04:00:00"                                          # Tiempo estimado en HH:MM:SS que el proceso durara como maximo' >> $EXPDIR/experimento.plgin
echo '### FINCRTM' >> $EXPDIR/experimento.plgin


