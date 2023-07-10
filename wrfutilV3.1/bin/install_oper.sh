#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF

        Ud. deberia usar este escript de la siguiente manera:
                $0 < entorno WPS > < entorno WRF > < entorno CRTM > < entorno LETKF > 
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${1?"$USO"}  ${2?"$USO"}  ${3?"$USO"} ${4?"$USO"} ${WRFUTILDIR:?"$USO"}

WPS=$1
WRF=$2
CRTM=$3
LETKF=$4

###
# CalibGFS
###
NOMBRE="CalibGFS"
$WRFUTILDIR/bin/crearExpFCST.sh $NOMBRE $WPS $WRF 
$WRFUTILDIR/bin/installCalibracion.sh $NOMBRE
$WRFUTILDIR/bin/installPostDF.sh $NOMBRE

###
# CalibGEFS
###
NOMBRE="CalibGEFS"
$WRFUTILDIR/bin/crearExpFCST.sh $NOMBRE $WPS $WRF 
$WRFUTILDIR/bin/installCalibracion.sh $NOMBRE
$WRFUTILDIR/bin/installPostDF.sh $NOMBRE

###
# deterministico
###
NOMBRE="deterministico"
$WRFUTILDIR/bin/crearExpFCST.sh $NOMBRE $WPS $WRF 
$WRFUTILDIR/bin/installCRTM.sh $NOMBRE $CRTM
$WRFUTILDIR/bin/installCalibracion.sh $NOMBRE
$WRFUTILDIR/bin/installPostDF.sh $NOMBRE


###
# ensamble
###
NOMBRE="ensamble"
$WRFUTILDIR/bin/crearExpFCST.sh $NOMBRE $WPS $WRF 
$WRFUTILDIR/bin/installCalibracion.sh $NOMBRE
$WRFUTILDIR/bin/installPostDF.sh $NOMBRE


###
# asimialcion
###
NOMBRE="asimilacion"
$WRFUTILDIR/bin/crearExpASIM.sh $NOMBRE $WPS $WRF $LETKF

###
# Configuracion de los Experimentos
###
cd $WRFUTILDIR/RUNs/
git clone https://gitlab.smn.gov.ar/msacco/wrfutil_oper_conf.git
cp -R wrfutil_oper_conf/* .
rm -r wrfutil_oper_conf
