#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF
        Ud. deberia usar este escript de la siguiente manera:
                $0 <nuevo nombre entonro CRTM > </path/al/directorio/a/replicar/CRTM >  </path/al/directorio/a/replicar/ISS > <entorno.env>
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF

: ${1:?"$USO"} ${2:?"$USO"} ${3:?"$USO"} ${4:?"$USO"} ${WRFUTILDIR:?"$USO"}

[ ! -f $4 ] && echo "archivo de configuracion de entorno inexistente" && exit 1 

### PARAMETROS

CRTMENV=$1
CRTMPATH=$2
ISSPATH=$3
ENTORNO=$4

### CONFIGURACION

source $WRFUTILDIR/lib/errores.env
CONFIG=$WRFUTILDIR/config.env
#[ -f $WRFUTILDIR/RUNs/$NOMBRE/config.env ] && CONFIG=$WRFUTILDIR/RUNs/$NOMBRE/config.env 
#CONFIG=$(readlink -e ${CONFIG} 2> /dev/null)
[ -z "$CONFIG" ] && dispararError 4 "config.env"
source $CONFIG

##### FIN INICIALIZACION ######


## Verifica que no exista la replica
[ -f "$WRFUTILDIR/CRTMs/${CRTMENV}.tar" ] && dispararError 3 $CRTMENV
## Verifica que exista el original
[ ! -d "$CRTMPATH" ] && dispararError 4 $CRTMPATH

### CRTM 

mkdir -p $TMPDIR/CRTM
[ $? -ne 0 ] && dispararError 2 "$TMPDIR/CRTM"
cd $TMPDIR/CRTM

fixpath="${CRTMPATH}/fix"

ln -sf ${fixpath}/CloudCoeff/Little_Endian/CloudCoeff.bin CloudCoeff.bin
ln -sf ${fixpath}/AerosolCoeff/Little_Endian/AerosolCoeff.bin AerosolCoeff.bin
ln -sf ${fixpath}/EmisCoeff/IR_Land/SEcategory/Little_Endian/IGBP.IRland.EmisCoeff.bin IGBP.IRland.EmisCoeff.bin
ln -sf ${fixpath}/EmisCoeff/IR_Snow/SEcategory/Little_Endian/NPOESS.IRsnow.EmisCoeff.bin NPOESS.IRsnow.EmisCoeff.bin
ln -sf ${fixpath}/EmisCoeff/IR_Ice/SEcategory/Little_Endian/NPOESS.IRice.EmisCoeff.bin NPOESS.IRice.EmisCoeff.bin
ln -sf ${fixpath}/EmisCoeff/IR_Water/Little_Endian/Nalli.IRwater.EmisCoeff.bin Nalli.IRwater.EmisCoeff.bin
ln -sf ${fixpath}/SpcCoeff/Little_Endian/abi_g16.SpcCoeff.bin abi_g16.SpcCoeff.bin
ln -sf ${fixpath}/TauCoeff/ODPS/Little_Endian/abi_g16.TauCoeff.bin abi_g16.TauCoeff.bin
ln -sf $ISSPATH/run_CRTM
ln -sf $ISSPATH/grafica_CRTM.py
ln -sf $ISSPATH/LOGO_SMN.png
ln -sf $ISSPATH/SVGAWVX_TEMP.cpt
ln -sf $ISSPATH/cpt_convert.py
ln -sf $ISSPATH/smn_topes.cpt
cp $ENTORNO .

tar cf ../${CRTMENV}.tar  .
mkdir -p $WRFUTILDIR/CRTMs/
mv ../${CRTMENV}.tar $WRFUTILDIR/CRTMs/
cd ..
rm -r $TMPDIR/CRTM

