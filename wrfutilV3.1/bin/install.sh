#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF

        Ud. deberia usar este escript de la siguiente manera:
                $0 <path absoluto a la base de datos del goegrid> <path absoluto al repositorio de GFS> <path absoluto al repositorio de GFS para calibracion>
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${1:?"$USO"} ${2:?"$USO"} ${3:?"$USO"} ${WRFUTILDIR:?"$USO"}

CONDADIR=$(which conda 2> /dev/null)
CONDADIR=${CONDADIR%/*}
[[ -z $CONDADIR ]] && echo "NO se encuentra el comando conda. NO se puede continuar con la instalacion" && exit 1

mkdir $WRFUTILDIR/RUNs
mkdir $WRFUTILDIR/tmp
mkdir $WRFUTILDIR/RUNs/HIST
ln -s $2  $WRFUTILDIR/RUNs/HIST/GFS
ln -s $2  $WRFUTILDIR/RUNs/GFS
ln -s $3  $WRFUTILDIR/RUNs/HIST/GFS_CAL
mkdir $WRFUTILDIR/RUNs/HIST/LOGS
mkdir $WRFUTILDIR/RUNs/HIST/OBS

mkdir -p $WRFUTILDIR/RUNs/HIST/WPS/last
mkdir $WRFUTILDIR/RUNs/PROCS
mkdir $WRFUTILDIR/WPSs
mkdir $WRFUTILDIR/WRFs
mkdir $WRFUTILDIR/LETKFs
mkdir $WRFUTILDIR/CRTMs
ln -s $1 $WRFUTILDIR/RUNs/geog
$CONDADIR/conda env create -f $WRFUTILDIR/templates/wrfutil.yml
cp $WRFUTILDIR/templates/config.env $WRFUTILDIR/
cp $WRFUTILDIR/templates/limpiar.conf  $WRFUTILDIR/RUNs/HIST/
sed -i -e "s|__NOMBRE__|\$1|g" $WRFUTILDIR/config.env
sed -i -e "s|__GEOG__|$1|g" $WRFUTILDIR/config.env
sed -i -e "s|__CONDADIR__|$CONDADIR|g" $WRFUTILDIR/config.env
sed -i -e "s|__WRFUTILDIR__|$WRFUTILDIR|g" $WRFUTILDIR/config.env
cd $WRFUTILDIR/bin/python/OBS/airs_deco/
./compile_HM.sh
cd $WRFUTILDIR/bin/python/OBS/radar_so/utils
./make_gfortran.sh
