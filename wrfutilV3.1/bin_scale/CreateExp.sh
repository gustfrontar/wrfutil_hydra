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

Use: createExp.sh  "TEMPLATE_NAME" "BASEDIR"

EOF
: ${2?"$USO"}

WRFUTILDIR=$(pwd)/../

### PARAMETROS

export TEMPLATE_NAME=$1    #Which template will be used to create the experiment (see folders in /conf )
export BASEDIR=$2          #Which is the root directory of the experiment.

##### FIN INICIALIZACION ######

EXPDIR=$BASEDIR
[ -d "$BASEDIR" ]  && "$BASEDIR" exists
mkdir -p $BASEDIR 


### Copio la configuracion y genero los archivos de configuracion para el experimento.
mkdir $BASEDIR/conf
cp $WRFUTILDIR/conf/${TEMPLATE_NAME}/* $BASEDIR/conf/
sed -i -e "s|__BASEDIR__|$BASEDIR|g"   $BASEDIR/conf/config.env
source $BASEDIR/conf/config.env


###  Creando Entorno 
mkdir -p $LOGDIR
mkdir -p $PROCSDIR
mkdir -p $SCALEDIR
mkdir -p $SCALEPPDIR
mkdir -p $LETKFDIR
mkdir -p $OBSDIR
cp -r $WRFUTILDIR/bin_scale $EXPDIR/bin
cp -r $WRFUTILDIR/lib $EXPDIR
mkdir -p $NAMELISTDIR
cp    $WRFUTILDIR/namelists/$MODEL_CONF_SET/* $NAMELISTDIR/
mkdir -p $HISTDIR
mkdir -p $PERTDIR

### Link the executables 

ln -s $SCALEPATH/bin/scale-rm_pp${SCALE_opt}   $SCALEPPDIR/
ln -s $SCALEPATH/bin/scale-rm_init${SCALE_opt} $SCALEPPDIR/
ln -s $SCALEPATH/bin/sno${SCALE_opt}           $SCALEPPDIR/
ln -s $SCALEPATH/bin/scale-rm${SCALE_opt}      $SCALEDIR/
ln -s $SCALEPATH/bin/sno${SCALE_opt}           $SCALEDIR/

ln -s $SCALELETKFPATH $LETKFDIR/letkf.exe
ln -s $SCALEPERTPATH $PERTDIR/pert_init_bdy.exe
