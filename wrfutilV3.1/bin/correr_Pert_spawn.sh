#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

### CONFIGURACION
source $BASEDIR/conf/config.env  #BASEDIR tiene que estar seteado como variable de entorno.

[ ! -f $BASEDIR/lib/errores.env ] && exit 1
source $BASEDIR/lib/errores.env
[ ! -f "$BASEDIR/conf/$EXPMACH" ] && dispararError 4 "$BASEDIR/conf/$EXPMACH"
source $BASEDIR/conf/$EXPMACH
[ ! -f "$BASEDIR/conf/$EXPCONF" ] && dispararError 4 "$BASEDIR/conf/$EXPCONF"
source $BASEDIR/conf/$EXPCONF
[ ! -f "$BASEDIR/conf/$EXPMODELCONF" ] && dispararError 4 "$BASEDIR/conf/$EXPMODELCONF"
source $BASEDIR/conf/$EXPMODELCONF

##### FIN INICIALIZACION ######
cd $WPSDIR 
cp $BASEDIR/bin/python/main_perturb_met_em_par.py $WPSDIR/main_perturb_met_em.py
cp $BASEDIR/bin/python/module_pert_met_em.py  $WPSDIR/module_pert_met_em.py

FECHA_INI_PASO=$(date -u -d "$FECHA_INI UTC +$(($WPS_INI_FREQ*$PASO)) seconds" +"%Y-%m-%d %T")
FECHA_INI_BDY=$(date_floor "$FECHA_INI_PASO" $INTERVALO_BDY )
FECHA_INI_BDY=$(date -u -d "$FECHA_INI_BDY" +"%Y%m%d%H%M%S")

sed -i -e "s|__BASE_PATH_ORI__|'$HISTDIR/WPS/met_em_ori/${FECHA_INI_BDY}/'|g"   $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__BASE_PATH_OUT__|'$HISTDIR/WPS/met_em/${FECHA_INI_BDY}'|g"        $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__TARGET_ENSEMBLE_SIZE__|$TARGET_ENSEMBLE_SIZE|g"                  $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__ACTUAL_ENSEMBLE_SIZE__|$ACTUAL_ENSEMBLE_SIZE|g"                  $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__PERT_AMP__|$PERT_AMP|g"                                          $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__PERT_TYPE__|'$PERT_TYPE'|g"                                      $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__VAR_LIST__|$VAR_LIST|g"                                          $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__THREADS_NUM__|$PERTTHREADS|g"                                    $WPSDIR/main_perturb_met_em.py

#read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
cd $WPSDIR
python -u ./main_perturb_met_em.py > main_perturb_met_em.out 
EC=$?
[[ $EC -ne 0 ]] && dispararError 9 "main_perturb_met_em.py"
#EOF

#QPROC_NAME=PERT_$PASO
#QNODE=$PERTNODE
#QPROC=$PERTPROC
#QTHREAD=$PERTTHREAD
#QMIEM=01
#QWALLTIME=$WPSWALLTIME
#queue 00 00
#check_proc 00 00
