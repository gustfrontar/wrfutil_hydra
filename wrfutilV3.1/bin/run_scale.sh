#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

### CONFIGURACION
BASEDIR=$(pwd)/../
source $BASEDIR/../setup_spack.sh
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/${EXPTYPE}.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source ${BASEDIR}/lib/encolar${QUEUESYS}.sh                     # Selecciona el metodo de encolado segun el systema QUEUESYS elegido
##### FIN INICIALIZACION ######

cd $SCALEDIR
#Seteamos las fechas de inicio y final de los forecasts.
FECHA_FORECAST_INI=$(date -u -d "$FECHA_INI UTC +$(($FORECAST_INI_FREQ*$PASO)) seconds" +"%Y-%m-%d %T")
FECHA_FORECAST_END=$(date -u -d "$FECHA_INI UTC +$((($FORECAST_INI_FREQ*$PASO)+$FORECAST_LEAD_TIME )) seconds" +"%Y-%m-%d %T")

echo "=========================="
echo "Forecast for STEP " $PASO
echo "Starting at " $FECHA_FORECAST_INI
echo "Ending   at " $FECHA_FORECAST_END

#Desglozamos las fechas
read -r IY IM ID IH Im Is  <<< $(date -u -d "$FECHA_FORECAST_INI UTC" +"%Y %m %d %H %M %S")
read -r FY FM FD FH Fm Fs  <<< $(date -u -d "$FECHA_FORECAST_END UTC" +"%Y %m %d %H %M %S")

# Editamos el namelist.
cp $NAMELISTDIR/namelist.scale $SCALEDIR/
cp $NAMELISTDIR/namelist.sno_fcst $SCALEDIR/

#Edit the namelist from the template

if [ $MAP_PROJ == "lambert" ];then
  MAP_PROJ_SCALE="LC"
fi

SCALEDATA="$SCALEPATH/data"
INI_SCALE="$(date -ud "$FECHA_INI" +'%Y,%m,%d,%H,%M,%S' )"
E_WE_LOC=$((E_WE / NPROC_X))
E_SN_LOC=$((E_SN / NPROC_Y))
SCALEDATA="$SCALEPATH/data"

sed -i -e "s|__INI_SCALE__|$INI_SCALE|g"      $SCALEDIR/namelist.scale
sed -i -e "s|__LEN__|$FORECAST_LEAD_TIME|g"   $SCALEDIR/namelist.scale
sed -i -e "s|__INT_R_SCALE__|$INT_R_SCALE|g"  $SCALEDIR/namelist.scale
sed -i -e "s|__INT_F_SCALE__|$INT_F_SCALE|g"  $SCALEDIR/namelist.scale
sed -i -e "s|__E_WE_LOC__|$E_WE_LOC|g"        $SCALEDIR/namelist.scale
sed -i -e "s|__E_SN_LOC__|$E_SN_LOC|g"        $SCALEDIR/namelist.scale
sed -i -e "s|__NPROC_X__|$NPROC_X|g"          $SCALEDIR/namelist.scale
sed -i -e "s|__NPROC_Y__|$NPROC_Y|g"          $SCALEDIR/namelist.scale
sed -i -e "s|__DX__|$DX|g"                    $SCALEDIR/namelist.scale
sed -i -e "s|__DY__|$DY|g"                    $SCALEDIR/namelist.scale
sed -i -e "s|__MAP_PROJ__|$MAP_PROJ_SCALE|g"  $SCALEDIR/namelist.scale
sed -i -e "s|__REF_LAT__|$REF_LAT|g"          $SCALEDIR/namelist.scale
sed -i -e "s|__REF_LON__|$REF_LON|g"          $SCALEDIR/namelist.scale
sed -i -e "s|__TRUELAT1__|$TRUELAT1|g"        $SCALEDIR/namelist.scale
sed -i -e "s|__TRUELAT2__|$TRUELAT2|g"        $SCALEDIR/namelist.scale
sed -i -e "s|__SCALEDATA__|$SCALEDATA|g"      $SCALEDIR/namelist.scale
sed -i -e "s|__HISTDIR__|fcst|g"              $SCALEDIR/namelist.scale
sed -i -e "s|__DT__|$DT|g"                    $SCALEDIR/namelist.scale
sed -i -e "s|__RADT__|$RADT|g"                $SCALEDIR/namelist.scale
sed -i -e "s|__DYNDT__|$DYNDT|g"              $SCALEDIR/namelist.scale
sed -i -e "s|__LNDDT__|$DT2|g"                $SCALEDIR/namelist.scale
sed -i -e "s|__OCNDT__|$DT2|g"                $SCALEDIR/namelist.scale
sed -i -e "s|__URBDT__|$DT2|g"                $SCALEDIR/namelist.scale

###
# Create the script that will run scale for each ensemble member.
###
echo "Corriendo el scale para el paso " $PASO

read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
source $BASEDIR/conf/model.conf
echo "Procesando Miembro $MIEM" 

mkdir -p $SCALEDIR/$MIEM
cd $SCALEDIR/$MIEM

ln -sf $SCALEDIR/scale-rm${SCALE_opt} ./scale-rm
ln -sf $SCALEDIR/sno${SCALE_opt} ./sno
cp $SCALEDIR/namelist.scale . 
cp $SCALEDIR/namelist.sno_fcst . 

ln -sf  $HISTDIR/const/topo .
ln -sf  $HISTDIR/const/landuse .
mkdir -p bdy
mkdir -p init
mkdir -p fcst
mkdir -p log/scale
mkdir -p log/sno
DIR_FECHA_INI=$(date -u -d "$FECHA_FORECAST_INI UTC" +"%Y%m%d%H%M%S")
ln -sf $HISTDIR/init/${DIR_FECHA_INI}/$MIEM/*.nc ./init/
ln -sf $HISTDIR/bdy/${DIR_FECHA_INI}/$MIEM/*.nc ./bdy/

export FORT90L=${SCALE_RUNTIME_FLAGS}

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SCALE_NETCDF_C}/lib:${SCALE_NETCDF_F}/lib:${SCALE_PNETCDF}/lib:${SCALE_HDF}/lib:$LD_LIBRARY_PATH

echo "Running scale for member $MIEM"
$MPIEXE ./scale-rm namelist.scale
ERROR=$(( $ERROR + $? ))

echo "Running SNO for member $MIEM"
$MPIEXESERIAL ./sno namelist.sno_fcst
ERROR=$(( $ERROR + $? ))

EOF

#Node / core distribution parameters
QNODE=$SCALENODE
QPROC=$SCALEPROC
TPROC=$SCALEPROC
QTHREAD=$ITHREAD
QWALLTIME=$SCALEWALLTIME
QPROC_NAME=FCST_${PASO}
QCONF=${EXPTYPE}.conf
QWORKPATH=$SCALEDIR

#Execute the job 
echo "Time taken by scale"
queue $MIEMBRO_INI $MIEMBRO_FIN 
time check_proc $MIEMBRO_INI $MIEMBRO_FIN

#Copy wrfout files to its final destionation.
for QMIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
   OUTPUTPATH="$HISTDIR/FCST/${DIR_FECHA_INI}/${QMIEM}/"
   mkdir -p $OUTPUTPATH
   echo "Copying file $SCALEDIR/$QMIEM/fcst/history" 
   mv $SCALEDIR/$QMIEM/fcst/*.nc        $OUTPUTPATH
   mv $SCALEDIR/$QMIEM/log/scale        $OUTPUTPATH
   mv $SCALEDIR/$QMIEM/namelist*        $OUTPUTPATH
done


