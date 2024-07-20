#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

### CONFIGURACION
#Load experiment configuration
BASEDIR=$(pwd)/../
source $BASEDIR/../setup_spack.sh
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/${EXPTYPE}.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source ${BASEDIR}/lib/encolar${QUEUESYS}.sh                     # Selecciona el metodo de encolado segun el systema QUEUESYS elegido

##### FIN INICIALIZACION ######
if [ ! -z "${PJM_SHAREDTMP}" -a  ${USETMPDIR} -eq 1 ] ; then 
   echo "Using Fugaku's shared dir to run SCALE"
   SCALEDIR=${PJM_SHAREDTMP}/SCALE
fi

# Editamos el namelist.
cp $NAMELISTDIR/namelist.scale $SCALEDIR/
cp $NAMELISTDIR/namelist.sno_fcst $SCALEDIR/

cd $SCALEDIR
#Seteamos las fechas de inicio y final de los forecasts.
if [ $PASO -eq 0 ] ; then 
   FECHA_FORECAST_INI=$FECHA_INI 
   FECHA_FORECAST_END=$(date -u -d "$FECHA_INI UTC +$SPIN_UP_LENGTH seconds" +"%Y-%m-%d %T")
   FECHA_ANALISIS=$FECHA_FORECAST_END
   INTERVALO_WRF=$SPIN_UP_LENGTH
   echo " Running spin up "
   echo " Spin up will start at $FECHA_FORECAST_INI"
   echo " Spin up will end   at $FECHA_FORECAST_END"
   #Optimize the boundary data frequency for the spin-up.
   FORECAST_BDY_FREQ=$( maximum_common_divisor $SPIN_UP_LENGTH $INTERVALO_BDY  )   
else 
   echo "This is a regular data assimilation step"
   FECHA_FORECAST_INI=$(date -u -d "$FECHA_INI UTC +$(($ANALISIS_FREC*($PASO-1)+$SPIN_UP_LENGTH)) seconds" +"%Y-%m-%d %T")
   FECHA_FORECAST_END=$(date -u -d "$FECHA_FORECAST_INI UTC +$ANALISIS_WIN_FIN seconds" +"%Y-%m-%d %T")
   FECHA_ANALISIS=$(date -u -d "$FECHA_FORECAST_INI UTC +$ANALISIS_FREC seconds" +"%Y-%m-%d %T")
fi

echo "=========================="
echo "Forecast for STEP " $PASO
echo "Starting at " $FECHA_FORECAST_INI
echo "Ending   at " $FECHA_FORECAST_END

#Desglozamos las fechas
read -r IY IM ID IH Im Is  <<< $(date -u -d "$FECHA_FORECAST_INI UTC" +"%Y %m %d %H %M %S")
read -r FY FM FD FH Fm Fs  <<< $(date -u -d "$FECHA_FORECAST_END UTC" +"%Y %m %d %H %M %S")


# Editamos el namelist.
if [ $MAP_PROJ == "lambert" ];then
  MAP_PROJ_SCALE="LC"
fi

SCALEDATA="$SCALEPATH/data"
INI_SCALE="$(date -ud "$FECHA_FORECAST_INI" +'%Y,%m,%d,%H,%M,%S' )"
E_WE_LOC=$((E_WE / NPROC_X))
E_SN_LOC=$((E_SN / NPROC_Y))
SCALEDATA="$SCALEPATH/data"

sed -i -e "s|__INI_SCALE__|$INI_SCALE|g"      $SCALEDIR/namelist.scale
if [ $PASO == 0 ] ;then
sed -i -e "s|__LEN__|$SPIN_UP_LENGTH|g"     $SCALEDIR/namelist.scale
sed -i -e "s|__INT_R_SCALE__|$SPIN_UP_LENGTH|g" $SCALEDIR/namelist.scale
sed -i -e "s|__INT_F_SCALE__|$SPIN_UP_LENGTH|g" $SCALEDIR/namelist.scale
else
sed -i -e "s|__LEN__|$ANALISIS_WIN_FIN|g"     $SCALEDIR/namelist.scale
sed -i -e "s|__INT_R_SCALE__|$ANALISIS_FREC|g" $SCALEDIR/namelist.scale
sed -i -e "s|__INT_F_SCALE__|$ANALISIS_WIN_STEP|g" $SCALEDIR/namelist.scale
fi 
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
sed -i -e "s|__HISTDIR__|hist|g"              $SCALEDIR/namelist.scale
sed -i -e "s|__DT__|$DT|g"                    $SCALEDIR/namelist.scale
sed -i -e "s|__RADT__|$RADT|g"                $SCALEDIR/namelist.scale
sed -i -e "s|__DYNDT__|$DYNDT|g"              $SCALEDIR/namelist.scale
sed -i -e "s|__LNDDT__|$DT2|g"                $SCALEDIR/namelist.scale
sed -i -e "s|__OCNDT__|$DT2|g"                $SCALEDIR/namelist.scale
sed -i -e "s|__URBDT__|$DT2|g"                $SCALEDIR/namelist.scale

#Build the script to run SCALE
read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
source $BASEDIR/conf/model.conf
echo "Processing member $MIEM" 

if [ ! -z ${PJM_SHAREDTMP} -a  ${USETMPDIR} -eq 1 ] ; then
   echo "Using Fugaku's shared dir to run SCALE"
   SCALEDIR=${PJM_SHAREDTMP}/SCALE/    #SCALEDIR is redefined here
fi

cp $SCALEDIR/namelist.scale $SCALEDIR/$MIEM/namelist.scale

ln -sf $SCALEDIR/scale-rm${SCALE_opt} $SCALEDIR/$MIEM/scale-rm
ln -sf $SCALEDIR/sno${SCALE_opt} $SCALEDIR/$MIEM/sno

ln -sf  $HISTDIR/const/topo .
ln -sf  $HISTDIR/const/landuse .
mkdir -p bdy
mkdir -p init
mkdir -p hist
mkdir -p log/scale
mkdir -p log/sno

if [ $PASO -eq 0 ] ; then
   DATE_FORECAST_INI=$FECHA_INI 
   DATE_FORECAST_END=$(date -u -d "$DATE_FORECAST_INI UTC +$SPIN_UP_LENGTH seconds" +"%Y-%m-%d %T")
else 
   DATE_FORECAST_INI=$(date -u -d "$FECHA_INI UTC +$(($ANALISIS_FREC*($PASO-1)+$SPIN_UP_LENGTH)) seconds" +"%Y-%m-%d %T")
   DATE_FORECAST_END=$(date -u -d "$DATE_FORECAST_INI UTC +$ANALISIS_WIN_FIN seconds" +"%Y-%m-%d %T")
fi

if [ $WPS_CYCLE -eq 1 ] ; then
   INI_BDY_DATE=$(date_floor "$DATE_FORECAST_INI" $INTERVALO_INI_BDY )
   INI_BDY_DATE=$(date -u -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")
else
   INI_BDY_DATE=$(date_floor "$FECHA_INI" $INTERVALO_BDY )
   INI_BDY_DATE=$(date -u -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")
fi
DATE_FORECAST_INI=$(date -u -d "$DATE_FORECAST_INI" +"%Y%m%d%H%M%S")

if [ ${PASO} -eq 0 ] ; then #Assume spin-up period.
  if [ ! -z ${PJM_SHAREDTMP} -a  ${USETMPDIR} -eq 1 ] ; then
     BDY_DIR=${PJM_SHAREDTMP}/HIST/bdy/${INI_BDY_DATE}/$MIEM/    
     INI_DIR=${PJM_SHAREDTMP}/HIST/init/${INI_BDY_DATE}/$MIEM/    
  else 
     BDY_DIR=$HISTDIR/bdy/${INI_BDY_DATE}/$MIEM/
     INI_DIR=$HISTDIR/init/${INI_BDY_DATE}/$MIEM/
  fi
else
  if [ ! -z ${PJM_SHAREDTMP} -a  ${USETMPDIR} -eq 1 ] ; then
     BDY_DIR=${PJM_SHAREDTMP}/HIST/bdy/${INI_BDY_DATE}/$MIEM/    
     INI_DIR=${PJM_SHAREDTMP}/HIST/ANAL/${DATE_FORECAST_INI}/$MIEM/    
  else 
     BDY_DIR=$HISTDIR/bdy/${INI_BDY_DATE}/$MIEM/
     INI_DIR=$HISTDIR/ANAL/${DATE_FORECAST_INI}/$MIEM/
  fi
fi

#if [ ${PASO} -eq 0 ] ; then #Assume spin-up period.
#   #Optimize the boundary data frequency for the spin-up.
#   FORECAST_BDY_FREQ=$( maximum_common_divisor $SPIN_UP_LENGTH $INTERVALO_BDY )
#fi

ln -sf $BDY_DIR/*.nc ./bdy/
ln -sf $INI_DIR/*.nc ./init/

export FORT90L=${SCALE_RUNTIME_FLAGS}

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SCALE_NETCDF_C}/lib:${SCALE_NETCDF_F}/lib:${SCALE_PNETCDF}/lib:${SCALE_HDF}/lib:$LD_LIBRARY_PATH

echo "Running scale for member $MIEM"
$MPIEXE ./scale-rm namelist.scale
ERROR=$(( $ERROR + $? ))

#test=$(tail -n1 $WRFDIR/$MIEM/rsl.error.0000 | grep SUCCESS ) && res="OK"
#mv rsl.error.0000 ./wrf_${PASO}_${MIEM}.log

EOF

#Node / core distribution parameters
QNODE=$SCALENODE
QPROC=$SCALEPROC
TPROC=$SCALEPROC
QTHREAD=$ITHREAD
QWALLTIME=$SCALEWALLTIME
QPROC_NAME=GUESS_${PASO}
QCONF=${EXPTYPE}.conf
QWORKPATH=$SCALEDIR

#Execute the job 
echo "Tiempo en correr el scale"
queue $MIEMBRO_INI $MIEMBRO_FIN 
time check_proc $MIEMBRO_INI $MIEMBRO_FIN

if [ $PASO -eq 0  ] ; then  #Copy the spin up output as the analysis for the next cycle.
   for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
       OUTPUTPATH="$HISTDIR/ANAL/$(date -u -d "$FECHA_ANALISIS UTC" +"%Y%m%d%H%M%S")/$MIEM"
       mkdir -p $OUTPUTPATH
       echo "Copying file $SCALEDIR/$MIEM/init/init_$(date -u -d "$FECHA_ANALISIS UTC" +"%Y%m%d-%H%M%S.000" )"
       mv $SCALEDIR/$MIEM/init/init_$(date -u -d "$FECHA_ANALISIS UTC" +"%Y%m%d-%H%M%S.000")* $OUTPUTPATH/
       mv $SCALEDIR/$MIEM/log/scale/*                                                    $OUTPUTPATH/
       mv $SCALEDIR/$MIEM/namelist*                                                      $OUTPUTPATH/
   done

else 
  #Copy the guess files corresponding to the analysis time.
  if [[ ! -z "$GUARDOGUESS" ]] && [[ $GUARDOGUESS -eq 1 ]] ; then
     for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
       OUTPUTPATH="$HISTDIR/GUES/$(date -u -d "$FECHA_ANALISIS UTC" +"%Y%m%d%H%M%S")/$MIEM"
       mkdir -p $OUTPUTPATH
       echo "Copying file $SCALEDIR/$MIEM/init/init_$(date -u -d "$FECHA_ANALISIS UTC" +"%Y%m%d-%H%M%S.000" )"
       cp $SCALEDIR/$MIEM/init/init_$(date -u -d "$FECHA_ANALISIS UTC" +"%Y%m%d-%H%M%S.000")* $OUTPUTPATH/
       mv $SCALEDIR/$MIEM/log/scale/*                                                   $OUTPUTPATH/
       mv $SCALEDIR/$MIEM/namelist*                                                     $OUTPUTPATH/
     done
  fi
fi







