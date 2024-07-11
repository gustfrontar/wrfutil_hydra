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
cp $NAMELISTDIR/namelist.scale_init $SCALEDIR/
cp $NAMELISTDIR/namelist.scale $SCALEDIR/
cp $NAMELISTDIR/namelist.ncinput    $SCALEDIR/

DATE_INI_BDY_f=$(date_floor "$FECHA_FORECAST_INI" $INTERVALO_BDY )   #Get the closest prior date in which BDY data is available.
DATE_INI_BDY=$(date -u -d "$DATE_INI_BDY_f UTC" +"%Y-%m-%d_%T" ) #Cambio al formato WPS
DATE_END_BDY_f=$(date_ceil  "$FECHA_FORECAST_END" $INTERVALO_BDY )   #Get the closest posterior date in which BDY data is available. 
DATE_END_BDY=$(date -u -d "$DATE_END_BDY_f UTC" +"%Y-%m-%d_%T" ) #Cambio al formato WPS

#Edit the namelist from the template

E_WE_LOC=$((E_WE / NPROC_X))
E_SN_LOC=$((E_SN / NPROC_Y))
DATE_INI_SCALE="$(date -ud "$DATE_INI_BDY_f" +'%Y,%m,%d,%H,%M,%S' )"
if [ $MAP_PROJ == "lambert" ];then
  MAP_PROJ_SCALE="LC"
fi

SCALEDATA="$SCALEPATH/data"
LEN_BDY=$(( $(date -ud "$DATE_END_BDY_f" +%s) - $(date -ud "$DATE_INI_BDY_f" +%s) ))
NFILES_BDY=$(( LEN_BDY / INTERVALO_BDY + 1 ))

sed -i -e "s|__FECHA_INI__|$DATE_INI_SCALE|g" $SCALEDIR/namelist.scale_init
sed -i -e "s|__INTERVALO__|$INTERVALO_BDY|g"  $SCALEDIR/namelist.scale_init
sed -i -e "s|__E_WE_LOC__|$E_WE_LOC|g"        $SCALEDIR/namelist.scale_init
sed -i -e "s|__E_SN_LOC__|$E_SN_LOC|g"        $SCALEDIR/namelist.scale_init
sed -i -e "s|__NPROC_X__|$NPROC_X|g"          $SCALEDIR/namelist.scale_init
sed -i -e "s|__NPROC_Y__|$NPROC_Y|g"          $SCALEDIR/namelist.scale_init
sed -i -e "s|__DX__|$DX|g"                    $SCALEDIR/namelist.scale_init
sed -i -e "s|__DY__|$DY|g"                    $SCALEDIR/namelist.scale_init
sed -i -e "s|__MAP_PROJ__|$MAP_PROJ_SCALE|g"  $SCALEDIR/namelist.scale_init
sed -i -e "s|__REF_LAT__|$REF_LAT|g"          $SCALEDIR/namelist.scale_init
sed -i -e "s|__REF_LON__|$REF_LON|g"          $SCALEDIR/namelist.scale_init
sed -i -e "s|__TRUELAT1__|$TRUELAT1|g"        $SCALEDIR/namelist.scale_init
sed -i -e "s|__TRUELAT2__|$TRUELAT2|g"        $SCALEDIR/namelist.scale_init

sed -i -e "s|__SCALEDATA__|$SCALEDATA|g"      $SCALEDIR/namelist.scale_init
sed -i -e "s|__NFILES_BDY__|$NFILES_BDY|g"    $SCALEDIR/namelist.scale_init

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


###
# Create the script that will run scale_init and scale for each ensemble member.
###
echo "Corriendo el scale init and scale para el paso " $PASO

read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
source $BASEDIR/conf/model.conf
echo "Procesando Miembro $MIEM" 

mkdir -p $SCALEDIR/$MIEM
cd $SCALEDIR/$MIEM

ln -sf $SCALEDIR/scale-rm_init${SCALE_opt}   $SCALEDIR/$MIEM/scale-rm_init
ln -sf $SCALEDIR/scale-rm${SCALE_opt} ./scale-rm
cp $SCALEDIR/namelist.scale_init $SCALEDIR/$MIEM/
cp $SCALEDIR/namelist.ncinput    $SCALEDIR/$MIEM/
cp $SCALEDIR/namelist.scale . 

ln -sf  $HISTDIR/const/topo .
ln -sf  $HISTDIR/const/landuse .
mkdir -p bdyorg
mkdir -p bdy
mkdir -p init
mkdir -p fcst
mkdir -p log/scale_init
mkdir -p log/scale

if [ $WPS_CYCLE -eq 1 ] ; then
   FECHA_INI_PASO=$(date -u -d "$FECHA_INI UTC +$(($WPS_INI_FREQ*$PASO)) seconds" +"%Y-%m-%d %T")
   FECHA_END_PASO=$(date -u -d "$FECHA_INI UTC +$((($WPS_INI_FREQ*$PASO)+$WPS_LEAD_TIME )) seconds" +"%Y-%m-%d %T")
else 
   FECHA_INI_PASO=$FECHA_INI
   FECHA_END_PASO=$FECHA_FIN
fi
DATE_INI_BDY=$(date_floor "$FECHA_INI_PASO" $INTERVALO_INI_BDY )  #Get the closest prior date in which BDY data is available.
DATE_END_BDY=$(date_ceil  "$FECHA_END_PASO" $INTERVALO_BDY )      #Get the closest posterior date in which BDY data is available.

export FORT90L=${SCALE_RUNTIME_FLAGS}

if [ $WPS_DATA_SOURCE == 'GFS' ] ; then 

   echo "currently not supported"
   exit 1

   BDYBASE=$BDYPATH/gefs.$(date -d "$DATE_INI_BDY" +"%Y%m%d")/$(date -d "$DATE_INI_BDY" +"%H")/$BDYPREFIX/$MIEM/
   echo "Selected data source is GFS, we will run ungrib to decode the data"
   echo "I'm lloking for the BDY files in the folder  $BDYBASE"
   ln -sf $WPSDIR/$BDYVTABLE ./Vtable
   echo "Decoding gribs in : $WPSDIR/$MIEM"
   #Linkeo los gribs
   cd $WPSDIR/$MIEM
   ./link_grib.csh $BDYBASE/$BDYSEARCHSTRING  
   ERROR=$(( $ERROR + $? ))

   #Corro el Ungrib 
   OMP_NUM_THREADS=1
   $MPIEXESERIAL ./ungrib.exe > ungrib.log
   ERROR=$(( $ERROR + $? ))

elif [ $WPS_DATA_SOURCE == 'WRF' ]  ; then

   BDYBASE=$BDYPATH/$(date -d "$DATE_INI_BDY" +"%Y%m%d%H%M%S")/$MIEM/
   echo "Selected data source is WRF, we will use wrfout directly in scale-rm_init"

   #Need to loop over wrfout files.
   CDATE=$DATE_INI_BDY 
   cnt=0
   while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $(date -d "$DATE_END_BDY" +"%Y%m%d%H%M%S") ] ; do 
      WRFFILE=$BDYBASE/wrfout_d01_$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%T" )
      BDYORGPATH=$SCALEDIR/$MIEM/bdyorg/bdyorg_$(printf %05g $cnt)
      ln -sf $WRFFILE $BDYORGPATH
      #Update CDATE
      CDATE=$(date -u -d "$CDATE UTC + $INTERVALO_BDY seconds" +"%Y-%m-%d %T")
      cnt=$((cnt+1))
   done
  
else 

  echo "Not recognized WPS_DATA_SOURCE option: "$WPS_DATA_SOURCE
  echo "I can do nothing"
  ERROR=$(( $ERROR + 1 ))

fi

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SPACK_NETCDF_C}/lib:${SPACK_NETCDF_F}/lib:${SPACK_PNETCDF}/lib:${SPACK_HDF}/lib:$LD_LIBRARY_PATH

echo "Running scale init for member $MIEM"
$MPIEXE ./scale-rm_init namelist.scale_init
ERROR=$(( $ERROR + $? ))

echo "Running scale for member $MIEM"
$MPIEXE ./scale-rm namelist.scale
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
echo "Time taken by scale init and scale"
queue $MIEMBRO_INI $MIEMBRO_FIN 
time check_proc $MIEMBRO_INI $MIEMBRO_FIN

#Copy wrfout files to its final destionation.
for QMIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
   DIR_FECHA_INI=$(date -u -d "$FECHA_FORECAST_INI UTC" +"%Y%m%d%H%M%S")
   OUTPUTPATH="$HISTDIR/init/${DIR_FECHA_INI}/$QMIEM/"
   mkdir -p $OUTPUTPATH
   mv $SCALEDIR/$QMIEM/init/*.nc $OUTPUTPATH
   OUTPUTPATH="$HISTDIR/bdy/${DIR_FECHA_INI}/$QMIEM/"
   mkdir -p $OUTPUTPATH
   mv $SCALEDIR/$QMIEM/bdy/*.nc $OUTPUTPATH
   OUTPUTPATH="$HISTDIR/FCST/${DIR_FECHA_INI}/${QMIEM}/"
   mkdir -p $OUTPUTPATH
   echo "Copying file $SCALEDIR/$QMIEM/fcst/history" 
   mv $SCALEDIR/$QMIEM/fcst/*.nc        $OUTPUTPATH
#    mv $SCALEDIR/$QMIEM/*.log            $OUTPUTPATH
#    mv $SCALEDIR/$QMIEM/namelist*        $OUTPUTPATH
done


