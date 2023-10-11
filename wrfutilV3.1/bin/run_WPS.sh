#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

#Load experiment configuration
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/${EXPTYPE}.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source ${BASEDIR}/lib/encolar${QUEUESYS}.sh                    

##### FIN INICIALIZACION ######

#Decompress tar files
if [ ! -e $WPSDIR/code/geogrid.exe ] ; then
   echo "Descomprimiendo WPS"
   mkdir -p $WPSDIR/code
   tar -xf $WPSDIR/wps.tar -C $WPSDIR/code
   #Remove namelist.wps (it will be created later from the templates)
   if [ -e $WPSDIR/code/namelist.wps ] ; then
      rm -f $WPSDIR/code/namelist.wps 
   fi
fi

if [ ! -e $WPSDIR/code/wrf_to_wps.exe ] ; then
   echo "Descomprimiendo WRF_TO_WPS"
   mkdir -p $WPSDIR/code
   tar -xf $WPSDIR/wrf_to_wps.tar -C $WPSDIR/code
fi


mkdir -p $WPSDIR
cp $NAMELISTDIR/namelist.wps $WPSDIR/

#We found the closests dates in the boundary data to the initial and final times of the forecast. 

if [ $WPS_CYCLE -eq 1 ] ; then 
   #WPS_CYCLE = 1 means that boundary data will be taken from forecasts initialized at different times. 
   #so there is a bdy forecast initialization frequency and a bdy data output frequency. 
   #Find the prior bdy data initialization date closest to the current step date 
   FECHA_INI_PASO=$(date -u -d "$FECHA_INI UTC +$(($WPS_INI_FREQ*$PASO)) seconds" +"%Y-%m-%d %T")
   #Find the posterior bdy data initialization date closest to the forecast end date
   FECHA_END_PASO=$(date -u -d "$FECHA_INI UTC +$((($WPS_INI_FREQ*$PASO)+$WPS_LEAD_TIME )) seconds" +"%Y-%m-%d %T")
else 
   #Asume that bdy data for the entire experiment is available in the same folder (eg using analysis as bdy or for short experiments that use a single forecasts as bdy)	
   #In this case WPS is run at PASO=0
   if [ $PASO -ge 1 ] ; then   
      echo "WPS with WPS_CYCLE=0 is run only at PASO=0. And currently PASO=",$PASO	   
      return 0
   fi
   FECHA_INI_PASO=$FECHA_INI
   FECHA_END_PASO=$FECHA_FIN
fi

DATE_INI_BDY=$(date_floor "$FECHA_INI_PASO" $INTERVALO_BDY )   #Get the closest prior date in which BDY data is available.
DATE_INI_BDY=$(date -u -d "$DATE_INI_BDY UTC" +"%Y-%m-%d_%T" ) #Cambio al formato WPS
DATE_END_BDY=$(date_ceil  "$FECHA_END_PASO" $INTERVALO_BDY )   #Get the closest posterior date in which BDY data is available. 
DATE_END_BDY=$(date -u -d "$DATE_END_BDY UTC" +"%Y-%m-%d_%T" ) #Cambio al formato WPS

#Edit the namelist from the template
sed -i -e "s|__FECHA_INI__|$DATE_INI_BDY|g"   $WPSDIR/namelist.wps
sed -i -e "s|__FECHA_FIN__|$DATE_END_BDY|g"   $WPSDIR/namelist.wps
sed -i -e "s|__INTERVALO__|$INTERVALO_WPS|g"  $WPSDIR/namelist.wps
sed -i -e "s|__E_WE__|$E_WE|g"                $WPSDIR/namelist.wps
sed -i -e "s|__E_SN__|$E_SN|g"                $WPSDIR/namelist.wps
sed -i -e "s|__DX__|$DX|g"                    $WPSDIR/namelist.wps
sed -i -e "s|__DY__|$DY|g"                    $WPSDIR/namelist.wps
sed -i -e "s|__MAP_PROJ__|$MAP_PROJ|g"        $WPSDIR/namelist.wps
sed -i -e "s|__REF_LAT__|$REF_LAT|g"          $WPSDIR/namelist.wps
sed -i -e "s|__REF_LON__|$REF_LON|g"          $WPSDIR/namelist.wps
sed -i -e "s|__TRUELAT1__|$TRUELAT1|g"        $WPSDIR/namelist.wps
sed -i -e "s|__TRUELAT2__|$TRUELAT2|g"        $WPSDIR/namelist.wps
sed -i -e "s|__STAND_LON__|$STAND_LON|g"      $WPSDIR/namelist.wps
sed -i -e "s|__GEOG__|$GEOG|g"                $WPSDIR/namelist.wps
sed -i -e "s|__IOTYPE__|$IOTYPE|g"            $WPSDIR/namelist.wps


####
## Corremos geogrid (si es que no fue corrido)
###
if [ ! -e $WPSDIR/geogrid/geo_em.d01.nc ] ; then
   echo "Ejecutando geogrid"
   mkdir -p $WPSDIR/geogrid
   ln -sf $WPSDIR/code/*   $WPSDIR/geogrid/
   cp $WPSDIR/namelist.wps $WPSDIR/geogrid/ 
read -r -d '' QSCRIPTCMD << "EOF"
        cd $WPSDIR/geogrid
	ulimit -s unlimited
        $MPIEXE ./geogrid.exe $WPS_RUNTIME_FLAGS
EOF
        QPROC_NAME=GEOG_$PASO
	QPROC=$WPSPROC
	QNODE=$WPSNODE
	QTHREADS=$WPSTHREAD
	QMIEM=00
	QWALLTIME=$WPSWALLTIME
        QCONF=${EXPTYPE}.conf
	queue 00 00 
	check_proc 00 00
fi

###
# Create the script that will run WPS for each ensemble member.
###
echo "Corriendo el WPS para el paso " $PASO

read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
echo "Processing member $MIEM"
cd $WPSDIR
rm -r $WPSDIR/$MIEM
mkdir -p $WPSDIR/$MIEM
ln -sf $WPSDIR/code/* $WPSDIR/$MIEM
cp $WPSDIR/namelist.wps $WPSDIR/$MIEM/
cd $WPSDIR/$MIEM
ln -sf $WPSDIR/geogrid/geo_em* .
if [ $WPS_CYCLE -eq 1 ] ; then
   FECHA_INI_PASO=$(date -u -d "$FECHA_INI UTC +$(($WPS_INI_FREQ*$PASO)) seconds" +"%Y-%m-%d %T")
else 
   FECHA_INI_PASO=$FECHA_INI
fi
DATE_INI_BDY=$(date_floor "$FECHA_INI_PASO" $INTERVALO_INI_BDY )  #Get the closest prior date in which BDY data is available.
DATE_END_BDY=$(date_ceil  "$FECHA_END_PASO" $INTERVALO_BDY )      #Get the closest posterior date in which BDY data is available.


if [ $WPS_DATA_SOURCE == 'GFS' ] ; then 
   BDYBASE=$BDYPATH/gefs.$(date -d "$DATE_INI_BDY" +"%Y%m%d")/$(date -d "$DATE_INI_BDY" +"%H")/$BDYPREFIX/$MIEM/
   echo "Selected data source is GFS, we will run ungrib to decode the data"
   echo "I'm lloking for the BDY files in the folder  $BDYBASE"
   ln -sf $WPSDIR/$BDYVTABLE ./Vtable
   echo "Decoding gribs in : $WPSDIR/$MIEM"
   #Linkeo los gribs
   cd $WPSDIR/$MIEM
   ./link_grib.csh $BDYBASE/$BDYSEARCHSTRING  

   [[ $? -ne 0 ]] && dispararError 9 "link_grib.csh "

   #Corro el Ungrib 
   OMP_NUM_THREADS=1
   $MPIEXESERIAL ./ungrib.exe $WPS_RUNTIME_FLAGS > ungrib.log
   [[ $? -ne 0 ]] && dispararError 9 "ungrib.exe"

elif[ $WPS_DATA_SOURCE == 'WRF' ] ; then 

   BDYBASE=$BDYPATH/$(date -d "$DATE_INI_BDY" +"%Y%m%d%H%M%S")/$MIEM/
   echo "Selected data source is WRF, we will use wrf_to_wps tool to decode the data"
   echo "I'm lloking for the BDY files in the folder $BDYBASE"
   echo "Decoding wrfouts in : $WPSDIR/$MIEM" 
   #Need to loop over wrfout files.
   CDATE=$DATE_INI_BDY 
   while [ $CDATE -le $DATE_END_BDY ] ; do 
      WRFFILE=BDYBASE/wrfout_d01_$(date -u -d "$FECHA_ANALISIS UTC" +"%Y-%m-%d_%T" )
      WPSFILE=$WPSDIR/$MIEM/FILE:$(date -u -d "$FECHA_ANALISIS UTC" +"%Y-%m-%d_%H" )
      $MPIEXESERIAL ./wrf_to_wps.exe $WRFFILE $WPSFILE
      #Update CDATE
      CDATE=$(date -u -d "$CDATE + $INTERVALO_BDY seconds" +"%Y-%m-%d %T")
   done
else 
  echo "Not recognized WPS_DATA_SOURCE option: "$WPS_DATA_SOURCE
  echo "I can do nothing"
  disparaError 9 "run_wps.sh"

fi

#Corro el Metgrid en paralelo
$MPIEXE ./metgrid.exe $WPS_RUNTIME_FLAGS
EC=$?
[[ $EC -ne 0 ]] && dispararError 9 "metgrid.exe"

EOF

# Parametros de encolamiento
## Calculo cuantos miembros hay que correr
QNODE=$WPSNODE
QPROC=$WPSPROC
QTHREAD=$WPSTHREAD
QWALLTIME=$WPSWALLTIME
QPROC_NAME=WPS_${PASO}
QCONF=${EXPTYPE}.conf

# Encolar
queue $BDY_MIEMBRO_INI $BDY_MIEMBRO_FIN
check_proc $BDY_MIEMBRO_INI $BDY_MIEMBRO_FIN

FECHA_INI_BDY=$(date_floor "$FECHA_INI_PASO" $INTERVALO_BDY )    #Get the closest prior date in which BDY data is available.
FECHA_INI_BDY=$(date -u -d "$FECHA_INI_BDY UTC" +"%Y%m%d%H%M%S" ) #Cambio al formato WPS
#Copiamos los archivos del directorio 
for QMIEM in $(seq -w $BDY_MIEMBRO_INI $BDY_MIEMBRO_FIN) ; do
   OUTPUTPATH="$HISTDIR/WPS/met_em_ori/${FECHA_INI_BDY}/$QMIEM/"
   mkdir -p $OUTPUTPATH
   mv $WPSDIR/$QMIEM/met_em* $OUTPUTPATH
done

echo "Termine de correr el WPS"

