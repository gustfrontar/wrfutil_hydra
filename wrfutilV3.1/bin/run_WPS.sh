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

if [ $PASO -ge 1 && $WPS_CYCLE -eq 0 ] ; then
   echo "WPS with WPS_CYCLE=0 is run only at PASO=0. And currently PASO=",$PASO         
   return 0
fi

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

#Edit the namelist from the template
sed -i -e "s|__FECHA_INI__|"2000-01-01_00:00:00"|g"   $WPSDIR/namelist.wps  #We use some random date for runing geogrid.
sed -i -e "s|__FECHA_FIN__|"2000-01-01_00:00:00"|g"   $WPSDIR/namelist.wps
sed -i -e "s|__INTERVALO__|$INTERVALO_BDY|g"  $WPSDIR/namelist.wps
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
   #ln -sf $WPSDIR/code/*   $WPSDIR/geogrid/
   #cp $WPSDIR/namelist.wps $WPSDIR/geogrid/ 
read -r -d '' QSCRIPTCMD << "EOF"
        cd $WPSDIR/geogrid
        ln -sf $WPSDIR/code/*   ./
        cp $WPSDIR/namelist.wps ./
	ulimit -s unlimited
        $MPIEXE ./geogrid.exe 
        ERROR=$(( $ERROR + $? ))
        cp geo_em.d01.nc $WPSDIR/geogrid/
EOF
        QPROC_NAME=GEOG_$PASO
	QPROC=$WPSPROC
	QNODE=$WPSNODE
	QTHREADS=$WPSTHREAD
	QMIEM=00
	QWALLTIME=$WPSWALLTIME
        QCONF=${EXPTYPE}.conf
	QWORKPATH=$WPSDIR
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
ln -sf $WPSDIR/code/* $WPSDIR/$MIEM
cd $WPSDIR/$MIEM
cp $WPSDIR/geogrid/geo_em* .

#Define some variables that depends if we are running in assimilation
#mode or in forecast mode.
if [ "$EXPTYPE" = "assimilation" ] ; then
   CYCLE_FREQ=$ANALISIS_FREC
   LEAD=$ANALISIS_WIN_FIN
   SPIN_UP=$SPIN_UP_LENGTH
elif [ $EXPTYPE -eq "forecast" ] ; then
   CYCLE_FREQ=$FORECAST_INI_FREQ
   LEAD=$FORECAST_LEAD_TIME
   SPIN_UP=0
fi

#Define the start date and enddate of the WPS files to be generated.
#If WPS_CYCLE is 0 then we will process all the files from FECHA_INI to FECHA_FIN
if [ $WPS_CYCLE -eq 1 ] ; then 
   if [ $PASO -eq 0 ] ; then
      INI_DATE_STEP=$FECHA_INI
      END_DATE_STEP=$(date -u -d "$FECHA_INI UTC +$SPIN_UP seconds" +"%Y-%m-%d %T")
   else
      INI_DATE_STEP=$(date -u -d "$FECHA_INI UTC +$(($CYCLE_FREQ*($PASO-1)+$SPIN_UP)) seconds" +"%Y-%m-%d %T")
      END_DATE_STEP=$(date -u -d "$INI_DATE_STEP UTC +$LEAD seconds" +"%Y-%m-%d %T")
   fi
else  # WPS_CYCLE -eq 0 
   INI_DATE_STEP=$FECHA_INI
   END_DATE_STEP=$FECHA_FIN
fi

DATE_INI_BDY=$(date_floor "$INI_DATE_STEP" $INTERVALO_INI_BDY )   #get the closest initialization of the boudndary data in which BDY data is available.
INI_DATE_STEP=$(date_floor "$INI_DATE_STEP" $INTERVALO_BDY )      #Get the closest prior date in which BDY data is available.
END_DATE_STEP=$(date_ceil  "$END_DATE_STEP" $INTERVALO_BDY )      #Get the closest posterior date in which BDY data is available.

cp $NAMELISTDIR/namelist.wps $WPSDIR/$MIEM/
sed -i -e "s|__FECHA_INI__|$INI_DATE_STEP|g"  $WPSDIR/$MIEM/namelist.wps  #We use some random date for runing geogrid.
sed -i -e "s|__FECHA_FIN__|$END_DATE_STEP|g"  $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__INTERVALO__|$INTERVALO_BDY|g"  $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__E_WE__|$E_WE|g"                $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__E_SN__|$E_SN|g"                $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__DX__|$DX|g"                    $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__DY__|$DY|g"                    $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__MAP_PROJ__|$MAP_PROJ|g"        $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__REF_LAT__|$REF_LAT|g"          $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__REF_LON__|$REF_LON|g"          $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__TRUELAT1__|$TRUELAT1|g"        $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__TRUELAT2__|$TRUELAT2|g"        $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__STAND_LON__|$STAND_LON|g"      $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__GEOG__|$GEOG|g"                $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__IOTYPE__|$IOTYPE|g"            $WPSDIR/$MIEM/namelist.wps

if [ $WPS_DATA_SOURCE == 'GFS' ] ; then 

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
   echo "Selected data source is WRF, we will use wrf_to_wps tool to decode the data"
   echo "I'm lloking for the BDY files in the folder $BDYBASE"
   echo "Decoding wrfouts in : $WPSDIR/$MIEM" 
   #Need to loop over wrfout files.

   #Set the WPS_FILE_FORMAT [this depends on the file frequency]
   if [ $INTERVALO_BDY -lt 60 ] ; then
      WPS_FILE_DATE_FORMAT="%Y-%m-%d_%H:%M:%S"
   elif [ $INTERVALO_BDY -lt 3600 ] ; then
      WPS_FILE_DATE_FORMAT="%Y-%m-%d_%H:%M"
   else
      WPS_FILE_DATE_FORMAT="%Y-%m-%d_%H"
   fi

   CDATE=$INI_DATE_STEP
   while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $(date -d "$END_DATE_STEP" +"%Y%m%d%H%M%S") ] ; do 
      WRFFILE=$BDYBASE/wrfout_d01_$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%T" )
      WPSFILE=$WPSDIR/$MIEM/FILE:$(date -u -d  "$CDATE UTC" +"$WPS_FILE_DATE_FORMAT" )
      $MPIEXESERIAL ./wrf_to_wps.exe $WRFFILE $WPSFILE 
      ERROR=$(( $ERROR + $? ))
      #Update CDATE
      CDATE=$(date -u -d "$CDATE UTC + $INTERVALO_BDY seconds" +"%Y-%m-%d %T")
   done
  
else 

  echo "Not recognized WPS_DATA_SOURCE option: "$WPS_DATA_SOURCE
  echo "I can do nothing"
  ERROR=$(( $ERROR + 1 ))

fi

echo "Running METGRID for member $MIEM"
$MPIEXE ./metgrid.exe 
ERROR=$(( $ERROR + $? ))

#Copy data to the met_em_ori directory. 
echo "Copy data"
OUTPUTPATH="$HISTDIR/WPS/met_em_ori/${DATE_INI_BDY}/$MIEM/"
mkdir -p $OUTPUTPATH
mv $WPSDIR/$MIEM/met_em* $OUTPUTPATH

EOF

# Parametros de encolamiento
## Calculo cuantos miembros hay que correr
QNODE=$WPSNODE
QPROC=$WPSPROC
QTHREAD=$WPSTHREAD
QWALLTIME=$WPSWALLTIME
QPROC_NAME=WPS_${PASO}
QCONF=${EXPTYPE}.conf
QWORKPATH=$WPSDIR

# Encolar
queue $BDY_MIEMBRO_INI $BDY_MIEMBRO_FIN
check_proc $BDY_MIEMBRO_INI $BDY_MIEMBRO_FIN

echo "Termine de correr el WPS"

