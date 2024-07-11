#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

#Load experiment configuration
BASEDIR=$(pwd)/../
source $BASEDIR/../setup_spack.sh
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/${EXPTYPE}.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source $BASEDIR/lib/encolar${QUEUESYS}.sh                    

##### FIN INICIALIZACION ######

mkdir -p $SCALEPPDIR
cp $NAMELISTDIR/namelist.scale_pp   $SCALEPPDIR/
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

DATE_INI_BDY_f=$(date_floor "$FECHA_INI_PASO" $INTERVALO_BDY )   #Get the closest prior date in which BDY data is available.
DATE_INI_BDY=$(date -u -d "$DATE_INI_BDY_f UTC" +"%Y-%m-%d_%T" ) #Cambio al formato WPS
DATE_END_BDY_f=$(date_ceil  "$FECHA_END_PASO" $INTERVALO_BDY )   #Get the closest posterior date in which BDY data is available. 
DATE_END_BDY=$(date -u -d "$DATE_END_BDY_f UTC" +"%Y-%m-%d_%T" ) #Cambio al formato WPS

#Edit the namelist from the template

E_WE_LOC=$((E_WE / NPROC_X))
E_SN_LOC=$((E_SN / NPROC_Y))
DATE_INI_SCALE="$(date -ud "$DATE_INI_BDY_f" +'%Y,%m,%d,%H,%M,%S' )"
if [ $MAP_PROJ == "lambert" ];then
  MAP_PROJ_SCALE="LC"
fi
GTOPO30=$TOPO_LAND_ORG_SCALE/topo/GTOPO30/Products
GLCCv2=$TOPO_LAND_ORG_SCALE/landuse/GLCCv2/Products

sed -i -e "s|__FECHA_INI__|$DATE_INI_SCALE|g" $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__E_WE_LOC__|$E_WE_LOC|g"        $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__E_SN_LOC__|$E_SN_LOC|g"        $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__NPROC_X__|$NPROC_X|g"          $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__NPROC_Y__|$NPROC_Y|g"          $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__DX__|$DX|g"                    $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__DY__|$DY|g"                    $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__MAP_PROJ__|$MAP_PROJ_SCALE|g"  $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__REF_LAT__|$REF_LAT|g"          $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__REF_LON__|$REF_LON|g"          $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__TRUELAT1__|$TRUELAT1|g"        $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__TRUELAT2__|$TRUELAT2|g"        $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__GTOPO30__|\"$GTOPO30\"|g"          $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__GLCCv2__|\"$GLCCv2\"|g"            $SCALEPPDIR/namelist.scale_pp

####
## Corremos scale pp (si es que no fue corrido)
###
if [ ! -e $SCALEPPDIR/mktopo/topo/topo.pe000000.nc ] ; then
if [ -e $HISTDIR/const/topo/topo.pe000000.nc ] ; then
   echo "copiando topo and landuse"
   mkdir -p $SCALEPPDIR/mktopo/topo
   mkdir -p $SCALEPPDIR/mktopo/landuse
   cp $HISTDIR/const/topo/*.nc $SCALEPPDIR/mktopo/topo/
   cp $HISTDIR/const/landuse/*.nc $SCALEPPDIR/mktopo/landuse/
else
   echo "Ejecutando scale pp"
   mkdir -p $SCALEPPDIR/mktopo/topo
   mkdir -p $SCALEPPDIR/mktopo/landuse
   mkdir -p $SCALEPPDIR/mktopo/log/scale_pp
   #ln -sf $WPSDIR/code/*   $WPSDIR/geogrid/
   #cp $SCALEPPDIR/namelist.scale_pp $WPSDIR/geogrid/ 
read -r -d '' QSCRIPTCMD << "EOF"
        cd $SCALEPPDIR/mktopo
        ln -sf $SCALEPPDIR/scale-rm_pp${SCALE_opt} ./scale-rm_pp
        cp $SCALEPPDIR/namelist.scale_pp ./
	ulimit -s unlimited
        export FORT90L=${SCALE_RUNTIME_FLAGS}
        export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SPACK_NETCDF_C}/lib:${SPACK_NETCDF_F}/lib:${SPACK_PNETCDF}/lib:${SPACK_HDF}/lib:$LD_LIBRARY_PATH
        $MPIEXE ./scale-rm_pp namelist.scale_pp
        ERROR=$(( $ERROR + $? ))
EOF
        QPROC_NAME=scale_pp_$PASO
	QPROC=$SCALEPROC
	QNODE=$SCALENODE
	QTHREADS=$SCALETHREAD
	QMIEM=00
	QWALLTIME="00:10:00"
        QCONF=${EXPTYPE}.conf
	QWORKPATH=$SCALEPPDIR
	queue 00 00 
	check_proc 00 00

#Copiamos los archivos del directorio 
OUTPUTPATH="$HISTDIR/const/"
mkdir -p $OUTPUTPATH
cp -r $SCALEPPDIR/mktopo/topo    $OUTPUTPATH
cp -r $SCALEPPDIR/mktopo/landuse $OUTPUTPATH

fi
fi

echo "Termine de correr el SCALE pp"

