#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

#Load experiment configuration
source $BASEDIR/conf/config.env
source $BASEDIR/conf/step.conf
source $BASEDIR/conf/exp.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source ${BASEDIR}/lib/encolar${QUEUESYS}.sh                    

##### FIN INICIALIZACION ######

if [ $STEP -ge 1 ] && [ $WPS_CYCLE -eq 0 ] ; then
   echo "WPS with WPS_CYCLE=0 is run only at STEP=0. And currently STEP=",$STEP        
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
sed -i -e "s|__INTERVALO__|$BDY_FREQ|g"       $WPSDIR/namelist.wps
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
        QPROC_NAME=GEOG_$STEP
	QPROC=$WPSPROC
	QNODE=$WPSNODE
	QTHREADS=$WPSTHREAD
	QMIEM=00
	QWALLTIME=$WPSWALLTIME
	QWORKPATH=$WPSDIR
	queue 00 00 
	check_proc 00 00
        if [ $? -ne 0 ] ; then
          echo "Error: Some members do not finish OK"
          echo "Aborting this step"
          exit 1 
        fi
fi

###
# Create the script that will run WPS for each ensemble member.
###
echo "Corriendo el WPS para el paso " $STEP

read -r -d '' QSCRIPTCMD << "EOF"
source $BASEDIR/conf/model.conf
source $BASEDIR/conf/exp.conf
source $BASEDIR/conf/step.conf

ulimit -s unlimited
echo "Processing member $MIEM"
ln -sf $WPSDIR/code/* $WPSDIR/$MIEM
cd $WPSDIR/$MIEM
cp $WPSDIR/geogrid/geo_em* .

cp $NAMELISTDIR/namelist.wps $WPSDIR/$MIEM/
INI_DATE_NML=$(date -u -d "$BDY_INI_DATE UTC" +"%Y-%m-%d_%T")
END_DATE_NML=$(date -u -d "$BDY_END_DATE UTC" +"%Y-%m-%d_%T")

sed -i -e "s|__FECHA_INI__|$INI_DATE_NML|g"   $WPSDIR/$MIEM/namelist.wps  #We use some random date for runing geogrid.
sed -i -e "s|__FECHA_FIN__|$END_DATE_NML|g"   $WPSDIR/$MIEM/namelist.wps
sed -i -e "s|__INTERVALO__|$BDY_FREQ|g"       $WPSDIR/$MIEM/namelist.wps
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


#Loop over the expected output files. If all the files are present then skip WPS step.
CDATE=$BDY_INI_DATE
FILE_COUNTER=0
while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $(date -d "$BDY_END_DATE" +"%Y%m%d%H%M%S") ] ; do

   MYPATH=$HISTDIR/WPS/met_em_ori/$(date -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")/$MIEM/
   MYFILE=$MYPATH/met_em.d01.$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%H:%M:%S" ).nc
   if ! [[ -e $MYFILE ]] ; then #My file do not exist
      echo "Not found file: " $MYFILE
      FILE_COUNTER=$(( $FILE_COUNTER + 1 ))
   fi
   CDATE=$(date -u -d "$CDATE UTC + $BDY_FREQ seconds" +"%Y-%m-%d %T")
done

if [ $FILE_COUNTER -eq 0 ] ; then 
   #All the required files are already there. 
   echo "We found all the required met_em files for member $MIEM, skiping ungrib and metgrid."
   ERROR=0

else #Runing ungrib and metgrid 
   echo $FILE_COUNTER " files are missing. We need to run UNGRIB and METGRID"
   #We need to run ungrib/wrftowps and metgrid to generate the required met_em files
   if [ $WPS_DATA_SOURCE == 'GFS' ] ; then 
 
      BDYBASE=$BDYPATH/gefs.$(date -d "$INI_BDY_DATE" +"%Y%m%d")/$(date -d "$INI_BDY_DATE" +"%H")/$BDYPREFIX/$MIEM/
      echo "Selected data source is GFS, we will run ungrib to decode the data"
      echo "I'm lloking for the BDY files in the folder  $BDYBASE"
      ln -sf $WPSDIR/$BDYVTABLE ./Vtable
      echo "Decoding gribs in : $WPSDIR/$MIEM"
      #Linking grib
      cd $WPSDIR/$MIEM
      ./link_grib.csh $BDYBASE/$BDYSEARCHSTRING  
      ERROR=$(( $ERROR + $? ))

      #Corro el Ungrib 
      OMP_NUM_THREADS=1
      rm ungrib.log
      $MPIEXESERIAL ./ungrib.exe > ungrib.log
      ERROR=$(( $ERROR + $? ))

   elif [ $WPS_DATA_SOURCE == 'WRF' ]  ; then

      BDYBASE=$BDYPATH/$(date -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")/$MIEM/
      echo "Selected data source is WRF, we will use wrf_to_wps tool to decode the data"
      echo "I'm lloking for the BDY files in the folder $BDYBASE"
      echo "Decoding wrfouts in : $WPSDIR/$MIEM" 
      #Need to loop over wrfout files.

      #Set the WPS_FILE_FORMAT [this depends on the file frequency]
      if [ $BDY_FREQ -lt 60 ] ; then
         WPS_FILE_DATE_FORMAT="%Y-%m-%d_%H:%M:%S"
      elif [ $BDY_FREQ -lt 3600 ] ; then
         WPS_FILE_DATE_FORMAT="%Y-%m-%d_%H:%M"
      else
         WPS_FILE_DATE_FORMAT="%Y-%m-%d_%H"
      fi

      CDATE=$BDY_INI_DATE
      while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $(date -d "$BDY_END_DATE" +"%Y%m%d%H%M%S") ] ; do 
         WRFFILE=$BDYBASE/wrfout_d01_$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%T" )
         WPSFILE=$WPSDIR/$MIEM/FILE:$(date -u -d  "$CDATE UTC" +"$WPS_FILE_DATE_FORMAT" )
         $MPIEXESERIAL ./wrf_to_wps.exe $WRFFILE $WPSFILE 
         ERROR=$(( $ERROR + $? ))
         #Update CDATE
         CDATE=$(date -u -d "$CDATE UTC + $BDY_FREQ seconds" +"%Y-%m-%d %T")
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
   OUTPUTPATH="$HISTDIR/WPS/met_em_ori/$(date -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")/$MIEM/"
   mkdir -p $OUTPUTPATH
   mv $WPSDIR/$MIEM/met_em* $OUTPUTPATH

fi

EOF

# Parametros de encolamiento
## Calculo cuantos miembros hay que correr
QNODE=$WPSNODE
QPROC=$WPSPROC
QTHREAD=$WPSTHREAD
QWALLTIME=$WPSWALLTIME
QPROC_NAME=WPS_${STEP}
QWORKPATH=$WPSDIR
QSKIP=$WPSSKIP

# Encolar
queue $BDY_MEM_INI $BDY_MEM_END
check_proc $BDY_MEM_INI $BDY_MEM_END
if [ $? -ne 0 ] ; then
   echo "Error: Some members do not finish OK"
   echo "Aborting this step"
   exit 1 
fi



echo "Termine de correr el WPS"

