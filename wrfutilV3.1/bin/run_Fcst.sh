#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

### CONFIGURACION
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/${EXPTYPE}.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source ${BASEDIR}/lib/encolar${QUEUESYS}.sh                     # Selecciona el metodo de encolado segun el systema QUEUESYS elegido
##### FIN INICIALIZACION ######

cd $WRFDIR
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
if [ $MULTIMODEL -eq 0 ] ; then  
   MULTIMODEL_CONF_INI=$MODEL_CONF;MULTIMODEL_CONF_FIN=$MODEL_CONF
   echo "We will run a single configuration ensemble"
else 
   echo "We will run a multimodel configuration ensemble with $NCONF configurations"
   MULTIMODEL_CONF_INI=1;MULTIMODEL_CONF_FIN=$NCONF
fi

for ICONF in $(seq -w $MULTIMODEL_CONF_INI $MULTIMODEL_CONF_FIN) ; do
   WRFCONF=$(printf '%03d' $((10#$ICONF)))
   echo "Model configuration is " $WRFCONF
   cp $NAMELISTDIR/namelist.input.${WRFCONF} $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__START_YEAR__|$IY|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__START_MONTH__|$IM|g"                              $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__START_DAY__|$ID|g"                                $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__START_HOUR__|$IH|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__START_MINUTE__|$Im|g"                             $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__START_SECOND__|$Is|g"                             $WRFDIR/namelist.input.${WRFCONF} 
   sed -i -e "s|__END_YEAR__|$FY|g"                                 $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__END_MONTH__|$FM|g"                                $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__END_DAY__|$FD|g"                                  $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__END_HOUR__|$FH|g"                                 $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__END_MINUTE__|$Fm|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__END_SECOND__|$Fs|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__INTERVALO_WPS__|$FORECAST_BDY_FREQ|g"             $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__INTERVALO_WRF__|$(($INTERVALO_WRF/60))|g"         $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__E_WE__|$E_WE|g"                                   $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__E_SN__|$E_SN|g"                                   $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__DX__|$DX|g"                                       $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__DY__|$DY|g"                                       $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__DT__|$DT|g"                                       $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__IOTYPE__|$IOTYPE|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__METLEV__|$METLEV|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__NUMTILE__|$NUMTILE|g"                             $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__NIOT__|$NIOT|g"                                   $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__NIOG__|$NIOG|g"                                   $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__E_VERT__|$E_VERT|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__P_TOP__|$P_TOP|g"                                 $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__RADT__|$RADT|g"                                   $WRFDIR/namelist.input.${WRFCONF}
done

#Descomprimimos el archivo .tar (si es que no fue descomprimido)
if [ ! -e $WRFDIR/code/real.exe ] ; then
   echo "Decompressing executables ..."
   mkdir -p $WRFDIR/code
   tar -xf wrf.tar -C $WRFDIR/code

   tar -xf pert_met_em.tar -C $WRFDIR/code
   #Si existe el namelist.input lo borramos para que no interfiera
   #con los que crea el sistema de asimilacion
   if [ -e $WRFDIR/code/namelist.input ] ; then
      rm -f $WRFDIR/code/namelist.input
   fi
fi

#Build the script to run REAL/DA-UPDATE-BC/WRF
read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
source $BASEDIR/conf/model.conf
echo "Procesando Miembro $MIEM" 

if [ $MULTIMODEL -eq 0 ] ; then
   NLCONF=$(printf '%03d' ${MODEL_CONF} )
else 
   NLCONF=$(printf '%03d' $(( ( ( 10#$MIEM - 1 ) % 10#$NCONF ) + 1 )) )
fi
cp $WRFDIR/namelist.input.${NLCONF} $WRFDIR/$MIEM/namelist.input

ln -sf $WRFDIR/code/* . 

DATE_FORECAST_INI=$(date -u -d "$FECHA_INI UTC +$(($FORECAST_INI_FREQ*$PASO)) seconds" +"%Y-%m-%d %T")
DATE_FORECAST_END=$(date -u -d "$FECHA_INI UTC +$((($FORECAST_INI_FREQ*$PASO)+$FORECAST_LEAD_TIME )) seconds" +"%Y-%m-%d %T")

if [ $WPS_CYCLE -eq 1 ] ; then
   INI_STEP_DATE=$(date -u -d "$FECHA_INI UTC +$(($FORECAST_INI_FREQ*$PASO)) seconds" +"%Y-%m-%d %T")
   INI_BDY_DATE=$(date_floor "$INI_STEP_DATE" $INTERVALO_INI_BDY )
   INI_BDY_DATE=$(date -u -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")
else
   INI_BDY_DATE=$(date_floor "$FECHA_INI" $INTERVALO_BDY )
   INI_BDY_DATE=$(date -u -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")
fi

MET_EM_DIR=$HISTDIR/WPS/met_em/${INI_BDY_DATE}/$MIEM/
if [ $FORECAST_BDY_FREQ -eq $INTERVALO_BDY ] ; then
  ln -sf $MET_EM_DIR/met_em* $WRFDIR/$MIEM/
else    
   #We will conduct interpolation of the met_em files.
   echo "Interpolating files in time to reach $FORECAST_BDY_FREQ time frequency."
   CDATE=$DATE_FORECAST_INI
   WPS_FILE_DATE_FORMAT="%Y-%m-%d_%H:%M:%S"

   while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $(date -d "$DATE_FORECAST_END" +"%Y%m%d%H%M%S") ] ; do
     FILE_TAR=met_em.d01.$(date -u -d "$CDATE UTC" +"$WPS_FILE_DATE_FORMAT" ).nc

     if [ -e $MET_EM_DIR/$FILE_TAR ] ; then #Target file exists. 
        ln -sf $MET_EM_DIR/$FILE_TAR  ./
     else 
        DATE_INI=$(date_floor "$CDATE" $INTERVALO_BDY )
        DATE_END=$(date -u -d "$DATE_INI UTC + $INTERVALO_BDY seconds" +"%Y-%m-%d %T")
        FILE_INI=met_em.d01.$(date -u -d "$DATE_INI UTC" +"$WPS_FILE_DATE_FORMAT" ).nc
        FILE_END=met_em.d01.$(date -u -d "$DATE_END UTC" +"$WPS_FILE_DATE_FORMAT" ).nc
        #File does not exist. Interpolate data in time to create it.  
        echo "&interp                          "  > ./pertmetem.namelist
        echo "time_ini=0                       " >> ./pertmetem.namelist
        echo "time_end=$INTERVALO_BDY          " >> ./pertmetem.namelist
        echo "time_tar=$((($(date -d "$CDATE" +%s) - $(date -d "$DATE_INI" +%s))))    " >> ./pertmetem.namelist
        echo "date_tar='$(date -u -d "$CDATE UTC" +"$WPS_FILE_DATE_FORMAT" )'         " >> ./pertmetem.namelist
        echo "file_ini='$MET_EM_DIR/$FILE_INI'   " >> ./pertmetem.namelist
        echo "file_end='$MET_EM_DIR/$FILE_END'   " >> ./pertmetem.namelist
        echo "file_tar='$WRFDIR/$MIEM/$FILE_TAR' " >> ./pertmetem.namelist
        echo "/                                " >> ./pertmetem.namelist
        echo "Running INTERP MET EM for member $MIEM"
        $MPIEXESERIAL ./interp_met_em.exe > ./interp_met_em.log  
        ERROR=$(( $ERROR + $? )) 
        #ln -sf $WPSDIR/$MIEM/$FILE_TAR  ./
        #TODO: Check if would be better to create the new files in the WPS_MET_EM_DIR so
        #Interpolated files can be used again in the next cycle.
     fi
     #Update CDATE
     CDATE=$(date -u -d "$CDATE UTC + $FORECAST_BDY_FREQ seconds" +"%Y-%m-%d %T")
   done
fi

echo "Running REAL for member $MIEM"
$MPIEXE $WRFDIR/$MIEM/real.exe 
ERROR=$(( $ERROR + $? ))
mv rsl.error.0000 ./real_${PASO}_${MIEM}.log
EXCOD=$?
[[ $EXCOD -ne 0 ]] && dispararError 9 "real.exe"

echo "Running WRF for member $MIEM"
$MPIEXE $WRFDIR/$MIEM/wrf.exe 
ERROR=$(( $ERROR + $? ))
mv rsl.error.0000 ./wrf_${PASO}_${MIEM}.log

EOF

#Node / core distribution parameters
QNODE=$WRFNODE
QPROC=$WRFPROC
TPROC=$WRFTPROC
QTHREAD=$WRFTHREAD
QWALLTIME=$WRFWALLTIME
QPROC_NAME=FCST_${PASO}
QCONF=${EXPTYPE}.conf
QWORKPATH=$WRFDIR

#Execute the job 
echo "Time taken by real.exe and wrf.exe"
queue $MIEMBRO_INI $MIEMBRO_FIN 
time check_proc $MIEMBRO_INI $MIEMBRO_FIN

#Copy wrfout files to its final destionation.
for QMIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
    OUTPUTPATH="$HISTDIR/FCST/$(date -u -d "$FECHA_FORECAST_INI UTC" +"%Y%m%d%H%M%S")/${QMIEM}/"
    mkdir -p $OUTPUTPATH
    echo "Copying file $WRFDIR/$QMIEM/wrfout_d01_$(date -u -d "$FECHA_FORECAST_INI UTC" +"%Y-%m-%d_%T" )"
    mv $WRFDIR/$QMIEM/wrfout_d01_*     $OUTPUTPATH
    mv $WRFDIR/$QMIEM/*.log            $OUTPUTPATH
    mv $WRFDIR/$QMIEM/namelist*        $OUTPUTPATH
done


