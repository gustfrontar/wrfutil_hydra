#!/bin/bash
#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

### CONFIGURACION
#Load experiment configuration
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/assimilation.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source $BASEDIR/lib/spawn_utils.sh
#Set some environmental parameters
eval "$ENVSET"
##### FIN INICIALIZACION ######

cd $WRFDIR
#Seteamos las fechas de inicio y final de los forecasts.
FECHA_FORECAST_INI=$(date -u -d "$FECHA_INI UTC +$(($ANALISIS_FREC*$PASO)) seconds" +"%Y-%m-%d %T")
FECHA_FORECAST_END=$(date -u -d "$FECHA_INI UTC +$((($ANALISIS_FREC*$PASO)+$ANALISIS_WIN_FIN)) seconds" +"%Y-%m-%d %T")
FECHA_ANALISIS=$(date -u -d "$FECHA_INI UTC +$((($ANALISIS_FREC*$PASO)+$ANALISIS_FREC)) seconds" +"%Y-%m-%d %T")

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
   cp $NAMELISTDIR/namelist.input.asim.${WRFCONF} $WRFDIR/namelist.input.${WRFCONF}
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
   sed -i -e "s|__INTERVALO_WPS__|$INTERVALO_WPS|g"                 $WRFDIR/namelist.input.${WRFCONF}
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
done

#Descomprimimos el archivo .tar (si es que no fue descomprimido)
if [ ! -e $WRFDIR/code/real.exe ] ; then
   echo "Descomprimiendo WRF"
   mkdir -p $WRFDIR/code
   cd $WRFDIR
   tar -xf wrf.tar -C $WRFDIR/code
   #Si existe el namelist.input lo borramos para que no interfiera
   #con los que crea el sistema de asimilacion
   if [ -e $WRFDIR/code/namelist.input ] ; then
      rm -f $WRFDIR/code/namelist.input
   fi
fi
#Descomprimimos el wrfda.tar (si es que no fue descomprimido)
if [ ! -e $WRFDIR/code/da_update_bc.exe ] ; then
   echo "Descomprimiendo DA_UPDATE_BC"
   mkdir -p $WRFDIR/code
   tar -xf wrfda.tar -C $WRFDIR/code
fi
#Descomprimimos el spawn.tar (si es que no fue descomprimido)
if [ ! -e $WRFDIR/code/spawn.exe ] ; then
   echo "Descomprimiendo SPAWN"
   mkdir -p $WRFDIR/code
   cd $WRFDIR
   tar -xf spawn.tar -C $WRFDIR/code
fi


ulimit -s unlimited
OMP_NUM_THREADS=$REALTHREADS
OMP_STACKSIZE=512M

#Get the date corresponding to the met_em folder.
if [ $WPS_CYCLE -eq 1 ] ; then
   #Met_em source is updated periodicaly (as in a long DA cycling experiment)
   INI_STEP_DATE=$(date -u -d "$FECHA_INI UTC +$(($FORECAST_INI_FREQ*$PASO)) seconds" +"%Y-%m-%d %T")
   INI_BDY_DATE=$(date_floor "$INI_STEP_DATE" $INTERVALO_INI_BDY )
   INI_BDY_DATE=$(date -u -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")
else
   #Met_em source is fixed (e.g. a single GFS or WRF forecast). Usually for a short DA cycling experiment. 
   INI_BDY_DATE=$(date_floor "$FECHA_INI" $INTERVALO_INI_BDY )
   INI_BDY_DATE=$(date -u -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")
fi

for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
   cd $WRFDIR
   [[ -d $WRFDIR/$MIEM ]] && rm -r $WRFDIR/$MIEM
   mkdir $WRFDIR/$MIEM
   cd $WRFDIR/$MIEM
   if [ $MULTIMODEL -eq 0 ] ; then
      NLCONF=$(printf '%03d' $((10#$MODEL_CONF)) )
   else
      NLCONF=$(printf '%03d' $(( ( (10#$MIEM-1) % 10#$NCONF) + 1 )))
   fi
   cp $WRFDIR/namelist.input.${NLCONF} $WRFDIR/$MIEM/namelist.input
   ln -sf $WRFDIR/code/* .
   ln -sf $HISTDIR/WPS/met_em/$INI_BDY_DATE/$MIEM/met_em* $WRFDIR/$MIEM/
done

#Run real.exe using the spawn for parallel programs
IS_SERIAL=0
cd $WRFDIR/code/
echo "Running real.exe"
spawn ./real.exe $WRFDIR $MIEMBRO_INI $MIEMBRO_FIN $WRFNODE $WRFPROC $WRF_RUNTIME_FLAGS


#RUN DA_UPDATE_BC.EXE
if [ $PASO -gt 0 ] ; then
   echo "Running da_update_bc"
   for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
      cd $WRFDIR/$MIEM
      cp $NAMELISTDIR/parame* .
      mv $WRFDIR/$MIEM/wrfinput_d01 $WRFDIR/$MIEM/wrfinput_d01.org
      cp $HISTDIR/ANAL/$(date -u -d "$FECHA_INI UTC +$(($ANALISIS_FREC*$PASO)) seconds" +"%Y%m%d%H%M%S")/anal$(printf %05d $((10#$MIEM))) $WRFDIR/$MIEM/wrfinput_d01
      ln -sf $WRFDIR/code/da_update_bc.exe $WRFDIR/$MIEM/da_update_bc.exe
      #time $WRFDIR/$MIEM/da_update_bc.exe > ./da_update_bc_${PASO}_${MIEM}.log  &
   done
   #Run da_update_bc.exe using the spawn for serial programs
   cd $WRFDIR/code/
   spawn ./mpi_updatebc_wrapper.exe $WRFDIR $MIEMBRO_INI $MIEMBRO_FIN $WRFNODE $WRFPROC $WRF_RUNTIME_FLAGS
fi

export OMP_NUM_THREADS=1
#Run wrf.exe using the spawn for parallel programs
cd $WRFDIR/code
IS_SERIAL=0
echo "Running wrf.exe" 
spawn ./wrf.exe $WRFDIR $MIEMBRO_INI $MIEMBRO_FIN $WRFNODE $WRFPROC $WRF_RUNTIME_FLAGS

#Copy wrfout files corresponding to the analysis time to the history folder.
if [[ ! -z "$GUARDOGUESS" ]] && [[ $GUARDOGUESS -eq 1 ]] ; then
   for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
       OUTPUTPATH="$HISTDIR/GUES/$(date -u -d "$FECHA_ANALISIS UTC" +"%Y%m%d%H%M%S")/"
       mkdir -p $OUTPUTPATH
       echo "Copying file $WRFDIR/$MIEM/wrfout_d01_$(date -u -d "$FECHA_ANALISIS UTC" +"%Y-%m-%d_%T" )"
       cp $WRFDIR/$MIEM/wrfout_d01_$(date -u -d "$FECHA_ANALISIS UTC" +"%Y-%m-%d_%T") $OUTPUTPATH/gues$(printf %05d $((10#$MIEM)))
       mv $WRFDIR/$MIEM/rsl.error.0000                                                $OUTPUTPATH/wrf.log
       mv $WRFDIR/$MIEM/namelist*                                                     $OUTPUTPATH
   done
fi


