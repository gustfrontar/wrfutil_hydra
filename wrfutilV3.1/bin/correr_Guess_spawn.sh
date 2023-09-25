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
source $BASEDIR/conf/forecast.conf
source $BASEDIR/conf/assimilation.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf

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

#Obtenemos el numero total de procesadores a usar por cada WRF
WRFTPROC=$(( $WRFNODE * $WRFPROC )) #Total number of cores to be used by each ensemble member.
TCORES=$(( $ICORE * INODE ))   #Total number of cores available to run the ensemble.
MAX_SIM_MEM=$(( $TCORES / $WRFTPROC ))  #Floor rounding (bash default) Maximumu number of simultaneous members
if [ $MAX_SIM_MEM -gt $MIEMBRO_FIN ] ; then
   MAX_SIM_MEM=$MIEMBRO_FIN
fi


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
   if [ -e namelist.input ] ; then
      rm -f namelist.input
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
   ln -sf $HISTDIR/WPS/met_em/$FECHA_INI_BDY/$MIEM/met_em* $WRFDIR/$MIEM/
done


#Run REAL.EXE
cd $WRFDIR
ini_mem=$(( $MIEMBRO_INI ))
end_mem=$(( $MIEMBRO_INI + $MAX_SIM_MEM - 1 ))
exe_group=1
#Run several instances of real.exe using the spawner.
while [ $ini_mem -le $MIEMBRO_FIN ] ; do
   echo "Executing group number " $exe_group "Ini member = "$ini_mem " End member = "$end_mem
   ./spawn.exe ./real.exe $WRFTPROC $WRFDIR $ini_mem $end_mem
   #Set ini_mem and end_mem for the next round.
   ini_mem=$(( $ini_mem + $MAX_SIM_MEM ))
   end_mem=$(( $end_mem + $MAX_SIM_MEM ))
   if [ $end_mem -gt $MIEMBRO_FIN ] ; then
      end_mem=$MIEMBRO_FIN
   fi
   exe_group=$(( $exe_group + 1 )) #This is just to count the number of cycles performed. 
done

#Check if real was successfull
for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
   if [ -f $WRFDIR/$MIEM/rsl.error.0000 ] ; then
      mv $WRFDIR/$MIEM/rsl.error.0000 $WRFDIR/$MIEM/real.log
   fi
   grep SUCCESSFUL $WRFDIR/$MIEM/real.log
   if [ $? -ne 0 ] ; then
      echo "$MIEM"
      dispararError 9 "real.exe"
   fi
done

#RUN DA_UPDATE_BC.EXE
#TODO: DA_UPDATE_BC is a serial program, but we can include a dummy call to MPI_INITIALIZE / FINALIZE to launche it with the spawn.exe utility.
#En los scripts de la K computer habia un dummy.f90 que solo llamaba a MPI_INIT y MPI_FINALIZE quiza se pueda armar un script que ejecute primero el dummy.f90 y luego el programa serial desado. 
echo "Running da_update_bc"
for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
   cd $WRFDIR/$MIEM
   if [ $PASO -gt 0 ] ; then
      cp $NAMELISTDIR/parame* .
      mv $WRFDIR/$MIEM/wrfinput_d01 $WRFDIR/$MIEM/wrfinput_d01.org
      cp $HISTDIR/ANAL/$(date -u -d "$FECHA_INI UTC +$(($ANALISIS_FREC*$PASO)) seconds" +"%Y%m%d%H%M%S")/anal$(printf %05d $((10#$MIEM))) $WRFDIR/$MIEM/wrfinput_d01
      ln -sf $WRFDIR/code/da_update_bc.exe $WRFDIR/$MIEM/da_update_bc.exe
      time $WRFDIR/$MIEM/da_update_bc.exe > ./da_update_bc_${PASO}_${MIEM}.log  &
   fi
done
time wait  #Wait for the different instances of da_update_bc.

export OMP_NUM_THREADS=1
#Run WRF.EXE
cd $WRFDIR
ini_mem=$(( $MIEMBRO_INI ))
end_mem=$(( $MIEMBRO_INI + $MAX_SIM_MEM - 1 ))
exe_group=1
#Run several instances of wrf.exe using the spawner.
while [ $ini_mem -le $MIEMBRO_FIN ] ; do
   echo "Executing group number " $exe_group "Ini member = "$ini_mem " End member = "$end_mem
   ./spawn.exe ./wrf.exe $WRFTPROC $WRFDIR $ini_mem $end_mem
   #Set ini_mem and end_mem for the next round.
   ini_mem=$(( $ini_mem + $MAX_SIM_MEM ))
   end_mem=$(( $end_mem + $MAX_SIM_MEM ))
   if [ $end_mem -gt $MIEMBRO_FIN ] ; then
      end_mem=$MIEMBRO_FIN
   fi
   exe_group=$(( $exe_group + 1 )) #This is just to count the number of cycles performed. 
done

#Check if real was successfull
for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
   if [ -f $WRFDIR/$MIEM/rsl.error.0000 ] ; then
      mv $WRFDIR/$MIEM/rsl.error.0000 $WRFDIR/$MIEM/wrf.log
   fi
   grep SUCCESSFUL $WRFDIR/$MIEM/wrf.log
   if [ $? -ne 0 ] ; then
      echo "$MIEM"
      dispararError 9 "wrf.exe"
   fi
done


#Copiamos los archivos del Guess correspondientes a la hora del analisis.
if [[ ! -z "$GUARDOGUESS" ]] && [[ $GUARDOGUESS -eq 1 ]] ; then
   for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
       OUTPUTPATH="$HISTDIR/GUES/$(date -u -d "$FECHA_ANALISIS UTC" +"%Y%m%d%H%M%S")/"
       mkdir -p $OUTPUTPATH
       echo "Copying file $WRFDIR/$MIEM/wrfout_d01_$(date -u -d "$FECHA_ANALISIS UTC" +"%Y-%m-%d_%T" )"
       cp $WRFDIR/$MIEM/wrfout_d01_$(date -u -d "$FECHA_ANALISIS UTC" +"%Y-%m-%d_%T") $OUTPUTPATH/gues$(printf %05d $((10#$MIEM)))
       mv $WRFDIR/$MIEM/*.log                                                         $OUTPUTPATH
       mv $WRFDIR/$MIEM/namelist*                                                     $OUTPUTPATH
   done
fi


