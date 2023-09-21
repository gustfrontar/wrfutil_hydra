#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

### CONFIGURACION
source $BASEDIR/conf/config.env  #BASEDIR tiene que estar seteado como variable de entorno.

[ ! -f $BASEDIR/lib/errores.env ] && exit 1
source $BASEDIR/lib/errores.env
[ ! -f "$BASEDIR/conf/$EXPMACH" ] && dispararError 4 "$BASEDIR/conf/$EXPMACH"
source $BASEDIR/conf/$EXPMACH
[ ! -f "$BASEDIR/conf/$EXPCONF" ] && dispararError 4 "$BASEDIR/conf/$EXPCONF"
source $BASEDIR/conf/$EXPCONF
[ ! -f "$BASEDIR/conf/$EXPMODELCONF" ] && dispararError 4 "$BASEDIR/conf/$EXPMODELCONF"
source $BASEDIR/conf/$EXPMODELCONF
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
   MULTIMODEL_CONF_INI=0;MULTIMODEL_CONF_FIN=0
   echo "We will run a single configuration ensemble"
else 
   echo "We will run a multimodel configuration ensemble with $NCONF configurations"
   MULTIMODEL_CONF_INI=1;MULTIMODEL_CONF_FIN=$NCONF
fi

for ICONF in $(seq -w $MULTIMODEL_CONF_INI $MULTIMODEL_CONF_FIN) ; do
   WRFCONF=$(printf '%03d' $((10#$ICONF % 10#$MULTIMODEL_CONF_FIN)))
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

ulimit -s unlimited
OMP_NUM_THREADS=$REALTHREADS
OMP_STACKSIZE=512M

for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
   cd $WRFDIR
   [[ -d $WRFDIR/$MIEM ]] && rm -r $WRFDIR/$MIEM
   cd $WRFDIR/$MIEM
   if [ $MULTIMODEL -eq 0 ] ; then
      NLCONF=000
   else
      NLCONF=$(printf '%03d' $((10#$MIEM % 10#$NCONF)))
   fi
   cp $WRFDIR/namelist.input.${NLCONF} $WRFDIR/$MIEM/namelist.input
   ln -sf $WRFDIR/code/* .
   ln -sf $HISTDIR/WPS/met_em/$MIEM/met_em* $WRFDIR/$MIEM/
done

$MPIEXE $WRFDIR/$MIEM/real.exe  ####TODO Aca hay que reemplazar por el spawner. 

do MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
   cd $WRFDIR/$MIEM
   mv rsl.error.0000 ./real_${PASO}_${MIEM}.log
   EXCOD=$?
   [[ $EXCOD -ne 0 ]] && dispararError 9 "real.exe"

   if [ $PASO -gt 0 ] ; then
      echo "Corremos el update_bc"
      cp $NAMELISTDIR/parame* .
      mv $WRFDIR/$MIEM/wrfinput_d01 $WRFDIR/$MIEM/wrfinput_d01.org
      cp $HISTDIR/ANAL/$(date -u -d "$FECHA_INI UTC +$(($ANALISIS_FREC*$PASO)) seconds" +"%Y%m%d%H%M%S")/anal$(printf %05d $((10#$MIEM))) $WRFDIR/$MIEM/wrfinput_d01
      ln -sf $WRFDIR/code/da_update_bc.exe $WRFDIR/$MIEM/da_update_bc.exe
      time $WRFDIR/$MIEM/da_update_bc.exe > ./da_update_bc_${PASO}_${MIEM}.log  &
   fi
done
time wait  #Wait for the different instances of da_update_bc.

export OMP_NUM_THREADS=1

$MPIEXE  $MPIARGS $WRFDIR/$MIEM/wrf.exe  ####TODO Aca hay que reemplazar por el spawner

do MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
   cd $WRFDIR/$MIEM
   if [! $(tail -n1 $WRFDIR/$MIEM/rsl.error.0000 | grep SUCCESS ) ] ; then
      dispararError 9 "wrf.exe"
   fi
   mv rsl.error.0000 ./wrf_${PASO}_${MIEM}.log
done


#cd $WRFDIR

# Parametros de encolamiento
## Calculo cuantos miembros hay que correr
#QNODE=$WRFNODE
#QPROC=$WRFPROC
#TPROC=$WRFTPROC
#QTHREAD=$WRFTHREAD
#QWALLTIME=$WRFWALLTIME
#QPROC_NAME=GUESS_${PASO}

## Encolar
#echo "Tiempo en correr el real, update bc y wrf"
#queue $MIEMBRO_INI $MIEMBRO_FIN 
#time check_proc $MIEMBRO_INI $MIEMBRO_FIN

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


