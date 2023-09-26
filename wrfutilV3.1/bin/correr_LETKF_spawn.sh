#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

### CONFIGURACION
#Load experiment configuration
BASEDIR=$(pwd)/../
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/assimilation.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
##### FIN INICIALIZACION ######


LETKFDIRRUN=$LETKFDIR/run/
#rm -fr $LETKFDIR
#mkdir -p $LETKFDIR


#Descomprimimos el archivo .tar (si es que no fue descomprimido)
if [ ! -e $LETKFDIR/code/letkf.exe ] ; then
   echo "Descomprimiendo LETKF"
   mkdir -p $LETKFDIR/code/
   tar -xf $LETKFDIR/letkf.tar -C $LETKFDIR/code/
   #Si existe el namelist.input lo borramos para que no interfiera
   #con los que crea el sistema de asimilacion
   if [ -e $LETKFDIR/code/letkf.namelist ] ; then
      rm -f $LETKFDIR/code/letkf.namelist
   fi
fi

ln -sf $LETKFDIR/code/* $LETKFDIRUN

echo "Editing namelist"
cp $NAMELISTDIR/letkf.namelist $LETKFDIRRUN/letkf.namelist
NAMELISTFILE=$LETKFDIRRUN/letkf.namelist
NSLOTS=$(( ($ANALISIS_WIN_FIN-$ANALISIS_WIN_INI)/$ANALISIS_WIN_STEP+1))
NBSLOT=$(( ($ANALISIS_FREC-$ANALISIS_WIN_INI)/$ANALISIS_WIN_STEP+1))
NBV=$((10#$MIEMBRO_FIN-10#$MIEMBRO_INI+1))

#N_RADAR=($(echo $SO_INSTRUMENT_LIST | sed 's/,/\n/g') )
N_RADAR=0 #${#N_RADAR[@]}

sed -i -e "s|__COV_INFL_MUL__|$COV_INFL_MUL|g"             $NAMELISTFILE
sed -i -e "s|__SP_INFL_ADD__|$SP_INFL_ADD|g"               $NAMELISTFILE
sed -i -e "s|__RELAX_ALPHA_SPREAD__|$RELAX_ALPHA_SPREAD|g" $NAMELISTFILE
sed -i -e "s|__RELAX_ALPHA__|$RELAX_ALPHA|g"               $NAMELISTFILE
sed -i -e "s|__NSLOTS__|$NSLOTS|g"                         $NAMELISTFILE
sed -i -e "s|__NBSLOT__|$NBSLOT|g"                         $NAMELISTFILE
sed -i -e "s|__NBV__|$NBV|g"                               $NAMELISTFILE
sed -i -e "s|__SIGMA_OBS__|$SIGMA_OBS|g"                   $NAMELISTFILE
sed -i -e "s|__SIGMA_OBSV__|$SIGMA_OBSV|g"                 $NAMELISTFILE
sed -i -e "s|__SIGMA_OBS_RADAR__|$SIGMA_OBS_RADAR|g"       $NAMELISTFILE
sed -i -e "s|__N_RADAR__|$N_RADAR|g"                       $NAMELISTFILE

echo "Linking files"
FECHA_WINDOW_INI=$(date -u -d "$FECHA_INI UTC +$((($ANALISIS_FREC*$PASO)+$ANALISIS_WIN_INI)) seconds" +"%Y-%m-%d %T")
FECHA_ANALISIS=$(date -u -d "$FECHA_INI UTC +$((($ANALISIS_FREC*$PASO)+$ANALISIS_FREC)) seconds" +"%Y-%m-%d_%T")
for IMIEM in $(seq -f "%02g" $MIEMBRO_INI $MIEMBRO_FIN ) ; do
   for ISLOT in $(seq -f "%02g" 0 $(($NSLOTS-1))) ; do
      FECHA_SLOT=$(date -u -d "$FECHA_WINDOW_INI UTC +$(($ISLOT*$ANALISIS_WIN_STEP)) seconds" +"%Y-%m-%d_%T")
      ln -sf ${WRFDIR}/$IMIEM/wrfout_d01_$FECHA_SLOT ${LETKFDIRRUN}/gs$(printf %02d $((10#$ISLOT+1)))$(printf %05d $((10#$IMIEM)))
   done	   
done
cp  $WRFDIR/$MIEMBRO_FIN/wrfout_d01_$FECHA_ANALISIS ${LETKFDIRRUN}/guesemean
cp  $WRFDIR/$MIEMBRO_FIN/wrfout_d01_$FECHA_ANALISIS ${LETKFDIRRUN}/analemean

#Linkeamos las observaciones.
#DIROBSRAD=$HISTDIR/OBS/${NOMBRE}/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN+$ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/
#DIROBS=$DIROBSASIM
#rm ${LETKFDIR}/rad*
#rm ${LETKFDIR}/obs*.dat
#for NUM in $(seq -f "%02g" 0 $(($NSLOTS-1)))
#do
#       SLOT_TIME=$(( ($OBSFREC)*10#$NUM ))
#       ln -sf $DIROBS/obs_$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN+$ANALISIS-$DESPLAZAMIENTO+$SLOT_TIME)) minutes" +"%Y%m%d%H%M").dat ${LETKFDIR}/obs$(printf %02d $((10#$NUM+1))).dat
#       echo "antes del IF radar     ASIM_RADAR=$ASIM_RADAR y NUM= ${NUM} y radar_slot=${RADAR_SLOTSNUM[@]}"
#        if [[ "$ASIM_RADAR" == "1" &&  " ${RADAR_SLOTSNUM[@]} " =~ " ${NUM} "  ]]
#        then
#               echo "Entramos el IF radar"
#               count=1
#               #for radar in ${RADAR_LIST[@]}
#                for radar in $(echo $SO_INSTRUMENT_LIST | sed 's/,/\n/g') 
#               do
#                       radarFile=$DIROBSRAD/${radar}_$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN+$ANALISIS-$DESPLAZAMIENTO+$SLOT_TIME)) minutes" +"%Y%m%d_%H%M%S").dat        
#                       echo "Radar Files: $radarFile"
#                       [[ -e $radarFile ]] && ln -sf $radarFile ${LETKFDIR}/rad$(printf %02d $((10#$NUM+1)))$(printf %02d $count).dat
#                       count=$(($count+1))
#               done
#        fi
#done

#Ejecutamos el LETKF
ulimit -s unlimited
ulimit -l unlimited
cd $LETKFDIRRUN
OMP_NUM_THREADS=$LETKFTHREADS
OMP_STACKSIZE=512M
TCORES=$(( $LETKFNODE * $LETKFPROC ))
mpiexec -np $TCORES  ./letkf.exe
if [ $? != 0 ] ; then
   dispararError 9 letkf.exe
fi

##
# Guaramos los analisis
##
echo "Copying analysis files"
echo "Guardando analisis $MIEMBRO_INI - $MIEMBRO_FIN"
DIRANAL=$HISTDIR/ANAL/$(date -u -d "$FECHA_INI UTC +$((($ANALISIS_FREC*$PASO)+$ANALISIS_FREC)) seconds" +"%Y%m%d%H%M%S")
for MIEM in $(seq -f "%02g" $MIEMBRO_INI $MIEMBRO_FIN ) ; do
        mkdir -p  ${DIRANAL}/
        mv  $WRFDIR/$MIEM/wrfout_d01_$FECHA_ANALISIS $DIRANAL/anal$(printf %05d $((10#$MIEM)))
	echo "ncatted -h -a 'START_DATE',global,m,c,$FECHA_ANALYSIS            $DIRANAL/anal$(printf %05d $((10#$MIEM)))"
        ncatted -h -a 'START_DATE',global,m,c,$FECHA_ANALYSIS            $DIRANAL/anal$(printf %05d $((10#$MIEM)))
        ncatted -h -a 'SIMULATION_START_DATE',global,m,c,$FECHA_ANALISIS $DIRANAL/anal$(printf %05d $((10#$MIEM)))
done

#mv $LETKFDIRRUN/*_sp $DIRANAL
mv $LETKFDIRRUN/guesemean $DIRANAL   #Guess mean
mv $LETKFDIRRUN/analemean $DIRANAL   #Analysis mean
ncatted -h -a 'START_DATE',global,m,c,$FECHA_ANALISIS              $DIRANAL/analemean
ncatted -h -a 'SIMULATION_START_DATE',global,m,c,$FECHA_ANALISIS   $DIRANAL/analemean

#Copiamos las observaciones asimiladas.
cp $LETKFDIRRUN/obs.dat $DIRANAL

# Guardamos un NOUT
mv $LETKFDIRRUN/NOUT-000 $DIRANAL

echo "We are done!"


