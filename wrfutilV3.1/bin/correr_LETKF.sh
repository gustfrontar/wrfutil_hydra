#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF

        Ud. deberia usar este escript de la siguiente manera:
                $0 <nombre entorno >
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${1:?"$USO"}  ${WRFUTILDIR:?"$USO"}


### PARAMETROS

export NOMBRE=$1

### CONFIGURACION

[ ! -f $WRFUTILDIR/lib/errores.env ] && exit 1
source $WRFUTILDIR/lib/errores.env
CONFIG=$WRFUTILDIR/config.env
[ -f $WRFUTILDIR/RUNs/$NOMBRE/config.env ] && CONFIG=$WRFUTILDIR/RUNs/$NOMBRE/config.env
CONFIG=$(readlink -e ${CONFIG} 2> /dev/null)
[ -z "$CONFIG" ] && dispararError 4 "config.env"
source $CONFIG  $NOMBRE
[ ! -d "$BASEDIR" ] && dispararError 7 "$NOMBRE"
[ ! -f "$BASEDIR/$EXPCONF" ] && dispararError 4 "$BASEDIR/$EXPCONF"
source $BASEDIR/$EXPCONF
[ ! -f "$BASEDIR/$EXPDEP" ] && dispararError 4 "$BASEDIR/$EXPDEP"
source $BASEDIR/$EXPDEP
[ ! -f "$BASEDIR/$EXPASIM" ] && dispararError 4 "$BASEDIR/$EXPASIM"
source $BASEDIR/$EXPASIM


######### FIN VERIFICACION 

cd $BASEDIR/LETKF
tar xf LETKF.tar
cp $WRFUTILDIR/namelists/letkf.namelist $LETKFDIR/letkf.namelist
NAMELISTFILE=$LETKFDIR/letkf.namelist
NSLOTS=$(($OBSWIN/$OBSFREC+1))
NBSLOT=$(( $NSLOTS-$(($ANALISIS-$DESPLAZAMIENTO))/60 ))
NBV=$((10#$MIEMBRO_FIN-10#$MIEMBRO_INI+1))
#N_RADAR=${#RADAR_LIST[@]}
#N_RADAR=${#SO_INSTRUMENT_LIST[@]}
N_RADAR=($(echo $SO_INSTRUMENT_LIST | sed 's/,/\n/g') )
N_RADAR=${#N_RADAR[@]}

sed -i -e "s|__COV_INFL_MUL__|$COV_INFL_MUL|g" $NAMELISTFILE
sed -i -e "s|__SP_INFL_ADD__|$SP_INFL_ADD|g" $NAMELISTFILE
sed -i -e "s|__RELAX_ALPHA_SPREAD__|$RELAX_ALPHA_SPREAD|g" $NAMELISTFILE
sed -i -e "s|__RELAX_ALPHA__|$RELAX_ALPHA|g" $NAMELISTFILE
sed -i -e "s|__NSLOTS__|$NSLOTS|g" $NAMELISTFILE
sed -i -e "s|__NBSLOT__|$NBSLOT|g" $NAMELISTFILE
sed -i -e "s|__NBV__|$NBV|g" $NAMELISTFILE
sed -i -e "s|__SIGMA_OBS__|$SIGMA_OBS|g" $NAMELISTFILE
sed -i -e "s|__SIGMA_OBSV__|$SIGMA_OBSV|g" $NAMELISTFILE
sed -i -e "s|__SIGMA_OBS_RADAR__|$SIGMA_OBS_RADAR|g" $NAMELISTFILE
sed -i -e "s|__N_RADAR__|$N_RADAR|g" $NAMELISTFILE

#script de ejecucion
read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
ulimit -l unlimited
source envvars.sh
source $WRFUTILDIR/lib/recupero_miem.sh
NSLOTS=$(($OBSWIN/$OBSFREC+1))
MAXFALLAS=5
nextMiem(){
   miem=$1
   count=$(($2+1))
   [[ $count -gt $MAXFALLAS ]]  && dispararError 8 "Fallaron demasiados WRFs"
   WRFCONT=$(ls $BASEDIR/WRF/$miem/WRFOUT/wrfout_d01*| wc -l)
   if [[ $WRFCONT -lt $NSLOTS ]]
   then
       newmiem=${RECUPEROMAP[$MIEM]}
       miem=$(nextMiem "$newmiem" "$count" )
   fi
   echo $miem
}

contador=0
for MIEM in $(seq -f "%02g" $MIEMBRO_INI $MIEMBRO_FIN )
do
	echo "Procesando Miembro $MIEM"
	cd "$BASEDIR/WRF/$MIEM"
	WRFCONT=$(ls WRFOUT/wrfout_d01*| wc -l)
        MIEMORI=$MIEM	
	if [[ $WRFCONT -lt $NSLOTS ]] 
	then  
		contador=$(($contador+1))
		MIEMORI=$(nextMiem "$newmiem" "$count" )
		echo "$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN+$ANALISIS)) minutes" +"%Y-%m-%d_%H_%M_%S") -> Reemplazamos el Miembro $MIEM por $MIEMORI" >> $BASEDIR/LOGS/recuperoAsim.log
	        mkdir -p $BASEDIR/WRF/$MIEM/WRFOUT	
	        #cp $BASEDIR/WRF/$MIEMORI/WRFOUT/* $BASEDIR/WRF/$MIEM/WRFOUT/
	fi 
	[[ $contador -gt $MAXFALLAS ]]  && dispararError 8 "Fallaron demasiados WRFs"
	for NUM in $(seq -f "%02g" 0 $(($NSLOTS-1)))
	do
		SLOT_TIME=$(( ($OBSFREC)*10#$NUM ))
		WRFOUT=WRFOUT/wrfout_d01_$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN+$ANALISIS-$DESPLAZAMIENTO+$SLOT_TIME)) minutes" +"%Y-%m-%d_%H_%M_%S")
		WRFOUT=$BASEDIR/WRF/$MIEM/$WRFOUT
		ln -sf $WRFOUT ${LETKFDIR}/gs$(printf %02d $((10#$NUM+1)))$(printf %05d $((10#$MIEM)))
	done
done
WRFOUT=WRFOUT/wrfout_d01_$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN+$ANALISIS)) minutes" +"%Y-%m-%d_%H_%M_%S")
cp  $BASEDIR/WRF/$MIEMBRO_FIN/$WRFOUT ${LETKFDIR}/gues$(printf %05d $((10#$MIEMBRO_FIN+1)))
cp  $BASEDIR/WRF/$MIEMBRO_FIN/$WRFOUT ${LETKFDIR}/anal$(printf %05d $((10#$MIEMBRO_FIN+1)))

###
# Linkeamos las observaciones que ya estan procesadas
###

DIROBSRAD=$BASEDIR/HIST/OBS/${NOMBRE}/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN+$ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/
DIROBS=$DIROBSASIM
rm ${LETKFDIR}/rad*
rm ${LETKFDIR}/obs*.dat
for NUM in $(seq -f "%02g" 0 $(($NSLOTS-1)))
do
	SLOT_TIME=$(( ($OBSFREC)*10#$NUM ))
	ln -sf $DIROBS/obs_$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN+$ANALISIS-$DESPLAZAMIENTO+$SLOT_TIME)) minutes" +"%Y%m%d%H%M").dat ${LETKFDIR}/obs$(printf %02d $((10#$NUM+1))).dat
	echo "antes del IF radar     ASIM_RADAR=$ASIM_RADAR y NUM= ${NUM} y radar_slot=${RADAR_SLOTSNUM[@]}"
        if [[ "$ASIM_RADAR" == "1" &&  " ${RADAR_SLOTSNUM[@]} " =~ " ${NUM} "  ]]
        then
		echo "Entramos el IF radar"
		count=1
		#for radar in ${RADAR_LIST[@]}
                for radar in $(echo $SO_INSTRUMENT_LIST | sed 's/,/\n/g') 
		do
			radarFile=$DIROBSRAD/${radar}_$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN+$ANALISIS-$DESPLAZAMIENTO+$SLOT_TIME)) minutes" +"%Y%m%d_%H%M%S").dat	
			echo "Radar Files: $radarFile"
			[[ -e $radarFile ]] && ln -sf $radarFile ${LETKFDIR}/rad$(printf %02d $((10#$NUM+1)))$(printf %02d $count).dat
			count=$(($count+1))
		done
        fi
done

###
# Ahora si corremos el LETKF
###
cd $BASEDIR/LETKF/
OMP_NUM_THREADS=$LETKFTHREADS
OMP_STACKSIZE=512M
cd $BASEDIR/LETKF
$MPIEXE -rr  $BASEDIR/LETKF/letkf.exe
LETKF_ERROR=$?
if [[ $LETKF_ERROR != 0 ]]
then
        echo $LETKF_ERROR > $BASEDIR/FATAL.ERROR
	echo "ERROR: El letkf termino con error, se interrumpe el el ciclo"
        exit 10
fi

##
# Guaramos Todos los analisis
##
echo "Corriendo Guadado de Analisis Ensamble $MIEMBRO_INI - $MIEMBRO_FIN"
WRFOUT=WRFOUT/wrfout_d01_$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y-%m-%d_%H_%M_%S")
FECHAANAL=$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y-%m-%d_%H:%M:%S")
FECHAOBS=$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y%m%d_%H%M%S")
DIRANAL=$BASEDIR/HIST/ANA/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours +$((10#$CICLO_MIN + $ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/
for MIEM in $(seq -f "%02g" $MIEMBRO_INI $MIEMBRO_FIN )
do
        echo "Procesando Miembro $MIEM"
        mkdir -p  "$DIRANAL"
        mv  $BASEDIR/WRF/$MIEM/$WRFOUT $DIRANAL/anal$(printf %05d $((10#$MIEM)))
        ncatted -h -a 'START_DATE',global,m,c,$FECHAANAL $DIRANAL/anal$(printf %05d $((10#$MIEM)))
        ncatted -h -a 'SIMULATION_START_DATE',global,m,c,$FECHAANAL $DIRANAL/anal$(printf %05d $((10#$MIEM)))
done

##
# Guardamos la media y el Spread del ensamble
##
mv $BASEDIR/LETKF/*_sp $DIRANAL
mv $BASEDIR/LETKF/gues$(printf %05d $((10#$MIEMBRO_FIN+1))) $DIRANAL
mv $BASEDIR/LETKF/anal$(printf %05d $((10#$MIEMBRO_FIN+1))) $DIRANAL
cp $BASEDIR/LETKF/obs.dat $BASEDIR/HIST/OBS/$NOMBRE/$FECHAOBS/obs_${FECHAOBS}_asimiladas.dat
cp $BASEDIR/LETKF/obs.dat $DIRANAL
ncatted -h -a 'START_DATE',global,m,c,$FECHAANAL $DIRANAL/anal$(printf %05d $((10#$MIEMBRO_FIN+1)))
ncatted -h -a 'SIMULATION_START_DATE',global,m,c,$FECHAANAL $DIRANAL/anal$(printf %05d $((10#$MIEMBRO_FIN+1)))

##
# Salvamos un NOUT
##
mv $BASEDIR/LETKF/NOUT-000 $DIRANAL
##
# Salvamos met.dir que es el directorio que se uso para los mets
##
cp $BASEDIR/WRF/mets.dir $DIRANAL
EOF

# Parametros de encolamiento
QDEPEND_NAME=$DEPLETKF
QDEPENDTYPE=$TYPELETKF
QPROC_NAME=LETKF_$1_$PASO
QPROC=$LETKFPROC
QTHREADS=$LETKFTHREADS
QARRAY=""
QWALLTIME=$LETKFWALLTIME
QEXCLU=1

# Encolar
queue



