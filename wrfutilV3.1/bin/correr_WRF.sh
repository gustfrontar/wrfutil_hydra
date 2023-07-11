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


##### FIN INICIALIZACION ######




###
# Creacion y en colamiento del script
###
cd $BASEDIR/WRF
#script de ejecucion
read -r -d '' QSCRIPTCMD << "EOF"
MIEM=$(printf "%02g" $((10#$MIEMBRO_FIN-10#$MIEMBRO_INI)))
test $ARRAYID && MIEM=$(printf "%02g" $ARRAYID)
ulimit -s unlimited
ulimit -l unlimited
echo "Procesando Miembro $MIEM"
read -r IY IM ID IH Im  <<< $(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN)) minutes" +"%Y %m %d %H %M")
read -r FY FM FD FH Fm  <<< $(date -u -d "$FECHA_INI UTC +$((10#$CICLO+$PLAZO)) hours+$((10#$CICLO_MIN+10#$PLAZO_MIN)) minutes" +"%Y %m %d %H %M")
cd $WRFDIR
CNL=`ls $WRFUTILDIR/namelists/namelist.input.${NAMELIST_TIPO}* |wc -w`
NAMELISTFILE="$WRFDIR/$MIEM/namelist.input"
SIHVIGILAFILE="$WRFDIR/$MIEM/SIHvigila_d01.txt"
PATH=$PATH:$CONDADIR
rm -r "$WRFDIR/$MIEM/WRFOUT"
mkdir -p "$WRFDIR/$MIEM/WRFOUT"
tar  xf WRF.tar -C "$WRFDIR/$MIEM"
NL=$(printf '%03d' $((10#$MIEM % 10#$CNL)))
TIPO=det
test  $NAMELIST_TIPO != "det"  && TIPO=${NAMELIST_TIPO}.$NL
cp $WRFUTILDIR/namelists/namelist.input.${TIPO} $NAMELISTFILE
if [[ $NOMBRE == "deterministico" ]]
then
	mkdir -p "$WRFDIR/$MIEM/SIHVIGILA"
	mkdir -p "$WRFDIR/$MIEM/UH"
	cp $WRFUTILDIR/namelists/SIHvigila_d01.txt $SIHVIGILAFILE
fi

IOTYPE=11
[[  $REALPROC -gt 1 ]] && IOTYPE=102

sed -i -e "s|__RUN_HOURS__|$PLAZO|g" $NAMELISTFILE
sed -i -e "s|__RUN_MINUTES__|$PLAZO_MIN|g" $NAMELISTFILE
sed -i -e "s|__START_YEAR__|$IY|g" $NAMELISTFILE
sed -i -e "s|__START_MONTH__|$IM|g" $NAMELISTFILE
sed -i -e "s|__START_DAY__|$ID|g" $NAMELISTFILE
sed -i -e "s|__START_HOUR__|$IH|g" $NAMELISTFILE
sed -i -e "s|__START_MINUTE__|$Im|g" $NAMELISTFILE
sed -i -e "s|__END_YEAR__|$FY|g" $NAMELISTFILE
sed -i -e "s|__END_MONTH__|$FM|g" $NAMELISTFILE
sed -i -e "s|__END_DAY__|$FD|g" $NAMELISTFILE
sed -i -e "s|__END_HOUR__|$FH|g" $NAMELISTFILE
sed -i -e "s|__END_MINUTE__|$Fm|g" $NAMELISTFILE
sed -i -e "s|__INTERVALO_WRF__|$INTERVALO_WRF|g" $NAMELISTFILE
sed -i -e "s|__INTERVALO_WPS__|$(($INTERVALO_WPS*60))|g" $NAMELISTFILE
sed -i -e "s|__E_WE__|$E_WE|g" $NAMELISTFILE
sed -i -e "s|__E_SN__|$E_SN|g" $NAMELISTFILE
sed -i -e "s|__DX__|$DX|g" $NAMELISTFILE
sed -i -e "s|__DY__|$DY|g" $NAMELISTFILE
sed -i -e "s|__IOTYPE__|$IOTYPE|g" $NAMELISTFILE
sed -i -e "s|__NUMTILE__|$NUMTILE|g" $NAMELISTFILE
sed -i -e "s|__NIOT__|$NIOT|g" $NAMELISTFILE
sed -i -e "s|__NIOG__|$NIOG|g" $NAMELISTFILE
sed -i -e "s|__METLEV__|$METLEV|g" $NAMELISTFILE
sed -i -e "s| !nproc_x| nproc_x                             = $NPROC_X,|g" $NAMELISTFILE
sed -i -e "s| !nproc_y| nproc_y                             = $NPROC_Y,|g" $NAMELISTFILE


cd $BASEDIR/WRF/$MIEM
[ -f envvars.sh ] && source envvars.sh
[ -f entorno.sh ] && source entorno.sh

###
# Si es asimilacion entonces hay que hacer el MERGE
##
if [[ $PASO -gt 0 && $ASIMILACION -ne 0 ]]
then
        echo "Asimilando miembro $MIEM"
	cd $BASEDIR/WRF/$MIEM/
	cp $WRFUTILDIR/namelists/parame* .
        mv $BASEDIR/WRF/$MIEM/wrfinput_d01 $BASEDIR/WRF/$MIEM/wrfinput_d01.org
        cp $BASEDIR/HIST/ANA/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/anal$(printf %05d $((10#$MIEM))) $BASEDIR/WRF/$MIEM/wrfinput_d01
        ln -sf $BASEDIR/LETKF/da_update_bc.exe $BASEDIR/WRF/$MIEM/da_update_bc.exe
        $BASEDIR/WRF/$MIEM/da_update_bc.exe
fi

echo "Corriendo WRF en Miembro $MIEM"
cd $BASEDIR/WRF/$MIEM
export OMP_NUM_THREADS=1
export I_MPI_OFI_PROVIDER=psm2
export I_MPI_DEBUG=0
export I_MPI_ADJUST_REDUCE=1
export I_MPI_PIN_DOMAIN="core"
export I_MPI_PIN_ORDER="scatter"
export KMP_AFFINITY='verbose,granularity=core,scatter'
export WRFIO_NCD_LARGE_FILE_SUPPORT=1


mkdir -p $BASEDIR/HIST/LOGS/ENSAMBLE/$MIEM
LOGFILE=$BASEDIR/HIST/LOGS/ENSAMBLE/$MIEM/historial.txt
[[ -f $LOGFILE ]] || echo "Miembro ; Comienzo ; Fin ; Duracion ; Exit Code ; Tamano wrfout total ; Ultima linea del rsl.error.0000; Lista de nodos" >> $LOGFILE

comienzo=$(date  +"%Y%m%d_%H%M%S")
comienzoT=$(date +"%s")
$MPIEXE  $MPIARGS $BASEDIR/WRF/$MIEM/wrf.exe
excod=$?
res="ERROR"
test=$(tail -n1 $BASEDIR/WRF/$MIEM/rsl.error.0000 | grep SUCCESS ) && res="OK"
fin=$(date  +"%Y%m%d_%H%M%S")
finT=$(date +"%s")
duracion=$(($finT-$comienzoT))
tamano=$(du -shc $BASEDIR/WRF/$MIEM/WRFOUT/wrfout*| tail -n1 | cut -f1)
LLRSL=$(tail -n1 $BASEDIR/WRF/$MIEM/rsl.error.0000)
cp $BASEDIR/WRF/$MIEM/namelist.input $BASEDIR/HIST/LOGS/ENSAMBLE/$MIEM/${comienzo}_namelist.input
#REGISTRO="${comienzo};${fin};${duracion};${res};${tamano};${LLRSL};$SLURM_NODEID"
REGISTRO="${MIEM};${comienzo};${fin};${duracion};${res};${tamano};${LLRSL};$SLURM_JOB_NODELIST"
echo $REGISTRO >> $LOGFILE
#[[ $excod -ne 0 ]] &&  echo "wrf.exe" > $BASEDIR/FATAL_ERROR
[[ $excod -ne 0 ]] && dispararError 9 "wrf.exe"
touch $WRFDIR/$MIEM/WRFOUT/wrfout_termine

if [[ $NOMBRE == "deterministico" ]]
then
	touch $WRFDIR/$MIEM/UH/wrfout_termine
	source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil
	python $WRFUTILDIR/bin/python/File_Reducer_SIHpp.py $WRFDIR/$MIEM/SIHVIGILA/
fi

if [[ ! -z "$GUARDOGUESS" ]] && [[ $GUARDOGUESS -eq 1 ]]
then 
	WRFOUT=WRFOUT/wrfout_d01_$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN+$ANALISIS)) minutes" +"%Y-%m-%d_%H_%M_%S")
	dirname=$BASEDIR/HIST/GUESS/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN+$ANALISIS)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/$MIEM/
	mkdir -p "$dirname"
	cp $BASEDIR/WRF/$MIEM/$WRFOUT  $dirname
fi 
EOF

# Parametros de encolamiento
QDEPEND_NAME=$DEPWRF
QDEPENDTYPE=$TYPEWRF
QPROC_NAME=WRF_$1_$PASO
QPROC=$WRFPROC
TPROC=$WRFTPROC
QTHREADS=$WRFTHREADS
QARRAY="$((10#$MIEMBRO_INI))-$((10#$MIEMBRO_FIN))"
QWALLTIME=$WRFWALLTIME
QEXCLU=1

# Encolar
queue

