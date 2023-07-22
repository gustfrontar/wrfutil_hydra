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

ulimit -s unlimited

cd $BASEDIR/WRF/
MIEM=$(printf "%02g" $((10#$MIEMBRO_FIN-10#$MIEMBRO_INI)))
read -r -d '' QSCRIPTCMD << "EOF"


### INICIO ###
# Set fcst hour
HORA=0
test ! -z $ARRAYID && HORA=$(($ARRAYID%$(($PLAZO+1))))

# Set ensemble member
MIEM=00
TYPE=DET
if [[ $((10#$ARRAYCNT-10#$PLAZO-1)) -gt 0 ]]; then
   test ! -z $ARRAYID && MIEM=$(printf  "%02g" $(($MIEMBRO_INI+$(($ARRAYID/$(($PLAZO+1)))))))
   TYPE=ENS
fi

HORA2=$(printf "%02g" $HORA)
export UPP_DIR=/data/
ulimit -s unlimited
ulimit -l unlimited
echo "Procesando Plazo $HORA Miembro $MIEM"
read -r IY IM ID IH Im  <<< $(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN)) minutes" +"%Y %m %d %H %M")
read -r FY FM FD FH Fm  <<< $(date -u -d "$FECHA_INI UTC +$((10#$CICLO+$PLAZO)) hours+$((10#$CICLO_MIN+10#$PLAZO_MIN)) minutes" +"%Y %m %d %H %M")
read -r PY PM PD PH Pm  <<< $(date -u -d "$FECHA_INI UTC +$((10#$CICLO+$HORA)) hours+$((10#$CICLO_MIN)) minutes" +"%Y %m %d %H %M")
[ -f $BASEDIR/WRF/$MIEM/envvars.sh ] && source $BASEDIR/WRF/$MIEM/envvars.sh
[ -f $BASEDIR/WRF/$MIEM/entorno.sh ] && source $BASEDIR/WRF/$MIEM/entorno.sh
export RUNDIR="$BASEDIR/WRF/$MIEM/"
export FECHA=`echo $FECHA_INI$CICLO| cut -d'/' -s -f1,2,3 --output-delimiter=''`
export PATHOUT=$BASEDIR/HIST/POST/${IY}${IM}${ID}_${IH}${Im}00/$NOMBRE/$MIEM/
export PATHOUT_UH=$BASEDIR/HIST/POST/${IY}${IM}${ID}_${IH}${Im}00/$NOMBRE/$MIEM/UH/
mkdir -p $PATHOUT
mkdir -p $PATHOUT_UH
#export MIEM
PATH=$PATH:$CONDADIR
source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil
export FILEIN=$RUNDIR/WRFOUT/wrfout_d01_$(date -d  "$FECHA_INI +$((10#$CICLO)) hours+$((10#$CICLO_MIN)) minutes + $HORA hours" +'%Y-%m-%d_%H_%M')_00
export FILEIN_UH=$RUNDIR/UH/wrfout_d01_$(date -d  "$FECHA_INI +$((10#$CICLO)) hours+$((10#$CICLO_MIN)) minutes + $HORA hours" +'%Y-%m-%d_%H_%M')_00



### UPP ###
if [[ $MIEM == "00" ]]; then 
	export startdate=$FECHA
	export incrementhr=1
	export DOMAINPATH=$UPPDIR/$HORA/
	export fhr=$HORA
	export lastfhr=$HORA
	rm -r ${DOMAINPATH}/postprd
	mkdir -p ${DOMAINPATH}/postprd

	export UNIPOST_HOME=${UPP_DIR}/UPPV4.1
	export POSTEXEC=${UNIPOST_HOME}/exec
	export SCRIPTS=${UNIPOST_HOME}/scripts
	export modelDataPath=${RUNDIR}/WRFOUT
	export paramFile=${DOMAINPATH}/wrf_cntrl.parm # grib1 (WRF only)
	[ -f $BASEDIR/WRF/$MIEM/envvars.sh ] && source $BASEDIR/WRF/$MIEM/envvars.sh
	[ -f $BASEDIR/WRF/$MIEM/entorno.sh ] && source $BASEDIR/WRF/$MIEM/entorno.sh
	if [[ -f "$BASEDIR/WRF/wrf_cntrl.parm" ]]; then
		cp $BASEDIR/WRF/wrf_cntrl.parm $paramFile
	else
		cp $WRFUTILDIR/templates/wrf_cntrl.parm $paramFile
	fi
	if [[ -f "$BASEDIR/WRF/run_unipost" ]]; then
		cp $BASEDIR/WRF/run_unipost $DOMAINPATH
	else
		cp $WRFUTILDIR/templates/run_unipost $DOMAINPATH
	fi
	if [[ -f "$BASEDIR/WRF/alpha.npy" ]]; then
		cp $BASEDIR/WRF/alpha.npy $DOMAINPATH
	else
		cp $WRFUTILDIR/templates/alpha_$NOMBRE.npy $DOMAINPATH/alpha.npy
	fi


	cd $DOMAINPATH
	run_unipost
	
	source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil
	python $WRFUTILDIR/bin/python/rot_V.py $DOMAINPATH/alpha.npy $DOMAINPATH/postprd/WRFPRS_d01.$HORA2 $PATHOUT/WRFPRS_d01.$HORA2


fi
EOF

# Parametros de encolamiento
QDEPEND_NAME=$DEPPOST
QDEPENDTYPE=$TYPEPOST
QPROC_NAME=POST_$1_$PASO
QPROC=1
QTHREADS=16
QARRAY=0-$(($(($(($PLAZO+1))*$(($MIEMBRO_FIN-$MIEMBRO_INI+1))))-1))%480
#QARRAY=0-$(($(($PLAZO+1))*$(($MIEMBRO_FIN-$MIEMBRO_INI+1))))
QWALLTIME=$POSTWALLTIME
QEXCLU=0

# Encolar
queue
