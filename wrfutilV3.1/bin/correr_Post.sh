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

#####
# extraerIndices ( hora, wrfprs , variables)
#       hora            -> hora del pronostico a la que se desea extraer las variables
#       wrfprs          -> archivo de donde se extraeran la variables
#       variables       ->
# Devuelve una un Expresion Regular para usar con el grep y extraer inventoryFile
extraerIndices(){
	IFSOLD=$IFS
	IFS=:
	while read -r  var nivel acum   
	do
		if [[ $acum -eq 3 ]]
		then 
			P1=0
		else
			P1=$((10#$1-10#$acum))
			[[ $P1 -lt 0 ]] && P1=0
		fi
		indice=$(wgrib $2 |grep -E ":${var}"| grep -E "P1=$P1" |grep -E "${nivel//\"/}" |cut -f1 -d":")
		
		IFS=$IFSOLD
		indice2=($indice)
		ord=0
		## Lineas muy dificiles !!! Resuleve ambiguedades para hora 0 y 1 de variables acumuladas
		[[  10#$acum -eq 3 ]] && ord=$(( ${#indice2[@]} > 1 ? $((${#indice2[@]}-1)) : 0 ))
		[[ 10#$1 -le 1 ]]  && [[ 10#$acum -gt 0 ]] && indice=${indice2[$ord]}
		for ind in $indice      # la ambiguedad de las descripcion de variables.
		do
			variable+="|$ind"
		done
		IFS=:
	done < $3
	IFS=$IFSOLD
	todasER=${variable//\|/\:\|}
	echo "^(${todasER:2})"
}
#####
# copiarVariables ( indiceER , wrfprsIN , wrfprsOUT)
# indiceER      -> salida de la funcion extraerIndices
# wrfprsIN      -> archivo con todas las variables originales
# wrfprsOUT     -> archivo donde seran copiadas las variables
copiarVariables(){
	#wgrib $2 -i -dwdgrib -o $3 <<< "$(wgrib $2 |grep -E  $1)"
	wgrib $2 -i -grib -h -H -o $3 <<< "$(wgrib $2 |grep -E  $1)"
}

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

if [[ $MIEM == "00" ]]
then
	# para el deterministico, 
	python $WRFUTILDIR/bin/python/write_post_smn_det_oper_sfc.py &
	python $WRFUTILDIR/bin/python/write_post_smn_det_oper_lev.py &
	wait
	python $WRFUTILDIR/bin/python/write_post_smn_det_oper_uh.py 

	#echo "Generando recorte HISTorico"
	#cd $PATHOUT
	#HORA3=$(printf "%03g" $HORA)
	#fileSFC=$PATHOUT/model.WRF_DET_4km.${IY}${IM}${ID}_${IH}${Im}00.${HORA3}.OPERSFC.nc
	#ncks -v XLAT,XLONG,PSFC,T2,Q2,PP,REFL1KM,REFL4KM,MDBZ,MCAPE,CIN,Umet10,Vmet10,SNOWNC,SLP,W $fileSFC ${fileSFC/OPER/HIST}
	#if [[ $(($HORA%6)) -eq 0 ]]
	#then
	#    fileLEV=$PATHOUT/model.WRF_DET_4km.${IY}${IM}${ID}_${IH}${Im}00.${HORA3}.OPERLEV.nc
	#    ncks -v XLAT,XLONG,T,Q,GEOPT,Umet,Vmet $fileLEV ${fileLEV/OPER/HIST}
	#fi

else
	# Para el ensamble
	export MIEM
	python $WRFUTILDIR/bin/python/write_post_smn_ens_oper_sfc.py  &
	python $WRFUTILDIR/bin/python/write_post_smn_ens_oper_lev.py  &
	wait

	echo "Generando recorte HISTorico"
	cd $PATHOUT
	#for file in $(ls $(eval echo $PATHOUT/*.{000..$PLAZO..6}.*LEV.nc))
	#do
	#	ncks -v XLAT,XLONG,T,Q,GEOPT,Umet,Vmet $file ${file/OPER/HIST}
	#done
	#for file in $(ls $PATHOUT/*SFC.nc)
	#do
	#	ncks -v XLAT,XLONG,PSFC,T2,Q2,PP,REFL1KM,REFL4KM,MDBZ,MCAPE,CIN,Umet10,Vmet10,SNOWNC,SLP,W $file ${file/OPER/HIST}
	#done

fi

### RECORTE HISTORICO
echo "Generando recorte HISTorico"
cd $PATHOUT
HORA3=$(printf "%03g" $HORA)
if [[ $TYPE == 'DET' ]]; then
	fileSFC=$PATHOUT/model.WRF_${TYPE}_4km.${IY}${IM}${ID}_${IH}${Im}00.${HORA3}.OPERSFC.nc
        fileLEV=$PATHOUT/model.WRF_${TYPE}_4km.${IY}${IM}${ID}_${IH}${Im}00.${HORA3}.OPERLEV.nc
else
        fileSFC=$PATHOUT/model.WRF_${TYPE}_4km.${IY}${IM}${ID}_${IH}${Im}00.${HORA3}.M${MIEM}_OPERSFC.nc
        fileLEV=$PATHOUT/model.WRF_${TYPE}_4km.${IY}${IM}${ID}_${IH}${Im}00.${HORA3}.M${MIEM}_OPERLEV.nc
fi
ncks -v XLAT,XLONG,PSFC,T2,Q2,PP,REFL1KM,REFL4KM,MDBZ,MCAPE,CIN,Umet10,Vmet10,SNOWNC,SLP,W $fileSFC ${fileSFC/OPER/HIST}
if [[ $(($HORA%6)) -eq 0 ]]; then
     ncks -v XLAT,XLONG,T,Q,GEOPT,Umet,Vmet $fileLEV ${fileLEV/OPER/HIST}
fi

### GUARDADO DE WRFOUT
if [[ $GUARDOWRFOUT -eq 1 ]]
then
	echo "Salvando las salidas del WRF"
	dirname=$BASEDIR/HIST/FCST/$(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN)) minutes" +"%Y%m%d_%H%M%S")/$NOMBRE/$MIEM/
	mkdir -p "$dirname"
	cp $BASEDIR/WRF/$MIEM/WRFOUT/wrfout_d01_${PY}-${PM}-${PD}_${PH}_${Pm}_00 $dirname/
fi
echo "Terminando"

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

	FECHA=$(echo ${FECHA_INI}_${CICLO}| cut -d'/' -s -f1,2,3 --output-delimiter='')
	cd $UPPDIR
	for productoCMD in $(ls guardo*.cmd)
	do 
		producto=$(echo ${productoCMD/guardo_/}| cut -d . -f1)
		echo "Corriendo $producto "
		patternER=$(extraerIndices $HORA2 ${PATHOUT}/WRFPRS_d01.$HORA2 $productoCMD)
		mkdir -p ${PATHOUT}/${producto}
		if [[ "$producto" == "meteofactory"  ]] 
		then
			copiarVariables $patternER ${PATHOUT}/WRFPRS_d01.${HORA2} ${PATHOUT}/${producto}/wrfprs_d01.${HORA2}
		else
			copiarVariables $patternER ${PATHOUT}/WRFPRS_d01.${HORA2} ${PATHOUT}/${producto}/wrfprs_${FECHA}00.0${HORA2}
		fi
	done
	

	#Edito el ctl del grib grande que se va a enviar para Wintem
	cp $WRFUTILDIR/templates/WRFPRS_d01.ctl ./ 
	sed -i -e "s|__FECHA__|$(date --date $(echo ${FECHA_INI} |cut -d '_' -f1) +"${CICLO}Z%d%b%Y")|g"  WRFPRS_d01.ctl
	cp WRFPRS_d01.ctl ${PATHOUT}

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
