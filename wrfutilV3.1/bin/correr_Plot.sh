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
source $CONFIG  $NOMBRE
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

cd $BASEDIR
cp $WRFUTILDIR/templates/PInteres.dat /data/salida/$AMBIENTE/intra/

DIRPLOT=$WRFUTILDIR/RUNs/HIST/PLOT/
STARTDATE="${FECHA_INI//\/}_${CICLO}0000"
mkdir -p ${DIRPLOT}/${STARTDATE}/${NOMBRE}/
###
# Creacion y en colamiento del script
###
read -r -d '' QSCRIPTCMD << "EOF"
export DIRPOST=$WRFUTILDIR/RUNs/HIST/POST/
export DIRWRF=$WRFUTILDIR/RUNs/HIST/FCST/
export DIRPLOT=$WRFUTILDIR/RUNs/HIST/PLOT/
export DIRSALIDA=/data/salida/$AMBIENTE/intra/
export STARTDATE="${FECHA_INI//\/}"
export OMP_NUM_THREAD=$PLOTTHREADS
PATH=$PATH:$CONDADIR
export BASEPLOTDIR=$WRFUTILDIR/bin/python/ploteos/
export PERFFREC=6
export PERFILES=$BASEPLOTDIR/perfiles.txt
#export CONCURR=$(($ICORE-4-$(wc -l < $BASEPLOTDIR/regionalizacion/zonaGraficadoWrf.csv)))
export CONCURR=$(($ICORE-4-$(wc -l < $BASEPLOTDIR/sectores_regiones/zonaGraficadoWrf.csv)))

source $CONDADIR/activate wrfutil2  || $CONDADIR/conda activate wrfutil2 
cd $BASEPLOTDIR
if [[ $ARRAYID -gt 0 ]]
then 
	python -W ignore horariosProno.py $STARTDATE $CICLO $ARRAYID  &
fi

#cd $BASEPLOTDIR/regionalizacion
cd $BASEPLOTDIR/sectores_regiones
IFSOLD=$IFS
IFS=,
if [[ $ARRAYID -gt 0 ]]
then 
	while IFS=' ' read -r name latN latS lonO lonE 
	do 
#		python -W ignore horariosPronoReg.py $STARTDATE $CICLO $ARRAYID $name $latN $latS $lonO $lonE &      
		python -W ignore horariosPronoSec.py $STARTDATE $CICLO $ARRAYID $name $latN $latS $lonO $lonE &
#	done < $BASEPLOTDIR/regionalizacion/zonasNCS.csv
	done < $BASEPLOTDIR/sectores_regiones/zonasNCS.csv
fi

cd $BASEPLOTDIR
perfiles=( )
parnum=7
while IFS=' ' read -r name idest lat lon 
do 
   for hora in $(seq -s "," 0 $PERFFREC $PLAZO)
   do 
     perfiles+=( $STARTDATE $CICLO $hora $name $idest $lat $lon )	
   done
done <<< $(python $BASEPLOTDIR/leerTablas.py $PERFILES ${STARTDATE}_${CICLO}0000 $NOMBRE)
IFS=$IFSOLD
len_perfiles=$((${#perfiles[@]}/$parnum))
ppn=$(($parnum*($len_perfiles/$ARRAYCNT) ))
[[ $(($len_perfiles%$ARRAYCNT)) -gt 0 ]] && ppn=$(($ppn+$parnum))
if [[ $ppn -eq 0 ]]
then
	if  [[ len_perfiles -lt $ARRAYID ]]
	then 
		python -W ignore $BASEPLOTDIR/perfilesProno.py ${perfiles[@]:$ARRAYID*$parnum:$parnum}
	fi
else
	echo ${perfiles[@]:$(($ARRAYID*$ppn)):$ppn} | xargs -P$(($CONCURR)) -n $parnum python -W ignore $BASEPLOTDIR/perfilesProno.py
fi

wait
echo "Terminando Ploteos: Hora:$ARRAYID y perfiles:$ppn"

EOF


# Parametros de encolamiento
QDEPEND_NAME=$DEPPLOT
QDEPENDTYPE=$TYPEPLOT
QPROC_NAME="PLOTH_${NOMBRE}_${PASO}"
QPROC=32
QTHREADS=1
QARRAY="0-$PLAZO"
QWALLTIME="20:00"
QEXCLU=1

# Encolar
queue


###
# Creacion y en colamiento del script
###
read -r -d '' QSCRIPTCMD << "EOF"
export DIRPOST=$WRFUTILDIR/RUNs/HIST/POST/
export DIRWRF=$WRFUTILDIR/RUNs/HIST/FCST/
export DIRPLOT=$WRFUTILDIR/RUNs/HIST/PLOT/
export DIRSALIDA=/data/salida/$AMBIENTE/intra/
export STARTDATE="${FECHA_INI//\/}"
export OMP_NUM_THREAD=$PLOTTHREADS
PATH=$PATH:$CONDADIR
export BASEPLOTDIR=$WRFUTILDIR/bin/python/ploteos/

source $CONDADIR/activate wrfutil  || $CONDADIR/conda activate wrfutil 
#cd $BASEPLOTDIR/regionalizacion
cd $BASEPLOTDIR/sectores_regiones
IFSOLD=$IFS
IFS=,
while IFS=' ' read name latN latS lonO lonE 
do 
#	python -W ignore acumPronoReg.py $STARTDATE $CICLO 01 $PLAZO $name $latN $latS $lonO $lonE &     
	python -W ignore acumPronoSec.py $STARTDATE $CICLO 01 $PLAZO $name $latN $latS $lonO $lonE &
#done < $BASEPLOTDIR/regionalizacion/zonasNCS.csv
done < $BASEPLOTDIR/sectores_regiones/zonasNCS.csv
wait
IFS=$IFSOLD
EOF

# Parametros de encolamiento
QDEPEND_NAME=$DEPPLOT
QDEPENDTYPE=$TYPEPLOT
QPROC_NAME="PLOTA1_${NOMBRE}_${PASO}"
QPROC=1
QTHREADS=1
QARRAY="0-0"
QWALLTIME="20:00"
QEXCLU=1

# Encolar
queue

###
# Creacion y en colamiento del script
###
read -r -d '' QSCRIPTCMD << "EOF"
export DIRPOST=$WRFUTILDIR/RUNs/HIST/POST/
export DIRWRF=$WRFUTILDIR/RUNs/HIST/FCST/
export DIRPLOT=$WRFUTILDIR/RUNs/HIST/PLOT/
export DIRSALIDA=/data/salida/$AMBIENTE/intra/
export STARTDATE="${FECHA_INI//\/}"
export OMP_NUM_THREAD=$PLOTTHREADS
PATH=$PATH:$CONDADIR
export BASEPLOTDIR=$WRFUTILDIR/bin/python/ploteos/
export CONCURR=10

source $CONDADIR/activate wrfutil  || $CONDADIR/conda activate wrfutil 
cd $BASEPLOTDIR
python -W ignore meteogramasDET.py $STARTDATE $CICLO  
EOF


# Parametros de encolamiento
QDEPEND_NAME=$DEPPLOT
QDEPENDTYPE=$TYPEPLOT
QPROC_NAME="PLOTA2_${NOMBRE}_${PASO}"
QPROC=1
QTHREADS=1
QARRAY="0-0"
QWALLTIME="20:00"
QEXCLU=1

# Encolar
queue

read -r -d '' QSCRIPTCMD << "EOF"
export DIRPOST=$WRFUTILDIR/RUNs/HIST/POST/
export DIRWRF=$WRFUTILDIR/RUNs/HIST/FCST/
export DIRPLOT=$WRFUTILDIR/RUNs/HIST/PLOT/
export DIRSALIDA=/data/salida/$AMBIENTE/intra/
export STARTDATE="${FECHA_INI//\/}"
export OMP_NUM_THREAD=$PLOTTHREADS
PATH=$PATH:$CONDADIR
export BASEPLOTDIR=$WRFUTILDIR/bin/python/ploteos/
export CONCURR=10
#cantReg=$(cat $BASEPLOTDIR/regionalizacion/zonaGraficadoWrf.csv | wc -l)
cantReg=$(cat $BASEPLOTDIR/sectores_regiones/zonaGraficadoWrf.csv | wc -l)
source $CONDADIR/activate wrfutil  || $CONDADIR/conda activate wrfutil 
#cd $BASEPLOTDIR/regionalizacion
cd $BASEPLOTDIR/sectores_regiones
IFSOLD=$IFS
IFS=,
while IFS=' ' read -r name latN latS lonO lonE 
do 
      regacum+=( $STARTDATE $CICLO 01 $PLAZO $name  $latN $latS $lonO $lonE )	
#done < $BASEPLOTDIR/regionalizacion/zonaGraficadoWrf.csv
done < $BASEPLOTDIR/sectores_regiones/zonaGraficadoWrf.csv
IFS=$IFSOLD
parnum=9
echo ${regacum[@]:$(($ARRAYID*5*$parnum)):5*$parnum} | xargs -P5 -n $parnum python -W ignore regnacSup.py


wait
echo "Terminando Ploteos: Hora:$ARRAYID y perfiles:$ppn"


EOF
#cantReg=$(cat $WRFUTILDIR/bin/python/ploteos/regionalizacion/zonaGraficadoWrf.csv | wc -l)
cantReg=$(cat $WRFUTILDIR/bin/python/ploteos/sectores_regiones/zonaGraficadoWrf.csv | wc -l)
# Parametros de encolamiento
QDEPEND_NAME=$DEPPLOT
QDEPENDTYPE=$TYPEPLOT
QPROC_NAME="PLOTRA_${NOMBRE}_${PASO}"
QPROC=1
QTHREADS=1
QARRAY="0-$((($cantReg-1)/5))"
QWALLTIME="20:00"
QEXCLU=1

# Encolar
queue
read -r -d '' QSCRIPTCMD << "EOF"
export DIRPOST=$WRFUTILDIR/RUNs/HIST/POST/
export DIRWRF=$WRFUTILDIR/RUNs/HIST/FCST/
export DIRPLOT=$WRFUTILDIR/RUNs/HIST/PLOT/
export DIRSALIDA=/data/salida/$AMBIENTE/intra/
export STARTDATE="${FECHA_INI//\/}"
export OMP_NUM_THREAD=$PLOTTHREADS
PATH=$PATH:$CONDADIR
export BASEPLOTDIR=$WRFUTILDIR/bin/python/ploteos/
export CONCURR=10

source $CONDADIR/activate wrfutil  || $CONDADIR/conda activate wrfutil 
cd $BASEPLOTDIR/
python -W ignore acumProno.py $STARTDATE $CICLO 01 $PLAZO

EOF

# Parametros de encolamiento
QDEPEND_NAME=$DEPPLOT
QDEPENDTYPE=$TYPEPLOT
QPROC_NAME="PLOTA3_${NOMBRE}_${PASO}"
QPROC=1
QTHREADS=1
QARRAY="0-0"
QWALLTIME="20:00"
QEXCLU=1

# Encolar
queue
