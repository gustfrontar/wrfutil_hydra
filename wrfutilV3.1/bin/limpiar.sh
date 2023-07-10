#!/bin/bash 

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${WRFUTILDIR:?"$USO"}


### PARAMETROS

TESTFILE=$1

### CONFIGURACION

[ ! -f $WRFUTILDIR/lib/errores.env ] && exit 1
source $WRFUTILDIR/lib/errores.env
[ ! -f $WRFUTILDIR/RUNs/HIST/limpiar.conf ] && dispararError 4 "limpiar.conf"
source $WRFUTILDIR/RUNs/HIST/limpiar.conf

##### FIN INICIALIZACION ######

exec 2> $LOGFILE
exec 1> $LOGFILE

function limpiar(){
if [[ $TIMEPOGUARDO -gt 0 ]]
then
	for iter in $(seq 1 $VENTANA) 
	do 
		finddir=$WRFUTILDIR/RUNs/HIST/$LIMPBASEDIR
		fecha=$( date -d "today  -$iter days" +"%Y%m%d_")
		for fechadir in $(ls -d $finddir/${fecha}* 2> /dev/null)
		do
			find $fechadir/$NOMBRE -name "$MATCHFILE" -mtime +$TIMEPOGUARDO -exec bash -c $EXOP "file={}; $ECHO $PROCESAR "  \;
			#find $fechadir/$NOMBRE -type d -empty -delete
		done
		find $finddir -type d -empty -delete
	done
fi
}
EXOP=""
[[ $DEBUG == 1 ]] && EXOP=" -x "
ECHO=""
[[ $TESTRUN == 1 ]] && ECHO="echo"

SECONDS=0
if [[ -f $TESTFILE ]]
then 
	source $TESTFILE
else 
	for exp in ${ENTORNOS[@]}  
	do 
		echo "experimento -> $exp"  
		export NOMBRE=$exp
		source $WRFUTILDIR/RUNs/$exp/experimento.limp
	done
fi

duration=$SECONDS
echo "TIEMPO de EJECUCION: $(($duration / 60)) minutos y $(($duration % 60)) segundos."
