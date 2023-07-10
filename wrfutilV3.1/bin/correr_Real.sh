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
read -r IY IM ID IH Im  <<< $(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN)) minutes" +"%Y %m %d %H %M")

#Elijo  de que esperimento voy a sacar los met_ems*
[[ -z $METIN ]] && METIN=$NOMBRE

#Elijo la clase de mets que voy a usar, (originales o perturbados)
METCLASS="originales"                                    
[[ $PERTURBAR -eq 1 ]] &&  METCLASS="perturbados"

# Calculos el ciclo de inicializacio de METS que necesito en funcion de mi pronostico
## Ej: si FRECDIARIA=4, Para los IH del CICLO (experimento.conf) [0-5] da CICLOMET=0; Para los IH del CICLO (experimento.conf) [6-11] da CICLOMET=6; etc
## por default usa este METDIR:
PERDIARIO=$((24/$FRECDIARIA))
CICLOMET=$(( (10#$IH/$PERDIARIO)*$PERDIARIO ))
necesario=${IY}${IM}${ID}_$(printf "%02g" $CICLOMET)0000
METDIR=$BASEDIR/HIST/WPS/$necesario/$METIN/$METCLASS

## Esto es por si quiero elejir los ultimos disponibles (pisa el METDIR):
ultimo=""
[[ $METLAST -ne 0 ]] && ultimo=$(cat $BASEDIR/HIST/WPS/last/${METIN}_${METCLASS}.txt ) 
[[ ! -z $ultimo ]] && METDIR=$BASEDIR/HIST/WPS/$ultimo/$METIN/$METCLASS/
[[ $METLAST -ne 0 ]] && [[  -z  $ultimo ]] && dispararError 4 "$BASEDIR/HIST/WPS/last/${METIN}_${METCLASS}.txt"

## Esto es para cuando recuperamos y el last es mucho mas nuevo de lo que necesitamos (pisa el METDIR): 
[[ $METLAST -ne 0 ]] && [[ $necesario  < $ultimo ]] && METDIR=$BASEDIR/HIST/WPS/$necesario/$METIN/$METCLASS

## Si es asimilacion y voy a hacer pronostico, hay que Usar los mets que se usron para hacer el analisis (pisa el METDIR):  
[[ $ASIMILACION -ne 0 ]] && [[ $PASO -gt 0 ]] && [[ -f $BASEDIR/HIST/ANA/${IY}-${IM}-${ID}_${IH}_${Im}_00/$NOMBRE/mets.dir ]] && METDIR=$(cat $BASEDIR/HIST/ANA/${IY}-${IM}-${ID}_${IH}_${Im}_00/$NOMBRE/mets.dir)


# Dejamos registro que que mets se van a utilizar
echo "usando los METGRID  de $METDIR"
echo $METDIR >$BASEDIR/WRF/mets.dir

#script de ejecucion
#
read -r -d '' QSCRIPTCMD << "EOF"
MIEM=$(printf "%02g" $((10#$MIEMBRO_FIN-10#$MIEMBRO_INI)))
test $ARRAYID && MIEM=$(printf "%02g" $ARRAYID)
ulimit -s unlimited
echo "Procesando Miembro $MIEM" 
read -r IY IM ID IH Im  <<< $(date -u -d "$FECHA_INI UTC +$((10#$CICLO)) hours+$((10#$CICLO_MIN)) minutes" +"%Y %m %d %H %M")
read -r FY FM FD FH Fm  <<< $(date -u -d "$FECHA_INI UTC +$((10#$CICLO+$PLAZO)) hours+$((10#$CICLO_MIN+10#$PLAZO_MIN)) minutes" +"%Y %m %d %H %M")
cd $WRFDIR
CNL=`ls $WRFUTILDIR/namelists/namelist.input.${NAMELIST_TIPO}* |wc -w`
NAMELISTFILE="$WRFDIR/$MIEM/namelist.input"
[[ -d $WRFDIR/$MIEM ]] && rm -r $WRFDIR/$MIEM
mkdir -p "$WRFDIR/$MIEM"
tar  xf WRF.tar -C "$WRFDIR/$MIEM"
NL=$(printf '%03d' $((10#$MIEM % 10#$CNL)))
TIPO=det
test  $NAMELIST_TIPO != "det"  && TIPO=${NAMELIST_TIPO}.$NL
cp $WRFUTILDIR/namelists/namelist.input.${TIPO} $NAMELISTFILE
echo "cp $BASEDIR/namelists/namelist.input.${TIPO} $NAMELISTFILE"
IOTYPE=2
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
sed -i -e "s|__INTERVALO_WPS__|$(($INTERVALO_WPS*60))|g" $NAMELISTFILE
sed -i -e "s|__INTERVALO_WRF__|$INTERVALO_WRF|g" $NAMELISTFILE
sed -i -e "s|__E_WE__|$E_WE|g" $NAMELISTFILE
sed -i -e "s|__E_SN__|$E_SN|g" $NAMELISTFILE
sed -i -e "s|__DX__|$DX|g" $NAMELISTFILE
sed -i -e "s|__DY__|$DY|g" $NAMELISTFILE
sed -i -e "s|__IOTYPE__|$IOTYPE|g" $NAMELISTFILE
sed -i -e "s|__METLEV__|$METLEV|g" $NAMELISTFILE
sed -i -e "s|__NUMTILE__|$NUMTILE|g" $NAMELISTFILE
sed -i -e "s|__NIOT__|$NIOT|g" $NAMELISTFILE
sed -i -e "s|__NIOG__|$NIOG|g" $NAMELISTFILE


cd $BASEDIR/WRF/$MIEM
[ -f envvars.sh ] && source envvars.sh
[ -f entorno.sh ] && source entorno.sh

# Esto los hacemos para calcular que met tenemos que usar de los dispoibles
# rotamos estilo round robin
METDIR=$(cat $BASEDIR/WRF/mets.dir)
cantMiem=$(ls $METDIR |grep -v 00 |wc -w)    # Cantidad de miebros en los que se preproceso mets
#minMiem=( $(ls $METDIR |sort -n) ) # El minimio miembro . . . puede ser distnto de 0
#minMiem=${minMiem[0]} # ahora si el minimo
#METMIEM=$(printf '%02d' $(( (10#$MIEM % $cantMiem) + 10#$minMiem - 10#$MIEMBRO_INI)) ) # elijo el met que me corrsponde segun el miembro que soy
if [[ $((10#$MIEM)) -le $(($cantMiem)) ]]
then 
	METMIEM=$MIEM
else 
	METMIEM=$(printf '%02d' $((10#$MIEMBRO_FIN - 10#$MIEM +1)))
fi
echo "Soy el miembro $MIEM y uso los metgrid del miembro $METMIEM"
ln -sf $METDIR/$METMIEM/met_* .
OMP_NUM_THREADS=$REALTHREADS
OMP_STACKSIZE=512M
$MPIEXE  $BASEDIR/WRF/$MIEM/real.exe
EXCOD=$?
#[[ $EXCOD -ne 0 ]] && echo "real.exe" > $BASEDIR/FATAL.ERROR 
[[ $EXCOD -ne 0 ]] && dispararError 9 "real.exe"
EOF


# Parametros de encolamiento
QDEPEND_NAME=$DEPREAL
QDEPENDTYPE=$TYPEREAL
QPROC_NAME=REAL_$1_$PASO
QPROC=$REALPROC
QTHREADS=$REALTHREADS
QARRAY="$((10#$MIEMBRO_INI))-$((10#$MIEMBRO_FIN))"
QWALLTIME=$REALWALLTIME
QEXCLU=1


# Encolar
queue
