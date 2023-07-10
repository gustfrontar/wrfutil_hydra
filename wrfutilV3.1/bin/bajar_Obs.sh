#!/bin/bash 
#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

###  Leyenda de USO

read -r -d '' USO << EOF

        Ud. deberia usar este escript de la siguiente manera:
                $0 <nombre entonro > 
        Recuerde exportar la variable de entorno WRFUTILDIR=<directorio de instalacion del wrfutil>
EOF
: ${1:?"$USO"}  ${WRFUTILDIR:?"$USO"}


### PARAMETROS

NOMBRE=$1

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
[ ! -f "$BASEDIR/experimento.asim" ] && dispararError 4 "$BASEDIR/$EXPCONF"
source $BASEDIR/experimento.asim




mkdir -p $DIROBSREPO 

echo " DEscargando Datos en : $DIROBSREPO"
cd $DIROBSREPO



# Fecha del analisis 
read -r anio mes dia hora min  <<< $(date -u -d "$FECHA_OBS +$CICLO_OBS hours +$CICLO_OBS_MIN minutes" +"%Y %m %d %H %M")  

# fecha fin de la ventana de observaciones
read -r Fanio Fmes Fdia Fhora Fmin  <<< $(date -u -d "$FECHA_OBS +$CICLO_OBS hours +$((10#$CICLO_OBS_MIN-10#$DESPLAZAMIENTO+10#$OBSWIN)) minutes" +"%Y %m %d %H %M")  

#fecha de inicio de la ventana de observaciones
read -r Ianio Imes Idia Ihora Imin  <<< $(date -u -d "$FECHA_OBS +$CICLO_OBS hours +$((10#$CICLO_OBS_MIN-10#$DESPLAZAMIENTO)) minutes" +"%Y %m %d %H %M")


fin=$(date -u -d "${Fanio}/${Fmes}/${Fdia} ${Fhora}:${Fmin} "  +"%Y-%m-%d %H:%M %z")
ini=$(date -u -d "${Ianio}/${Imes}/${Idia} ${Ihora}:${Imin} "  +"%Y-%m-%d %H:%M %z")

inc=$OBSFREC     ## Este incremento es en minutos                               


#######
# Datos de Oracle
###
# BOYAS , SYNOP , TEMP , SHIP
#######
#cp /data/OBS/desarrollo/asm_*_??.lst $OBSREPO/
cp /data/OBS/desarrollo/asm_*_${Ianio}${Imes}${Idia}${Ihora}*.lst $DIROBSREPO
cp /data/OBS/desarrollo/asm_*_${Fanio}${Fmes}${Fdia}${Fhora}*.lst $DIROBSREPO

# compiando los MOtionVectors
cp /data/OBS/AMV/* $DIROBSREPO/



#######
## Datos de Carlos 
#######
rm amdar_asim.tar
wget -v -t0 http://amdar.smn.gob.ar/amdar/asim/amdar_asim.tar
tar xvf amdar_asim.tar

#######
## Datos de Aerolineas Argentinas
#######

mkdir -p $DIROBSREPO/AA_AMDAR/
lftp  -u aerolineas2,'#ftp_aaa!14dnpt.' ftp.smn.gob.ar -e "set net:max-retries 1; mget -O $DIROBSREPO/AA_AMDAR/ entrada/*_${Ianio}${Imes}${Idia}${Ihora}*txt ; exit"
lftp  -u aerolineas2,'#ftp_aaa!14dnpt.' ftp.smn.gob.ar -e "set net:max-retries 1; mget -O $DIROBSREPO/AA_AMDAR/ entrada/*_${Fanio}${Fmes}${Fdia}${Fhora}*txt ; exit"



#######
## Estaciones Automaticas 
#######
obs_ini=$(date -u -d "$ini - $(($inc/2)) minutes" +"%Y-%m-%dT%H:%M:%S")
obs_fin=$(date -u -d "$fin + $(($inc/2)) minutes" +"%Y-%m-%dT%H:%M:%S")

PATH=$CONDADIR:$PATH
source $CONDADIR/activate wrfutil || $CONDADIR/conda activate wrfutil

mkdir -p $DIROBSREPO/EMAs/
python ${WRFUTILDIR}/bin/python/OBS/baja_emas.py $obs_ini $obs_fin

source $CONDADIR/deactivate || $CONDADIR/conda deactivate

cd $DIROBSREPO
export Fjuliano=$(date -d "${Fanio}/${Fmes}/${Fdia}" +"%j")
export Ijuliano=$(date -d "${Ianio}/${Imes}/${Idia}" +"%j")


#####
## Obtenemos las observaciones del AIRSs
#####
descargarAIRs(){
	oldpwd=$(pwd)
	mkdir -p $DIROBSREPO/AIRS
	cd $DIROBSREPO/AIRS
	rm  index.html
	wget --http-user=smn_mdillon --http-password=Dorrego19  "https://discnrt1.gesdisc.eosdis.nasa.gov/data/Aqua_AIRS_NRT/AIRS2RET_NRT.006/$Panio/$Pjuliano/"
	perl -lane 'm/(AIRS\..*?\.hdf)\"*/ &&print $1' index.html |sort -u > disponibles.airs
	ls AIRS.${Panio}.${Pmes}.${Pdia}.* | sort -u > descargados.airs
	comm -3 disponibles.airs descargados.airs > descargar.airs

	for file in $(cat descargar.airs)
	do
       		wget --http-user=smn_mdillon --http-password=Dorrego19 "https://discnrt1.gesdisc.eosdis.nasa.gov/data/Aqua_AIRS_NRT/AIRS2RET_NRT.006/$Panio/$Pjuliano/${file}"
	done
	cd $oldpwd
}

export Pjuliano=$Ijuliano
export Panio=$Ianio
export Pmes=$Imes
export Pdia=$Idia	
descargarAIRs
if [[ "$((10#$Idia))" -ne "$((10#$Fdia))" ]] 
then 
	export Pjuliano=$Fjuliano
	export Panio=$Fanio
	export Pmes=$Fmes
	export Pdia=$Fdia	
	descargarAIRs
fi
