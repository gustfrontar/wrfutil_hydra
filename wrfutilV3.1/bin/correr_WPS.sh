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

####
## Configuracion del namelist.wps
## 1era Parte
###

#Descomprimimos el archivo .tar (si es que no fue descomprimido)
if [ ! -e $WPSDIR/code/geogrid.exe ] ; then
   echo "Descomprimiendo WPS"
   mkdir -p $WPSDIR/code
   mv $WPSDIR/wps.tar $WPSDIR/code
   cd $WPSDIR/code 
   tar -xf wps.tar
   #Si existe el namelist.wps lo borramos para que no interfiera
   #con los que crea el sistema de asimilacion
   if [ -e namelist.wps ] ; then
      rm -f namelist.wps 
   fi
fi

#Editamos el namelist.wps [Esto es para un experimento asi que asumimos que ya esta toda la info de los bordes disponibles en formato grib disponible]
mkdir -p $WPSDIR
cp $NAMELISTDIR/namelist.wps $WPSDIR/
#Buscamos la fecha inmediatamente inferior al inicio del experimento en en base a la frecuencia de los archivos del ensamble de condiciones de borde.  
FECHAINIWPS=$(date_floor "$FECHA_INI" $INTERVALO_BDY )
FECHAINIWPS=$(date -u -d "$FECHAINIWPS UTC" +"%Y-%m-%d_%T" ) #Cambio al formato WPS
FECHAFINWPS=$(date_ceil  "$FECHA_FIN" $INTERVALO_BDY )
FECHAFINWPS=$(date -u -d "$FECHAFINWPS UTC" +"%Y-%m-%d_%T" ) #Cambio al formato WPS

sed -i -e "s|__FECHA_INI__|$FECHAINIWPS|g" $WPSDIR/namelist.wps
sed -i -e "s|__FECHA_FIN__|$FECHAFINWPS|g" $WPSDIR/namelist.wps
sed -i -e "s|__INTERVALO__|$INTERVALO_WPS|g"  $WPSDIR/namelist.wps
sed -i -e "s|__E_WE__|$E_WE|g" $WPSDIR/namelist.wps
sed -i -e "s|__E_SN__|$E_SN|g" $WPSDIR/namelist.wps
sed -i -e "s|__DX__|$DX|g" $WPSDIR/namelist.wps
sed -i -e "s|__DY__|$DY|g" $WPSDIR/namelist.wps
sed -i -e "s|__MAP_PROJ__|$MAP_PROJ|g" $WPSDIR/namelist.wps
sed -i -e "s|__REF_LAT__|$REF_LAT|g" $WPSDIR/namelist.wps
sed -i -e "s|__REF_LON__|$REF_LON|g" $WPSDIR/namelist.wps
sed -i -e "s|__TRUELAT1__|$TRUELAT1|g" $WPSDIR/namelist.wps
sed -i -e "s|__TRUELAT2__|$TRUELAT2|g" $WPSDIR/namelist.wps
sed -i -e "s|__STAND_LON__|$STAND_LON|g" $WPSDIR/namelist.wps
sed -i -e "s|__GEOG__|$GEOG|g" $WPSDIR/namelist.wps
sed -i -e "s|__IOTYPE__|$IOTYPE|g" $WPSDIR/namelist.wps


####
## Corremos geogrid (si es que no fue corrido)
###
cd $WPSDIR
if [ ! -e $WPSDIR/geogrid/geo_em.d01.nc ] ; then
   echo "Ejecutando geogrid"
   mkdir -p $WPSDIR/geogrid
   ln -sf $WPSDIR/code/*   $WPSDIR/geogrid/
   cp $WPSDIR/namelist.wps $WPSDIR/geogrid/ 
read -r -d '' QSCRIPTCMD << "EOF"
        cd $WPSDIR/geogrid
	ulimit -s unlimited
        $MPIEXE ./geogrid.exe
EOF
        QPROC_NAME=GEOG_$PASO
	QPROC=$WPSPROC
	TPROC=$WPSTPROC
	QTHREADS=$WPSTHREADS
	QMIEM=01
	QWALLTIME=$WPSWALLTIME
	queue 01 01 
	check_proc 01 01
fi

###
# Creacion y en colamiento del script
###
echo "Corriendo el WPS para el paso " $PASO
cd $WPSDIR
mkdir -p $HISTDIR/WPS

read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
echo "Procesando Miembro $MIEM"
cd $WPSDIR
rm -r $WPSDIR/$MIEM
mkdir -p $WPSDIR/$MIEM
ln -sf $WPSDIR/code/* $WPSDIR/$MIEM
cp $WPSDIR/namelist.wps $WPSDIR/$MIEM/
cd $WPSDIR/$MIEM
ln -sf $WPSDIR/geogrid/geo_em* .
FECHA_INI_BDY=$(date_floor "$FECHA_INI" $INTERVALO_INI_BDY )
BDYBASE=$BDYDIR/gefs.$(date -d "$FECHA_INI_BDY" +"%Y%m%d")/$(date -d "$FECHA_INI_BDY" +"%H")/$BDYPREFIX/$MIEM/
echo "Estoy buscando los archivos del BDY en la carpeta $BDYBASE"
## Si el miembro es 00 entonces es deterministico, sino Ensamble
ln -sf $WPSDIR/$BDYVTABLE ./Vtable
echo "Generando met_em a partir de BDYs en: $WPSDIR/$MIEM"
#Linkeo los gribs
cd $WPSDIR/$MIEM
./link_grib.csh $BDYBASE/$BDYSEARCHSTRING  

[[ $? -ne 0 ]] && dispararError 9 "link_grib.csh "

#Corro el Ungrib 
OMP_NUM_THREADS=1
./ungrib.exe > ungrib.log
[[ $? -ne 0 ]] && dispararError 9 "ungrib.exe"

#Corro el Metgrid en paralelo
$MPIEXE ./metgrid.exe
EC=$?
[[ $EC -ne 0 ]] && dispararError 9 "metgrid.exe"

EOF

# Parametros de encolamiento
## Calculo cuantos miembros hay que correr
QNODE=$WPSNODE
QPROC=$WPSPROC
TPROC=$WPSTPROC
QTHREADS=$WPSTHREADS
QWALLTIME=$WPSWALLTIME
QPROC_NAME=WPS_${PASO}

# Encolar
queue $BDY_MIEMBRO_INI $BDY_MIEMBRO_FIN
check_proc $BDY_MIEMBRO_INI $BDY_MIEMBRO_FIN

#Copiamos los archivos del directorio 
for QMIEM in $(seq -w $BDY_MIEMBRO_INI $BDY_MIEMBRO_FIN) ; do
   OUTPUTPATH="$HISTDIR/WPS/met_em_ori/$QMIEM/"
   mkdir -p $OUTPUTPATH
   mv $WPSDIR/$QMIEM/met_em* $OUTPUTPATH
done



echo "Termine de correr el WPS"

