#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

#Load experiment configuration
BASEDIR=$(pwd)/../
source $BASEDIR/conf/config.env
source $BASEDIR/conf/step.conf
source $BASEDIR/conf/exp.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source $BASEDIR/lib/encolar${QUEUESYS}.sh                    

if [ $MACHINE == "FUGAKU" ] ;then
  source $TOPDIR/setup_spack.sh
fi
##### FIN INICIALIZACION ######

cp $NAMELISTDIR/namelist.scale_pp   $SCALEPPDIR/
cp $NAMELISTDIR/namelist.sno_prep   $SCALEPPDIR/
cp $NAMELISTDIR/namelist.scale_init $SCALEPPDIR/
cp $NAMELISTDIR/namelist.ncinput    $SCALEPPDIR/
cp $NAMELISTDIR/namelist.gradsinput $SCALEPPDIR/
FZTEXT=$(cat $NAMELISTDIR/fz_$((E_VERT-1))lev.txt)
if [ $? -ne 0 ]; then
  echo "no vertical level file fz_$((E_VERT-1))lev.txt" ]
  return 1
fi 

#We found the closests dates in the boundary data to the initial and final times of the forecast. 

if [ $STEP -ge 1 ] && [ $WPS_CYCLE -eq 0 ] ; then
   echo "WPS with WPS_CYCLE=0 is run only at STEP=0. And currently STEP=",$STEP        
   return 0
fi

#Edit the namelist from the template

E_WE_LOC=$((E_WE / NPROC_X))
E_SN_LOC=$((E_SN / NPROC_Y))
DATE_INI_SCALE="$(date -ud "$BDY_INI_DATE" +'%Y,%m,%d,%H,%M,%S' )"
if [ $MAP_PROJ == "lambert" ];then
  MAP_PROJ_SCALE="LC"
fi
GTOPO30=$TOPO_LAND_ORG_SCALE/topo/GTOPO30/Products
GLCCv2=$TOPO_LAND_ORG_SCALE/landuse/GLCCv2/Products

sed -i -e "s|__FECHA_INI__|$DATE_INI_SCALE|g" $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__E_WE_LOC__|$E_WE_LOC|g"        $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__E_SN_LOC__|$E_SN_LOC|g"        $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__E_VERT__|$((E_VERT-1))|g"      $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__NPROC_X__|$NPROC_X|g"          $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__NPROC_Y__|$NPROC_Y|g"          $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__DX__|$DX|g"                    $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__DY__|$DY|g"                    $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__MAP_PROJ__|$MAP_PROJ_SCALE|g"  $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__REF_LAT__|$REF_LAT|g"          $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__REF_LON__|$REF_LON|g"          $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__TRUELAT1__|$TRUELAT1|g"        $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__TRUELAT2__|$TRUELAT2|g"        $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__GTOPO30__|\"$GTOPO30\"|g"          $SCALEPPDIR/namelist.scale_pp
sed -i -e "s|__GLCCv2__|\"$GLCCv2\"|g"            $SCALEPPDIR/namelist.scale_pp
sed -i -e "/!---FZ---/a $(echo $FZTEXT | sed -e 's/\./\\./g')" $SCALEPPDIR/namelist.scale_pp


E_WE_LOC=$((E_WE / NPROC_X))
E_SN_LOC=$((E_SN / NPROC_Y))
DATE_INI_SCALE="$(date -ud "$BDY_INI_DATE" +'%Y,%m,%d,%H,%M,%S' )"
if [ $MAP_PROJ == "lambert" ];then
  MAP_PROJ_SCALE="LC"
fi

SCALEDATA="$SCALEPATH/data"
LEN_BDY=$(( $(date -ud "$BDY_END_DATE" +%s) - $(date -ud "$BDY_INI_DATE" +%s) ))
NFILES_BDY=$(( LEN_BDY / BDY_FREQ + 1 ))

if [ ${WPS_DATA_SOURCE} == "WRF" ] ; then
  FTYPE_ORG="NetCDF"
  BNAME_ORG="bdyorg/bdyorg"
else
  FTYPE_ORG="GrADS"
  BNAME_ORG="namelist.gradsinput"
fi

sed -i -e "s|__INI_SCALE__|$DATE_INI_SCALE|g" $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__LEN__|$LEN_BDY|g"              $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__INTERVALO__|$BDY_FREQ|g"       $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__E_VERT__|$((E_VERT-1))|g"      $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__E_WE_LOC__|$E_WE_LOC|g"        $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__E_SN_LOC__|$E_SN_LOC|g"        $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__NPROC_X__|$NPROC_X|g"          $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__NPROC_Y__|$NPROC_Y|g"          $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__DX__|$DX|g"                    $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__DY__|$DY|g"                    $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__MAP_PROJ__|$MAP_PROJ_SCALE|g"  $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__REF_LAT__|$REF_LAT|g"          $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__REF_LON__|$REF_LON|g"          $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__TRUELAT1__|$TRUELAT1|g"        $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__TRUELAT2__|$TRUELAT2|g"        $SCALEPPDIR/namelist.scale_init

sed -i -e "s|__BASENAME_ORG__|$BNAME_ORG|g"   $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__FILETYPE_ORG__|$FTYPE_ORG|g"   $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__SCALEDATA__|$SCALEDATA|g"      $SCALEPPDIR/namelist.scale_init
sed -i -e "s|__NFILES_BDY__|$NFILES_BDY|g"    $SCALEPPDIR/namelist.scale_init

sed -i -e "/!---FZ---/a $(echo $FZTEXT | sed -e 's/\./\\./g')" $SCALEPPDIR/namelist.scale_init

####
## Corremos scale pp (si es que no fue corrido)
###
if [ ! -e $SCALEPPDIR/mktopo/topo/topo.pe000000.nc ] ; then
if [ -e $HISTDIR/const/topo/topo.pe000000.nc ] ; then
   echo "copiando topo and landuse"
   mkdir -p $SCALEPPDIR/mktopo/topo
   mkdir -p $SCALEPPDIR/mktopo/landuse
   cp $HISTDIR/const/topo/*.nc $SCALEPPDIR/mktopo/topo/
   cp $HISTDIR/const/landuse/*.nc $SCALEPPDIR/mktopo/landuse/
else
   echo "Ejecutando scale pp"
   mkdir -p $SCALEPPDIR/mktopo/topo
   mkdir -p $SCALEPPDIR/mktopo/landuse
   mkdir -p $SCALEPPDIR/mktopo/log/scale_pp
   mkdir -p $SCALEPPDIR/mktopo/log/sno
   #ln -sf $WPSDIR/code/*   $WPSDIR/geogrid/
   #cp $SCALEPPDIR/namelist.scale_pp $WPSDIR/geogrid/ 
read -r -d '' QSCRIPTCMD << "EOF"
        cd $SCALEPPDIR/mktopo
        ln -sf $SCALEPPDIR/scale-rm_pp${SCALE_opt} ./scale-rm_pp
        ln -sf $SCALEPPDIR/sno${SCALE_opt} ./sno
        cp $SCALEPPDIR/namelist.scale_pp ./
        cat $SCALEPPDIR/namelist.sno_prep | sed -e "s/__ITEM__/topo/g" > ./namelist.sno_topo
        cat $SCALEPPDIR/namelist.sno_prep | sed -e "s/__ITEM__/landuse/g" > ./namelist.sno_landuse
	ulimit -s unlimited
        export FORT90L=${SCALE_RUNTIME_FLAGS}
        export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SCALE_NETCDF_C}/lib:${SCALE_NETCDF_F}/lib:${SCALE_PNETCDF}/lib:${SCALE_HDF}/lib:$LD_LIBRARY_PATH
        $MPIEXE ./scale-rm_pp namelist.scale_pp
        ERROR=$(( $ERROR + $? ))
        $MPIEXESERIAL ./sno namelist.sno_topo
        ERROR=$(( $ERROR + $? ))
        $MPIEXESERIAL ./sno namelist.sno_landuse
        ERROR=$(( $ERROR + $? ))
EOF
        QPROC_NAME=scale_pp_$STEP
	QPROC=$SCALEPROC
	QNODE=$SCALENODE
	QSKIP=$SCALESKIP
        QOMP=$SCALEOMP
	QMIEM=00
	QWALLTIME="00:10:00"
#        QCONF=${EXPTYPE}.conf
	QWORKPATH=$SCALEPPDIR
	queue 00 00 
	check_proc 00 00

fi

#Copiamos los archivos del directorio 
OUTPUTPATH="$HISTDIR/const/"
mkdir -p $OUTPUTPATH/topo
mkdir -p $OUTPUTPATH/landuse
cp -r $SCALEPPDIR/mktopo/topo/*    $OUTPUTPATH/topo
cp -r $SCALEPPDIR/mktopo/landuse/* $OUTPUTPATH/landuse

fi

###
# Create the script that will run scale init for each ensemble member.
###

echo "Corriendo el scale init para el paso " $STEP

read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
source $BASEDIR/conf/model.conf
source $BASEDIR/conf/exp.conf
source $BASEDIR/conf/step.conf

echo "Processing member $MEM" 

mkdir -p $SCALEPPDIR/$MEM
cd $SCALEPPDIR/$MEM

ln -sf $SCALEPPDIR/scale-rm_init${SCALE_opt}   $SCALEPPDIR/$MEM/scale-rm_init
cp $SCALEPPDIR/namelist.scale_init $SCALEPPDIR/$MEM/
cp $SCALEPPDIR/namelist.ncinput    $SCALEPPDIR/$MEM/
cp $SCALEPPDIR/namelist.gradsinput    $SCALEPPDIR/$MEM/

ln -sf  $HISTDIR/const/topo .
ln -sf  $HISTDIR/const/landuse .
mkdir -p bdyorg
mkdir -p bdy
mkdir -p init
mkdir -p log/scale_init

export FORT90L=${SCALE_RUNTIME_FLAGS}

if [ $WPS_DATA_SOURCE == 'GFS' ] ; then 

   BDYBASE=$BDYPATH/gefs.$(date -d "$INI_BDY_DATE" +"%Y%m%d")/$(date -d "$INI_BDY_DATE" +"%H")/$BDYPREFIX/$MEM/
   echo "Selected data source is GFS in grads format in $BDYBASE"

   #Need to loop over wrfout files.
   CDATE=$BDY_INI_DATE 
   cnt=0
   while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $(date -d "$BDY_END_DATE" +"%Y%m%d%H%M%S") ] ; do 
      for item in atm sfc land ; do
         GFSFILE=$BDYBASE/${item}_$(date -u -d "$CDATE UTC" +"%Y%m%d%H%M%S").grd
         BDYORGPATH=$SCALEPPDIR/$MEM/bdyorg/bdy${item}_$(printf %05g $cnt).grd
         ln -sf $GFSFILE $BDYORGPATH
      done
      #Update CDATE
      CDATE=$(date -u -d "$CDATE UTC + $BDY_FREQ seconds" +"%Y-%m-%d %T")
      cnt=$((cnt+1))
   done

elif [ $WPS_DATA_SOURCE == 'WRF' ]  ; then

   BDYBASE=$BDYPATH/$(date -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")/$MEM/
   echo "Selected data source is WRF, we will use wrfout directly in scale-rm_init"

   #Need to loop over wrfout files.
   CDATE=$BDY_INI_DATE 
   cnt=0
   while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $(date -d "$BDY_END_DATE" +"%Y%m%d%H%M%S") ] ; do 
      WRFFILE=$BDYBASE/wrfout_d01_$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%T" )
      BDYORGPATH=$SCALEPPDIR/$MEM/bdyorg/bdyorg_$(printf %05g $cnt)
      ln -sf $WRFFILE $BDYORGPATH
      #Update CDATE
      CDATE=$(date -u -d "$CDATE UTC + $BDY_FREQ seconds" +"%Y-%m-%d %T")
      cnt=$((cnt+1))
   done
  
else 

  echo "Not recognized WPS_DATA_SOURCE option: "$WPS_DATA_SOURCE
  echo "I can do nothing"
  ERROR=$(( $ERROR + 1 ))

fi

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SCALE_NETCDF_C}/lib:${SCALE_NETCDF_F}/lib:${SCALE_PNETCDF}/lib:${SCALE_HDF}/lib:$LD_LIBRARY_PATH

echo "Running scale init for member $MEM"
$MPIEXE ./scale-rm_init namelist.scale_init
ERROR=$(( $ERROR + $? ))

if [ $ERROR -gt 0 ] ; then
   echo "Error: SCALE init step finished with errors"   
fi

EOF

#Node / core distribution parameters
QNODE=$SCALENODE
QPROC=$SCALEPROC
QSKIP=$SCALESKIP
QOMP=$SCALEOMP
QWALLTIME="00:10:00"
QPROC_NAME=scale_init_${STEP}
#QCONF=${EXPTYPE}.conf
QWORKPATH=$SCALEPPDIR

#Execute the job 
echo "Time taken by scale init and scale"
queue $BDY_MEM_INI $BDY_MEM_END 
time check_proc $BDY_MEM_INI $BDY_MEM_END
if [ $? -ne 0 ] ; then
   echo "Error: Some members do not finish OK"
   echo "Aborting this step"
   exit 1 
fi

#Copy init and bdy files to its final destionation.
for QMEM in $(seq -w $BDY_MEM_INI $BDY_MEM_END) ; do
   DIR_FECHA_INI=$(date -u -d "$BDY_INI_DATE UTC" +"%Y%m%d%H%M%S")
   OUTPUTPATH="$HISTDIR/init/${DIR_FECHA_INI}/$QMEM/"
   mkdir -p $OUTPUTPATH
   mv $SCALEPPDIR/$QMEM/init/*.nc $OUTPUTPATH
   OUTPUTPATH="$HISTDIR/bdy/${DIR_FECHA_INI}/$QMEM/"
   mkdir -p $OUTPUTPATH
   mv $SCALEPPDIR/$QMEM/bdy/*.nc $OUTPUTPATH
done


echo "Termine de correr el SCALE pp"

