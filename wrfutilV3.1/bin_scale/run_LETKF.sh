#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

### CONFIGURACION
#Load experiment configuration
BASEDIR=$(pwd)/../
source $BASEDIR/conf/config.env
source $BASEDIR/conf/step.conf
source $BASEDIR/conf/exp.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source $BASEDIR/conf/letkf.conf
source $BASEDIR/conf/obs.conf
source ${BASEDIR}/lib/encolar${QUEUESYS}.sh                     # Selecciona el metodo de encolado segun el systema QUEUESYS elegido
if [ $MACHINE == "FUGAKU" ] ;then
  source $TOPDIR/setup_spack.sh
fi

if [ ! -z ${PJM_SHAREDTMP} ] ; then
   if [ ${USETMPDIR} -eq 1 ] ; then
      echo "Using Fugaku's shared dir to run WRF"
      LETKFDIR=${PJM_SHAREDTMP}/LETKF
   fi
fi

LETKFDIRRUN=$LETKFDIR/00/   #We need to add 00 in order to be consistent with the paralelization lib.
rm -fr $LETKFDIRRUN
mkdir -p $LETKFDIRRUN

cp $LETKFDIR/letkf.exe $LETKFDIRRUN/

mkdir -p $LETKFDIRRUN/obsdep
mkdir -p $LETKFDIRRUN/log/scale
mkdir -p $LETKFDIRRUN/log/letkf
mkdir -p $LETKFDIRRUN/topo
mkdir -p $LETKFDIRRUN/landuse
cp $HISTDIR/const/topo/*.nc $LETKFDIRRUN/topo/
cp $HISTDIR/const/landuse/*.nc $LETKFDIRRUN/landuse/
#if [ ! -e $LETKFDIR/code/letkf.exe ] ; then
#   echo "Descomprimiendo LETKF"
#   mkdir -p $LETKFDIR/code/
#   tar -xf $LETKFDIR/letkf.tar -C $LETKFDIR/code/
#   #Si existe el namelist.input lo borramos para que no interfiera
#   #con los que crea el sistema de asimilacion
#   if [ -e $LETKFDIR/code/letkf.namelist ] ; then
#      rm -f $LETKFDIR/code/letkf.namelist
#   fi
#fi

#ln -sf $LETKFDIR/code/* $LETKFDIRRUN

echo "Editing namelist"
NAMELISTFILE=$LETKFDIRRUN/namelist.scale-letkf
cp $NAMELISTDIR/namelist.scale-letkf $NAMELISTFILE

NSLOTS=$(( ($ANALYSIS_WIN_END-$ANALYSIS_WIN_INI)/$ANALYSIS_WIN_STEP+1))
NBSLOT=$(( ($ANALYSIS_FREQ-$ANALYSIS_WIN_INI)/$ANALYSIS_WIN_STEP+1))
NBV=$(( 10#$MEM_END - 10#$MEM_INI + 1 ))

#Check if we are going to assimilate radar data and how many radar
#are being assimilated.
#if [[ ${OBS_LIST[@]} =~ "RADARC" ]] ; then 
#   N_RADAR=${#RADARC_LIST[@]}  #Get the number of radars. 
#else 
#   N_RADAR=0 
#fi

OBS_IN_NUM=${#OBS_LIST[@]}
OBS_IN_NAME=""
OBS_IN_FORMAT=""
OBSDA_RUN=""

echo "Linking model files"
#FECHA_WINDOW_INI=$(date -u -d "$FECHA_INI UTC +$((($ANALISIS_FREC*($PASO-1))+$SPIN_UP_LENGTH+$ANALISIS_WIN_INI)) seconds" +"%Y-%m-%d %T")
#ANALYSIS_DATE_WFMT=$(date -u -d "$FECHA_INI UTC +$((($ANALISIS_FREC*($PASO-1))+$SPIN_UP_LENGTH+$ANALISIS_FREC)) seconds" +"%Y-%m-%d_%T")   #WRF format
ANALYSIS_DATE_SFMT=$(date -u -d "$DA_INI_DATE UTC +$(( (ANALYSIS_FREQ*(STEP-1))+SPIN_UP_LENGTH+ANALYSIS_FREQ )) seconds" +"%Y%m%d-%H%M%S.000")   #WRF format
ANALYSIS_DATE_PFMT=$(date -u -d "$DA_INI_DATE UTC +$(( (ANALYSIS_FREQ*(STEP-1))+SPIN_UP_LENGTH+ANALYSIS_FREQ )) seconds" +"%Y%m%d%H%M%S")  #Path/folder format

#for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN ) ; do
#   for ISLOT in $(seq -f "%02g" 0 $(($NSLOTS-1))) ; do
#      FECHA_SLOT=$(date -u -d "$FECHA_WINDOW_INI UTC +$(($ISLOT*$ANALISIS_WIN_STEP)) seconds" +"%Y-%m-%d_%T")
#      ln -sf ${WRFDIR}/${MIEM}/wrfout_d01_$FECHA_SLOT ${LETKFDIRRUN}/gs$(printf %02d $((10#$ISLOT+1)))$(printf %05d $((10#$MIEM)))
#   done	   
#done

for MEM in $(seq -w $MEM_INI $MEM_END) ; do
MEMS=$(printf %04g $((10#$MEM)) ) 
mkdir -p $LETKFDIRRUN/gues/$MEMS
mkdir -p $LETKFDIRRUN/anal/$MEMS
#cp ${SCALEDIR}/${MIEM}/init/init_$ANALYSIS_DATE_SFMT* ${LETKFDIRRUN}/gues/$MIEM/
#ln -sf ${SCALEDIR}/${MIEM}/init/init_$ANALYSIS_DATE_SFMT* ${LETKFDIRRUN}/anal/$MIEM/
cp ${HISTDIR}/GUES/$ANALYSIS_DATE_PFMT/$MEM/init_$ANALYSIS_DATE_SFMT* ${LETKFDIRRUN}/gues/$MEMS/
cp ${HISTDIR}/GUES/$ANALYSIS_DATE_PFMT/$MEM/init_$ANALYSIS_DATE_SFMT* ${LETKFDIRRUN}/anal/$MEMS/

mkdir -p $LETKFDIRRUN/hist/$MEMS
ln -s ${SCALEDIR}/$MEM/fcst/history* ${LETKFDIRRUN}/hist/$MEMS/
done

mkdir -p $LETKFDIRRUN/gues/mean
mkdir -p $LETKFDIRRUN/anal/mean
cp ${HISTDIR}/GUES/$ANALYSIS_DATE_PFMT/$MEM_END/init_$ANALYSIS_DATE_SFMT* ${LETKFDIRRUN}/gues/mean/
cp ${HISTDIR}/GUES/$ANALYSIS_DATE_PFMT/$MEM_END/init_$ANALYSIS_DATE_SFMT* ${LETKFDIRRUN}/anal/mean/
mkdir -p $LETKFDIRRUN/gues/sprd
mkdir -p $LETKFDIRRUN/anal/sprd
cp ${HISTDIR}/GUES/$ANALYSIS_DATE_PFMT/$MEM_END/init_$ANALYSIS_DATE_SFMT* ${LETKFDIRRUN}/gues/sprd/
cp ${HISTDIR}/GUES/$ANALYSIS_DATE_PFMT/$MEM_END/init_$ANALYSIS_DATE_SFMT* ${LETKFDIRRUN}/anal/sprd/

#CTIME_OFMT=$(date -u -d "$CTIME UTC" +"%Y%m%d%H%M%S")
CTIME_OFMT=$ANALYSIS_DATE_PFMT
for j in $(seq $OBS_IN_NUM) ;do
  i=$((j-1))
  OBS_TYPE=${OBS_LIST[$i]}
  if [  $OBS_TYPE == "RADARC" ] ;then
    OBS_IN_FORMAT=${OBS_IN_FORMAT}"'RADAR', "
    OBSPATHT="'$OBSPATH/${OBS_TYPE}/WRF_${OBS_TYPE}_${RADARC_LIST[$i]}_${CTIME_OFMT}.dat', "
  elif [ $OBS_TYPE == "ADPAUT" ] || [ $OBS_TYPE == "ADPSFC" ] ;then
    OBS_IN_FORMAT=${OBS_IN_FORMAT}"'AWS', "
    OBSPATHT="'${OBSPATH}/${OBS_TYPE}/WRF_${OBS_TYPE}_${CTIME_OFMT}.dat', "
  else
    OBS_IN_FORMAT=${OBS_IN_FORMAT}"'PREPBUFR', "
    OBSPATHT="'${OBSPATH}/${OBS_TYPE}/WRF_${OBS_TYPE}_${CTIME_OFMT}.dat', "
  fi
  OBS_IN_NAME=${OBS_IN_NAME}${OBSPATHT}
  OBSDA_RUN=${OBSDA_RUN}".true., "
done

FNAME_RST="init_$ANALYSIS_DATE_SFMT"
NPROC=$((SCALEPROC/SCALESKIP))
PPN=$((ICORE/SCALESKIP))

sed -i -e "s|__COV_INFL_MUL__|$COV_INFL_MUL|g"             $NAMELISTFILE
#sed -i -e "s|__SP_INFL_ADD__|$SP_INFL_ADD|g"               $NAMELISTFILE
sed -i -e "s|__RELAX_ALPHA_SPREAD__|$RELAX_ALPHA_SPREAD|g" $NAMELISTFILE
sed -i -e "s|__RELAX_ALPHA__|$RELAX_ALPHA|g"               $NAMELISTFILE
sed -i -e "s|__NSLOTS__|$NSLOTS|g"                         $NAMELISTFILE
sed -i -e "s|__NBSLOT__|$NBSLOT|g"                         $NAMELISTFILE
sed -i -e "s|__NBV__|$NBV|g"                               $NAMELISTFILE
sed -i -e "s|__SIGMA_OBS__|$SIGMA_OBS|g"                   $NAMELISTFILE
sed -i -e "s|__SIGMA_OBSV__|$SIGMA_OBSV|g"                 $NAMELISTFILE
sed -i -e "s|__SIGMA_OBS_RADAR__|$SIGMA_OBS_RADAR|g"       $NAMELISTFILE
sed -i -e "s|__SIGMA_OBSZ__|$SIGMA_OBSZ|g"                 $NAMELISTFILE
#sed -i -e "s|__SIGMA_OBST__|$SIGMA_OBST|g"                 $NAMELISTFILE
#sed -i -e "s|__N_RADAR__|$N_RADAR|g"                       $NAMELISTFILE
sed -i -e "s|__LEV_UPDATE_Q__|$LEV_UPDATE_Q|g"             $NAMELISTFILE
#sed -i -e "s|__THRESHOLD_DZ__|$THRESHOLD_DZ|g"             $NAMELISTFILE
sed -i -e "s|__GROSS_ERROR__|$GROSS_ERROR|g"               $NAMELISTFILE
sed -i -e "s|__GROSS_ERROR_REFLECTIVITY__|$GROSS_ERROR_REFLECTIVITY|g"  $NAMELISTFILE

sed -i -e "s|__PPN__|$PPN|g"                               $NAMELISTFILE
sed -i -e "s|__MEM_NODES__|$SCALENODE|g"                   $NAMELISTFILE
sed -i -e "s|__NPROC__|$NPROC|g"                           $NAMELISTFILE
sed -i -e "s|__WIN_STEP__|$ANALYSIS_WIN_STEP|g"            $NAMELISTFILE
sed -i -e "s|__OBS_IN_NUM__|$OBS_IN_NUM|g"                 $NAMELISTFILE
sed -i -e "s|__OBS_IN_NAME__|$OBS_IN_NAME|g"               $NAMELISTFILE
sed -i -e "s|__OBS_IN_FORMAT__|$OBS_IN_FORMAT|g"           $NAMELISTFILE
sed -i -e "s|__OBSDA_RUN__|$OBSDA_RUN|g"                   $NAMELISTFILE
sed -i -e "s|__FNAME_RST__|$FNAME_RST|g"                   $NAMELISTFILE

cat $SCALEDIR/namelist.scale >> $NAMELISTFILE 

#script de ejecucion
read -r -d '' QSCRIPTCMD << "EOF"
  cd $LETKFDIR/00/

  export FORT90L=${SCALE_RUNTIME_FLAGS}

  export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SCALE_NETCDF_C}/lib:${SCALE_NETCDF_F}/lib:${SCALE_PNETCDF}/lib:${SCALE_HDF}/lib:$LD_LIBRARY_PATH

  $MPIEXE ./letkf.exe namelist.scale-letkf ./log/letkf/NOUT
  ERROR=$(( $ERROR + $? ))
EOF

# Parametros de encolamiento
QPROC_NAME=LETKF_$STEP
QNODE=$LETKFNODE
QPROC=$LETKFPROC
QWALLTIME=$LETKFWALLTIME
QSKIP=$LETKFSKIP
QOMP=$LETKFOMP
QWALLTIME=$LETKFWALLTIME
QWORKPATH=$LETKFDIR

echo "Sending letkf script to the queue"
cd $LETKFDIR
queue 00 00
check_proc 00 00
if [ $? -ne 0 ] ; then
   echo "Error: Some members do not finish OK"
   echo "Aborting this step"
   exit 1
fi

##
# Save the analysis
##
echo "Copying analysis files"
DIRGUES=$HISTDIR/GUES/$ANALYSIS_DATE_PFMT
DIRANAL=$HISTDIR/ANAL/$ANALYSIS_DATE_PFMT
for MEM in $(seq -w $MEM_INI $MEM_END ) ; do
        MEMS=$(printf %04g $((10#$MEM)) )
        mkdir -p  ${DIRANAL}/$MEM
#	echo "Updating the date in the analysis file " $WRFDIR/$MIEM/wrfout_d01_$ANALYSIS_DATE_WFMT $ANALYSIS_DATE_WFMT
#	$LETKFDIR/code/update_wrf_time.exe $WRFDIR/$MIEM/wrfout_d01_$ANALYSIS_DATE_WFMT $ANALYSIS_DATE_WFMT
#        mv  $WRFDIR/$MIEM/wrfout_d01_$ANALYSIS_DATE_WFMT $DIRANAL/anal$(printf %05d $((10#$MIEM)))
        mv  $LETKFDIRRUN/anal/$MEMS/init_$ANALYSIS_DATE_SFMT* $DIRANAL/$MEM/
done

#mv $LETKFDIRRUN/*_sp $DIRANAL
#$LETKFDIR/code/update_wrf_time.exe $LETKFDIRRUN/analemean $ANALYSIS_DATE_WFMT

mkdir -p  ${DIRGUES}/mean
mkdir -p  ${DIRANAL}/mean

cp $LETKFDIRRUN/gues/mean/init_$ANALYSIS_DATE_SFMT* $DIRGUES/mean/
cp $LETKFDIRRUN/anal/mean/init_$ANALYSIS_DATE_SFMT* $DIRANAL/mean/

#Copiamos las observaciones asimiladas.
###cp $LETKFDIRRUN/obs.dat $DIRANAL

# Guardamos un NOUT
###mv $LETKFDIRRUN/NOUT-00000 $DIRANAL

echo "Successfully finished running LETKF"


