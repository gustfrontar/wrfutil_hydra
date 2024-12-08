#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

### CONFIGURACION
#Load experiment configuration
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/${EXPTYPE}.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source $BASEDIR/conf/letkf.conf
source $BASEDIR/conf/obs.conf
source ${BASEDIR}/lib/encolar${QUEUESYS}.sh                     # Selecciona el metodo de encolado segun el systema QUEUESYS elegido

if [ ! -z ${PJM_SHAREDTMP} ] ; then
   if [ ${USETMPDIR} -eq 1 ] ; then
      echo "Using Fugaku's shared dir to run WRF"
      LETKFDIR=${PJM_SHAREDTMP}/LETKF
   fi
fi

LETKFDIRRUN=$LETKFDIR/00/   #We need to add 00 in order to be consistent with the paralelization lib.
rm -fr $LETKFDIRRUN
mkdir -p $LETKFDIRRUN

#Descomprimimos el archivo .tar (si es que no fue descomprimido)
if [ ! -e $LETKFDIR/code/letkf.exe ] ; then
   echo "Descomprimiendo LETKF"
   mkdir -p $LETKFDIR/code/
   tar -xf $LETKFDIR/letkf.tar -C $LETKFDIR/code/
   #Si existe el namelist.input lo borramos para que no interfiera
   #con los que crea el sistema de asimilacion
   if [ -e $LETKFDIR/code/letkf.namelist ] ; then
      rm -f $LETKFDIR/code/letkf.namelist
   fi
fi

ln -sf $LETKFDIR/code/* $LETKFDIRRUN

echo "Editing namelist"
cp $NAMELISTDIR/letkf.namelist $LETKFDIRRUN/letkf.namelist
NAMELISTFILE=$LETKFDIRRUN/letkf.namelist
NSLOTS=$(( ($ANALISIS_WIN_FIN-$ANALISIS_WIN_INI)/$ANALISIS_WIN_STEP+1))
NBSLOT=$(( ($ANALISIS_FREC-$ANALISIS_WIN_INI)/$ANALISIS_WIN_STEP+1))
NBV=$(( 10#$MIEMBRO_FIN - 10#$MIEMBRO_INI + 1 ))

#Check if we are going to assimilate radar data and how many radar
#are being assimilated.
if [[ ${OBS_LIST[@]} =~ "RADARC" ]] ; then 
   N_RADAR=${#RADARC_LIST[@]}  #Get the number of radars. 
else 
   N_RADAR=0 
fi

sed -i -e "s|__COV_INFL_MUL__|$COV_INFL_MUL|g"             $NAMELISTFILE
sed -i -e "s|__SP_INFL_ADD__|$SP_INFL_ADD|g"               $NAMELISTFILE
sed -i -e "s|__RELAX_ALPHA_SPREAD__|$RELAX_ALPHA_SPREAD|g" $NAMELISTFILE
sed -i -e "s|__RELAX_ALPHA__|$RELAX_ALPHA|g"               $NAMELISTFILE
sed -i -e "s|__NSLOTS__|$NSLOTS|g"                         $NAMELISTFILE
sed -i -e "s|__NBSLOT__|$NBSLOT|g"                         $NAMELISTFILE
sed -i -e "s|__NBV__|$NBV|g"                               $NAMELISTFILE
sed -i -e "s|__SIGMA_OBS__|$SIGMA_OBS|g"                   $NAMELISTFILE
sed -i -e "s|__SIGMA_OBSV__|$SIGMA_OBSV|g"                 $NAMELISTFILE
sed -i -e "s|__SIGMA_OBS_RADAR__|$SIGMA_OBS_RADAR|g"       $NAMELISTFILE
sed -i -e "s|__SIGMA_OBSZ__|$SIGMA_OBSZ|g"                 $NAMELISTFILE
sed -i -e "s|__SIGMA_OBST__|$SIGMA_OBST|g"                 $NAMELISTFILE
sed -i -e "s|__N_RADAR__|$N_RADAR|g"                       $NAMELISTFILE
sed -i -e "s|__LEV_UPDATE_Q__|$LEV_UPDATE_Q|g"             $NAMELISTFILE
sed -i -e "s|__THRESHOLD_DZ__|$THRESHOLD_DZ|g"             $NAMELISTFILE
sed -i -e "s|__GROSS_ERROR__|$GROSS_ERROR|g"               $NAMELISTFILE
sed -i -e "s|__GROSS_ERROR_REFLECTIVITY__|$GROSS_ERROR_REFLECTIVITY|g"  $NAMELISTFILE

echo "Linking model files"
FECHA_WINDOW_INI=$(date -u -d "$FECHA_INI UTC +$((($ANALISIS_FREC*($PASO-1))+$SPIN_UP_LENGTH+$ANALISIS_WIN_INI)) seconds" +"%Y-%m-%d %T")
ANALYSIS_DATE_WFMT=$(date -u -d "$FECHA_INI UTC +$((($ANALISIS_FREC*($PASO-1))+$SPIN_UP_LENGTH+$ANALISIS_FREC)) seconds" +"%Y-%m-%d_%T")   #WRF format
ANALYSIS_DATE_PFMT=$(date -u -d "$FECHA_INI UTC +$((($ANALISIS_FREC*($PASO-1))+$SPIN_UP_LENGTH+$ANALISIS_FREC)) seconds" +"%Y%m%d%H%M%S")  #Path/folder format

for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN ) ; do
   for ISLOT in $(seq 0 $(($NSLOTS-1))) ; do
      FECHA_SLOT=$(date -u -d "$FECHA_WINDOW_INI UTC +$(( 10#$ISLOT*$ANALISIS_WIN_STEP)) seconds" +"%Y-%m-%d_%T")
      ln -sf ${WRFDIR}/${MIEM}/wrfout_d01_$FECHA_SLOT ${LETKFDIRRUN}/gs$(printf %02d $((10#$ISLOT+1)))$(printf %05d $((10#$MIEM)))
      echo ${WRFDIR}/${MIEM}/wrfout_d01_$FECHA_SLOT ${LETKFDIRRUN}/gs$(printf %02d $((10#$ISLOT+1)))$(printf %05d $((10#$MIEM)))
   done	   
done

exit

cp  ${WRFDIR}/${MIEMBRO_FIN}/wrfout_d01_$ANALYSIS_DATE_WFMT ${LETKFDIRRUN}/guesemean
cp  ${WRFDIR}/${MIEMBRO_FIN}/wrfout_d01_$ANALYSIS_DATE_WFMT ${LETKFDIRRUN}/analemean

echo "Linking obs files" 
CTIME=$FECHA_WINDOW_INI
for SLOT in $(seq -f "%02g" 1 $(($NSLOTS))) ; do  #Loop over time slot within the assimilation window.
     CTIME_OFMT=$(date -u -d "$CTIME UTC" +"%Y%m%d%H%M%S")
     NOBS=1 #Counter for the concatenation of conventional obs. files.
     for OBS_TYPE in "${OBS_LIST[@]}" ; do 
        echo "Linking obs type : $OBS_TYPE"
        if [ $OBS_TYPE == "RADARC" ] ; then
          NRAD=1
          for RADARC in "${RADARC_LIST[@]}" ; do
            OBSFILE=${OBSPATH}/${OBS_TYPE}/WRF_${OBS_TYPE}_${RADARC}_${CTIME_OFMT}.dat 
            [[ -e ${OBSFILE} ]] && echo "Linking $OBSFILE"
            [[ -e ${OBSFILE} ]] && ln -sf ${OBSFILE} ${LETKFDIRRUN}/rad$(printf %02d $((10#$SLOT)))$(printf %02d $((10#$NRAD))).dat
            NRAD=$(( $NRAD + 1 ))
          done
        else
          #Concatenate conventional sources of observations into one single file. 
          OBSFILE=${OBSPATH}/${OBS_TYPE}/WRF_${OBS_TYPE}_${CTIME_OFMT}.dat
          echo $OBSFILE
          if [ $NOBS -eq 1 ] ; then 
             [[ -e ${OBSFILE} ]] &&  echo "Linking ${OBSFILE}"
             [[ -e ${OBSFILE} ]] &&  cp  ${OBSFILE} ${LETKFDIRRUN}/obs$(printf %02d $((10#$SLOT))).dat
          else 
             [[ -e ${OBSFILE} ]] &&  echo "Linking ${OBSFILE}"
             [[ -e ${OBSFILE} ]] &&  cat ${OBSFILE} ${LETKFDIRRUN}/obs$(printf %02d $((10#$SLOT))).dat >> ${LETKFDIRRUN}/obs$(printf %02d $((10#$SLOT))).dat
          fi
          NOBS=$(( $NOBS + 1 ))
        fi
     done
   CTIME=$(date -u -d "$CTIME UTC +$ANALYSIS_WIN_STEP seconds" +"%Y-%m-%d %T")
done

#script de ejecucion
read -r -d '' QSCRIPTCMD << "EOF"
  cd $LETKFDIR/00/
  export FORT90L=""
  time $MPIEXE  ./letkf.exe 
  ERROR=$(( $ERROR + $? ))
EOF

# Parametros de encolamiento
QPROC_NAME=LETKF_$PASO
QNODE=$LETKFNODE
QPROC=$LETKFPROC
QTHREAD=$LETKFTHREAD
QWALLTIME=$LETKFWALLTIME
QCONF=${EXPTYPE}.conf
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
echo "Guardando analisis $MIEMBRO_INI - $MIEMBRO_FIN"
DIRANAL=$HISTDIR/ANAL/$ANALYSIS_DATE_PFMT
for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN ) ; do
        mkdir -p  ${DIRANAL}/
	echo "Updating the date in the analysis file " $WRFDIR/$MIEM/wrfout_d01_$ANALYSIS_DATE_WFMT $ANALYSIS_DATE_WFMT
	$LETKFDIR/code/update_wrf_time.exe $WRFDIR/$MIEM/wrfout_d01_$ANALYSIS_DATE_WFMT $ANALYSIS_DATE_WFMT
        mv  $WRFDIR/$MIEM/wrfout_d01_$ANALYSIS_DATE_WFMT $DIRANAL/anal$(printf %05d $((10#$MIEM)))
done

#mv $LETKFDIRRUN/*_sp $DIRANAL
$LETKFDIR/code/update_wrf_time.exe $LETKFDIRRUN/analemean $ANALYSIS_DATE_WFMT
mv $LETKFDIRRUN/guesemean $DIRANAL   #Guess mean
mv $LETKFDIRRUN/analemean $DIRANAL   #Analysis mean

#Copiamos las observaciones asimiladas.
cp $LETKFDIRRUN/obs.dat $DIRANAL

# Guardamos un NOUT
mv $LETKFDIRRUN/NOUT-00000 $DIRANAL

echo "Successfully finished running WRF-LETKF"


