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
source $BASEDIR/conf/exp.conf
source $BASEDIR/conf/step.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source $BASEDIR/lib/encolar${QUEUESYS}.sh       
if [ $MACHINE == "FUGAKU" ] ;then
  source $TOPDIR/setup_spack.sh
fi


##### To use Fugaku's shared filesystem  ######
if [ ! -z ${PJM_SHAREDTMP} ] && [  ${USETMPDIR} -eq 1 ] ; then 
   echo "Using Fugaku's shared dir to run WRF"
   WRFDIR=${PJM_SHAREDTMP}/WRF
fi

cd $SCALEDIR
echo "=========================="
echo "Forecast for STEP " $STEP
echo "Starting at " $INI_DATE_FCST
echo "Ending   at " $END_DATE_FCST

#Desglozamos las fechas
read -r IY IM ID IH Im Is  <<< $(date -u -d "$INI_DATE_FCST UTC" +"%Y %m %d %H %M %S")
read -r FY FM FD FH Fm Fs  <<< $(date -u -d "$END_DATE_FCST UTC" +"%Y %m %d %H %M %S")

#Prepare the namelists
if [ $MULTIMODEL -eq 0 ] ; then  
   MULTIMODEL_CONF_INI=$MODEL_CONF;MULTIMODEL_CONF_FIN=$MODEL_CONF
   echo "We will run a single configuration ensemble"
else 
   echo "SCALE does not support multimodel configuration. MULTIMODEL is set to 0."
   MULTIMODEL=0
   MULTIMODEL_CONF_INI=$MODEL_CONF;MULTIMODEL_CONF_FIN=$MODEL_CONF
   echo "We will run a single configuration ensemble"
fi

# Edit the namelist.
cp $NAMELISTDIR/namelist.scale $SCALEDIR/
cp $NAMELISTDIR/namelist.sno_fcst $SCALEDIR/

#Edit the namelist from the template

if [ $MAP_PROJ == "lambert" ];then
  MAP_PROJ_SCALE="LC"
fi

SCALEDATA="$SCALEPATH/data"
INI_SCALE="$(date -ud "$INI_DATE_FCST" +'%Y,%m,%d,%H,%M,%S' )"
E_WE_LOC=$((E_WE / NPROC_X))
E_SN_LOC=$((E_SN / NPROC_Y))
FCSTLEN=$(( $(date -ud "$END_DATE_FCST" +%s) - $(date -ud "$INI_DATE_FCST" +%s) ))
INT_R_SCALE=$FCSTLEN
if [ $EXPTYPE == "DACYCLE" ] ;then
   INT_F_SCALE=$ANALYSIS_WIN_STEP
else 
   INT_F_SCALE=$FCST_OFREQ
fi

sed -i -e "s|__INI_SCALE__|$INI_SCALE|g"      $SCALEDIR/namelist.scale
sed -i -e "s|__LEN__|$FCSTLEN|g"              $SCALEDIR/namelist.scale
sed -i -e "s|__INT_R_SCALE__|$INT_R_SCALE|g"   $SCALEDIR/namelist.scale
sed -i -e "s|__INT_F_SCALE__|$INT_F_SCALE|g"   $SCALEDIR/namelist.scale
sed -i -e "s|__E_WE_LOC__|$E_WE_LOC|g"        $SCALEDIR/namelist.scale
sed -i -e "s|__E_SN_LOC__|$E_SN_LOC|g"        $SCALEDIR/namelist.scale
sed -i -e "s|__NPROC_X__|$NPROC_X|g"          $SCALEDIR/namelist.scale
sed -i -e "s|__NPROC_Y__|$NPROC_Y|g"          $SCALEDIR/namelist.scale
sed -i -e "s|__DX__|$DX|g"                    $SCALEDIR/namelist.scale
sed -i -e "s|__DY__|$DY|g"                    $SCALEDIR/namelist.scale
sed -i -e "s|__MAP_PROJ__|$MAP_PROJ_SCALE|g"  $SCALEDIR/namelist.scale
sed -i -e "s|__REF_LAT__|$REF_LAT|g"          $SCALEDIR/namelist.scale
sed -i -e "s|__REF_LON__|$REF_LON|g"          $SCALEDIR/namelist.scale
sed -i -e "s|__TRUELAT1__|$TRUELAT1|g"        $SCALEDIR/namelist.scale
sed -i -e "s|__TRUELAT2__|$TRUELAT2|g"        $SCALEDIR/namelist.scale
sed -i -e "s|__SCALEDATA__|$SCALEDATA|g"      $SCALEDIR/namelist.scale
sed -i -e "s|__HISTDIR__|fcst|g"              $SCALEDIR/namelist.scale
sed -i -e "s|__DT__|$DT|g"                    $SCALEDIR/namelist.scale
sed -i -e "s|__RADT__|$RADT|g"                $SCALEDIR/namelist.scale
sed -i -e "s|__DYNDT__|$DYNDT|g"              $SCALEDIR/namelist.scale
sed -i -e "s|__LNDDT__|$DT2|g"                $SCALEDIR/namelist.scale
sed -i -e "s|__OCNDT__|$DT2|g"                $SCALEDIR/namelist.scale
sed -i -e "s|__URBDT__|$DT2|g"                $SCALEDIR/namelist.scale

###
# Create the script that will run scale for each ensemble member.
###
echo "Running the scale for the step" $STEP

read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
source $BASEDIR/conf/model.conf
source $BASEDIR/conf/step.conf
source $BASEDIR/conf/exp.conf
echo "Processing member $MEM"

mkdir -p $SCALEDIR/$MEM
cd $SCALEDIR/$MEM

ln -sf $SCALEDIR/scale-rm${SCALE_opt} ./scale-rm
ln -sf $SCALEDIR/sno${SCALE_opt} ./sno
cp $SCALEDIR/namelist.scale . 
cp $SCALEDIR/namelist.sno_fcst . 

ln -sf  $HISTDIR/const/topo .
ln -sf  $HISTDIR/const/landuse .
mkdir -p bdy
mkdir -p init
mkdir -p fcst
mkdir -p log/scale
mkdir -p log/sno
DIR_DATE_INI=$(date -u -d "$INI_DATE_FCST UTC" +"%Y%m%d%H%M%S")
if [ "$EXPTYPE" == "DAFCST" ]; then
   ln -sf $HISTDIR/ANAL/${DIR_DATE_INI}/$MEM/*.nc ./init/
elif [ "$EXPTYPE" == "DACYCLE" ] && [ $STEP -gt 0 ]; then 
   ln -sf $HISTDIR/ANAL/${DIR_DATE_INI}/$MEM/*.nc ./init/
else
   ln -sf $HISTDIR/init/${DIR_DATE_INI}/$MEM/*.nc ./init/
fi
ln -sf $HISTDIR/bdy/${DIR_DATE_INI}/$MEM/*.nc ./bdy/

export FORT90L=${SCALE_RUNTIME_FLAGS}

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SCALE_NETCDF_C}/lib:${SCALE_NETCDF_F}/lib:${SCALE_PNETCDF}/lib:${SCALE_HDF}/lib:$LD_LIBRARY_PATH

echo "Running scale for member $MEM"
$MPIEXE ./scale-rm namelist.scale
ERROR=$(( $ERROR + $? ))

if [ $ERROR -gt 0 ] ; then
   echo "Error: SCALE step finished with errors"   
else
  echo "Running SNO for member $MEM"
  $MPIEXESERIAL ./sno namelist.sno_fcst
  if [ $? -gt 0 ] ; then
     echo "Warning: SNO finished with errors"
  fi
fi

#If there are no errors so far copy the files to their final destination
if [ $ERROR -eq 0 ] ; then
   ANALYSIS_DATE_PFMT=$(date -u -d "$ANALYSIS_DATE UTC" +"%Y%m%d%H%M%S")    
   ANALYSIS_DATE_SFMT=$(date -u -d "$ANALYSIS_DATE UTC" +"%Y%m%d-%H%M%S.000")     
   if [ $EXPTYPE == "DACYCLE" ] ; then
      #Copy the guess files corresponding to the analysis time.
      if [[ ! -z "$SAVEGUESS" ]] && [[ $SAVEGUESS -eq 1 ]] ; then
         OUTPUTPATH="$HISTDIR/GUES/$ANALYSIS_DATE_PFMT/$MEM"
         mkdir -p $OUTPUTPATH
         echo "Copying first guess files"
         cp $SCALEDIR/$MEM/init/init_${ANALYSIS_DATE_SFMT}*.nc $OUTPUTPATH/
         mv $SCALEDIR/$MEM/log/scale/* $OUTPUTPATH/
         mv $SCALEDIR/$MEM/namelist*   $OUTPUTPATH/
      fi
      if [ $STEP -eq 0  ] ; then  #Copy the spin up output as the analysis for the next cycle.
         OUTPUTPATH="$HISTDIR/ANAL/$ANALYSIS_DATE_PFMT/$MEM"
         mkdir -p $OUTPUTPATH
         echo "Copying spinup files"
         mv $SCALEDIR/$MEM/init/init_${ANALYSIS_DATE_SFMT}*.nc   $OUTPUTPATH/
         mv $SCALEDIR/$MEM/log/scale/* $OUTPUTPATH/
         mv $SCALEDIR/$MEM/namelist*   $OUTPUTPATH/
      fi
   elif [ $EXPTYPE == "DAFCST" ] || [ $EXPTYPE == "FCST" ] ; then
      #Copy history files to its final destionation.
      OUTPUTPATH="$HISTDIR/$EXPTYPE/${DIR_DATE_INI}/${MEM}/"
      mkdir -p $OUTPUTPATH
      mv $SCALEDIR/$MEM/fcst/*.nc   $OUTPUTPATH/
      mv $SCALEDIR/$MEM/log/scale/* $OUTPUTPATH/
      mv $SCALEDIR/$MEM/namelist*   $OUTPUTPATH/
   fi
fi


EOF

#Node / core distribution parameters
QNODE=$SCALENODE
QPROC=$SCALEPROC
QSKIP=$SCALESKIP
QOMP=$SCALEOMP
QWALLTIME=$SCALEWALLTIME
QPROC_NAME=ENSFCST_${STEP}
QWORKPATH=$SCALEDIR

#Execute the job 
echo "Time taken by scale"
echo "queue $MEM_INI $MEM_END "
queue $MEM_INI $MEM_END 
time check_proc $MEM_INI $MEM_END



