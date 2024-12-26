#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
# Ported to hydra and Fugaku
# Date: 11/2023
#############

### CONFIGURACION
#Load experiment configuration
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/dates.conf
source $BASEDIR/conf/${EXPTYPE}.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source ${BASEDIR}/lib/encolar${QUEUESYS}.sh       

##### To use Fugaku's shared filesystem  ######
if [ ! -z ${PJM_SHAREDTMP} ] && [  ${USETMPDIR} -eq 1 ] ; then 
   echo "Using Fugaku's shared dir to run WRF"
   WRFDIR=${PJM_SHAREDTMP}/WRF
fi

cd $WRFDIR

#Desglozamos las fechas
read -r IY IM ID IH Im Is  <<< $(date -u -d "$INI_DATE_FCST UTC" +"%Y %m %d %H %M %S")
read -r FY FM FD FH Fm Fs  <<< $(date -u -d "$END_DATE_FCST UTC" +"%Y %m %d %H %M %S")


#Prepare the namelists
if [ $MULTIMODEL -eq 0 ] ; then  
   MULTIMODEL_CONF_INI=$MODEL_CONF;MULTIMODEL_CONF_FIN=$MODEL_CONF
   echo "We will run a single configuration ensemble"
else 
   echo "We will run a multimodel configuration ensemble with $NCONF configurations"
   MULTIMODEL_CONF_INI=1;MULTIMODEL_CONF_FIN=$NCONF
fi

for ICONF in $(seq -w $MULTIMODEL_CONF_INI $MULTIMODEL_CONF_FIN) ; do
   WRFCONF=$(printf '%03d' $((10#$ICONF)))
   cp $NAMELISTDIR/namelist.input.${WRFCONF} $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__START_YEAR__|$IY|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__START_MONTH__|$IM|g"                              $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__START_DAY__|$ID|g"                                $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__START_HOUR__|$IH|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__START_MINUTE__|$Im|g"                             $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__START_SECOND__|$Is|g"                             $WRFDIR/namelist.input.${WRFCONF} 
   sed -i -e "s|__END_YEAR__|$FY|g"                                 $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__END_MONTH__|$FM|g"                                $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__END_DAY__|$FD|g"                                  $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__END_HOUR__|$FH|g"                                 $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__END_MINUTE__|$Fm|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__END_SECOND__|$Fs|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__INTERVALO_WPS__|$FCST_BDY_FREQ|g"                 $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__INTERVALO_WRF__|$((10#$FCST_OFREQ/60))|g"         $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__E_WE__|$E_WE|g"                                   $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__E_SN__|$E_SN|g"                                   $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__DX__|$DX|g"                                       $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__DY__|$DY|g"                                       $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__DT__|$DT|g"                                       $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__IOTYPE__|$IOTYPE|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__METLEV__|$METLEV|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__NUMTILE__|$NUMTILE|g"                             $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__NIOT__|$NIOT|g"                                   $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__NIOG__|$NIOG|g"                                   $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__E_VERT__|$E_VERT|g"                               $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__P_TOP__|$P_TOP|g"                                 $WRFDIR/namelist.input.${WRFCONF}
   sed -i -e "s|__RADT__|$RADT|g"                                   $WRFDIR/namelist.input.${WRFCONF}
done

#Untar executables
if [ ! -e $WRFDIR/code/real.exe ] ; then
   echo "Decompressing executables ..."
   mkdir -p $WRFDIR/code
   cd $WRFDIR
   tar -xf wrf.tar -C $WRFDIR/code
   tar -xf wrfda.tar -C $WRFDIR/code
   tar -xf pert_met_em.tar -C $WRFDIR/code
   if [ -e $WRFDIR/code/namelist.input ] ; then
      rm -f $WRFDIR/code/namelist.input
   fi
fi


#Build the script to run REAL/DA-UPDATE-BC/WRF
read -r -d '' QSCRIPTCMD << "EOF"
ulimit -s unlimited
source $BASEDIR/conf/model.conf
source $BASEDIR/conf/dates.conf
echo "Processing member $MIEM" 

if [ ! -z ${PJM_SHAREDTMP} ] && [  ${USETMPDIR} -eq 1 ] ; then
   echo "Using Fugaku's shared dir to run WRF"
   WRFDIR=${PJM_SHAREDTMP}/WRF/    #WRFDIR is redefined here
fi

if [ $MULTIMODEL -eq 0 ] ; then
   NLCONF=$(printf '%03d' ${MODEL_CONF} )
else 
   NLCONF=$(printf '%03d' $(( ( ( 10#$MIEM - 1 ) % 10#$NCONF ) + 1 )) )
fi
cp $WRFDIR/namelist.input.${NLCONF} $WRFDIR/$MIEM/namelist.input

ln -sf $WRFDIR/code/* . 

if [ ! -z ${PJM_SHAREDTMP} ] && [  ${USETMPDIR} -eq 1 ] ; then
   MET_EM_DIR=${PJM_SHAREDTMP}/HIST/WPS/met_em/$(date -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")/$MIEM/    #MET_EM_DIR is redefined here
else 
   MET_EM_DIR=$HISTDIR/WPS/met_em/$(date -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")/$MIEM/
fi

#Loop over the expected output files. If all the files are present then skip WRF step.
CDATE=$INI_DATE_FCST
FILE_COUNTER=0
while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $(date -d "$END_DATE_FCST" +"%Y%m%d%H%M%S") ] ; do
   MYPATH=./
   MYFILE=$MYPATH/wrfout_d01_$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%H:%M:%S" )
   if ! [[ -e $MYFILE ]] ; then #My file do not exist
      echo "Not found file: " $MYFILE
      FILE_COUNTER=$(( $FILE_COUNTER + 1 ))
   fi
   CDATE=$(date -u -d "$CDATE UTC + $FCST_OFREQ seconds" +"%Y-%m-%d %T")
done

if [ $FILE_COUNTER -eq 0 ] ; then 
   #All the required files are already there. 
   echo "We found all the required wrfout files for member $MIEM, skiping DA_UPDATE_BC, REAL and WRF."
   ERROR=0

else
   #We need to run DA_UPDATE_BC, REAL and WRF.
   if [ $FCST_BDY_FREQ -eq $BDY_FREQ ] ; then
      ln -sf $MET_EM_DIR/met_em* $WRFDIR/$MIEM/
   else    
      #We will conduct interpolation of the met_em files.
      echo "Interpolating files in time to reach $FCST_BDY_FREQ time frequency."
      CDATE=$INI_DATE_FCST
      WPS_FILE_DATE_FORMAT="%Y-%m-%d_%H:%M:%S"

      while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $(date -d "$END_DATE_FCST" +"%Y%m%d%H%M%S") ] ; do
        FILE_TAR=met_em.d01.$(date -u -d "$CDATE UTC" +"$WPS_FILE_DATE_FORMAT" ).nc

        if [ -e $MET_EM_DIR/$FILE_TAR ] ; then #Target file exists. 
           ln -sf $MET_EM_DIR/$FILE_TAR  ./
        else 
           DATE_INI=$(date_floor "$CDATE" $BDY_FREQ )
           DATE_END=$(date -u -d "$DATE_INI UTC + $BDY_FREQ seconds" +"%Y-%m-%d %T")
           FILE_INI=met_em.d01.$(date -u -d "$DATE_INI UTC" +"$WPS_FILE_DATE_FORMAT" ).nc
           FILE_END=met_em.d01.$(date -u -d "$DATE_END UTC" +"$WPS_FILE_DATE_FORMAT" ).nc
           #File does not exist. Interpolate data in time to create it.  
           echo "&general                         "  > ./pertmetem.namelist
           echo "/                                " >> ./pertmetem.namelist
           echo "&interp                          " >> ./pertmetem.namelist
           echo "time_ini=0                       " >> ./pertmetem.namelist
           echo "time_end=$BDY_FREQ               " >> ./pertmetem.namelist
           echo "time_tar=$((($(date -d "$CDATE" +%s) - $(date -d "$DATE_INI" +%s))))    " >> ./pertmetem.namelist
           echo "date_tar='$(date -u -d "$CDATE UTC" +"$WPS_FILE_DATE_FORMAT" )'         " >> ./pertmetem.namelist
           echo "file_ini='$MET_EM_DIR/$FILE_INI'   " >> ./pertmetem.namelist
           echo "file_end='$MET_EM_DIR/$FILE_END'   " >> ./pertmetem.namelist
           echo "file_tar='$WRFDIR/$MIEM/$FILE_TAR' " >> ./pertmetem.namelist
           echo "/                                " >> ./pertmetem.namelist
           echo "Running INTERP_MET_EM for member $MIEM"
           $MPIEXESERIAL ./interp_met_em.exe > interp_met_em.log 
           ERROR=$(( $ERROR + $? ))
           #ln -sf $WPSDIR/$MIEM/$FILE_TAR  ./
           #TODO: Check if would be better to create the new files in the WPS_MET_EM_DIR so
           #Interpolated files can be used again in the next cycle.
        fi
        #Update CDATE
        CDATE=$(date -u -d "$CDATE UTC + $FCST_BDY_FREQ seconds" +"%Y-%m-%d %T")
     done
     if [ $ERROR -gt 0 ] ; then
        echo "Error: INTERP MET EM step finished with errors"   
     fi
   fi

   #If there are no errors so far proceed with the REAL
   if [ $ERROR -eq 0 ] ; then 
      echo "Running REAL for member $MIEM"
      time $MPIEXE $WRFDIR/$MIEM/real.exe 
      ERROR=$(( $ERROR + $? ))
      mv rsl.error.0000 ./real_${STEP}_${MIEM}.log
      if [ $ERROR -gt 0 ] ; then
         echo "Error: REAL step finished with errors"   
      fi
   fi

   #If there are no errors so far proceed with the DA_UPDATE_BC
   if [ $ERROR -eq 0 ] ; then
      echo "Running da_update_bc"
      cp $NAMELISTDIR/parame* .
      mv $WRFDIR/$MIEM/wrfinput_d01 $WRFDIR/$MIEM/wrfinput_d01.org
      cp $HISTDIR/ANAL/$(date -u -d "$INI_DATE_FCST" +"%Y%m%d%H%M%S")/anal$(printf %05d $((10#$MIEM))) $WRFDIR/$MIEM/wrfinput_d01
      echo "Running DA_UPDATE_BC for member $MIEM"
      ln -sf $WRFDIR/code/da_update_bc.exe $WRFDIR/$MIEM/da_update_bc.exe
      time $MPIEXESERIAL $WRFDIR/$MIEM/da_update_bc.exe > ./da_update_bc_${STEP}_${MIEM}.log
      ERROR=$(( $ERROR + $? ))
   fi
   if [ $ERROR -gt 0 ] ; then
      echo "Error: DA_UPDATE_BC step finished with errors"   
   fi  

   #If there are no errors so far proceed with WRF
   if [ $ERROR -eq 0 ] ; then
      echo "Running WRF for member $MIEM"
      time $MPIEXE $WRFDIR/$MIEM/wrf.exe 
      ERROR=$(( $ERROR + $? ))
      mv rsl.error.0000 ./wrf_${STEP}_${MIEM}.log
   fi
   if [ $ERROR -gt 0 ] ; then
      echo "Error: WRF step finished with errors"   
   else
      #Copy the guess files corresponding to the analysis time.
      OUTPUTPATH="$HISTDIR/DAFCST/$(date -u -d "$INI_DATE_FCST UTC" +"%Y%m%d%H%M%S")/$MIEM/"
      mkdir -p $OUTPUTPATH
      mv $WRFDIR/$MIEM/wrfout_d01_* $OUTPUTPATH/
      mv $WRFDIR/$MIEM/*.log*       $OUTPUTPATH/
      mv $WRFDIR/$MIEM/namelist*    $OUTPUTPATH/
   fi
fi

EOF

#Node / core distribution parameters
QNODE=$WRFNODE
QPROC=$WRFPROC
TPROC=$WRFTPROC
QTHREAD=$WRFTHREAD
QWALLTIME=$WRFWALLTIME
QPROC_NAME=DAFCST_${STEP}
QCONF=${EXPTYPE}.conf
QWORKPATH=$WRFDIR
QSKIP=$WRFSKIP

#Execute the job 
queue $FORECAST_MEMBER_INI $FORECAST_MEMBER_END
time check_proc $FORECAST_MEMBER_INI $FORECAST_MEMBER_END
if [ $? -ne 0 ] ; then
   echo "Error: Some members do not finish OK"
   echo "Aborting this step"
   exit 1
fi

