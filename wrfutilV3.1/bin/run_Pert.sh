#!/bin/bash

#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

### CONFIGURACION
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/${EXPTYPE}.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source ${BASEDIR}/lib/encolar${QUEUESYS}.sh                     # Selecciona el metodo de encolado segun el systema QUEUESYS elegido

##### FIN INICIALIZACION ######
if [ $PASO -ge 1 ] && [ $WPS_CYCLE -eq 0 ] ; then
   echo "Pert with WPS_CYCLE=0 is run only at PASO=0. And currently PASO=",$PASO         
   return 0
fi

#Decompress tar files
if [ ! -e $PERTDIR/code/pert_met_em.exe ] ; then
   echo "Decompressing Pert Met Em"
   mkdir -p $PERTDIR/code/
   tar -xf $PERTDIR/pert_met_em.tar -C $PERTDIR/code
   #Remove pertmetem.namelist (it will be created later from the templates)
   if [ -e $PERTDIR/code/pertmetem.namelist ] ; then
      rm -f $PERTDIR/code/pertmetem.namelist
   fi
fi

cd $PERTDIR

#script de ejecucion
read -r -d '' QSCRIPTCMD << "EOF"


#Define some variables that depends if we are running in assimilation
#mode or in forecast mode.
if [ "$EXPTYPE" = "assimilation" ] ; then
   CYCLE_FREQ=$ANALISIS_FREC
   LEAD=$ANALISIS_WIN_FIN
   SPIN_UP=$SPIN_UP_LENGTH
elif [ $EXPTYPE -eq "forecast" ] ; then
   CYCLE_FREQ=$FORECAST_INI_FREQ
   LEAD=$FORECAST_LEAD_TIME
   SPIN_UP=0
fi

#Define the start date and enddate of the WPS files to be generated.
#If WPS_CYCLE is 0 then we will process all the files from FECHA_INI to FECHA_FIN
if [ $WPS_CYCLE -eq 1 ] ; then
   if [ $PASO -eq 0 ] ; then
      INI_DATE_STEP=$FECHA_INI
      END_DATE_STEP=$(date -u -d "$FECHA_INI UTC +$SPIN_UP seconds" +"%Y-%m-%d %T")
   else
      INI_DATE_STEP=$(date -u -d "$FECHA_INI UTC +$(($CYCLE_FREQ*($PASO-1)+$SPIN_UP)) seconds" +"%Y-%m-%d %T")
      END_DATE_STEP=$(date -u -d "$INI_DATE_STEP UTC +$LEAD seconds" +"%Y-%m-%d %T")
   fi
else  # WPS_CYCLE -eq 0 
   INI_DATE_STEP=$FECHA_INI
   END_DATE_STEP=$FECHA_FIN
fi

INI_DATE_BDY=$(date_floor "$INI_DATE_STEP" $INTERVALO_INI_BDY )   #get the closest initialization of the boudndary data in which BDY data is available.
INI_DATE_STEP=$(date_floor "$INI_DATE_STEP" $INTERVALO_BDY )      #Get the closest prior date in which BDY data is available.
END_DATE_STEP=$(date_ceil  "$END_DATE_STEP" $INTERVALO_BDY )      #Get the closest posterior date in which BDY data is available.

INI_DATE_STEP_INT=$(date -d "$INI_DATE_STEP" +"%Y%m%d%H%M%S")
END_DATE_STEP_INT=$(date -d "$END_DATE_STEP" +"%Y%m%d%H%M%S")
INI_DATE_BDY_INT=$(date -d "$INI_DATE_BDY" +"%Y%m%d%H%M%S")

#Loop over the expected output files. If all the files are present then skip PERT step.
CDATE=$INI_DATE_STEP
FILE_COUNTER=0
while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $(date -d "$END_DATE_STEP" +"%Y%m%d%H%M%S") ] ; do
   MYFILE=met_em.d01.$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%H:%M:%S" ).nc
   for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
      MYPATH=$HISTDIR/WPS/met_em/$(date -d "$INI_DATE_BDY" +"%Y%m%d%H%M%S")/$MIEM/
      if ! [[ -e $MYPATH/$MYFILE ]] ; then #My file do not exist
         echo $MYPATH/$MYFILE " does not exist"
         FILE_COUNTER=$(( $FILE_COUNTER + 1 ))
      fi
   done
   CDATE=$(date -u -d "$CDATE UTC + $INTERVALO_BDY seconds" +"%Y-%m-%d %T")
done


if [ $FILE_COUNTER -eq 0 ] ; then
   echo "We found all the required perturbed met_em files for all members, skiping met_em perturbation"
   exit 0
else
   echo $FILE_COUNTER "perturbed met_ems are missing, we need to run pert_met_em."
fi

#Copy met em files to run pert metem
mkdir -p $HISTDIR/WPS/met_em/
ICP=1                                                                      #Counter for the number of cp commands.
MAX_SIM_CP=$(( 10#$BDY_MIEMBRO_FIN - 10#$BDY_MIEMBRO_INI + 1 ))                  #Maximum number of simultaneous cp commands.

#Creating directories
echo "Making a copy of met_em files to create the expanded ensemble"
for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
   #member1=$(printf '%02d' $((10#$MIEM)))
   mkdir -p $HISTDIR/WPS/met_em/$INI_DATE_BDY_INT/$MIEM/
done

#Copying the files
CDATE=$INI_DATE_STEP
while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $END_DATE_STEP_INT ] ; do
   echo "Copying met ems for date "$CDATE
   for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
       CDATE_WPS=$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%T" )
       cp $HISTDIR/WPS/met_em_ori/$INI_DATE_BDY_INT/$(printf '%0'${#MIEMBRO_INI}'d' $ICP)/met_em.d01.$CDATE_WPS.nc $HISTDIR/WPS/met_em/$INI_DATE_BDY_INT/$MIEM/met_em.d01.$CDATE_WPS.nc &
       ICP=$(( $ICP + 1 ))
       if [ $ICP -gt $MAX_SIM_CP ] ; then
          time wait
          ICP=1
       fi
   done
   CDATE=$(date -u -d "$CDATE UTC + $INTERVALO_BDY seconds" +"%Y-%m-%d %T")
done
time wait


echo "Creating links to the pert_met_em.exe application"
#Create links for the met_em_files in the run directory.
CDATE=$INI_DATE_STEP
ITIME=1
while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $END_DATE_STEP_INT ] ; do
 met_em_time=$(printf '%02d' $((10#$ITIME)))
  for MIEM in $(seq -w $MIEMBRO_INI $MIEMBRO_FIN) ; do
    CDATE_WPS=$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%T" )
    met_em_file_name=$HISTDIR/WPS/met_em/$INI_DATE_BDY_INT/$MIEM/met_em.d01.$CDATE_WPS.nc
    ln -sf $met_em_file_name $PERTDIR/00/ep${met_em_time}$(printf '%05d' $((10#$MIEM))) &
   done
   time wait
   for MIEM in $(seq -w $BDY_MIEMBRO_INI $BDY_MIEMBRO_FIN ) ; do
     CDATE_WPS=$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%T" )
     met_em_file_name=$HISTDIR/WPS/met_em_ori/$INI_DATE_BDY_INT/$MIEM/met_em.d01.$CDATE_WPS.nc
     ln -sf $met_em_file_name $PERTDIR/00/eo${met_em_time}$(printf '%05d' $((10#$MIEM))) &
   done
   time wait
   ITIME=$(( $ITIME + 1 ))
   CDATE=$(date -u -d "$CDATE UTC + $INTERVALO_BDY seconds" +"%Y-%m-%d %T")
done

echo "Generating the namelist"
#Prepare the namelist.
NBV_ORI=$(( 10#$BDY_MIEMBRO_FIN - 10#$BDY_MIEMBRO_INI + 1 ))
NBV_TAR=$(( 10#$MIEMBRO_FIN     - 10#$MIEMBRO_INI     + 1 ))
NTIMES=$(( 10#$ITIME - 1 ))

cp $NAMELISTDIR/pertmetem.namelist $PERTDIR/00/pertmetem.namelist

#Edit the namelist from the template
sed -i -e "s|__NBV_ORI__|$NBV_ORI|g"   $PERTDIR/00/pertmetem.namelist
sed -i -e "s|__NBV_TAR__|$NBV_TAR|g"   $PERTDIR/00/pertmetem.namelist
sed -i -e "s|__NTIMES__|$NTIMES|g"     $PERTDIR/00/pertmetem.namelist
sed -i -e "s|__METHOD__|$PERT_TYPE|g"  $PERTDIR/00/pertmetem.namelist
sed -i -e "s|__SIGMA__|$PERT_AMP|g"    $PERTDIR/00/pertmetem.namelist
sed -i -e "s|__NITER__|$NITER|g"       $PERTDIR/00/pertmetem.namelist

ln -sf $PERTDIR/code/pert_met_em.exe $PERTDIR/00/

time $MPIEXE  ./pert_met_em.exe
ERROR=$(( $ERROR + $? ))

cp NOUT-000 $LOGDIR/NOUT-000-PertMetEm_$PASO

EOF

# Parametros de encolamiento
QPROC_NAME=PERTMETEM_$PASO
QNODE=$PERTNODE
QPROC=$PERTPROC
QTHREAD=$PERTTHREAD
QWALLTIME=$PERTWALLTIME
QCONF=${EXPTYPE}.conf
QWORKPATH=$PERTDIR

echo "Sending pert met em script to the queue"
# Encolar
queue 00 00
check_proc 00 00

echo "Sucessfully finished perturbation of met_em files"





