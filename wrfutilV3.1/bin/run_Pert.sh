#!/bin/bash
#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

### CONFIGURACION
source $BASEDIR/conf/config.env
source $BASEDIR/conf/step.conf
source $BASEDIR/conf/exp.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source ${BASEDIR}/lib/encolar${QUEUESYS}.sh                     # Selecciona el metodo de encolado segun el systema QUEUESYS elegido

##### FIN INICIALIZACION ######
if [ $STEP -ge 1 ] && [ $WPS_CYCLE -eq 0 ] ; then
   echo "Pert with WPS_CYCLE=0 is run only at STEP=0. And currently STEP=",$STEP        
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
source $BASEDIR/conf/step.conf
source $BASEDIR/conf/exp.conf

BDY_INI_DATE_INT=$(date -d "$BDY_INI_DATE" +"%Y%m%d%H%M%S")
BDY_END_DATE_INT=$(date -d "$BDY_END_DATE" +"%Y%m%d%H%M%S")
INI_BDY_DATE_INT=$(date -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")

#Loop over the expected output files. If all the files are present then skip PERT step.
CDATE=$BDY_INI_DATE
FILE_COUNTER=0
while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $(date -d "$BDY_END_DATE" +"%Y%m%d%H%M%S") ] ; do
   MYFILE=met_em.d01.$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%H:%M:%S" ).nc
   for MEM in $(seq -w $MEM_INI $MEM_END) ; do
      MYPATH=$HISTDIR/WPS/met_em/$(date -d "$INI_BDY_DATE" +"%Y%m%d%H%M%S")/$MEM/
      if ! [[ -e $MYPATH/$MYFILE ]] ; then #My file do not exist
         echo $MYPATH/$MYFILE " does not exist"
         FILE_COUNTER=$(( $FILE_COUNTER + 1 ))
      fi
   done
   CDATE=$(date -u -d "$CDATE UTC + $BDY_FREQ seconds" +"%Y-%m-%d %T")
done


if [ $FILE_COUNTER -eq 0 ] ; then
   echo "We found all the required perturbed met_em files for all members, skiping met_em perturbation"
   ERROR=0
else
   echo $FILE_COUNTER "perturbed met_ems are missing, we need to run pert_met_em."

   #Copy met em files to run pert metem
   mkdir -p $HISTDIR/WPS/met_em/
   ICP=1                                                                      #Counter for the number of cp commands.
   MAX_SIM_CP=$(( 10#$BDY_MEM_END - 10#$BDY_MEM_INI + 1 ))                  #Maximum number of simultaneous cp commands.

   #Creating directories
   echo "Making a copy of met_em files to create the expanded ensemble"
   for MEM in $(seq -w $MEM_INI $MEM_END) ; do
      #member1=$(printf '%02d' $((10#$MEM)))
      mkdir -p $HISTDIR/WPS/met_em/$INI_BDY_DATE_INT/$MEM/
   done

   #Copying the files
   CDATE=$BDY_INI_DATE
   while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $BDY_END_DATE_INT ] ; do
      echo "Copying met ems for date "$CDATE
      for MEM in $(seq -w $MEM_INI $MEM_END) ; do
         CDATE_WPS=$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%T" )
         cp $HISTDIR/WPS/met_em_ori/$INI_BDY_DATE_INT/$(printf '%0'${#MEM_INI}'d' $ICP)/met_em.d01.$CDATE_WPS.nc $HISTDIR/WPS/met_em/$INI_BDY_DATE_INT/$MEM/met_em.d01.$CDATE_WPS.nc &
         ICP=$(( $ICP + 1 ))
         if [ $ICP -gt $MAX_SIM_CP ] ; then
            time wait
            ICP=1
         fi
      done
      CDATE=$(date -u -d "$CDATE UTC + $BDY_FREQ seconds" +"%Y-%m-%d %T")
   done
   time wait


   echo "Creating links to the pert_met_em.exe application"
   #Create links for the met_em_files in the run directory.
   CDATE=$BDY_INI_DATE
   ITIME=1
   while [ $(date -d "$CDATE" +"%Y%m%d%H%M%S") -le $BDY_END_DATE_INT ] ; do
      met_em_time=$(printf '%02d' $((10#$ITIME)))
      for MEM in $(seq -w $MEM_INI $MEM_END) ; do
         CDATE_WPS=$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%T" )
         met_em_file_name=$HISTDIR/WPS/met_em/$INI_BDY_DATE_INT/$MEM/met_em.d01.$CDATE_WPS.nc
         ln -sf $met_em_file_name $PERTDIR/00/ep${met_em_time}$(printf '%05d' $((10#$MEM))) &
      done
      time wait
      for MEM in $(seq -w $BDY_MEM_INI $BDY_MEM_END ) ; do
        CDATE_WPS=$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%T" )
        met_em_file_name=$HISTDIR/WPS/met_em_ori/$INI_BDY_DATE_INT/$MEM/met_em.d01.$CDATE_WPS.nc
        ln -sf $met_em_file_name $PERTDIR/00/eo${met_em_time}$(printf '%05d' $((10#$MEM))) &
      done
      time wait
      ITIME=$(( $ITIME + 1 ))
      CDATE=$(date -u -d "$CDATE UTC + $BDY_FREQ seconds" +"%Y-%m-%d %T")
   done

   echo "Generating the namelist"
   #Prepare the namelist.
   NBV_ORI=$(( 10#$BDY_MEM_END - 10#$BDY_MEM_INI + 1 ))
   NBV_TAR=$(( 10#$MEM_END     - 10#$MEM_INI     + 1 ))
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

   time $MPIEXE ./pert_met_em.exe
   ERROR=$(( $ERROR + $? ))

   cp NOUT-000 $LOGDIR/NOUT-000-PertMetEm_$STEP
fi

EOF

# Parametros de encolamiento
QPROC_NAME=PERTMETEM_$STEP
QNODE=$PERTNODE
QPROC=$PERTPROC
QSKIP=$PERTSKIP
QWALLTIME=$PERTWALLTIME
QWORKPATH=$PERTDIR

echo "Sending pert met em script to the queue"
# Encolar
queue 00 00
check_proc 00 00
if [ $? -ne 0 ] ; then
   echo "Error: Some members do not finish OK"
   echo "Aborting this step"
   exit 1 
fi

echo "Sucessfully finished perturbation of met_em files"





