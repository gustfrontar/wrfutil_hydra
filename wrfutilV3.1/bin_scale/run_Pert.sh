#!/bin/bash
#############
# Servicio Meteorologico Nacional
# Autores: Maximiliano A. Sacco; Yanina Garcia Skabar; Maria Eugenia Dillon; Cynthia Mariana Matsudo
# Fecha: 03/2019
#############

BASEDIR=$(pwd)/../

### CONFIGURACION
source $BASEDIR/conf/config.env
source $BASEDIR/conf/step.conf
source $BASEDIR/conf/exp.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source ${BASEDIR}/lib/encolar${QUEUESYS}.sh                     # Selecciona el metodo de encolado segun el systema QUEUESYS elegido

if [ $MACHINE == "FUGAKU" ] ;then
  source $TOPDIR/setup_spack.sh
  source $BASEDIR/conf/machine.conf
fi

##### FIN INICIALIZACION ######
if [ $STEP -ge 1 ] && [ $WPS_CYCLE -eq 0 ] ; then
   echo "Pert with WPS_CYCLE=0 is run only at STEP=0. And currently STEP=",$STEP        
   return 0
fi

cd $PERTDIR

PERTDIRRUN=$PERTDIR/00/   #We need to add 00 in order to be consistent with the paralelization lib.
rm -fr $PERTDIRRUN
mkdir -p $PERTDIRRUN
mkdir -p $PERTDIRRUN/log/scale_init
mkdir -p $PERTDIRRUN/log/pert
ln -sf $HISTDIR/const/topo $PERTDIRRUN/
ln -sf $HISTDIR/const/landuse  $PERTDIRRUN/

ln -s $PERTDIR/pert_init_bdy.exe $PERTDIRRUN/

BDY_INI_DATE_INT=$(date -d "$BDY_INI_DATE" +"%Y%m%d%H%M%S")
BDY_INI_DATE_SCALE=$(date -d "$BDY_INI_DATE" +"%Y%m%d-%H%M%S.000")
NPROC=$((SCALEPROC/SCALESKIP))
PPN=$((ICORE/SCALESKIP))

#Loop over the expected output files. If all the files are present then skip PERT step.
CDATE=$BDY_INI_DATE
FILE_COUNTER=0
petest=pe$(printf %06f $((NPROC-1)))
MYFILE=init_${BDY_INI_DATE_SCALE}.${petest}.nc
for IMEM in $(seq -w $MEM_INI $MEM_END) ; do
   MYPATH=$HISTDIR/init/${BDY_INI_DATE_INT}/${IMEM}_pert/
   if ! [[ -e $MYPATH/$MYFILE ]] ; then #My file do not exist
#      echo $MYPATH/$MYFILE " does not exist"
      FILE_COUNTER=$(( $FILE_COUNTER + 1 ))
   fi
done
MYFILE=boundary.${petest}.nc
for IMEM in $(seq -w $MEM_INI $MEM_END) ; do
   MYPATH=$HISTDIR/bdy/${BDY_INI_DATE_INT}/${IMEM}_pert/
   if ! [[ -e $MYPATH/$MYFILE ]] ; then #My file do not exist
#      echo $MYPATH/$MYFILE " does not exist"
      FILE_COUNTER=$(( $FILE_COUNTER + 1 ))
   fi
done

if [ $FILE_COUNTER -eq 0 ] ; then
   echo "We found all the required perturbed files for all members, skiping init and bdy perturbation"
   ERROR=0
else

   #Copy met em files to run pert metem
   ICP=1                                                                      #Counter for the number of cp commands.
   MAX_SIM_CP=$(( 10#$BDY_MEM_END - 10#$BDY_MEM_INI + 1 ))                  #Maximum number of simultaneous cp commands.

   #Creating directories
   echo "Making a copy of init and bdy files to create the expanded ensemble"
   for IMEM in $(seq -w $MEM_INI $MEM_END) ; do
      mkdir -p $HISTDIR/init/$BDY_INI_DATE_INT/${IMEM}_pert/
      mkdir -p $HISTDIR/bdy/$BDY_INI_DATE_INT/${IMEM}_pert/
   done

   #Copying the files
   for pe in $(seq -f %06g 0 $((NPROC-1))) ; do
      MYFILE=init_${BDY_INI_DATE_SCALE}.pe${pe}.nc
      for IMEM in $(seq -w $MEM_INI $MEM_END) ; do
         CDATE_WPS=$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%T" )
         cp $HISTDIR/init/$BDY_INI_DATE_INT/$(printf '%0'${#MEM_INI}'d' $ICP)/$MYFILE $HISTDIR/init/$BDY_INI_DATE_INT/${IMEM}_pert/ &
         ICP=$(( $ICP + 1 ))
         if [ $ICP -gt $MAX_SIM_CP ] ; then
            wait
            ICP=1
         fi
      done
   done
   wait
   ICP=1
   for pe in $(seq -f %06g 0 $((NPROC-1))) ; do
      MYFILE=boundary.pe${pe}.nc
      for IMEM in $(seq -w $MEM_INI $MEM_END) ; do
         CDATE_WPS=$(date -u -d "$CDATE UTC" +"%Y-%m-%d_%T" )
         cp $HISTDIR/bdy/$BDY_INI_DATE_INT/$(printf '%0'${#MEM_INI}'d' $ICP)/$MYFILE $HISTDIR/bdy/$BDY_INI_DATE_INT/${IMEM}_pert/ &
         ICP=$(( $ICP + 1 ))
         if [ $ICP -gt $MAX_SIM_CP ] ; then
            wait
            ICP=1
         fi
      done
   done
   wait

   echo "Generating the namelist"
   #Prepare the namelist.
   NAMELISTFILE=$PERTDIR/namelist.pert_init_bdy
   NBV=$(( 10#$BDY_MEM_END - 10#$BDY_MEM_INI + 1 ))
   NBV_TAR=$(( 10#$MEM_END     - 10#$MEM_INI     + 1 ))

   FILES_IN_RESTART="$HISTDIR/init/$BDY_INI_DATE_INT/<member>/init_${BDY_INI_DATE_SCALE}"
   FILES_OUT_RESTART="$HISTDIR/init/$BDY_INI_DATE_INT/<member>_pert/init_${BDY_INI_DATE_SCALE}"
   FILES_IN_BOUNDARY="$HISTDIR/bdy/$BDY_INI_DATE_INT/<member>/boundary"
   FILES_OUT_BOUNDARY="$HISTDIR/bdy/$BDY_INI_DATE_INT/<member>_pert/boundary"

   INI_SCALE="$(date -ud "$INI_DATE_FCST" +'%Y,%m,%d,%H,%M,%S' )"
   FCSTLEN=$(( $(date -ud "$END_DATE_FCST" +%s) - $(date -ud "$INI_DATE_FCST" +%s) ))
   E_WE_LOC=$((E_WE / NPROC_X))
   E_SN_LOC=$((E_SN / NPROC_Y))

   cp $NAMELISTDIR/namelist.pert_init_bdy $NAMELISTFILE

   #Edit the namelist from the template
   sed -i -e "s|__FILES_IN_BOUNDARY__|\'$FILES_IN_BOUNDARY\'|g"   $NAMELISTFILE
   sed -i -e "s|__FILES_OUT_BOUNDARY__|\'$FILES_OUT_BOUNDARY\'|g" $NAMELISTFILE
   sed -i -e "s|__FILES_IN_RESTART__|\'$FILES_IN_RESTART\'|g"     $NAMELISTFILE
   sed -i -e "s|__FILES_OUT_RESTART__|\'$FILES_OUT_RESTART\'|g"   $NAMELISTFILE
   sed -i -e "s|__NBV_TAR__|$NBV_TAR|g"                       $NAMELISTFILE
   sed -i -e "s|__METHOD__|$PERT_TYPE|g"                      $NAMELISTFILE
   sed -i -e "s|__SIGMA__|$PERT_AMP|g"                        $NAMELISTFILE

   sed -i -e "s|__NBV__|$NBV|g"                               $NAMELISTFILE
   sed -i -e "s|__PPN__|$PPN|g"                               $NAMELISTFILE
   sed -i -e "s|__MEM_NODES__|$SCALENODE|g"                   $NAMELISTFILE
   sed -i -e "s|__NPROC__|$NPROC|g"                           $NAMELISTFILE

   sed -i -e "s|__INI_SCALE__|$INI_SCALE|g"                   $NAMELISTFILE
   sed -i -e "s|__DT__|$DT|g"                                 $NAMELISTFILE
   sed -i -e "s|__LEN__|$FCSTLEN|g"                           $NAMELISTFILE
   sed -i -e "s|__E_VERT__|$((E_VERT-1))|g"                   $NAMELISTFILE
   sed -i -e "s|__E_WE_LOC__|$E_WE_LOC|g"                      $NAMELISTFILE
   sed -i -e "s|__E_SN_LOC__|$E_SN_LOC|g"                      $NAMELISTFILE
   sed -i -e "s|__NPROC_X__|$NPROC_X|g"          $NAMELISTFILE
   sed -i -e "s|__NPROC_Y__|$NPROC_Y|g"          $NAMELISTFILE
   sed -i -e "s|__DX__|$DX|g"                    $NAMELISTFILE
   sed -i -e "s|__DY__|$DY|g"                    $NAMELISTFILE

FZTEXT=$(cat $NAMELISTDIR/fz_$((E_VERT-1))lev.txt)
if [ $? -ne 0 ]; then
  echo "no vertical level file fz_$((E_VERT-1))lev.txt" ]
  return 1
fi
   sed -i -e "/!---FZ---/a $(echo $FZTEXT | sed -e 's/\./\\./g')" $NAMELISTFILE


#   cat $SCALEPPDIR/namelist.scale_init >> $NAMELISTFILE 
#
#   sed -i -e "|&PARAM_TIME|a TIME_DT=12.0D0,|"                         $NAMELISTFILE

   ln -sf $PERTDIR/pert_init_bdy.exe $PERTDIR/00/
   ln -sf $NAMELISTFILE              $PERTDIR/00/

#script de ejecucion
read -r -d '' QSCRIPTCMD << "EOF"
   source $BASEDIR/conf/step.conf
   source $BASEDIR/conf/exp.conf

   time $MPIEXE ./pert_init_bdy.exe ./namelist.pert_init_bdy ./log/pert/NOUT 2>&1
   ERROR=$(( $ERROR + $? ))
   cp log/pert/NOUT-000000 $LOGDIR/NOUT-000-Pert_$STEP

EOF

# Parametros de encolamiento
QPROC_NAME=PERT_$STEP
QNODE=$PERTNODE
QPROC=$PERTPROC
QSKIP=$PERTSKIP
QOMP=$PERTOMP
QWALLTIME=$PERTWALLTIME
QWORKPATH=$PERTDIR

echo "Sending pert_init_bdy script to the queue"
# Encolar
queue 00 00
check_proc 00 00
if [ $? -ne 0 ] ; then
   echo "Error: Some members do not finish OK"
   echo "Aborting this step"
   exit 1 
fi

fi

echo "Sucessfully finished perturbation of init and bdy files"


