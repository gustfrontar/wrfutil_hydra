verification_against_obs () {

#This function uses the offline observation operator to 
#generate verfication files for the forecast.

local CDATEL=$1    #Forecast start date
local EDATEL=$2    #Forecast end date  (if analysis or gues then should be equal)
local WORKDIR=$TMPDIR/verification/


rm -fr $TMPDIR/*.txt
rm -fr $TMPDIR/*.grd

#Prepare namelist for observation operator
cp $NAMELISTOBSOPE $TMPDIR/verification/obsope.namelist

if [ $RUN_ONLY_MEAN -eq 1 ] ; then
  #In this case we force the ensemble size to be 1.
  sed -i 's/@NBV@/'1'/g' $TMPDIR/verification/obsope.namelist
fi

edit_namelist_obsope $TMPDIR/verification/obsope.namelist

if [ $RUN_ONLY_MEAN -eq 1 ] ; then
   INIMEMBER=$MEANMEMBER
   ENDMEMBER=$MEANMEMBER
else
   INIMEMBER=1
   ENDMEMBER=$MEMBER
fi


if [ $FORECAST -eq 1  ] ; then

my_domain=1
while [ $my_domain -le $MAX_DOM ] ; do


 if [ $my_domain -lt 10 ] ; then
   my_domain=0$my_domain
 fi

 local output_dir=${RESULTDIRG}/obsver_d$my_domain
 mkdir -p $output_dir

 local my_date=$CDATEL
 local S=1
 while [ $my_date -le $EDATEL ] ; do

  local my_file_name=`wrfout_file_name ${my_date} $my_domain `
  local SLOT=`add_zeros $S 2 `

  if [ $RUN_ONLY_MEAN -eq 1 ] ; then
    local MEM=`ens_member $MEANMEMBER `
    local MEM1=`ens_member 1 `
    ln -sf ${RESULTDIRG}/${MEM}/${my_file_name}         ${WORKDIR}/gs${SLOT}${MEM1}
  else
    local M=$INIMEMBER
    while [ $M -le $ENDMEMBER ] ; do
      local MEM=`ens_member $M `
      ln -sf ${RESULTDIRG}/${MEM}/${my_file_name}       ${WORKDIR}/gs${SLOT}${MEM}              
      M=`expr $M + 1 `
    done
  fi

  #TODO: Add radar observations here.
  ln -sf ${OBSDIR}/${my_date}.dat                       ${WORKDIR}/obs${SLOT}.dat      

  my_date=`date_edit2 $my_date $WINDOW_FREC `
  S=`expr $S + 1 `
 done

  #Create run script
  echo "#!/bin/bash                                                                         " >  ${WORKDIR}/tmp.sh
  echo "export LD_LIBRARY_PATH=$RUNTIMELIBS:\$LD_LIBRARY_PATH                               " >> ${WORKDIR}/tmp.sh
  if [ $SYSTEM -eq  1 ] ; then
     echo " ulimit -s unlimited                                                             " >> ${WORKDIR}/tmp.sh
  fi
  echo "cd ${WORKDIR}/                                                                      " >> ${WORKDIR}/tmp.sh
  echo "$MPIBIN -np ${MAX_RUNNING} ./obsope.exe   > obsop.log                               " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/*.grd       ${output_dir}                                             " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/*.txt       ${output_dir}                                             " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/obsope????? ${output_dir}                                             " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/OBSO-???    ${output_dir}                                             " >> ${WORKDIR}/tmp.sh
  chmod 755 ${WORKDIR}/tmp.sh

  ssh $PPSSERVER ${WORKDIR}/tmp.sh

  my_domain=`expr $my_domain + 1 `
done


fi #[Script for forecast verification (supports muliple nests)]

if [ $ANALYSIS -eq 1  ] ; then


 #Run for the first guess
 local my_date=$CDATEL
 local S=1

  local SLOT=`add_zeros $S 2`

  local M=$INIMEMBER
  while [ $M -le $ENDMEMBER ] ; do
   local MEM=`ens_member $M `
   ln -sf ${RESULTDIRG}/gues${MEM}                  ${WORKDIR}/gs${SLOT}${MEM}
   M=`expr $M + 1 `
  done
  ln -sf ${OBSDIR}/${my_date}.dat                   ${WORKDIR}/obs${SLOT}.dat

  #Create run script
  echo "#!/bin/bash                                                                         " >  ${WORKDIR}/tmp.sh
  echo "export LD_LIBRARY_PATH=$RUNTIMELIBS:\$LD_LIBRARY_PATH                               " >> ${WORKDIR}/tmp.sh
  if [ $SYSTEM -eq  1 ] ; then
     echo " ulimit -s unlimited                                                             " >> ${WORKDIR}/tmp.sh
  fi
  echo "cd ${WORKDIR}/                                                                      " >> ${WORKDIR}/tmp.sh
  echo "$MPIBIN -np ${MAX_RUNNING} ./obsope.exe           > obsop.log                       " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/obsmean01_01.grd    ${RESULTDIRG}/obsmeangues.grd                     " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/obsmerr01_01.grd    ${RESULTDIRG}/obsmerrgues.grd                     " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/obssprd01_01.grd    ${RESULTDIRG}/obssprdgues.grd                     " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/obsonum01_01.grd    ${RESULTDIRG}/obsnumgues.grd                      " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/oarea01_01.txt      ${RESULTDIRG}/oareagues.txt                       " >> ${WORKDIR}/tmp.sh
  
  echo "mv ${WORKDIR}/obsope????? ${RESULTDIRG}                                             " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/OBSO-???    ${RESULTDIRG}                                             " >> ${WORKDIR}/tmp.sh
  chmod 755 ${WORKDIR}/tmp.sh

  ssh $PPSSERVER ${WORKDIR}/tmp.sh

 #Run for the analysis
 local my_date=$CDATEL
 local S=1

  local SLOT=`add_zeros $S 2`

  local M=$INIMEMBER
  while [ $M -le $ENDMEMBER ] ; do
   local MEM=`ens_member $M `
   ln -sf ${RESULTDIRA}/anal${MEM}                  ${WORKDIR}/gs${SLOT}${MEM}
   M=`expr $M + 1 `
  done
  ln -sf ${OBSDIR}/${my_date}.dat                   ${WORKDIR}/obs${SLOT}.dat

  #Create run script
  echo "#!/bin/bash                                                                         " >  ${WORKDIR}/tmp.sh
  echo "export LD_LIBRARY_PATH=$RUNTIMELIBS:\$LD_LIBRARY_PATH                               " >> ${WORKDIR}/tmp.sh
  if [ $SYSTEM -eq  1 ] ; then
     echo " ulimit -s unlimited                                                             " >> ${WORKDIR}/tmp.sh
  fi
  echo "cd ${WORKDIR}/                                                                      " >> ${WORKDIR}/tmp.sh
  echo "$MPIBIN -np ${MAX_RUNNING} ./obsope.exe         > obsop.log                         " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/obsmean01_01.grd    ${RESULTDIRA}/obsmeananal.grd                     " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/obsmerr01_01.grd    ${RESULTDIRA}/obsmerranal.grd                     " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/obssprd01_01.grd    ${RESULTDIRA}/obssprdanal.grd                     " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/obsonum01_01.grd    ${RESULTDIRA}/obsnumanal.grd                      " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/oarea01_01.txt    ${RESULTDIRA}/oareaanal.txt                        " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/obsope????? ${RESULTDIRA}                                             " >> ${WORKDIR}/tmp.sh
  echo "mv ${WORKDIR}/OBSO-???    ${RESULTDIRA}                                             " >> ${WORKDIR}/tmp.sh
  chmod 755 ${WORKDIR}/tmp.sh

  ssh $PPSSERVER ${WORKDIR}/tmp.sh

fi #[Script for analysis and gues verification ]

