COMMAND=$1  #The script to be submitted to the queue.

BASEDIR=$(pwd)/../
source $BASEDIR/conf/machine.conf  #Load information about the requested number of nodes 
                                   #and the requested number of processors per node. 

#PBSSSH - the main script is submitted to the queue.
if [ "$QUEUESYS" = "PBSSSH" ] ; then
   if [ ! -z ${INTERACTIVE_JOB} ] && [ ${INTERACTIVE_JOB} -eq 1 ] ; then 
      echo "Submitting an interactive PBS JOB"
      qsub -I -l nodes=${INODE}:ppn=${ICORE} -l walltime=${TOTAL_TIME_LIMIT} -q ${QUEUE} ${COMMAND}   #Submit as an interactive JOB
   else 
      echo "Submitting a regular PBS JOB"
      qsub -l nodes=${INODE}:ppn=${ICORE} -l walltime=${TOTAL_TIME_LIMIT} -q ${QUEUE} ${COMMAND}      #Submit as a regular JOB
   fi 
fi

#SLURM
if [ "$QUEUESYS" = "SLURMSSH" ] ; then
   if [ ! -z ${INTERACTIVE_JOB} ] && [ ${INTERACTIVE_JOB} -eq 1 ] ; then
      echo "Submitting an interactive SLURM JOB"
      srun --time=${TOTAL_TIME_LIMIT} -N $INODE -c $ICORE -e error.log -p ${QUEUE} ${COMMAND}
   else
      echo "Submitting an regular SLURM JOB"
      sbatch --time=${TOTAL_TIME_LIMIT} -N $INODE -c $ICORE -e error.log -p ${QUEUE} ${COMMAND}
   fi
fi

#PJM  FUGAKU / FX100 / FX10 
if [ "$QUEUESYS" = "PJM" ] ; then
   if [ ${INODE} -gt 384 ] ; then
      queue="large"
   else 
      queue="small"
   fi
   if [ ! -z ${INTERACTIVE_JOB} ] && [ ${INTERACTIVE_JOB} -eq 1 ] ; then
      echo "Submitting an interactive PJM JOB"	   
      pjsub --interact -g ${FUGAKU_GROUP} -L "node=${INODE}" --mpi "max-proc-per-node=${ICORE}" -x PJM_LLIO_GFSCACHE=${LLIO_VOL} --llio "sharedtmp-size=${SHAREDCACHESIZE}" -L "elapse=${TOTAL_TIME_LIMIT}" -j -s --sparam wait-time=1200  ${COMMAND}
   else
      echo "Submitting a regular PJM JOB"
      pjsub -g ${FUGAKU_GROUP} -L "node=${INODE}" -L "rscgrp=${queue}" --mpi "max-proc-per-node=${ICORE}" -x PJM_LLIO_GFSCACHE=${LLIO_VOL} --llio "sharedtmp-size=${SHAREDCACHESIZE}" -L "elapse=${TOTAL_TIME_LIMIT}" -j -s ${COMMAND}
   fi

fi 

#STAND ALONE SERVER WITH NO QUEUE SYSTEM
if [ "$QUEUESYS" = "SINGLENODE" ] ; then
   echo "Submitting an interactive JOB on a SINGLENODE server"
      nohup ./${COMMAND} > ${COMMAND}.log &
fi



