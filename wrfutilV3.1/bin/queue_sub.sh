COMMAND=$1  #The script to be submitted to the queue.

BASEDIR=$(pwd)/../
source $BASEDIR/conf/machine.conf  #Load information about the requested number of nodes 
                                   #and the requested number of processors per node. 

#PBSSSH - the main script is submitted to the queue.
if [ "$QUEUESYS" = "PBSSSH" ] ; then
   qsub -l nodes=${INODE}:ppn=${ICORE} -l walltime=${TOTAL_TIME_LIMIT} -q ${QUEUE} ${COMMAND}
fi

#SLURM
if [ "$QUEUESYS" = "SLURMSSH" ] ; then
   sbatch -N $INODE -c $ICORE -e error.log -p mem1 ${COMMAND}
fi

#PJM  FUGAKU / FX100 / FX10 
if [ "$QUEUESYS" = "PJM" ] ; then
   pjsub --interact -g ${FUGAKU_GROUP} -L "node=${INODE}" --mpi "max-proc-per-node=${ICORE}" -x PJM_LLIO_GFSCACHE=${LLIO_VOL} -L "elapse=${TOTAL_TIME_LIMIT}" -j -s --no-stging --sparam wait-time=1200  ${COMMAND}
fi 



