COMMAND=$1  #The script to be submitted to the queue.

BASEDIR=$(pwd)/../
source $BASEDIR/conf/machine.conf  #Load information about the requested number of nodes 
                              #and the requested number of processors per node. 

#PBS and PBSpro
if [ "$QUEUESYS" = "PBS" ] ; then

   qsub -l nodes=${INODE}:ppn=${ICORE} -q ${QUEUE} ${COMMAND}

fi

#SLURM
if [ "$QUEUESYS" = "SLURM" ] ; then
   #NOTA: I could not make this run in batch mode in Fugaku PPS. 
   #It only runs in interactive mode (with srun) but it hangs at
   #the mpirun call in the non interactive mode. 
   #Is this something related to the use of MPI_SPAWN?
   #To do this less interactive this script can be called with nohup. 
   srun --export=ALL -c $ICORE --mem 30G -e error.log -p mem1 ${COMMAND}
   #sbatch -N 1 --cpus-per-task 4 -t 00:10:00 -p mem1  ${COMMAND}

fi

#PJM  FUGAKU / FX100 / FX10 
if [ "$QUEUESYS" = "PJM" ] ; then

   pjsub --interact -g ${FUGAKU_GROUP} -L "node=${INODE}" --mpi "max-proc-per-node=${ICORE}" -x PJM_LLIO_GFSCACHE=${LLIO_VOL} -L "elapse=${TOTAL_TIME_LIMIT}" -j -s --no-stging --sparam wait-time=1200  ${COMMAND}

fi 



