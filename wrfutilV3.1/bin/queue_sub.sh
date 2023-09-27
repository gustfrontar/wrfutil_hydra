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

   sbatch -N $INODE -n $ICORE -p mem1 ${COMMAND}

fi






