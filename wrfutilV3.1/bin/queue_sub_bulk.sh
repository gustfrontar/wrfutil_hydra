memstart=$1 ;shift
memend=$1 ;shift
time_limit=$1 ;shift
COMMAND="$*"  #The script to be submitted to the queue.

BASEDIR=$(pwd)/../
source $BASEDIR/conf/machine.conf  #Load information about the requested number of nodes 
                                   #and the requested number of processors per node. 

#PJM  FUGAKU / FX100 / FX10 
if [ "$QUEUESYS" = "PJM" ] ; then
   if [ ! -z ${INTERACTIVE_JOB} ] && [ ${INTERACTIVE_JOB} -eq 1 ] ; then
      echo "bulk job not supported in interactive node"	   
      exit 1
   else
      echo "Submitting a bulk PJM JOB"
      source /vol0301/data/hp150019/u10335/test_wrf/wrfutil_hydra/set_path.sh   
      echo "export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" > ${COMMAND}.tmp
      echo "" >> ${COMMAND}.tmp
      cat ${COMMAND} >> ${COMMAND}.tmp
      mv ${COMMAND}.tmp $COMMAND
      ret=$(pjsub --bulk --sparam "${memstart}-${memend}" -g ${FUGAKU_GROUP} -L "node=${INODE}" --mpi "max-proc-per-node=${ICORE}" -x PJM_LLIO_GFSCACHE=${LLIO_VOL} -L "elapse=${time_limit}" -j -s --no-stging ${COMMAND})
      if [ "$(echo $ret | cut -d " " -f 7)" == "submitted." ] ;then 
        jobid=$(echo $ret | cut -d " " -f 6)
        echo $jobid > $BASEDIR/WPS/running
        echo "a bulk PJM job $jobid submitted" 
      else
        echo "error submitting a bulk PJM job. " 
        echo $ret
        exit 1
      fi
   fi

fi 

