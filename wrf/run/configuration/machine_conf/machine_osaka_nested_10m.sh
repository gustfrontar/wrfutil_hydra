#MACHINE CONFIGURATION
GROUP=hp150019
SYSTEM="0"          #0 - K computer , 1 - qsub cluster
PROC_PER_NODE=8     #K computer number of procs per node.
NODES_PER_MEMBER=20 #Number of nodes per ensemble member.
PPSSERVER=pps1      #Hostname of pps server (for perturbation generation and post processing)
MAX_RUNNING=10      #Maximum number of simultaneous processes running in PPS servers.
ELAPSE="00:45:00"   #MAXIMUM ELAPSE TIME (MODIFY ACCORDING TO THE SIZE OF THE DOMAIN AND THE RESOLUTION)
MAX_BACKGROUND_JOBS=128
LD_LIBRARY_PATH_ADD="/opt/intel/lib/intel64"
PATH_ADD="/home/ra000015/a03094/libintel/bin/:/home/ra000015/a03094/grads-2.0.1.oga.1/Contents/:/opt/intel/bin/"

TOTAL_NODES_FORECAST=400
TOTAL_NODES_LETKF=40

#These options control job split (in case of big jobs)
#Job split is performed authomaticaly if the number of ensemble members is
#larger than MAX_MEMBER_PER_JOB.
MAX_MEMBER_PER_JOB=2    #Number of members included on each job.
MAX_SUBMITT_JOB=10      #Maximum number of jobs that can be submitted to the queue.

