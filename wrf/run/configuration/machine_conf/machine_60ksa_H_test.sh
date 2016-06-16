#MACHINE CONFIGURATION
SYSTEM="1"  # 0 - K computer , 1 - qsub cluster
PROC_PER_NODE=24    #K computer number of procs per node.
PROC_PER_MEMBER=8   #Number  of procs per ensemble members (torque)
NODES_PER_MEMBER=1  #Number of nodes per ensemble member (K)
PPSSERVER=tormenta  #Hostname of pps server (for perturbation generation and post processing)
MAX_RUNNING=2       #Maximum number of simultaneous processes running in PPS servers.
ELAPSE="00:10:00"   #MAXIMUM ELAPSE TIME (MODIFY ACCORDING TO THE SIZE OF THE DOMAIN AND THE RESOLUTION)
MAX_BACKGROUND_JOBS=128
LD_LIBRARY_PATH_ADD=""
PATH_ADD=""

TOTAL_NODES_FORECAST=4
TOTAL_NODES_LETKF=1

#These options control job split (in case of big jobs)
#Job split is performed authomaticaly if the number of ensemble members is
#larger than MAX_MEMBER_PER_JOB.
MAX_MEMBER_PER_JOB=1024    #Number of members included on each job.
MAX_SUBMITT_JOB=10         #Maximum number of jobs that can be submitted to the queue.

