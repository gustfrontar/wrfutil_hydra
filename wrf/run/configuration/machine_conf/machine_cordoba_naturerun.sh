#MACHINE CONFIGURATION
SYSTEM="1"          # 0 - K computer , 1 - qsub cluster
PROC_PER_NODE=40    #Hakushu has 40, Hydra has 24
PROC_PER_MEMBER=40  #Number of nodes per ensemble member.
PPSSERVER=hakushu   #Hostname of pps server (for perturbation generation and post processing)
MAX_RUNNING=5       #Maximum number of simultaneous processes running in PPS servers.
ELAPSE=             #MAXIMUM ELAPSE TIME (MODIFY ACCORDING TO THE SIZE OF THE DOMAIN AND THE RESOLUTION)
MAX_BACKGROUND_JOBS=10000
LD_LIBRARY_PATH_ADD="/opt/netcdf-fortran/4.2/lib/"
PATH_ADD=""

TOTAL_NODES_FORECAST=1
TOTAL_NODES_LETKF=1

#These options control job split (in case of big jobs)
#Job split is performed authomaticaly if the number of ensemble members is
#larger than MAX_MEMBER_PER_JOB.
MAX_MEMBER_PER_JOB=30      #Number of members included on each job.
MAX_SUBMITT_JOB=1          #Maximum number of jobs that can be submitted to the queue.

