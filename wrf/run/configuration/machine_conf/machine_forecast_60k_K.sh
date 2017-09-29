#MACHINE CONFIGURATION
GROUP=hp150019
SYSTEM="0"  # 0 - K computer , 1 - qsub cluster
PROC_PER_NODE=8     #K computer number of procs per node.
NODES_PER_MEMBER=1  #Number of nodes per ensemble member.
PPSSERVER=pps1      #Hostname of pps server (for perturbation generation and post processing)
MAX_RUNNING=1       #Maximum number of simultaneous processes running in PPS servers.
ELAPSE="00:30:00"   #MAXIMUM ELAPSE TIME (MODIFY ACCORDING TO THE SIZE OF THE DOMAIN AND THE RESOLUTION)
MAX_BACKGROUND_JOBS=128

TOTAL_NODES_FORECAST=20
TOTAL_NODES_LETKF=20

#These options control job split (in case of big jobs)
#Job split is performed authomaticaly if the number of ensemble members is
#larger than MAX_MEMBER_PER_JOB.
MAX_MEMBER_PER_JOB=1024    #Number of members included on each job.
MAX_SUBMITT_JOB=10         #Maximum number of jobs that can be submitted to the queue.

