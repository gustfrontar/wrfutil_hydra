#MACHINE CONFIGURATION
GROUP=""
SYSTEM="1"  # 0 - K computer , 1 - qsub cluster
PROC_PER_NODE=40    #K computer number of procs per node.
PROC_PER_MEMBER=10  #Number  of procs per ensemble members (torque)
NODES_PER_MEMBER=1  #Number of nodes per ensemble member.
PPSSERVER=hakushu   #Hostname of pps server (for perturbation generation and post processing)
MAX_RUNNING=5       #Maximum number of simultaneous processes running in PPS servers.
ELAPSE="00:10:00"   #MAXIMUM ELAPSE TIME (MODIFY ACCORDING TO THE SIZE OF THE DOMAIN AND THE RESOLUTION)
MAX_BACKGROUND_JOBS=128
LD_LIBRARY_PATH_ADD="/home/jruiz/mpich_intel/lib/:/opt/intel/lib/intel64:/apps/SLES11/opt/netcdf-fortran/4.4.1_intel/lib://apps/SLES11/opt/hdf5/1.8.14_intel/lib/"
PATH_ADD="/home/jruiz/mpich_intel/bin/"

TOTAL_NODES_FORECAST=2
TOTAL_NODES_LETKF=2

#These options control job split (in case of big jobs)
#Job split is performed authomaticaly if the number of ensemble members is
#larger than MAX_MEMBER_PER_JOB.
MAX_MEMBER_PER_JOB=60      #Number of members included on each job.
MAX_SUBMITT_JOB=1          #Maximum number of jobs that can be submitted to the queue.

