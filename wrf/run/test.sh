#PBS -l nodes=2:ppn=24
#PBS -S /bin/bash
#PBS -l walltime=100:00:00

#=======================================================================
# This script runs multiple data assimilation cycles.
# The driver script sets the TMP directory with all the required data
# and this scripts executes multiple data assimilation cycle and the output
# goes to the TMP directory.
#=======================================================================
#set -x
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
#Get root dir
CDIR=`pwd`

cd /home/jruiz/share/LETKF_WRF/wrf/run

source util.sh

cat $PBS_NODEFILE
