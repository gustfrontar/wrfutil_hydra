#!/bin/bash
. /data/opt/Env/Env_intel
source /etc/profile

module load common_intel/python/3.4.3

export PYTHONPATH=./smooth_2d/:./common_functions/:$PYTHONPATH

python3 $1
#export PYTHONPATH=/data/gylien/scripts/python3  
#python3 $1

#Warning. This scripts do not work properly in blanton due to a bug related to the module that is loaded.

#After runing this script run python3


