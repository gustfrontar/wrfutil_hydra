#!/bin/bash

ulimit -s unlimited

. $HOME/set_env_anaconda3.sh

export OMP_NUM_THREADS=8

export COMMON_PYTHON=$HOME/share/LETKF_WRF/common_python/

export PYTHONPATH=$COMMON_PYTHON/common_functions/:$COMMON_PYTHON/scale_modules/:$COMMON_PYTHON/smooth_2d/:$PYTHONPATH

python3 $1


