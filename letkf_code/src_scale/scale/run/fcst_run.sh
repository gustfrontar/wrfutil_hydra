#!/bin/bash
#===============================================================================
#
#  Wrap fcst.sh in FUGAKU/Linux and run it.
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst_run.sh [..]
#
#===============================================================================

cd "$(dirname "$0")"
myname="$(basename "$0")"
job='fcst'
if [ "$GROUP" == "fugaku" ] ; then
  echo 'specify group name $GROUP in which you want to submit the job'
  exit 1
fi

#===============================================================================
# Configuration

. ./config.main || exit $?
. ./config.${job} || exit $?

. src/func_datetime.sh || exit $?
. src/func_distribute.sh || exit $?
. src/func_util.sh || exit $?

. src/func_common_static.sh || exit $?
. src/func_${job}_static.sh || exit $?

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Start $myname $@"

setting "$@" || exit $?


echo
print_setting || exit $?
echo

#===============================================================================
# Create and clean the temporary directory

echo "[$(datetime_now)] Create and clean the temporary directory"
 
#if [ -e "${TMP}" ]; then
#  echo "[Error] $0: \$TMP will be completely removed." >&2
#  echo "        \$TMP = '$TMP'" >&2
#  exit 1
#fi
safe_init_tmpdir $TMP || exit $?

#===============================================================================
# Determine the distibution schemes

echo "[$(datetime_now)] Determine the distibution schemes"

safe_init_tmpdir $NODEFILE_DIR || exit $?
distribute_fcst "$MEMBERS" $CYCLE "(0)" $NODEFILE_DIR || exit $?

if ((CYCLE == 0)); then
  CYCLE=$cycle_auto
fi

#===============================================================================
# Determine the staging list

echo "[$(datetime_now)] Determine the staging list"

cat $SCRP_DIR/config.main | \
    sed -e "/\(^DIR=\| DIR=\)/c DIR=\"$DIR\"" \
    > $TMP/config.main

echo "SCRP_DIR=\"\$TMPROOT\"" >> $TMP/config.main
echo "NODEFILE_DIR=\"\$TMPROOT/node\"" >> $TMPS/config.main
echo "RUN_LEVEL=4" >> $TMP/config.main

echo "PARENT_REF_TIME=$PARENT_REF_TIME" >> $TMP/config.main

safe_init_tmpdir $STAGING_DIR || exit $?
staging_list_static || exit $?
config_file_list $TMPS/config || exit $?

NNODES_USE=$(( fmember * ( SCALE_NP / PPN ) ))
echo "NNODES=$NNODES_USE" >> $TMP/config.main

#-------------------------------------------------------------------------------
# Add shell scripts and node distribution files into the staging list

cp ${SCRP_DIR}/config.rc $TMP/config.rc
cp ${SCRP_DIR}/config.${job} $TMP/config.${job}
cp ${SCRP_DIR}/src/${job}.sh $TMP/${job}.sh
cp -r ${SCRP_DIR}/src $TMP/src

#===============================================================================
# Stage in

echo "[$(datetime_now)] Initialization (stage in)"

stage_in server || exit $?

#===============================================================================
# Creat a job script

jobscrp="$TMP/${job}_job.sh"

echo "[$(datetime_now)] Create a job script '$jobscrp'"

# FUGAKU
if [ "$PRESET" = 'FUGAKU' ]; then

  if (( NNODES_USE > 384 )) ; then
    RSCGRP="large"
  else
    RSCGRP="small"
  fi

  TPROC=$((NNODES_USE*PPN))
  if [ "$RSCGRP_LARGE" == "T" ] ;then
    RSCGRP="large"
    TPROC=$((385*PPN))
  fi  

  VOLUMES="/"$(readlink /data/${GROUP} | cut -d "/" -f 2)
  if [ "$VOLUMES" != "/vol0004" ] ;then
    VOLUMES="${VOLUMES}:/vol0004" # spack
  fi

cat > $jobscrp << EOF
#!/bin/sh 
#
#PJM -g ${GROUP} 
#PJM -x PJM_LLIO_GFSCACHE=${VOLUMES}
#PJM -L "freq=2200,eco_state=2"
#PJM -L "rscgrp=${RSCGRP}"
#PJM -L "node=$(((TPROC+PPN-1)/PPN))"
#PJM -L "elapse=${TIME_LIMIT}"
#PJM --mpi "max-proc-per-node=${PPN}"
#PJM -j
#PJM -s
#
#
export PARALLEL=${THREADS}
export OMP_NUM_THREADS=\${PARALLEL}
export FORT90L=-Wl,-T
export PLE_MPI_STD_EMPTYFILE=off
export OMP_WAIT_POLICY=active
export FLIB_BARRIER=HARD

EOF

  if (( USE_LLIO_BIN == 1 )); then
    for i in $(seq $nsteps) ; do
      echo "llio_transfer ${stepexecbin[$i]}" >> $jobscrp 
    done
    echo "" >> $jobscrp
  fi

  if (( USE_LLIO_DAT == 1 )); then
    echo "/home/system/tool/dir_transfer -l ./ ${TMPROOT}/dat" >> $jobscrp
    echo "" >> $jobscrp
  fi

  if (( USE_SPACK > 0 )); then

cat << EOF >>  $jobscrp 
SPACK_FJVER=${SPACK_FJVER}
. /vol0004/apps/oss/spack/share/spack/setup-env.sh
spack load netcdf-c%fj@\${SPACK_FJVER}
spack load netcdf-fortran%fj@\${SPACK_FJVER}
spack load parallel-netcdf%fj@\${SPACK_FJVER}

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:\$LD_LIBRARY_PATH

EOF

  else

    if [ -z "$SCALE_NETCDF_C" ] || [ -z "$SCALE_NETCDF_F" ] || [ -z "$SCALE_PNETCDF" ] || [ -z "$SCALE_HDF" ] ; then
      echo "[Error] Export SCALE environmental parameters (e.g., SCALE_NETCDF_C)"
      exit 1
    fi

cat << EOF >>  $jobscrp 
export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SCALE_NETCDF_C}/lib:${SCALE_NETCDF_F}/lib:${SCALE_PNETCDF}/lib:${SCALE_HDF}/lib:\$LD_LIBRARY_PATH

EOF

  fi

cat << EOF >>  $jobscrp 
./${job}.sh "$STIME" "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP" "$CONF_MODE" || exit \$?
EOF

  if (( USE_LLIO_BIN == 1 )); then
    for i in $(seq $nsteps) ; do
      echo "llio_transfer --purge ${stepexecbin[$i]}" >> $jobscrp 
    done
  fi

  if (( USE_LLIO_DAT == 1 )); then
    echo "/home/system/tool/dir_transfer -p -l ./ ${TMPROOT}/dat" >> $jobscrp
    echo "" >> $jobscrp
  fi

  echo "[$(datetime_now)] Run ${job} job on PJM"
  echo

  job_submit_PJM $jobscrp
  echo
  
  job_end_check_PJM $jobid
  res=$?

# qsub
elif [ "$PRESET" = 'Linux_torque' ]; then

if [ $NNODES_USE -lt 4 ] ; then
  RSCGRP=s
elif [ $NNODES_USE -le 16 ] ; then
  RSCGRP=m
elif [ $NNODES_USE -le 24 ] ; then
  RSCGRP=l
else
  echo "too many nodes required. " $NNODES_USE " > 24"
  exit 1
fi

cat > $jobscrp << EOF
#!/bin/sh
#PBS -N $job
#PBS -q ${RSCGRP}
#PBS -l nodes=${NNODES_USE}:ppn=${PPN}
#PBS -l walltime=${TIME_LIMIT}
#
#

cd \${PBS_O_WORKDIR}

export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

EOF

if [ "$SCALE_SYS" == "Linux64-gnu-ompi" ] ; then

cat >> $jobscrp << EOF

source /etc/profile.d/modules.sh
module unload mpt/2.12
module unload intelcompiler/16.0.1.150
module unload intelmpi/5.1.2.150
module unload hdf5/1.8.16-intel
module unload netcdf4/4.3.3.1-intel
module unload netcdf4/fortran-4.4.2-intel
module load gcc/4.7.2
module load openmpi/2.0.4-gcc
module load hdf5/1.8.16
module load netcdf4/4.3.3.1
module load netcdf4/fortran-4.4.2
module load lapack/3.6.0

export OMP_NUM_THREADS=1
export KMP_AFFINITY=compact

EOF

else

cat >> $jobscrp << EOF

source /etc/profile.d/modules.sh 
module unload mpt/2.12
module load intelmpi/5.1.2.150

export OMP_NUM_THREADS=${THREADS}
export KMP_AFFINITY=compact

export LD_LIBRARY_PATH="/home/seiya/lib:$LD_LIBRARY_PATH"
EOF


fi 

cat >> $jobscrp << EOF

ulimit -s unlimited

./${job}.sh "$STIME" "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP" "$CONF_MODE" &> run_progress || exit \$?
EOF

  echo "[$(datetime_now)] Run ${job} job on PJM"
  echo

  job_submit_torque $jobscrp
  echo
  
  job_end_check_torque $jobid
  res=$?

# SMN Slurm
elif [ "$PRESET" = 'Linux_slurm' ]; then

cat >> $jobscrp <<EOF
#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name="scale-letkf_${job}"
#SBATCH --mail-type=NONE
#SBATCH -t ${TIME_LIMIT}
#SBATCH --exclusive
#SBATCH -o scale-letkf_${job}
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=${PPN}
#SBATCH -n $((PPN*NNODES_USE))
#SBATCH -N ${NNODES_USE}
ulimit -s unlimited

./${job}.sh "$STIME" "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP" "$CONF_MODE" &> run_progress || exit $?
EOF

  echo "[$(datetime_now)] Run ${job} job on PJM"
  echo

  job_submit_slurm $jobscrp
  echo
  
  job_end_check_slurm $jobid
  res=$?

# direct
elif [ "$PRESET" = 'Linux' ]; then

  echo "[$(datetime_now)] Run ${job} job on PJM"
  echo

  cd $TMPS

./${job}.sh "$STIME" "$ETIME" "$MEMBERS" "$CYCLE" "$CYCLE_SKIP" "$IF_VERF" "$IF_EFSO" "$ISTEP" "$FSTEP" "$CONF_MODE" &> run_progress || exit $?

else
  echo "PRESET '$PRESET' is not supported."
  exit 1
fi

#===============================================================================
# Stage out

echo "[$(datetime_now)] Finalization (stage out)"

stage_out server || exit $?

#===============================================================================
# Finalization

echo "[$(datetime_now)] Finalization"
echo

backup_exp_setting $job $TMP $jobid ${job}_job.sh 'o e'

archive_log

if ((CLEAR_TMP == 1)); then
  safe_rm_tmpdir $TMP
fi

#===============================================================================

echo "[$(datetime_now)] Finish $myname $@"

exit $res
