#!/bin/bash 
#######################################################
# This script is only for the FUGAKU computer.
# This script can be sent directly to the queue.
# sbatch main_Pert.sh 
# Note: This script is intensive in memory use. 
# memory scales with the number of ensemble members in
# the original ensemble and with the size of the domain
# in the experiment. 
# In the future we can translate this into fortran and 
# use MPI for more efficient memory administration and
# faster performance. 
#######################################################

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -p mem2      # Specifying a queue
#SBATCH -n 10        # Specifying the number of CPUs
#SBATCH -N 1         # Specifying the number of nodes (this script can only use one node)
#SBATCH --mem 100G   # Specifying memory usage [MB]
export OMP_NUM_THREADS=10

#Get the working directory
if [ ! -z ${PBS_O_WORKDIR}    ]; then cd ${PBS_O_WORKDIR}   ;fi
if [ ! -z ${PJM_O_WORKDIR}    ]; then cd ${PJM_O_WORKDIR}   ;fi
if [ ! -z ${SLURM_SUBMIT_DIR} ]; then cd ${SLURM_SUBMIT_DIR};fi

#Load experiment configuration
BASEDIR=$(pwd)/../
source $BASEDIR/lib/errores.env
source $BASEDIR/conf/config.env
source $BASEDIR/conf/machine.conf 
source $BASEDIR/conf/${EXPTYPE}.conf

##### FIN INICIALIZACION ######
cd $WPSDIR
cp $BASEDIR/bin/python/main_perturb_met_em_par.py $WPSDIR/main_perturb_met_em.py
cp $BASEDIR/bin/python/module_pert_met_em.py  $WPSDIR/module_pert_met_em.py

FECHA_INI_PASO=$(date -u -d "$FECHA_INI UTC +$(($WPS_INI_FREQ*$PASO)) seconds" +"%Y-%m-%d %T")
FECHA_INI_BDY=$(date_floor "$FECHA_INI_PASO" $INTERVALO_BDY )
FECHA_INI_BDY=$(date -u -d "$FECHA_INI_BDY" +"%Y%m%d%H%M%S")

sed -i -e "s|__BASE_PATH_ORI__|'$HISTDIR/WPS/met_em_ori/${FECHA_INI_BDY}/'|g"   $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__BASE_PATH_OUT__|'$HISTDIR/WPS/met_em/${FECHA_INI_BDY}'|g"        $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__TARGET_ENSEMBLE_SIZE__|$TARGET_ENSEMBLE_SIZE|g"                  $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__ACTUAL_ENSEMBLE_SIZE__|$ACTUAL_ENSEMBLE_SIZE|g"                  $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__PERT_AMP__|$PERT_AMP|g"                                          $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__PERT_TYPE__|'$PERT_TYPE'|g"                                      $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__VAR_LIST__|$VAR_LIST|g"                                          $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__THREADS_NUM__|$PERTTHREADS|g"                                    $WPSDIR/main_perturb_met_em.py
sed -i -e "s|__NETCDF4__|$NETCDF4|g"                                            $WPSDIR/main_perturb_met_em.py

ulimit -s unlimited
cd $WPSDIR
time $PYTHON -u ./main_perturb_met_em.py > $LOGDIR/log_pert_${PASO}.log
EC=$?
[[ $EC -ne 0 ]] && dispararError 9 "main_perturb_met_em.py"

echo "Successfully finish perturbing met_em files"
echo "We finished @"$(date )
exit 0


