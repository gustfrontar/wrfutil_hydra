#!/bin/sh 
#
#PJM -g hp150019 
#PJM -x PJM_LLIO_GFSCACHE=/vol0003:/vol0004
#PJM -L "rscgrp=small"
#PJM -L "node=1"
#PJM -L "elapse=00:10:00"
#PJM --mpi "max-proc-per-node=4"
#PJM -j
#PJM -s
#
#
export PARALLEL=12
export OMP_NUM_THREADS=${PARALLEL}
export FORT90L=-Wl,-T
export PLE_MPI_STD_EMPTYFILE=off
export OMP_WAIT_POLICY=active
export FLIB_BARRIER=HARD

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/netcdf-c-4.9.0-g462kcd2ivou7ewax6wddywoyrbz2oib/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/netcdf-fortran-4.6.0-mmdtg5243y4mwqsl3gcu3m2kh27raq5n/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/parallel-netcdf-1.12.3-avpnzm4pwv2tuu2mv73lacb4vhcwlnds/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/hdf5-1.12.2-kb4msz2kuwzsmqsshhpryqebui6tqcfs/lib:$LD_LIBRARY_PATH

cmem=$(printf %04d ${PJM_BULKNUM})
imem=${PJM_BULKNUM}

##### RUN
echo "exec sno"
mpiexec -n 1 -std-proc log/NOUT_sno_${cmem} ./sno  conf/sno_${cmem}.conf
./scale_post_fcst conf/sno_${cmem}.conf > log/NOUT_ref_${cmem}
echo "done."

