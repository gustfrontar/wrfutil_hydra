#!/bin/bash 

# source this before compile and run letkf on Fugaku

echo "hello"
#exit  0 

. /vol0004/apps/oss/spack/share/spack/setup-env.sh

NC_HASH=`spack find -lx netcdf-c%fj | grep netcdf-c | awk '{print $1}'`
echo NC_HASH=$NC_HASH
NF_HASH=`spack find -lx netcdf-fortran%fj | grep netcdf-fortran | awk '{print $1}'`
PN_HASH=`spack find -lx parallel-netcdf%fj | grep parallel-netcdf | awk '{print $1}'`
HDF_HASH=`spack find -l --deps /${NC_HASH} | grep hdf5 | awk '{print $1}'`
export SPACK_HDF=`spack location --install-dir /${HDF_HASH}`
echo SPACK_NETCDF_C=`spack location --install-dir /${NC_HASH}`
export SPACK_NETCDF_C=`spack location --install-dir /${NC_HASH}`
export SPACK_NETCDF_F=`spack location --install-dir /${NF_HASH}`
export SPACK_PNETCDF=`spack location --install-dir /${PN_HASH}`



