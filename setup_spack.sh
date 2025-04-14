#!/bin/bash 

# source this before compile and run letkf on Fugaku

. /vol0004/apps/oss/spack/share/spack/setup-env.sh

NC_HASH=$(spack find -lx netcdf-c%fj | grep netcdf-c | tail -n 1 | awk '{print $1}')
NF_HASH=$(spack find -lx netcdf-fortran%fj | grep netcdf-fortran | tail -n 1 | awk '{print $1}')
PN_HASH=$(spack find -lx parallel-netcdf%fj | grep parallel-netcdf | tail -n 1 | awk '{print $1}')
HDF_HASH=$(spack find -l --deps /${NC_HASH} | grep hdf5 | tail -n 1 | awk '{print $1}')
export SCALE_HDF=$(spack location --install-dir /${HDF_HASH})
export SCALE_NETCDF_C=$(spack location --install-dir /${NC_HASH})
export SCALE_NETCDF_F=$(spack location --install-dir /${NF_HASH})
export SCALE_PNETCDF=$(spack location --install-dir /${PN_HASH})

export SCALE_NETCDF_INCLUDE="-I${SCALE_NETCDF_C}/include -I${SCALE_NETCDF_F}/include"
export SCALE_NETCDF_LIBS="-L${SCALE_NETCDF_C}/lib -L${SCALE_NETCDF_F}/lib -L${SCALE_HDF}/lib -L${SCALE_PNETCDF}/lib -lpnetcdf -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lfjprofmpi -lmpi_cxx"


# for WRF 
export NETCDF=/data/ra000007/u10335/PREVENIR/wrfutil_hydra/LIBRARIES/netcdf
export NETCDF_classic=1
export GRIB2=/data/ra000007/u10335/PREVENIR/wrfutil_hydra/LIBRARIES/grib2

# for hdf_fortran
#export HDF2=/data/hp150019/u10335/test_wrf/wrfutil_hydra_newest/LIBRARIES/netcdf

export PATH=$SCALE_NETCDF_C/bin:$PATH
export LD_LIBRARY_PATH=$GRIB2/lib:$SCALE_NETCDF_C/lib:$SCALE_NETCDF_F/lib:$SCALE_PNETCDF/lib:$SCALE_HDF/lib:$LD_LIBRARY_PATH


export MACHINE=FUGAKU
export SCALE_SYS=FUGAKU
export SCALE_DB=/data/ra000007/u10335/scale_database
export SCALE_ENABLE_PNETCDF=F
export SCALE_USE_SINGLEFP=T

export JASPERLIB=$GRIB2/lib
export JASPERINC=$GRIB2/include

export FC=frtpx
export F90=frtpx
export CC=fccpx
export CXX=fccpx

