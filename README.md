# LETKF-WRF



# LETKF-SCALE

## Source codes 

The origins of the codes used in this branch are as follows: 

SCALE: A branch [develop](https://github.com/scale-met/scale-dev/tree/develop) (commit: a8ab69c) of a private repository [scale-dev](https://github.com/scale-met/scale-dev).
SCALE-LETKF: A private repository [scale-letkf](https://github.com/SCALE-LETKF-SMN/scale-letkf) (commit: 157564d) of a version for the PREVENIR project. 

[SCALE User's guide](https://scale.riken.jp/archives/scale_users_guide_En.v5.5.4.pdf) is available online.

## Compilation

### SCALE

SCALE-RM needs to be compiled before LETKF, as it is dependent on SCALE libraries. In compilation, you need the following environmental variables. 
Basic parameters: 
```
export SCALE_SYS=FUGAKU
export SCALE_ENABLE_PNETCDF=F
export SCALE_DISABLE_SDM=T
export SCALE_DISABLE_AMPS=T
export SCALE_ENABLE_JMAPPLIB=F
export SCALE_ENABLE_DA=F     ## SCALE built-in DA libraries are disabled   
export SCALE_USE_SINGLEFP=F  ## T if you use single precision

```
Parameters regarding library paths (use spack_scale.sh):  
```
. /vol0004/apps/oss/spack/share/spack/setup-env.sh

NC_HASH=`spack find -lx netcdf-c%fj | grep netcdf-c | awk '{print $1}'`
NF_HASH=`spack find -lx netcdf-fortran%fj | grep netcdf-fortran | awk '{print $1}'`
PN_HASH=`spack find -lx parallel-netcdf%fj | grep parallel-netcdf | awk '{print $1}'`
HDF_HASH=`spack find -l --deps /${NC_HASH} | grep hdf5 | awk '{print $1}'`
export SCALE_HDF=`spack location --install-dir /${HDF_HASH}`
export SCALE_NETCDF_C=`spack location --install-dir /${NC_HASH}`
export SCALE_NETCDF_F=`spack location --install-dir /${NF_HASH}`
export SCALE_PNETCDF=`spack location --install-dir /${PN_HASH}`

export SCALE_NETCDF_INCLUDE="-I${SCALE_NETCDF_C}/include -I${SCALE_NETCDF_F}/include"
export SCALE_NETCDF_LIBS="-L${SCALE_NETCDF_C}/lib -L${SCALE_NETCDF_F}/lib -L${SCALE_HDF}/lib -L${SCALE_PNETCDF}/lib -lpnetcdf -lnetcdff -lnetcdf -lhdf5_hl -lhdf5"

# for hdf_fortran
export HDF2=(your netcdf library)

export PATH=$SCALE_NETCDF_C/bin:$PATH
```

You can build SCALE-RM and SNO with GNU make. 
```
cd scale_code/scale-rm/src
make -j 
cd ../util/sno
make
```

### SCALE-LETKF

Once you compile SCALE-RM, you have made binary executables (in `bin`) a set of libraries (in `scale-rm/srm/.libs` or `.libs_single`). Now you can compile SCALE-LETKF.
```
cd letkf_code/src_scale/scale
make -j
```
And you will find a binary executable letkf/letkf created. 

## Structure

There are several structural differences between WRF and SCALE forecast procedures. Their similarities and differences are summarized in the figure below.   

<img src="workflow_scale_wrf.png" width=500px>  

Both WRF and SCALE start from generating topography data using external data and the model grid information. SCALE also generates landuse data in this phase. Next, WRF generates met_em files from the external atmospheric data, and those met_em files are used by real.exe to generate initial condition wrfinput file, whereas SCALE creates initial and boundary condition data directly from the external data by scale-rm_init. In the wrfutilV3.1 script, both WRF and SCALE have temporary runtime directories corresponding to pre-process and integration phase, with slight difference to each other. 

## File format 

SCALE input and output files have specific format based on NetCDF. They are all separated into partial files corresponding to subdomains for MPI parallelization. For example, when SCALE-RM runs with 32 processes (= subdomains), 32 history files from history.pe000000.nc to history.pe000031.nc are generated.   

#### Surface boundary condition  

- topography files: files containing 2-d constant surface topography height variable  
- landuse files: files containing 2-d constant landuse information variables (land-sea contrast, urban fraction, etc.)

#### Atmospheric variables

- restart files: files containing 3-d and 2-d instantaneous variables for restart model integration 
- boundary files: files containing 3-d and 2-d instantaneous variables to provide boundary data for model integration  
- history files: files containing time-series of 3-d and 2-d variables   

The list of variables of restart and boundary files are prescribed inside SCALE codes, whereas the variables in history files are specified by namelist parameters.  
Variables in restart files are not very useful for visualization, as they are MOMX(momentum in x) instead of U and RHOT(density times potential temperature) instead of T.  
You need to use history files or manually transform variables from restart files.  

## Post-processing of output files  

SNO (SCALE NetCDF Operator) is a tool for post-processing of SCALE output files. It combines separated NetCDF files into one NetCDF or GrADS binary file. It also supports temporal and spatial regridding. 


