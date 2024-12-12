##### COMPILATION FLAGS FOR DIFFERENT ARCHITECTURES #######
export COMPILATION_NAME="INTEL_HYDRA"   # FUJITSU_FUGAKU , INTEL_FUGAKU , INTEL_HYDRA
MAKE_LETKF="TRUE"
MAKE_PERT_MET_EM="FALSE"
MAKE_WRF_TO_WPS="FALSE"
DEBUG="FALSE"

#### FUJITSU COMPILER FUGAKU COMPUTER
if [ $COMPILATION_NAME == "FUJITSU_FUGAKU" ] ; then
  if [ $DEBUG == "TRUE" ] ; then
     DEBUG_FLAGS=" -g -traceback "
  fi     
  export PARF90="mpifrtpx"
  export F90="frtpx"
  export FFLAGS="-O3 -Nalloc_assign ${DEBUG_FLAGS}"
  export OMP="-Kopenmp"
  export NETCDF="/home/ra000007/a04037/data/comp_libs_fujitsu/netcdf/"
  export LIB_NETCDF="-L$NETCDF/lib/ -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5_hl -lhdf5  -lfjprofmpi -lmpi_cxx"
  export INC_NETCDF="-I$NETCDF/include/"
  export LBLAS="-SSL2BLAMP"
  export LSCALAPACK="-SCALAPACK"
fi
#####################################

#### INTEL   COMPILER FUGAKU COMPUTER
if [ $COMPILATION_NAME == "INTEL_FUGAKU" ] ; then
  . /opt/intel/oneapi/setvars.sh intel64  
  if [ $DEBUG == "TRUE" ] ; then
     DEBUG_FLAGS=" -g -traceback "        
  fi 
  export PARF90="mpiifort"
  export F90="ifort"
  export FFLAGS="-O3 -xHost -convert big_endian -FR ${DEBUG_FLAGS}"
  export OMP="-qopenmp"
  export NETCDF="/home/ra000007/a04037/data/comp_libs/netcdf4/"
  export LIB_NETCDF="-L$NETCDF/lib/ -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5_hl -lhdf5"
  export INC_NETCDF="-I$NETCDF/include/"
  export LBLAS="-qmkl"
  export LSCALAPACK=""
fi
#####################################

#### INTEL COMPILER HYDRA COMPUTER
if [ $COMPILATION_NAME == "INTEL_HYDRA" ] ; then
  source /opt/load-libs.sh 1             
  if [ $DEBUG == "TRUE" ] ; then
     DEBUG_FLAGS=" -g -traceback "        
  fi 
  export PARF90="mpiifort"
  export F90="ifort"
  export FFLAGS="-O3 -FR ${DEBUG_FLAGS}"
  export OMP="-qopenmp"
  export NETCDF="/opt/netcdf/netcdf-4/intel/2021.4.0/"
  export HDF5="/opt/hdf5/hdf5-1.10.5/intel/2021.4.0/"
  export LIB_NETCDF="-L$NETCDF/lib/ -lnetcdff -lnetcdf  -L$HDF5/lib -lhdf5_fortran -lhdf5_hl -lhdf5"
  export INC_NETCDF="-I$NETCDF/include/"
  export LBLAS="-qmkl"
  export LSCALAPACK=""
fi
#####################################


export current_dir=$(pwd)
if [ $MAKE_LETKF == "TRUE" ] ; then
echo "Compiling LETKF"
cd ${current_dir}/src/wrf/letkf/ 
./make.sh
fi
if [ $MAKE_PERT_MET_EM == "TRUE" ] ; then
echo "Compiling PERT MET EM"
cd ${current_dir}/src/wrf/pert_met_em/
./make.sh
fi
if [ $MAKE_WRF_TO_WPS == "TRUE" ] ; then
echo "Compiling WRF TO WPS"
cd ${current_dir}/src/wrf/wrf_to_wps/
./make.sh
fi













