##### COMPILATION FLAGS FOR DIFFERENT ARCHITECTURES #######
export COMPILATION_NAME="FUJITSU_FUGAKU"   # FUJITSU_FUGAKU , INTEL_FUGAKU , INTEL_HYDRA
MAKE_LETKF="TRUE"
MAKE_PERT_MET_EM="TRUE"
MAKE_WRF_TO_WPS="TRUE"
DEBUG="FALSE"

#### FUJITSU COMPILER FUGAKU COMPUTER
if [ $COMPILATION_NAME == "FUJITSU_FUGAKU" ] ; then
  if [ $DEBUG == "TRUE" ] ; then
     DEBUG_FLAGS="-O0 -Nlist=i,lst=t -X03 -v03s -v03d -v03o -Ncompdisp -Koptmsg=1 -Cpp -Ec -Eg -Ha -He -Hf -Ho -Hs -Hu -Hx -Ncheck_global"
#     DEBUG_FLAGS=" -g -traceback "
  fi     
  export PARF90="mpifrtpx"
  export F90="frtpx"
  export FFLAGS="-O3 -Nalloc_assign ${DEBUG_FLAGS}"
  export OMP="-Kopenmp"
  export LIB_NETCDF="-L$SPACK_NETCDF_F/lib/ -lnetcdff -L$SPACK_NETCDF_C/lib/ -lnetcdf -L$SPACK_HDF/lib/ -lhdf5_hl -lhdf5  -lfjprofmpi -lmpi_cxx"
  export INC_NETCDF="-I$SPACK_NETCDF_F/include/ -I$SPACK_NETCDF_C/include -I$SPACK_HDF/include"
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

#### INTEL COMPILER FUGAKU HYDRA
if [ $COMPILATION_NAME == "INTEL_HYDRA" ] ; then
#  source /opt/load-libs.sh 1             
  if [ $DEBUG == "TRUE" ] ; then
     DEBUG_FLAGS=" -g -traceback "        
  fi 
  export PARF90="mpiifort"
  export F90="ifort"
  export FFLAGS="-O3 -xHost -convert big_endian -FR -O3 ${DEBUG_FLAGS}"
  export OMP="-qopenmp"
  export NETCDF="/data/share/bin/WRF-4.0_I19/LIBRARIES/netcdf-4.4.1_I19/"
  export HDF5="/data/share/bin/WRF-4.0_I19/LIBRARIES/hdf5-1.10.1_I19/"
  export LIB_NETCDF="-L$NETCDF/lib/ -lnetcdff -lnetcdf  -L$HDF5/lib -lhdf5_fortran -lhdf5_hl -lhdf5"
  export INC_NETCDF="-I$NETCDF/include/"
  export LBLAS="-mkl"
  export LSCALAPACK=""
fi
#####################################

mkdir -p bin
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













