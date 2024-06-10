
clear

echo "
      --------------------------------------------------------------------------

       Entorno de variables para la compilacion de las dependencias de WRF Util

      --------------------------------------------------------------------------


"
export RAIZ='/data/share/bin/'

export LC_ALL="C"

export CVER="I19"
### ADAPTAR ESTE PATH
export SOURCEDIR="$RAIZ/sources"
export APPDIR="$RAIZ/WRF-4.0_${CVER}/"
###
export INTEL_BASE="/home/opt/intel/"
export COMPILERVARS_ARCHITECTUR="intel64"
export COMPILERVARS_PLATFORM="linux"
export INTEL_COMPILER_TOPDIR="$INTEL_BASE/parallel_studio_xe_2019"
export IMPI_BASE="$INTEL_BASE/impi/2019.1.144/"
source $INTEL_COMPILER_TOPDIR/psxevars.sh "$COMPILERVARS_ARCHITECTUR"
source $IMPI_BASE/intel64/bin/mpivars.sh "$COMPILERVARS_ARCHITECTUR"


export DIRLIB=$APPDIR/LIBRARIES/
export ZLIB=$DIRLIB/zlib-1.2.11_${CVER}
export SZLIB=$DIRLIB/szlib-2.1.1_${CVER}
export YASM=$DIRLIB/yasm-1.3.0_${CVER}
export HDF5=$DIRLIB/hdf5-1.10.1_${CVER}
export HDF4=$DIRLIB/hdf-4.2.13_${CVER}
export PHDF5=$DIRLIB/hdf5-1.10.1_${CVER}
export HDFEOS5=$DIRLIB/hdfeos5.1.15_${CVER}
export HDFEOS=$DIRLIB/hdfeos2.19.1_${CVER}
export PNETCDF=$DIRLIB/pnetcdf-1.8.1_${CVER}
export NETCDF=$DIRLIB/netcdf-4.4.1_${CVER}
export JASPER=$DIRLIB/jasper-2.0.12_${CVER}
export JASPERLIB=$JASPER/lib64
export JASPERINC=$JASPER/include
export JPEGLIB=$JASPER/
export LIBPNG=$DIRLIB/libpng-1.6.32_${CVER}
export WRF_DIR=$APPDIR
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export GADDIR=$DIRLIB/grads-2.1.a3
export GAUDPT=$GADIR=$GADDIR/udpt
export WGRIB2=$DIRLIB/grib2/wgrib2/
export CMAKE=$DIRLIB/cmake

export CC=icc
export CXX=icpc
export F77=ifort
export FC=ifort
export F90=ifort
# export CFLAGS=' -xHost'
# export CXXFLAGS=' -xHost'
# export FFLAGS=' -xHost'
# export FCFLAGS=' -xHost'
# export CPP='icc -E'
# export CXXCPP='icpc -E'
# export CC=$IMPI_BASE/intel64/bin/mpicc
# export CXX=$IMPI_BASE/intel64/bin/mpicxx
# export F77=$IMPI_BASE/intel64/bin/mpif77
# export FC=$IMPI_BASE/intel64/bin/mpifc
# export F90=$IMPI_BASE/intel64/bin/mpif90
export CFLAGS='-xHost -g -fPIC'
export CXXFLAGS='-xHost -g -fPIC'
export FFLAGS='-xHost -g -fPIC'
export FCFLAGS='-xHost -g -fPIC'
export FLDFLAGS='-fPIC'
export F90LDFLAGS='-fPIC'
export LDFLAGS='-fPIC'
export CPP='icc -E'
export CXXCPP='icpc -E'
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_FC=ifort
export I_MPI_F90=ifort
# export I_OMPI_CC=icc
# export I_OMPI_CXX=icpc
# export I_OMPI_F77=ifort
# export I_OMPI_FC=ifort
# export I_OMPI_F90=ifort
export MPIF90=mpiifort
export MPIF77=mpiifort
export MPICC=mpiicc
export MPICXX=mpiicpc
export DM_FC=mpiifort
export DM_CC=mpiicc
#export MPIF90=$IMPI_BASE/intel64/bin/mpif90
#export MPIF77=$IMPI_BASE/intel64/bin/mpif77
#export MPIFC=$IMPI_BASE/intel64/bin/mpifc
#export MPICC=$IMPI_BASE/intel64/bin/mpicc
#export MPICXX=$IMPI_BASE/intel64/bin/mpicxx
#export DM_FC=$IMPI_BASE/intel64/bin/mpifc
#export DM_CC=$IMPI_BASE/intel64/bin/mpicc
export DM_FC=$IMPI_BASE/intel64/bin/mpiifort
export DM_CC=$IMPI_BASE/intel64/bin/mpicc

export PATH=$JPEGLIB/bin:$NETCDF/bin:$HDF5/bin:$HDF4/bin:$JASPER/bin:/$LIBPNG/bin:$GADDIR/bin:$WGRIB2:$PATH:$DIRLIB/yasm-1.3.0_I19/bin:/home/opt/intel/intelpython3/bin/
export WRFIO_NCD_LARGE_FILE_SUPPORT=1

export LD_LIBRARY_PATH=$JPEGLIB/lib:$HDF4/lib:$HDFEOS/lib:$HDFEOS5/lib:$ZLIB/lib:$HDF5/lib:$NETCDF/lib:/$PNETCDF/lib:$LIBPNG/lib:$JASPERLIB:$SZLIB/lib:$LD_LIBRARY_PATH
export LD_RUN_PATH=$LD_LIBRARY_PATH:$LD_RUN_PATH
export INCLUDE=$JPEGLIB/include:$HDF4/include:$PNETCDF/include:$HDF5/include:$NETCDF/include:$ZLIB/include:$JASPERINC:$LIBPNG/include:$INLCUDE:$SZLIB/include

echo "


"
