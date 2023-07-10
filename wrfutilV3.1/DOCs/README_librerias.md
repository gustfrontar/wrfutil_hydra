# Instalación de dependencias de WRF
Este tutorial lo guiará en la instalación de un conjunto de librerías necesarias para instalar el WRF 4.0
Las librerías y sus versiones como figuran en este tutorial garantizan una compilación exitosa y el correcto funcionamiento del WRF, pero otro conjunto podría ser tambien valido

>Para mas información de como compilar el WRF ver:
>   https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php

crear un directorio donde guardar todos los fuentes descargados de internet

    SOURCEDIR=<path to sources>
    mkdir $SOURCEDIR
    cd $SOURCEDIR

crear un archivo *envvars.sh* donde se definirán todas las variables de entorno necesarias para compilar los fuentes.


    cat << 'EOF' >  envvars.sh
    export LC_ALL="C"

    export SOURCEDIR=/home/yanina/sources
    export APPDIR=/home/yanina/WRF-4.0_${CVER}/

    export INTEL_BASE=/home/opt/intel/
    export CVER="I19"
    export COMPILERVARS_ARCHITECTUR=intel64
    export COMPILERVARS_PLATFORM=linux
    export INTEL_COMPILER_TOPDIR=$INTEL_BASE/parallel_studio_xe_2019
    export IMPI_BASE=$INTEL_BASE/impi/2019.1.144/
    source $INTEL_COMPILER_TOPDIR/psxevars.sh $COMPILERVARS_ARCHITECTUR
    source $IMPI_BASE/intel64/bin/mpivars.sh $COMPILERVARS_ARCHITECTUR


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
    \# export CFLAGS=' -xHost'
    \# export CXXFLAGS=' -xHost'
    \# export FFLAGS=' -xHost'
    \# export FCFLAGS=' -xHost'
    \# export CPP='icc -E'
    \# export CXXCPP='icpc -E'
    \# export CC=$IMPI_BASE/intel64/bin/mpicc
    \# export CXX=$IMPI_BASE/intel64/bin/mpicxx
    \# export F77=$IMPI_BASE/intel64/bin/mpif77
    \# export FC=$IMPI_BASE/intel64/bin/mpifc
    \# export F90=$IMPI_BASE/intel64/bin/mpif90
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
    \# export I_OMPI_CC=icc
    \# export I_OMPI_CXX=icpc
    \# export I_OMPI_F77=ifort
    \# export I_OMPI_FC=ifort
    \# export I_OMPI_F90=ifort
    export MPIF90=mpiifort
    export MPIF77=mpiifort
    export MPICC=mpiicc
    export MPICXX=mpiicpc
    export DM_FC=mpiifort
    export DM_CC=mpiicc
    \#export MPIF90=$IMPI_BASE/intel64/bin/mpif90
    \#export MPIF77=$IMPI_BASE/intel64/bin/mpif77
    \#export MPIFC=$IMPI_BASE/intel64/bin/mpifc
    \#export MPICC=$IMPI_BASE/intel64/bin/mpicc
    \#export MPICXX=$IMPI_BASE/intel64/bin/mpicxx
    \#export DM_FC=$IMPI_BASE/intel64/bin/mpifc
    \#export DM_CC=$IMPI_BASE/intel64/bin/mpicc
    export DM_FC=$IMPI_BASE/intel64/bin/mpiifort
    export DM_CC=$IMPI_BASE/intel64/bin/mpicc

    export PATH=$JPEGLIB/bin:$NETCDF/bin:$HDF5/bin:$HDF4/bin:$JASPER/bin:/$LIBPNG/bin:$GADDIR/bin:$WGRIB2:$PATH
    export WRFIO_NCD_LARGE_FILE_SUPPORT=1

    export LD_LIBRARY_PATH=$JPEGLIB/lib:$HDF4/lib:$HDFEOS/lib:$HDFEOS5/lib:$ZLIB/lib:$HDF5/lib:$NETCDF/lib:/$PNETCDF/lib:$LIBPNG/lib:$JASPERLIB:$SZLIB/lib:$LD_LIBRARY_PATH
    export LD_RUN_PATH=$LD_LIBRARY_PATH:$LD_RUN_PATH
    export INCLUDE=$JPEGLIB/include:$HDF4/include:$PNETCDF/include:$HDF5/include:$NETCDF/include:$ZLIB/include:$JASPERINC:$LIBPNG/include:$INLCUDE:$SZLIB/include
    EOF

Edite las siguientes variables (como mínimo) con los datos de su instalación

>* export SOURCEDIR=\<path to sources>
>* export APPDIR=\<directorio base de instalacion>
>* export INTEL_BASE=\<directorio base del compilador de Intel>
>* export CVER=\<codigo de compilador p/e "I19" para intel2019>

con todas las variables adaptas, cargue las variables de entorno y cree el directorio destino de la instalación

    source envvars.sh
    mkdir -p $APPDIR/LIBRARY

#### Descargar
descargar en $SOURCEDIR los siguentes archivos :
* szip-2.1.1.tar.gz
* zlib-1.2.11.tar.gz
* libpng-1.6.32.tar.gz
* yasm-1.3.0.tar.gz
* cmake-3.9.1.tar.gz
* libjpeg-turbo-1.5.2.tar.gz
* jasper-2.0.12.tar.gz  
* hdf-4.2.13.tar.gz
* hdf5-1.10.1.tar.gz
* HDF-EOS5.1.15.tar.Z
* HDF-EOS2.19v1.00.tar.Z
* parallel-netcdf-1.8.1.tar.gz
* netcdf-4.4.1.tar.gz
* netcdf-fortran-4.4.4.tar.gz

#### Compilación e instalación
 El siguiente procedimiento instalará todas las dependencia en el directorio $DIRLIB definido en el archivo envvars.sh

##### zlib-1.2.11.tar.gz
    cd $SOURCEDIR
    tar xvzf zlib-1.2.11.tar.gz
    cd zlib-1.2.11
    make distclean
    ./configure --prefix=$ZLIB
    make -j 40
    make check
    make install
    cd ..
    #rm -r zlib-1.2.11

##### szip-2.1.1.tar.gz
    cd $SOURCEDIR
    tar xvzf szip-2.1.1.tar.gz
    cd szip-2.1.1
    ./configure --prefix=$SZLIB
    make
    make check
    make install
    cd ..
    #rm -r szip-2.1.1

##### libpng-1.6.32.tar.gz
    cd $SOURCEDIR
    tar xvzf libpng-1.6.32.tar.gz
    cd libpng-1.6.32
    make distclean
    ./configure --prefix=$LIBPNG
    make
    make check
    make install
    cd ..
    #rm -r libpng-1.6.32
##### yasm-1.3.0.tar.gz
    cd $SOURCEDIR
    tar xvzf yasm-1.3.0.tar.gz
    cd yasm-1.3.0
    make distclean
    ./configure --prefix=$YASM
    make
    make check
    make install
    cd ..
    #rm -r yasm-1.3.0
##### cmake-3.9.1.tar.gz

    cd $SOURCEDIR  
    tar xvzf cmake-3.9.1.tar.gz   
    cd cmake-3.9.1    
    ./bootstrap --prefix=$CMAKE     
    gmake install    
    cd ..    
    #rm -r cmake-3.9.1    

##### libjpeg-turbo-1.5.2.tar.gz
    cd $SOURCEDIR
    tar xvzf libjpeg-turbo-1.5.2.tar.gz
    cd libjpeg-turbo-1.5.2
    ./configure --prefix=$JASPER --with-jpeg8
    make
    make check
    make install
    cd ..
    #rm -r libjpeg-turbo-1.5.2

##### jasper-2.0.12.tar.gz
Durante el make check, el test 3 da error, ignorar

    cd $SOURCEDIR
    tar xvzf jasper-2.0.12.tar.gz
    cd jasper-2.0.12
    mkdir build
    $CMAKE/bin/cmake -G "Unix Makefiles" -H$SOURCEDIR/jasper-2.0.12 -Bbuild -    DCMAKE_INSTALL_PREFIX=$JASPER -DJAS_ENABLE_OPENGL=false
    cd build
    make clean all
    make test 
    make install
    cd ..
    rm -r jasper-2.0.12

##### hdf-4.2.13.tar.gz
Hay algunas opciones del check que dan error.

    cd $SOURCEDIR
    tar xvzf hdf-4.2.13.tar.gz
    cd hdf-4.2.13 ### Apartir de aca se puede instalar una version paralela
    CC=mpicc FC=mpifc ./configure --prefix=$HDF4 --enable-fortran --with-zlib=$ZLIB --with-szlib=$SZLIB --enable-parallel --disable-netcdf --with-jpeg=$JASPER
    make -j 40
    make check
    make install
    cd ..
    #rm -r hdf-4.2.13

##### hdf5-1.10.1.tar.gz
Hay algunas opciones del check que dan error.

    cd $SOURCEDIR
    tar xvzf hdf5-1.10.1.tar.gz
    cd hdf5-1.10.1
    make distclean
    CC=mpicc FC=mpifc ./configure --prefix=$HDF5 --enable-fortran --with-zlib=$ZLIB --enable-parallel --enable-shared --with-szlib=$SZLIB
    make -j 40
    make check
    make install
    cd ..
    #rm -r hdf5-1.10.1

##### HDF-EOS5.1.15.tar.Z
Para el decodificador de airs solo se necesita las EOS2  
Quedan errores en el check

    cd $SOURCEDIR
    tar xvzf HDF-EOS5.1.15.tar.Z
    cd hdfeos5
    ./configure CC="$HDF5/bin/h5pcc -Df2cFortran" --prefix=$HDFEOS5 --enable-fortran --with-szlib=$SZLIB
    make
    make check
    make install
    cd ..
    #rm -rf hdfeos5

##### HDF-EOS2.19v1.00.tar.Z
Quedan errores en el check

    cd $SOURCEDIR
    tar xvzf HDF-EOS2.19v1.00.tar.Z
    cd hdfeos
    ./configure CC="$HDF4/bin/h4cc -Df2cFortran" --prefix=$HDFEOS --enable-fortran --with-szlib=$SZLIB
    make
    make check
    make install
    cd ..
    #rm -rf hdfeos

##### parallel-netcdf-1.8.1.tar.gz
Tiene que andar todo OK
Si descargo los archivos de testing paralelo, descomente la y ejecute el **make ptest**

    cd $SOURCEDIR
    tar xvzf parallel-netcdf-1.8.1.tar.gz
    cd parallel-netcdf-1.8.1
    make distclean
    ./configure --prefix=$PNETCDF --enable-fortran --enable-largefile #--with-szlib=$SZLIB
    make -j20
    make testing
    #make ptest
    make install
    cd ..
    #rm -r parallel-netcdf-1.8.1

##### netcdf-4.4.1.tar.gz

Reexportar las variables como se indica a continuacion


    export CFLAGS="-xHost -g -fPIC -I${HDF5}/include -I${ZLIB}/include -I$ {PNETCDF}/include "
    export CXXFLAGS="-xHost -g -fPIC -I${HDF5}/include -I${ZLIB}/include -I$ {PNETCDF}/include "
    export FFLAGS="-xHost -g -fPIC -I${HDF5}/include -I${ZLIB}/include -I$ {PNETCDF}/include "
    export FCFLAGS="-xHost -g -fPIC -I${HDF5}/include -I${ZLIB}/include -I$ {PNETCDF}/include "
    export FLDFLAGS="-fPIC -L${HDF5}/lib -L${ZLIB}/lib -L${PNETCDF}/lib"
    export F90LDFLAGS="-fPIC -L${HDF5}/lib -L${ZLIB}/lib -L${PNETCDF}/lib"
    export LDFLAGS="-fPIC -L${HDF5}/lib -L${ZLIB}/lib -L${PNETCDF}/lib"

    cd $SOURCEDIR
    tar xvzf netcdf-4.4.1.tar.gz
    cd netcdf-4.4.1
    make distclean
    CC=mpicc FC=mpifc ./configure --prefix=$NETCDF --enable-fortran --disable-static --enable-shared --with-pic --enable-parallel-tests --enable-pnetcdf --enable-largefile
    make
    make check
    make install
    ##  --- > para recordar nc-config --help
    cd ..
    #rm -r netcdf-4.4.1 

##### netcdf-fortran-4.4.4.tar.gz


    export CFLAGS="-xHost -g -fPIC -I${HDF5}/include -I${ZLIB}/include -I$ {NETCDF}/include"
    export CXXFLAGS="-xHost -g -fPIC -I${HDF5}/include -I${ZLIB}/include -I$ {NETCDF}/include"
    export FFLAGS="-xHost -g -fPIC -I${HDF5}/include -I${ZLIB}/include -I$ {NETCDF}/include"
    export FCFLAGS="-xHost -g -fPIC -I${HDF5}/include -I${ZLIB}/include -I$ {NETCDF}/include"
    export FLDFLAGS="-fPIC -L${HDF5}/lib -L${ZLIB}/lib -L${NETCDF}/lib "
    export F90LDFLAGS="-fPIC -L${HDF5}/lib -L${ZLIB}/lib -L${NETCDF}/lib "
    export LDFLAGS="-fPIC -L${HDF5}/lib -L${ZLIB}/lib -L${NETCDF}/lib "

    cd $SOURCEDIR
    tar xvzf netcdf-fortran-4.4.4.tar.gz
    cd netcdf-fortran-4.4.4
    make distclean
    CC=mpicc FC=mpifc F90=mpif90 F77=mpif77 ./configure --prefix=$NETCDF --disable-static --enable-shared --with-pic --enable-parallel-tests --enable-largefile
    make
    make check
    make install
    ##para recordar nc-config --help
    cd ..
    #rm -r netcdf-fortran-4.4.4


Todas las dependencias necesarias para instalar el WPS, WRF y LETKF deberían estar instaladas.
continúe con la instalación de los modelos.
