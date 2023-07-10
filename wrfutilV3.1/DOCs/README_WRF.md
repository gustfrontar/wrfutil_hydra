## Compilación del WRF

Este tutorial lo guiará durante todo el proceso de instalación del Modelo WRF y la herramienta de pre-proceso WPS. Antes de comenzar, asegurese de tener instaladas todas las dependencias necesarias. Ver [acá](./README_librerias.md) para una guía de instalación de dependencias.

Los pasos sugeridos en esta guía están basados en los procesos de Compilación y optimización sugeridos en la siguiente bibliografía:

 * **WRF Installation Best Prectices**, *HPC Advisory Council*
 * **WRF Build with Parallel NetCDF Support**, *ARC Centre of Excellence for Climate System Science*
 * **Influence of th compiler on multi-CPU performance of WRFv3**, *Geoscientific Model Development*
 * **WRF code Optimization for Meso-scale Process Studies**, *(WOMPS) dCSE Projec Report*  
 * **Performance Analysis and Optimization of th Weather Research and Forecasting Model (WRF) on Intel Multicore and Manycore Architecture**, *Samuel Elliott,  University of Colorado at Boulder: Department of Applied Mathematics*
 * **Performance Hints for WRF on Intel architecture**, *Intel Developer Zone*
 * **WRF Quilting and Decomposition Notes**, *Eric Kemp MArch 2015*
 * **Opportunities for WRF Model Acceleration**,*John Michalakes, Computatinal Sciences Center NREL*
 * **23 tips for performance tuning with the Intel MPI Library**, *Software & Services Group, Developer Products Division, Intel*  
 * **Tuning WRF Performance for HPC Cloud Sandy Bridge System**, *2014 Penguin Computing*
y en un estudio llevado a cabo en las oficinas de I+D del SMN


Primero se debera bajar los fuentes del WRF y WPS y Geog

    cd $SOURCEDIR
    wget http://www2.mmm.ucar.edu/wrf/src/WPSV4.0.TAR.gz 
    wget http://www2.mmm.ucar.edu/wrf/src/WRFV4.0.TAR.gz 
    wget http://www2.mmm.ucar.edu/wrf/src/wps_files/geog_high_res_mandatory.tar.gz 
    wget http://www2.mmm.ucar.edu/wrf/src/iowrf.f 
    # Para la version 4.0 hay que bajar este parche 
    wget http://www2.mmm.ucar.edu/wrf/src/fix/start_em.F_theta_m_pressure_fix.tar 


Si se bajo el tar ejecutar si se so el git saltear el tar

    tar xvzf $SOURCEDIR/WRFV4.0.TAR.gz -C $WRF_DIR
    tar xvf $SOURCEDIR/start_em.F_theta_m_pressure_fix.tar
    cp $SOURCEDIR/start_em.F $WRF_DIR/WRF/dyn_em/
    cd $WRF_DIR
    cp -rp $WRF_DIR/WRF $WRF_DIR/WRFDA
    cd $WRF_DIR/WRF
    ./configure

Durante el proceso de configuración se solicitará responder algunas preguntas, para el instalación e HM del SMN seleccionar:

  1. **21** ((dm+sm) INTEL (ifort/icc): Xeon (SNB with AVX mods))
  2. **default**

Luego, editar $WRF_DIR/WRF/configure.wrf y modificar:
 * Agregar a la definición de la variable **ARCHFLAGS**  
   * -DPNETCDF_QUILT
 * Reemplazar la opcion de *openmp* por *qopenmp* en la definicion de las variables **OMP** y **OMPCC** :
   * OMP = -qopenmp -fpp -auto
   * OMPCC = -qopenmp -fpp -auto

    ./compile em_real >& log.compile

Revisar el archivo *log.compile* y verificar que el proceso de compilación haya concluido exitosamente. De ser así, se deberían haber generado en el directorio *main* los siguiente archivos:

    ls -ls main/*.exe

  - wrf.exe (model executable)
  - real.exe (real data initialization)
  - ndown.exe (one-way nesting)
  - tc.exe (for tc bogusing--serial only)


## Compilación del WRFDA

    cd $WRF_DIR/WRFDA ./configure wrfda


Durante el proceso de configuración se solicitará responder algunas preguntas, para el instalación e HM del SMN seleccionar:

      1. **21** ((dm+sm) INTEL (ifort/icc): Xeon (SNB with AVX mods))
      2. **default**

Luego, editar $WRF_DIR/WRFDA/configure.wrf y modificar:

  * Agregar a la definición de la variable **ARCHFLAGS**  
    * -DPNETCDF_QUILT
  * Reemplazar la opcion de *openmp* por *qopenmp* en la definicion de las variables **OMP** y **OMPCC** :
    * OMP = -qopenmp -fpp -auto
    * OMPCC = -qopenmp -fpp -auto


    ./compile all_wrfvar >& log.compile


## Instalación del WPS
Si se bajo el tar ejecutar 

    tar xvzf $SOURCEDIR/WPSV4.0.TAR.gz -C $WRF_DIR
    cd $WRF_DIR/WPS
    ./clean -a
    ./configure

Durante el proceso de configuración se solicitará responder algunas preguntas, para el instalación e HM del SMN seleccionar:
  * **19** (Linux x86_64, Intel compiler (dmpar))  

Luego editar el $WRF_DIR/WPS/configure.wps y modificar:

 * Corroborar la definición de la variable  WRF_DIR

       WRF_LIB = -L$(WRF_DIR)/external/io_grib1 -lio_grib1 \
       -L$ (INTEL_COMPILER_TOPDIR)/compilers_and_libraries_2016.0.109/linux/compiler/lib/in tel64_lin -liomp5 \
       -L$(WRF_DIR)/external/io_grib_share -lio_grib_share \
       -L$(WRF_DIR)/external/io_int -lwrfio_int \
       -L$(WRF_DIR)/external/io_netcdf -lwrfio_nf \
       -L$(NETCDF)/lib -lnetcdff -lnetcdf

y modificar las siguientes lineas  para que queden así:

    COMPRESSION_LIBS = -L$(LIBPNG)/lib -L$(JASPERLIB) -ljasper -lpng -lz
    COMPRESSION_INC = -I$(JASPERINC) -I$(LIBPNG)/include

también debe editar el archivo **Makefile** en los directorios **metgrid** y **geogrid** para que la definición de la variable *FFLAG y *CFLAGS* quede así:  

  * FFLAGS="$(FFLAGS) -qopenmp -fpp -auto" \  
  * CFLAGS="$(CFLAGS) -qopenmp "  \

y luego compilar 

    ./compile >& log.compile

y verificar que se generen los siguiente archivos:

  * geogrid.exe -> geogrid/src/geogrid.exe
  * ungrib.exe -> ungrib/src/ungrib.exe
  * metgrid.exe -> metgrid/src/metgrid.exe

## Instalación de la base de datos de Superficie
Descomprimo los datos de superficie (geog)

    cd $WRF_DIR
    tar -zxvf $SOURCES/geog_high_res_mandatory.tar.gz
