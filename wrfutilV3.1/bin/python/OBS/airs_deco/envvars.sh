export JPEGLIB=$DIRLIB/jpeg-turbo-1.5.2
export PATH=$JPEGLIB/bin:$NETCDF/bin:$HDF5/bin:$HDF4/bin:$JASPER/bin:/$LIBPNG/bin:$GADDIR/bin:$WGRIB2:$NCO/bin:$PATH
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export LD_LIBRARY_PATH=$JPEGLIB/lib:$NCO/lib:$HDF4/lib:$HDFEOS/lib:$HDFEOS5/lib:$ZLIB/lib:$HDF5/lib:$NETCDF/lib:/$PNETCDF/lib:$LIBPNG/lib:$JASPERLIB:$JASPER/lib:$SZLIB/lib:$LD_LIBRARY_PATH
export LD_RUN_PATH=$LD_LIBRARY_PATH:$LD_RUN_PATH
export INCLUDE=$JPEGLIB/include:$HDF4/include:$PNETCDF/include:$HDF5/include:$NETCDF/include:$ZLIB/include:$JASPERINC:$LIBPNG/include:$INLCUDE:$SZLIB/include
