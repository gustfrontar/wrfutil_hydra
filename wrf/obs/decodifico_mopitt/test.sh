#/bin/sh
#=======================================================================
# [PURPOSE:]
#   download AIRS IR, AMSU, HSB retrievals -> LETKF format output
#
# [REQUIRED LIBRARIES:]
# 1. HDF-related: hdfeos, hdf5, hdf
# 2. Other standard library (usually in /usr/lib): libjpeg.a libsz.a, libz.a
#
# [PREPARING AIRS DATA URL:]
# 1. Go to http://disc.sci.gsfc.nasa.gov/AIRS/data-holdings/by-data-product
# 2. Go down to AIRS Level-2 Products, and find a row of "AIRX2RET"
# 3. Click a link of "Search" on the rightmost column "GES DISC Data Access"
# 4. Input Time Span and Location (optional) and click "Search GES-DISK"
# 5. "Add Selected Files To Cart"
# 6. "Continue to Shopping Cart" and "Checkout". Note: this is free.
# 7. Click "URL List (Data)"
# 8. Save the file list to "airslist.txt", and remove the final line (README)
#
# [AUTHOR:] T. Miyoshi  March 3, 2011, College Park, MD, USA
# Adapted for CO MOPITT data by Juan Ruiz 2015.
#=======================================================================
set -e
PROXY="proxy.fcen.uba.ar:8080"
F90=gfortran
FFLAGS="-fconvert=big-endian -fno-underscoring -g "
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
#
# DIR
#

TMPDIR=../../DATA/tmp
SOURCEDIR=`pwd`


mkdir -p $TMPDIR

cd $TMPDIR
rm -f *.dat

cp $SOURCEDIR/dec_mopitt.f90 .

echo " >> Making decoder.."
$F90 $FFLAGS -o decoder dec_mopitt.f90 -lGctp -lgeos -lmfhdf -ldf -lhe5_hdfeos  -ljpeg -lsz  -lz -lhdf5 -lhdf5_hl 
#-L$HDFEOS -lhdfeos -lGctp -L$HDF5 -lhdf5 -L$HDF -lmfhdf -ldf 

#
# Download and decode
#
./decoder  MOP02J-20100801-L2V16.2.3.he5
