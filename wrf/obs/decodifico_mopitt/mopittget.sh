#/bin/sh
#=======================================================================
# [PURPOSE:]
#   download MOPITT retrievals -> LETKF format for CO retrievals
#
# [REQUIRED LIBRARIES:]
# 1. HDF-related: hdfeos, hdf5, hdf
# 2. Other standard library (usually in /usr/lib): libjpeg.a libsz.a, libz.a
#
# [PREPARING AIRS DATA URL:]
# 1. Go to ftp://l5eil01.larc.nasa.gov//misrl2l3/MOPITT/MOP02J.006/2010.08.03/
# 2. Download data (1 file per day)
# 3. Decode the data
#
# [AUTHOR:] T. Miyoshi  March 3, 2011, College Park, MD, USA
# Adapted for MOPITT CO retrievals by Juan Ruiz
#=======================================================================

set -e
PROXY="proxy.fcen.uba.ar:8080"
F90=gfortran
FFLAGS="-fconvert=big-endian -fno-underscoring"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/

#
# DIR
#

CDIR=`pwd`
cd ../../
WRF=`pwd`
OBSDIR=$WRF/DATA/obs/
DLDIR=$OBSDIR/downloads/
SAVEDIR=$OBSDIR/mopitt
mkdir -p $DLDIR/mopitt
mkdir -p $SAVEDIR

echo $OBSDIR
echo $SAVEDIR

#
# WKDIR
#
WKDIR=$OBSDIR/tmp
rm -rf $WKDIR
mkdir -p $WKDIR
cd $WKDIR
#
# NETWORK
#
export http_proxy=$PROXY


#
# Make decoder
#

echo " >> Making decoder.."
cp $CDIR/dec_mopitt.f90 .
cp $CDIR/*.inc          .
$F90 $FFLAGS -o decoder dec_mopitt.f90 -lGctp -lgeos -lmfhdf -ldf -lhe5_hdfeos  -ljpeg -lsz  -lz -lhdf5 -lhdf5_hl

#
# Download and decode
#

cp $CDIR/mopittlist.txt ./list.txt
NFILES=`wc -l list.txt | awk '{print $1}'`
N=1

while test $N -le $NFILES
do

URL=`cat list.txt | head -n $N | tail -n 1`
FILE=`basename $URL`

if test -f $DLDIR/mopitt/$FILE
then
echo " >> $FILE exists. skip downloading."
ln -sf $DLDIR/mopitt/$FILE ./current.hdf
else
echo " >> downloading.. $FILE"
curl $URL -o ./current.hdf
cp ./current.hdf $DLDIR/mopitt/$FILE
fi
./decoder ./current.hdf

N=`expr $N + 1`
done


mv *.dat $SAVEDIR/

echo "NORMAL END"
