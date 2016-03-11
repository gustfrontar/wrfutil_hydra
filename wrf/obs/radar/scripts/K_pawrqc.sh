#!/bin/bash

TMPDIR=$HOME/data/TMP/RADARQC/
EXECPATH=$HOME/share/LETKF_WRF/wrf/obs/radar/radar_qc/
STATPATH=$HOME/share/LETKF_WRF/wrf/obs/radar/radar_qc/statistics_122014/
TOPOPATH=$HOME/share/LETKF_WRF/wrf/obs/radar/radar_qc/terrain_data/

DATADIR=$HOME/share/OBS/osaka_pawr/

mkdir -p $TMPDIR

DATE=20130713

export OMP_NUM_THREADS=10
cp $EXECPATH/radar_qc.exe $TMPDIR/

cd $TMPDIR

#Link the data
ln -sf $DATADIR/VR_* .
ln -sf $DATADIR/Z_*  .

#Link the statistics and the corresponding topgraphy file.
ln -sf $STATPATH/*.txt .
ln -sf $TOPOPATH/terrain_$DATE.bin ./terrain.bin 

echo "RUNNING THE PROGRAM"

ulimit -s unlimited

./radar_qc.exe




