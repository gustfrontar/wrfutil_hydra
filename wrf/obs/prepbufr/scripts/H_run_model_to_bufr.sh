#This script generates pseudo-radar observations from WRF outputs.
ulimit -s unlimited

EXPERIMENT_NAME=OSSE_BUFR_SA

MODELDATAPATH="/salidas/felix.carrasco/TRUERUN/total/"
OBSDATAPATH="$HOME/share/OBS/$EXPERIMENT_NAME"   
OBSSOURCE="$HOME/share/OBS/PREPBUFRSA/"   #If a realistic network is desired then real observation files has to be provided.

TMPDIR="$HOME/data/TMP/WRF_TO_BUFR/"
EXEC="$HOME/share/LETKF_WRF/wrf/obs/prepbufr/wrf_to_bufr/wrf_to_bufr.exe"

source ../../../run/util.sh #Date functions.

INIDATE=20100801000000
ENDDATE=20100930180000
FREQ=3600

NAMELIST_WRF2RADAR="$TMPDIR/w2b.namelist"

#Namelist parameters
add_obs_error=".true."      #If an observation error will be added.
u_error="1.0d0"
v_error="1.0d0"
t_error="1.0d0"
q_error="1.0d-3"
ps_error="1.0d2"
uobs=".true."
vobs=".true."
tobs=".true."
qobs=".true."
psobs=".true."
obstype="3"  
obsdensity="1.0d0"  
skip="1"
skipz="1"  

mkdir  -p $OBSDATAPATH
mkdir  -p $TMPDIR


#------- CREATE NAMELIST  ------------------

echo "&general                               " >  $NAMELIST_WRF2RADAR
echo "add_obs_error=$add_obs_error           " >> $NAMELIST_WRF2RADAR
echo "u_error=$u_error                       " >> $NAMELIST_WRF2RADAR
echo "v_error=$v_error                       " >> $NAMELIST_WRF2RADAR
echo "t_error=$t_error                       " >> $NAMELIST_WRF2RADAR
echo "q_error=$q_error                       " >> $NAMELIST_WRF2RADAR
echo "ps_error=$ps_error                     " >> $NAMELIST_WRF2RADAR
echo "uobs=$uobs                             " >> $NAMELIST_WRF2RADAR
echo "vobs=$vobs                             " >> $NAMELIST_WRF2RADAR
echo "tobs=$tobs                             " >> $NAMELIST_WRF2RADAR
echo "qobs=$qobs                             " >> $NAMELIST_WRF2RADAR
echo "psobs=$psobs                           " >> $NAMELIST_WRF2RADAR
echo "obstype=$obstype                       " >> $NAMELIST_WRF2RADAR
echo "obsdensity=$obsdensity                 " >> $NAMELIST_WRF2RADAR
echo "skip=$skip                             " >> $NAMELIST_WRF2RADAR
echo "skipz=$skipz                           " >> $NAMELIST_WRF2RADAR
echo "/                                      " >> $NAMELIST_WRF2RADAR

#---------END OF NAMELIST PARAMETERS

cd $TMPDIR

cdate=$INIDATE


while [ $cdate -le  $ENDDATE ]
do

MODELFILE=`wrfout_file_name $cdate 01`
ln -sf $MODELDATAPATH/$MODELFILE $TMPDIR/minput

echo "Processing file $MODELDATAPATH/$MODELFILE"
ln -sf $EXEC ./wrf_to_bufr.exe 
ln -sf $OBSDATAPATH/OBS${cdate}.dat ./out.dat      #Output file with OSSE observations for the current time.
ln -sf $OBSSOURCE/OBS${cdate}.dat ./obs.dat         #File containing obs distribution if obs distribution type is 3.

./wrf_to_bufr.exe

cdate=`date_edit2 $cdate $FREQ `
done






