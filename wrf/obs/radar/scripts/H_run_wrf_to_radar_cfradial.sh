#This script generates pseudo-radar observations from WRF outputs.
ulimit -s unlimited

MODELDATAPATH="$HOME/salidas/EXPERIMENTS/ANALYSIS_PARANA_2KM_control_paranacfsr_newobs_60m_radar_grib_Hydra/"
RADARDATAPATH="$HOME/share/DATA/OBS/OBS_REAL_PARANA_20091117_CFRADIAL/" 

TMPDIR="$HOME/TMP/WRF_TO_RADAR/"
EXEC="$HOME/share/LETKF_WRF/wrf/obs/radar/wrf_to_radar_cfradial/wrf_to_radar.exe"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH


source ../../../run/util.sh #Date functions.

MEMBER=00061 #Member number corresponding to the ensemble mean.

INIDATE=20091117180000
ENDDATE=20091117230000

DT_TIME_TRESH=300 #If no model data is found within DT_TIME_TRESH seconds, then radar file will not be processed.

NAMELIST_SCALE2RADAR="$TMPDIR/w2r.namelist"


add_obs_error=".FALSE."        
reflectivity_error="2.5d0"
radialwind_error="1.0d0"   
use_wt=".TRUE."                   #Include or not the effect of terminal velocity      
n_radar="1"                       #Number of radars                 
n_times="1"                       
method_ref_calc=2
input_type=1                      #1-restart files, 2-history files


mkdir  -p $RADARDATAPATH
mkdir  -p $TMPDIR

#------- CREATE NAMELIST PARAMETERS ------------------

echo "&general                               " >  $NAMELIST_SCALE2RADAR
echo "add_obs_error=$add_obs_error           " >> $NAMELIST_SCALE2RADAR
echo "reflectivity_error=$reflectivity_error " >> $NAMELIST_SCALE2RADAR
echo "radialwind_error=$radialwind_error     " >> $NAMELIST_SCALE2RADAR
echo "use_wt=$use_wt                         " >> $NAMELIST_SCALE2RADAR
echo "n_radar=$n_radar                       " >> $NAMELIST_SCALE2RADAR
echo "n_times=$n_times                       " >> $NAMELIST_SCALE2RADAR
echo "method_ref_calc=$method_ref_calc       " >> $NAMELIST_SCALE2RADAR
echo "/                                      " >> $NAMELIST_SCALE2RADAR

#---------END OF NAMELIST PARAMETERS

#Create folders and link files.

cdate=$INIDATE

itime=0001
irad=001

mkdir -p $MODELDATAPATH/wrf_to_radar/


#while [ $cdate -le  $ENDDATE ]
#do

cd $TMPDIR

#Delete links from previous round
rm -fr  $TMPDIR/radar* $TMPDIR/model*

for ifile in $( ls $RADARDATAPATH/cfrad* ) ; do

  myfile=$(basename $ifile)

  #cfrad.20100111_062345.000_to_20100111_062721.999_ANG120_v39_SUR.nc_dealiased

  #radardate1=`echo $myfile | cut -c7-14`
  #radardate2=`echo $myfile | cut -c16-21`
  radardate=`echo $myfile | cut -c7-14``echo $myfile | cut -c16-21`

  #Search for the closest model file. Loop over all availble model files to get the closest
  MINTDIST=999999

  find_date=0

  for imodeldate in $( ls $MODELDATAPATH/anal/ ) ; do

    datedist=`date_diff $imodeldate $radardate `

    datedist=`abs_val $datedist`

    if [ $datedist -lt $MINTDIST -a $datedist -lt $DT_TIME_TRESH ] ; then

     #We have a closer date.
     modeldate=$imodeldate
     MINTDIST=$datedist
     find_date=1

    fi

  done


  if [ $find_date -eq 1 ] ; then

    echo "Radar file for $cdate is $radarfile"

    ln -sf $EXEC ./wrf_to_radar.exe 

    #Model data won't be modified we can link it.  
    ln -sf $MODELDATAPATH/anal/${modeldate}/anal${MEMBER} ./model${itime}.nc

    #Radar data will be modified (additional fields will be written in the file)
    #So this file is copied to its final destionation and then linked to the TMPPATH
    cp $RADARDATAPATH/${myfile} $MODELDATAPATH/wrf_to_radar/$myfile
    ln -sf $MODELDATAPATH/wrf_to_radar/$myfile  ./radar${irad}_${itime}.nc

    #We execute the wrf_to_radar module.
    ./wrf_to_radar.exe

  else

    echo "Could not find a close analysis time for radar file $ifile "

  fi

done


