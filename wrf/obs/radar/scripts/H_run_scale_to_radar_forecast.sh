#This script generates pseudo-radar observations from WRF outputs.
ulimit -s unlimited

EXPERIMENT_NAME=SCALE_TO_RADAR_TEST

MODELDATAPATH="/home/paula.maldonado/scale_juan/argentina_anguil_2km_control/"
MODELTOPOPATH="/home/paula.maldonado/scale_juan/argentina_anguil_2km_control/const/topo"
RADARDATAPATH="/home/paula.maldonado/datosmate/RADAROBS/ANGUIL/20100111/120/dealiased_data/" 

TMPDIR="$HOME/TMP/SCALE_TO_RADAR/"
EXEC="$HOME/share/LETKF_WRF/wrf/obs/radar/scale_to_radar/scale_to_radar.exe"

export LD_LIBRARY_PATH=/usr/local/netcdf4.intel/lib/:$LD_LIBRARY_PATH

#DOMAIN="03"

source ../../../run/util.sh #Date functions.

INIDATE=20100111150000
ENDDATE=20100111150000

NAMELIST_SCALE2RADAR="$TMPDIR/s2r.namelist"

freq=1800                           #Window length in seconds.

model_output_freq=300               #Model output frequency

forecast_length=7200                #Forecast length

#Scale 2 radar parameters
add_obs_error=".FALSE."        
reflectivity_error="2.5d0"
radialwind_error="1.0d0"   
radarfile=""                
use_wt=".TRUE."                   #Include or not the effect of terminal velocity      
n_radar="1"                       #Number of radars                 
n_times="36"                      #Number of forecast outputs 
model_split_output=".false."
method_ref_calc=2
input_type=2                      #1-restart files , 2-history files

#Grid parameters 
MPRJ_basepoint_lon="296.02D0"
MPRJ_basepoint_lat="-36.5D0"
MPRJ_type="LC"
MPRJ_LC_lat1="-37.D0"
MPRJ_LC_lat2="-36.D0"

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
echo "input_type=$input_type                 " >> $NAMELIST_SCALE2RADAR
echo "/                                      " >> $NAMELIST_SCALE2RADAR

echo "&param_mapproj                         " >> $NAMELIST_SCALE2RADAR
echo "MPRJ_basepoint_lon=$MPRJ_basepoint_lon " >> $NAMELIST_SCALE2RADAR 
echo "MPRJ_basepoint_lat=$MPRJ_basepoint_lat " >> $NAMELIST_SCALE2RADAR
echo "MPRJ_type=$MPRJ_type                   " >> $NAMELIST_SCALE2RADAR
echo "MPRJ_LC_lat1=$MPRJ_LC_lat1             " >> $NAMELIST_SCALE2RADAR
echo "MPRJ_LC_lat2=$MPRJ_LC_lat2             " >> $NAMELIST_SCALE2RADAR
echo "/                                      " >> $NAMELIST_SCALE2RADAR

#---------END OF NAMELIST PARAMETERS

#Create folders and link files.

cdate=$INIDATE

while [ $cdate -le  $ENDDATE ] #Loop over initialization times
do
  
  cd $TMPDIR

  #Delete links from previous round
  rm -fr $TMPDIR/topo* $TMPDIR/model* $TMPDIR/radar*

  ln -sf $MODELTOPOPATH/topo*.nc .

  ln -sf $EXEC ./scale_to_radar.exe

  for ifile in $( ls $MODELDATAPATH/$cdate/fcst/mean/history* ) ; do

       myfile=$(basename $ifile)
       sufix=`echo $myfile | cut -c9-20`
       ln -sf $ifile ./model0001.${sufix}

  done


 
  forecast_start_date=$cdate
  forecast_end_date=`date_edit2 $cdate $forecast_length `

  fdate=$cdate

  itime=1

  irad=1

  while [ $fdate -le $forecast_end_date ] #Loop over forecast dates for the current forecast
  do

    #Search for the closest radar file. Loop over all availble radar files to get the closest
    MINTDIST=999999
    for ifile in $( ls $RADARDATAPATH/cfrad* ) ; do


      myfile=$(basename $ifile)

      filedate1=`echo $myfile | cut -c7-14`
      filedate2=`echo $myfile | cut -c16-21`
      filedate=${filedate1}${filedate2}

      datedist=`date_diff $filedate $cdate`

      datedist=`abs_val $datedist`

      if [ $datedist -lt $MINTDIST ] ; then

        #We have a closer date.
        radarfile=$ifile
        MINTDIST=$datedist

      fi

    done

    echo "Radar file for $cdate is $radarfile"
 
    irad=`add_zeros $irad 3 `
    itime=`add_zeros $itime 4 `

    ln -sf $radarfile ./radar${irad}_${itime}.nc

    fdate=`date_edit2 $fdate $forecast_out_freq `
  
    itime=`expr $itime + 1`

  done

./scale_to_radar.exe

cdate=`date_edit2 $cdate $freq `

done






