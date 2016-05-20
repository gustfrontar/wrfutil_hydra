#This script generates pseudo-radar observations from WRF outputs.
ulimit -s unlimited

EXPERIMENT_NAME=OSSE_OBS_WERROR

MODELDATAPATH="$HOME/datos/EXPERIMENTS/NATURERUN_CORDOBA_2K_cordoba_naturerun_1mtest/forecast/20140122120000/00001/"
RADARDATAPATH="$HOME/share/OBS/$EXPERIMENT_NAME"   

TMPDIR="$HOME/data/TMP/WRF_TO_RADAR/"
EXEC="$HOME/share/LETKF_WRF/wrf/obs/radar/wrf_to_radar/wrf_to_radar.exe"

RADARFILE=""

source ../../../run/util.sh #Date functions.

INIDATE=20140122120000
ENDDATE=20140122120100

NAMELIST_WRF2RADAR="$TMPDIR/w2r.namelist"

frec=30                           #Observation frequency in seconds.
add_obs_error=".FALSE."        
reflectivity_error="2.5d0"
radialwind_error="1.0d0"   
fake_radar=".TRUE."               #If false provide a radarfile. 
radarfile=""                
fradar_lon="-65.0d0"               
fradar_lat="-35.0d0"  
fradar_z="0.0d0"                  
radar_az_res="1.0d0"            
radar_r_res="500.0d0"              
radar_el_res="0.5d0"           
radar_min_az="1.0d0"           
radar_max_az="360.0d0"           
radar_min_r="500.0d0"   
radar_max_r="240.0d3"
radar_min_el="0.0d0"
radar_max_el="30.0d0"
use_wt=".TRUE."                   #Include or not the effect of terminal velocity      
radar_lambda="10.0d0"            
use_level_list=".FALSE."       
n_vertical_levels="2"  
level_list="1.0d0,2.0d0"
n_radar="1"                       #Number of radars                 


mkdir  -p $RADARDATAPATH
mkdir  -p $TMPDIR


#------- CREATE NAMELIST PARAMETERS ------------------

echo "&general                               " >  $NAMELIST_WRF2RADAR
echo "add_obs_error=$add_obs_error           " >> $NAMELIST_WRF2RADAR
echo "reflectivity_error=$reflectivity_error " >> $NAMELIST_WRF2RADAR
echo "radialwind_error=$radialwind_error     " >> $NAMELIST_WRF2RADAR
echo "fake_radar=$fake_radar                 " >> $NAMELIST_WRF2RADAR
echo "fradar_lon=$fradar_lon                 " >> $NAMELIST_WRF2RADAR
echo "fradar_lat=$fradar_lat                 " >> $NAMELIST_WRF2RADAR
echo "fradar_z=$fradar_z                     " >> $NAMELIST_WRF2RADAR
echo "radar_az_res=$radar_az_res             " >> $NAMELIST_WRF2RADAR
echo "radar_r_res=$radar_r_res               " >> $NAMELIST_WRF2RADAR
echo "radar_el_res=$radar_el_res             " >> $NAMELIST_WRF2RADAR
echo "radar_min_az=$radar_min_az             " >> $NAMELIST_WRF2RADAR
echo "radar_max_az=$radar_max_az             " >> $NAMELIST_WRF2RADAR
echo "radar_min_r=$radar_min_r               " >> $NAMELIST_WRF2RADAR
echo "radar_max_r=$radar_max_r               " >> $NAMELIST_WRF2RADAR
echo "radar_min_el=$radar_min_el             " >> $NAMELIST_WRF2RADAR
echo "radar_max_el=$radar_max_el             " >> $NAMELIST_WRF2RADAR
echo "use_wt=$use_wt                         " >> $NAMELIST_WRF2RADAR
echo "radar_lambda=$radar_lambda             " >> $NAMELIST_WRF2RADAR
echo "use_level_list=$use_level_list         " >> $NAMELIST_WRF2RADAR
echo "n_vertical_levels=$n_vertical_levels   " >> $NAMELIST_WRF2RADAR
echo "level_list=$level_list                 " >> $NAMELIST_WRF2RADAR
echo "n_radar=$n_radar                       " >> $NAMELIST_WRF2RADAR
echo "n_model=1                              " >> $NAMELIST_WRF2RADAR
echo "/                                      " >> $NAMELIST_WRF2RADAR

#---------END OF NAMELIST PARAMETERS

#Create folders and link files.

cdate=$INIDATE

itime=1

while [ $cdate -le  $ENDDATE ]
do

cd $TMPDIR


MODELFILE=`wrfout_file_name $cdate 01`
itime=`add_zeros $itime 4 `


ln -sf $MODELDATAPATH/$MODELFILE $TMPDIR/input_model${itime}.nc
ln -sf $EXEC ./wrf_to_radar.exe 

./wrf_to_radar.exe

  #Move the data.
  iradar=1
  while [ $iradar -le $n_radar ] 
  do
    iradar=`add_zeros $iradar 4`
    mv oradar${iradar}_0001.grd $RADARDATAPATH/RADAR_${cdate}.grd
    iradar=`expr $iradar + 1 `
  done


cdate=`date_edit2 $cdate $frec `
itime=`expr $itime + 1 `

done






