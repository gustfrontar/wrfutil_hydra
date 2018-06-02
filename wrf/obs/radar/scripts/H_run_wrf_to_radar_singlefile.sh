#This script generates pseudo-radar observations from WRF outputs.
ulimit -s unlimited

MODELFILE="$HOME/salidas/EXPERIMENTS/ANALYSIS_PARANA_2KM_control_paranacfsr_newobs_60m_radar_grib_Hydra/wrfout...."  #Completar con el nombre del file del modelo
RADARDATAPATH="$HOME/share/DATA/OBS/OBS_REAL_PARANA_20091117_CFRADIAL/cfrad..... "                                   #Completar con el nombre del file de radar 

TMPDIR="$HOME/TMP/WRF_TO_RADAR/"
EXEC="/home/jruiz/share/LETKF_WRF/wrf/obs/radar/wrf_to_radar_cfradial/wrf_to_radar.exe"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH


source /home/jruiz/share/LETKF_WRF/wrf/run/util.sh #Date functions.


#Creamos el namelist para el WRF to radar

NAMELIST_SCALE2RADAR="$TMPDIR/w2r.namelist"


add_obs_error=".FALSE."           #Permite agregar errores a las observaciones.
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



cd $TMPDIR

ln -sf $EXEC ./wrf_to_radar.exe 

#Model data won't be modified we can link it.  
ln -sf $MODELFILE ./model0001.nc

cp $RADARFILE  ./radar001_0001.nc

#We execute the wrf_to_radar module.
./wrf_to_radar.exe

	#Si el programa se ejecuto bien el cfradial deberia tener nuevas variables dBZ_model y V_model que 
#contienen la reflectividad y la velocidad radial calculados a partir del modelo.



