#KALMAN FILTER CONFIGURATION
MEMBER=40
HOMEDIR=${HOME}/share/
DATADIR=${HOME}/data/
INTERPANA=0 #1 - If input analysis has to be interpolated to the current domain. (If not we will assume that Analysis source and input domains have the same grid)
FORECAST=1  #Identify this as a forecast job.
ANALYSIS=0  #0 means this is not an analysis cycle run.
RUN_ONLY_MEAN=0 #1 - run only the forecast initialized with the ensemble member, 0 run the full ensemble.
MAX_DOM=1
USE_ANALYSIS_BC=0 #1 - use analysis as BC , 0 - use forecasts as bc (e.g. global gfs)
                  # if 1 then bc data will be taken from exp_met_em folder in the corresponding INPUT folder.
                  # if 0 then bc data will be taken from for_met_em folder in the corresponding INPUT folder.
                  # default is 1
USE_ANALYSIS_IC=0 #1 - use global analysis as IC, 0 use LETKF analysis as IC
                  #if 0 then profide a LETKF-analysis source (ANALYSIS_SOURCE)
                  #default is 0

ANALYSIS_SOURCE=${HOME}/data/EXPERIMENTS/ANALYSIS_SINLAKU_60K_control40m/

NVERTEXP=27  #used in forecast and da experiments.
NVERTDB=38   #used for verification.


#AUXILIARY VARIABLE FOR ENSEMBLE SIZE
MM=$MEMBER                      #Variable for iteration limits.
MEANMEMBER=`expr $MEMBER + 1 `  #This is the member ID corresponding to the ensemble mean.

ASSIMILATION_FREC=21600 #Forecast initialization frequency (seconds)
GUESFT=259200           #Forecast length (secons)

WINDOW=21600        #Forecast initialization frequency (seconds)
WINDOW_START=0      #Window start (seconds from forecast initialization)
WINDOW_END=$GUESFT  #Window end   (seconds from forecast initialization)
WINDOW_FREC=3600    #Output frequency for the forecast


#OBSERVATION OPERATOR CONFIGURATION FOR FORECAST VERIFICATION.
SIGMA_OBS="0.0d0"            #NOT USED
SIGMA_OBSV="0.0d0"           #NOT USED
SIGMA_OBSZ="0.0d0"           #NOT USED
SIGMA_OBST="0.0d0"           #NOT USED
COV_INFL_MUL="0.0d0"         #NOT USED
SP_INFL_ADD="0.d0"           #NOT USED
RELAX_ALPHA_SPREAD="0.0d0"   #NOT USED
RELAX_ALPHA="0.0d0"          #NOT USED

#DOMAIN AND BOUNDARY DATA

BOUNDARY_DATA_FREC=10800              #Boundary data frequency. (seconds)
BOUNDARY_DATA_PERTURBATION_FREQ=21600 #Frequency of data used to perturb boundary conditions (seconds)

#INITIAL AND BOUNDARY PERTURBATIONS
SCALE_FACTOR="0.1"         #Perturbation scale factor.(balanced random)
RANDOM_SCALE_FACTOR="0.0"  #Random perturbation scale factor.
PERTURB_BOUNDARY=1         #If boundary perturbations will be applied. 
PERTURB_BOUNDARY_TYPE=1    #DUMMY

#POSTPROC CONFIGURATION
OUTLEVS="1000.,975.,950.,925.,900.,850.,800.,700.,600.,500.,400.,300.,250.,200.,150.,100.,"
OUTVARS="'umet,vmet,W,QVAPOR,PSFC,RAINC,RAINNC,HGT,TSK,SST,geopt,tk,rh,rh2,u10m,v10m,slp,mcape'"
ARWPOST_FREC=21600   # Post processing frequency (seconds)
INPUT_ROOT_NAME='wrfout'
INTERP_METHOD=1


#Compute DUMMY increment for boundary data preparation.
rest=`expr  $GUESFT % $BOUNDARY_DATA_FREC `   
if [ $rest == 0 ] ; then
DINC=$GUESFT
else
DINC=`expr $GUESFT + $BOUNDARY_DATA_FREC - $rest `
fi


### LETKF setting
EXP=FORECAST_${DOMAINCONF}_${CONFIGURATION}      # name of experiment

### initial date setting
IDATE=20080820000000 
EDATE=20080927000000

#### DATA

TMPDIR=${DATADIR}/TMP/$EXP/                                                            # work directory
OUTPUTDIR=${DATADIR}/EXPERIMENTS/$EXP/                                                 # Where results should appear.
INPUTDIR=${HOMEDIR}/INPUT/$DOMAINCONF/
GEOG=${HOME}/share/GEOG/
OBS="/ucar_airs_th3x3_corrected/"                # Name of observation folder.
OBSDIR=${HOMEDIR}/OBS/$OBS/                                                               # observations

###### EXECUTABLES
WRF=${HOMEDIR}/LETKF_WRF/wrf/
LETKF=$WRF/letkf/letkf.exe                     # LETKF module
RUNTIMELIBS=${HOMEDIR}/libs_sparc64/lib/
WRF=${HOMEDIR}/LETKF_WRF/wrf/
UPDATEBC=$WRF/model/WRFDA/da_update_bc.exe
WRFMODEL=$WRF/model/WRFV3K/ 
WRFMODELPPS=$WRF/model/WRFV3INTEL/
WPS=$WRF/model/WPSINTEL/
ARWPOST=$WRF/model/ARWpostINTEL/
WRF_TO_WPS=$WRF/wrf_to_wps/
SPAWN=$WRF/spawn/
MPIBIN=mpiexec

#### SCRIPTS
UTIL=$WRF/run/util.sh

#### NAMELIST
NAMELISTWRF=$WRF/run/configuration/$DOMAINCONF/namelist.input
NAMELISTWPS=$WRF/run/configuration/$DOMAINCONF/namelist.wps
NAMELISTARWPOST=$WRF/run/configuration/$DOMAINCONF/namelist.ARWpost
NAMELISTLETKF=$WRF/run/configuration/letkf.namelist.$LETKFNAMELIST
NAMELISTOBSOPE=$WRF/run/configuration/obsope.namelist.$OBSOPENAMELIST
