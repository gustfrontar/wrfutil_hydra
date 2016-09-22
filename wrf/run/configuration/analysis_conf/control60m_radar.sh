#KALMAN FILTER CONFIGURATION
MEMBER=60        #Number of ensemble members.
MAX_DOM=1
HOMEDIR=${HOME}/share/
DATADIR=${HOME}/data/
ANALYSIS=1       #Identify this job as an analysis job.
FORECAST=0       #This is not a forecast job.
INTERPANA=0      #This is used in forecast jobs (but we need to define it here too)
RUN_ONLY_MEAN=0  #This is used in forecast jobs (but we need to define it here too)

USE_ANALYSIS_BC=1 #1 - use analysis as BC , 0 - use forecasts as bc (e.g. global gfs)
                  # if 1 then bc data will be taken from exp_met_em folder in the corresponding INPUT folder.
                  # if 0 then bc data will be taken from for_met_em folder in the corresponding INPUT folder.
                  # default is 1
USE_ANALYSIS_IC=0 #1 - use global analysis as IC, 0 use LETKF analysis as IC
                  #if 0 then profide a LETKF-analysis source (ANALYSIS_SOURC)
                  #default is 0

NVERTEXP=27  #used in forecast and da experiments.
NVERTDB=38   #used for verification.

#AUXILIARY VARIABLE FOR ENSEMBLE SIZE
MM=$MEMBER                      #Variable for iteration limits.
MEANMEMBER=`expr $MEMBER + 1 `  #This is the member ID corresponding to the ensemble mean.

WINDOW=300        #Assimilation frequency. (seconds)
WINDOW_START=300  #Window start (seconds from forecast initialization)     
WINDOW_END=300    #Window end   (seconds from forecast initialization)
WINDOW_FREC=300   #Output frequency within window (seconds) should be the same as the maximum observation frequency.
ASSIMILATION_FREC=300 #Assimilation frequency  (seconds)
NSLOTS=`expr $WINDOW_END \/ $WINDOW_FREC - $WINDOW_START \/ $WINDOW_FREC  + 1 `        #Number of time slots. 
NBSLOT=`expr $ASSIMILATION_FREC \/ $WINDOW_FREC - $WINDOW_START \/ $WINDOW_FREC + 1 `  #Time slot corresponding to the analysis.
if [ $NBSLOT -lt 10 ] ; then
   NBSLOT=0$NBSLOT
fi
SIGMA_OBS="2.0d3"
SIGMA_OBSV="0.2d0"
SIGMA_OBSZ="2.0d3" 
SIGMA_OBST="3.0d0"
GROSS_ERROR="15.0d0" 
COV_INFL_MUL="1.1d0"
SP_INFL_ADD="0.d0"  
RELAX_ALPHA_SPREAD="0.90d0"
RELAX_ALPHA="0.0d0" 
USE_ADAPTIVE_INFLATION=0  #1 turn on addaptive inflation (Miyoshi 2011), 0 Turn off adaptive inflation
				  #Note that for addaptive inflation to work COV_INFL_MUL should be < 0.0
GUESFT=$WINDOW_END  # First guess forecast length (seconds)

#DOMAIN AND BOUNDARY DATA

BOUNDARY_DATA_FREC=21600              #Boundary data frequency. (seconds)
BOUNDARY_DATA_PERTURBATION_FREQ=21600 #Frequency of data used to perturb boundary conditions (seconds)

#INITIAL AND BOUNDARY PERTURBATIONS
SCALE_FACTOR="0.05"         #Perturbation scale factor.
RANDOM_SCALE_FACTOR="0.005" #Random perturbation scale factor.
PERTURB_BOUNDARY=1          #DUMMY
PERTURB_BOUNDARY_TYPE=1     #DUMMY

#POSTPROC CONFIGURATION
OUTLEVS="0.1,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,"
OUTVARS="'umet,vmet,W,QVAPOR,QCLOUD,QRAIN,QICE,QSNOW,QGRAUP,RAINNC,tk,u10m,v10m,slp,mcape,dbz,max_dbz'"
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
OBS=""                # Name of observation folder.
RADAROBS="/OBS_0.5ms_2db_wt_noatt/"              # Name of radar observation folder.
EXP=ANALYSIS_${DOMAINCONF}_${CONFIGURATION}      # name of experiment

### initial date setting
IDATE=20140122173000
EDATE=20140122200000

#### DATA
OBSDIR=${HOMEDIR}/OBS/$OBS/                                                               # observations
NRADARS=1
RADAROBSDIR=${HOMEDIR}/OBS/$RADAROBS/
TMPDIR=${DATADIR}/TMP/$EXP/                                                               # work directory
OUTPUTDIR=${DATADIR}/EXPERIMENTS/$EXP/                                                    # Where results should appear.
INPUTDIR=${HOMEDIR}/INPUT/$DOMAINCONF/                                                    # This folder contains text files with the dates to compute perturbations.

#### EXECUTABLES
RUNTIMELIBS=${HOMEDIR}/libs_sparc64/lib/
WRF=${HOMEDIR}/LETKF_WRF/wrf/
LETKF=$WRF/letkf/letkf.exe                     # LETKF module
UPDATEBC=$WRF/model/WRFDA/da_update_bc.exe     
WRFMODEL=$WRF/model/WRFV3/                     # WRF model that run in computing nodes.
WRFMODELPPS=$WRF/model/WRFV3/                  # WRF model that runs in pps server 
ARWPOST=$WRF/model/ARWpost/
SPAWN=$WRF/spawn/
MPIBIN=mpiexec

#### SCRIPTS
UTIL=$WRF/run/util.sh


#### NAMELIST
NAMELISTWRF=$WRF/run/configuration/$DOMAINCONF/namelist.input
NAMELISTLETKF=$WRF/run/configuration/letkf.namelist.$LETKFNAMELIST
NAMELISTARWPOST=$WRF/run/configuration/$DOMAINCONF/namelist.ARWpost
NAMELISTOBSOPE=$WRF/run/configuration/obsope.namelist.$OBSOPENAMELIST

