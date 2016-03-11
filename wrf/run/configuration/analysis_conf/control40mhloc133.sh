#KALMAN FILTER CONFIGURATION
MEMBER=40           #Number of ensemble members.
MAX_DOM=1
HOMEDIR=${HOME}/share/
DATADIR=${HOME}/data/
ANALYSIS=1   #Identify this job as an analysis job.
FORECAST=0   #This is not a forecast job.
INTERPANA=0  #This is used in forecast jobs (but we need to define it here too)
RUN_ONLY_MEAN=0 #This is used in forecast jobs (but we need to define it here too)

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

WINDOW=21600        # Assimilation frequency. (seconds)
WINDOW_START=10800  #Window start (seconds from forecast initialization)     
WINDOW_END=32400    #Window end   (seconds from forecast initialization)
WINDOW_FREC=3600    #Output frequency within window (seconds) should be the same as the maximum observation frequency.
ASSIMILATION_FREC=21600 #Assimilation frequency  (seconds)
NSLOTS=`expr $WINDOW_END \/ $WINDOW_FREC - $WINDOW_START \/ $WINDOW_FREC  + 1 `        #Number of time slots. 
NBSLOT=`expr $ASSIMILATION_FREC \/ $WINDOW_FREC - $WINDOW_START \/ $WINDOW_FREC + 1 `  #Time slot corresponding to the analysis.
if [ $NBSLOT -lt 10 ] ; then
   NBSLOT=0$NBSLOT
fi
FIRSTSLOT=`expr $WINDOW_START \/ $WINDOW_FREC `  #Time slot corresponding to the analysis.
SIGMA_OBS="1.33d5"
SIGMA_OBSV="0.2d0"
SIGMA_OBSZ="6.0d3" 
SIGMA_OBST="3.0d0"
GROSS_ERROR="3.0d0" 
COV_INFL_MUL="1.0d0"
SP_INFL_ADD="0.d0"  
RELAX_ALPHA_SPREAD="0.8d0"
RELAX_ALPHA="0.0d0" 
GUESFT=$WINDOW_END  # First guess forecast length (seconds)

#DOMAIN AND BOUNDARY DATA

BOUNDARY_DATA_FREC=21600              #Boundary data frequency. (seconds)
BOUNDARY_DATA_PERTURBATION_FREQ=21600 #Frequency of data used to perturb boundary conditions (seconds)

#INITIAL AND BOUNDARY PERTURBATIONS
SCALE_FACTOR="0.1"         #Perturbation scale factor.
RANDOM_SCALE_FACTOR="0.0"  #Random perturbation scale factor.
PERTURB_BOUNDARY=1      #DUMMY
PERTURB_BOUNDARY_TYPE=1 #DUMMY

#POSTPROC CONFIGURATION
OUTLEVS="1000.,950.,900.,850.,800.,750.,700.,650.,600.,500.,400.,300.,200.,100.,"
OUTVARS="'umet,vmet,W,QVAPOR,q2,RAINC,RAINNC,SST,geopt,tk,u10m,v10m,slp,mcape'"
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
OBS="/ucar_airs_th3x3_corrected/"                # Name of observation folder.
EXP=ANALYSIS_${DOMAINCONF}_${CONFIGURATION}      # name of experiment

### initial date setting
IDATE=20080807000000
EDATE=20080930180000

#### DATA
OBSDIR=${HOMEDIR}/OBS/$OBS/                                                               # observations
TMPDIR=${DATADIR}/TMP/$EXP/                                                               # work directory
OUTPUTDIR=${DATADIR}/EXPERIMENTS/$EXP/                                                    # Where results should appear.
INPUTDIR=${HOMEDIR}/INPUT/$DOMAINCONF/                                                    # This folder contains text files with the dates to compute perturbations.

#### EXECUTABLES
RUNTIMELIBS=${HOMEDIR}/libs_sparc64/lib/
WRF=${HOMEDIR}/LETKF_WRF/wrf/
LETKF=$WRF/letkf/letkf_noeigenexa.exe                     # LETKF module
UPDATEBC=$WRF/model/WRFDA/da_update_bc.exe
WRFMODEL=$WRF/model/WRFV3K/ 
WRFMODELPPS=$WRF/model/WRFV3INTEL/                          # WRF model that runs in pps server
ARWPOST=$WRF/model/ARWpostINTEL/
SPAWN=$WRF/spawn/
MPIBIN=mpiexec

#### SCRIPTS
UTIL=$WRF/run/util.sh


#### NAMELIST
NAMELISTWRF=$WRF/run/configuration/$DOMAINCONF/namelist.input
NAMELISTLETKF=$WRF/run/configuration/letkf.namelist.$LETKFNAMELIST
NAMELISTARWPOST=$WRF/run/configuration/$DOMAINCONF/namelist.ARWpost
NAMELISTOBSOPE=$WRF/run/configuration/obsope.namelist.$OBSOPENAMELIST


