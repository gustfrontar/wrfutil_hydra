#KALMAN FILTER CONFIGURATION
DOMAINCONF=PARANA_10KM                  #Define a domain
LETKFNAMELIST=control                   #Define a letkf namelist template

MEMBER=60        #Number of ensemble members.
MAX_DOM=1        #Maximum number of WRF domains.
HOMEDIR=${HOME}/share/   
DATADIR=${HOME}/salidas/
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

NVERTEXP=27  #Number of vertical levels in initial and boundary conditions input grib data.
NVERTDB=38   #Number of vertical levels in initial and boundary conditions perturbation input grib data.

#AUXILIARY VARIABLE FOR ENSEMBLE SIZE
MM=$MEMBER                      #Variable for iteration limits.
MEANMEMBER=`expr $MEMBER + 1 `  #This is the member ID corresponding to the ensemble mean.

WINDOW=10800        #Assimilation frequency. (seconds)
WINDOW_START=3600   #Window start (seconds from forecast initialization)     
WINDOW_END=10800    #Window end   (seconds from forecast initialization)
WINDOW_FREC=3600    #Output frequency within window (seconds) should be the same as the maximum observation frequency.
ASSIMILATION_FREC=10800 #Assimilation frequency  (seconds)
NSLOTS=`expr $WINDOW_END \/ $WINDOW_FREC - $WINDOW_START \/ $WINDOW_FREC  + 1 `        #Number of time slots. 
NBSLOT=`expr $ASSIMILATION_FREC \/ $WINDOW_FREC - $WINDOW_START \/ $WINDOW_FREC + 1 `  #Time slot corresponding to the analysis.
if [ $NBSLOT -lt 10 ] ; then
   NBSLOT=0$NBSLOT
fi
SIGMA_OBS="70.0d3"
SIGMA_OBSV="0.2d0"
SIGMA_OBSZ="5.0d3" 
SIGMA_OBST="3.0d0"
GROSS_ERROR="4.0d0" 
COV_INFL_MUL="1.1d0"
SP_INFL_ADD="0.d0"  
RELAX_ALPHA_SPREAD="0.9d0"
RELAX_ALPHA="0.0d0" 
USE_ADAPTIVE_INFLATION=0  #1 turn on addaptive inflation (Miyoshi 2011), 0 Turn off adaptive inflation
				  #Note that for addaptive inflation to work COV_INFL_MUL should be < 0.0
GUESFT=$WINDOW_END  # First guess forecast length (seconds)

#DOMAIN AND BOUNDARY DATA

BOUNDARY_DATA_FREQ=21600              #Boundary data frequency. (seconds)
BOUNDARY_DATA_PERTURBATION_FREQ=21600 #Frequency of data used to perturb boundary conditions (seconds)

#POSTPROC CONFIGURATION
OUTLEVS="0.1,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,"      #Level list
OUTVARS="'umet,vmet,W,QVAPOR,QCLOUD,QRAIN,QICE,QSNOW,QGRAUP,RAINNC,tk,u10m,v10m,slp,mcape,dbz,max_dbz'"  #Variable list.
ARWPOST_FREC=21600   # Post processing frequency (seconds)
INPUT_ROOT_NAME='wrfout'
INTERP_METHOD=1
ENABLE_UPP=1           #1 - generate grib sigma files (for da downscalling) , 0 - do not generate grib sigma files.

### LETKF setting
OBS="PREPBUFR_AND_SURFACE"                                 # Name of conventional observations folder.
RADAROBS="OBS_REAL_PARANA_20091117_SO15KM_NEWQC"           # Name of radar observation folder.
EXP=ANALYSIS_${DOMAINCONF}_${CONFIGURATION}                # name of experiment

### initial date setting
IDATE=20091117000000
EDATE=20091118000000     #20091117230000

#### DATA
OBSDIR=${HOMEDIR}/DATA/OBS/$OBS/                                                          # Folder where conventional observations are.
NRADARS=1                                                                                 # Number of available radars.
RADAROBSDIR=${HOMEDIR}/DATA/OBS/$RADAROBS/                                                # Folder where radar observations are.
TMPDIR=${HOME}/TMP/$EXP/                                                                  # Temporal work directory (should be accessible for all computation nodes)
OUTPUTDIR=${DATADIR}/EXPERIMENTS/$EXP/                                                    # Where results will be stored.
GRIBDIR=${HOMEDIR}/DATA/GRIB/FNL/HIRES/SA/                                                # Folder where bdy and initial grib files are located.
GRIBTABLE="Vtable.GFS"                                                                    # Bdy and init data source Vtable name.
MEMBER_BDY=1                                                                              # Total number of boundary conditions ensemble members.
MEANMEMBER_BDY=1                                                                          # Boundary conditions ensemble member corresponding to the ensemble mean.

PERTGRIBDIR=${HOMEDIR}/DATA/GRIB/CFSR/HIRES/ARGENTINA/                                    # Folder where data for perturbing bdy are located.
PERTGRIBTABLE="Vtable.CFSR2_web"                                                          # Bdy perturbation source vtable name.
GEOG=${HOMEDIR}/LETKF_WRF/wrf/model/GEOG/                                                 # Folder where WPS GEOG dataset is located.

#INITIAL AND BOUNDARY PERTURBATIONS
AMP_FACTOR="0.1"             #Perturbation scale factor.
RANDOM_AMP_FACTOR="0.0"       #Random perturbation scale factor.
PERTURB_BOUNDARY=1            #Wether boundary conditions are going to be perturbed.
PERTURB_ATMOSPHERE=".true."   #Wether atmospheric conditions will be perturbed (boundary and first cycle)
PERTURB_SST=".true."          #Wether SST will be perturbed.
PERTURB_SOIL=".true."         #Wether SOIL conditions will be perturbed (soil moisture and soil temperature)
PERTURB_T=".true."            #Wether ATM temperature will be perturbed.
PERTURB_RH=".true."           #Wether ATM RH will be perturbed
PERTURB_WIND=".true."         #Wether ATM winds will be perturbed.
PERTURB_T_AMP="0.5d0"         #T random perturbation amplitude
PERTURB_RH_AMP="5.0d0"        #RH random perturbation amplitude 
PERTURB_WIND_AMP="0.5d0"      #WIND random perturbation amplitude
PERTURB_T_SCLH="40000d0"      #T random perturbation horizontal scale
PERTURB_RH_SCLH="40000d0"     #RH random perturbation horizontal scale
PERTURB_WIND_SCLH="40000d0"   #WIND random perturbation horizontal scale
PERTURB_T_SCLV="5000d0"       #T random perturbation vertical scale
PERTURB_RH_SCLV="5000d0"      #RH random perturbation vertical scale
PERTURB_WIND_SCLV="5000d0"    #WIND random perturbation vertical scale

#Random dates for boundary perturbations.
INIPERTDATE=20060101000000    #Initial date in grib database (used for perturbing initial and boundary conditions)
ENDPERTDATE=20091231180000    #Final date in grib database (used for perturbing initial and boundary conditions)
PERTREFDATE=20091117000000    #At this date the initial perturbation dates will be taken. This date is used to keep consisntency among the perturbations
                              #used in forecast and analysis experiments. This date must be previous or equal to IDATE.
INPUT_PERT_DATES_FROM_FILE=0  #0 - generate a new set of random dates, 1 - read random dates from a file. 
INI_PERT_DATE_FILE=${HOMEDIR}/DATA/INITIAL_RANDOM_DATES/initial_perturbation_dates_60m  #List of initial random dates.

#### EXECUTABLES
RUNTIMELIBS=${HOMEDIR}/libs_sparc64/lib/       #Libs that will be included in LD_LIBRARY_PATH in computing nodes.
WRF=${HOMEDIR}/LETKF_WRF/wrf/                  # WRF folder (for computing nodes)
LETKF=$WRF/letkf/letkf.exe                     # LETKF module (for computing nodes)
UPDATEBC=$WRF/model/WRFDA/da_update_bc.exe     # Update bc tool (WRFVAR) (for computing nodes)
WRFMODEL=$WRF/model/WRFV3.6/                   # WRF model that run in computing nodes.
WRFMODELPPS=$WRF/model/WRFV3.6/                # WRF model that runs in pps server  (usually the same as the one for the computing nodes)
WPS=$WRF/model/WPS3.6/                         # WRF model pre processing utilities (for pps server)
ARWPOST=$WRF/model/ARWpost/                    # WRF model post processing utilities that run in computing nodes.
UPP=$WRF/model/UPPV3.0/
SPAWN=$WRF/spawn/
MPIBIN=mpiexec

#### SCRIPTS
UTIL=$WRF/run/util.sh                          # Script containing bash functions that will be used during execution.

#### NAMELIST
NAMELISTWRF=$WRF/run/configuration/domain_conf/$DOMAINCONF/namelist.input           #Namelist for WRF model.
NAMELISTWPS=$WRF/run/configuration/domain_conf/$DOMAINCONF/namelist.wps             #Namelist for WRF pre processing tools
NAMELISTLETKF=$WRF/run/configuration/letkf_conf/letkf.namelist.$LETKFNAMELIST       #Namelist for LETKF
NAMELISTARWPOST=$WRF/run/configuration/domain_conf/$DOMAINCONF/namelist.ARWpost     #Namelist for post-processing tools.
NAMELISTUPP=$WRF/run/configuration/domain_conf/$DOMAINCONF/wrf_cntrl.parm           #Namelist for grib post-processing.
NAMELISTOBSOPE=$WRF/run/configuration/letkf_conf/obsope.namelist.$OBSOPENAMELIST    #Namelist for observation operator.
NAMELISTPERTMETEM=$WRF/run/configuration/letkf_conf/pertmetem.namelist.control      #Namelist for boundary conditions perturbation.

