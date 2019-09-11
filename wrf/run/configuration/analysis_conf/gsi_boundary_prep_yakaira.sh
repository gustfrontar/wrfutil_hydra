#KALMAN FILTER CONFIGURATION
DOMAINCONF=PARANA_10KM                #Define a domain
LETKFNAMELIST=control
MEMBER=30   #Number of ensemble members.
MAX_DOM=1
HOMEDIR=/home/jruiz/share/
DATADIR=/home/jruiz/share/
ANALYSIS=1       #Identify this job as an analysis job.
FORECAST=0       #This is not a forecast job.
INTERPANA=0      #This is used in forecast jobs (but we need to define it here too)
RUN_ONLY_MEAN=0  #This is used in forecast jobs (but we need to define it here too)
RUN_CHEM=0       #1 - run wrf chem , 0 run wrf


USE_ANALYSIS_BC=1 #1 - use analysis as BC , 0 - use forecasts as bc (e.g. global gfs)
                  # if 1 then bc data will be taken from exp_met_em folder in the corresponding INPUT folder.
                  # if 0 then bc data will be taken from for_met_em folder in the corresponding INPUT folder.
                  # default is 1
USE_ANALYSIS_IC=1 #1 - use global analysis as IC, 0 use LETKF analysis as IC
                  #if 0 then profide a LETKF-analysis source (ANALYSIS_SOURC)
                  #default is 0

NVERTEXP=27  #used in forecast and da experiments.
NVERTDB=38   #used for verification.

#AUXILIARY VARIABLE FOR ENSEMBLE SIZE
MM=$MEMBER                      #Variable for iteration limits.
MEANMEMBER=`expr $MEMBER + 1 `  #This is the member ID corresponding to the ensemble mean.

ASSIMILATION_FREC=21600  #Forecast initialization frequency (seconds)
GUESFT=21600             #Forecast length (seconds)

WINDOW=21600        #Forecast initialization frequency (seconds)
WINDOW_START=0      #Window start (seconds from forecast initialization)
WINDOW_END=$GUESFT  #Window end   (seconds from forecast initialization)
WINDOW_FREC=21600   #Output frequency for the forecast

#OBSERVATION OPERATOR CONFIGURATION FOR FORECAST VERIFICATION.
SIGMA_OBS="0.0d0"            #NOT USED
SIGMA_OBSV="0.0d0"           #NOT USED
SIGMA_OBSZ="0.0d0"           #NOT USED
SIGMA_OBST="0.0d0"           #NOT USED
COV_INFL_MUL="0.0d0"         #NOT USED
SP_INFL_ADD="0.0d0"          #NOT USED
RELAX_ALPHA_SPREAD="0.0d0"   #NOT USED
RELAX_ALPHA="0.0d0"          #NOT USED


#DOMAIN AND BOUNDARY DATA
MEMBER_BDY=1                          #Global boundary conditions ensemble size.
MEANMEMBER_BDY=1
BOUNDARY_DATA_FREQ=21600              #Boundary data frequency. (seconds)
BOUNDARY_DATA_PERTURBATION_FREQ=21600 #Frequency of data used to perturb boundary conditions (seconds)

#INITIAL AND BOUNDARY PERTURBATIONS
AMP_FACTOR="0.1"              #Balanced perturbation scale factor.
RANDOM_AMP_FACTOR="0.0"       #Umbalanced Random perturbation scale factor.
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

#POSTPROC CONFIGURATION
OUTLEVS="0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0"
OUTVARS="'umet,vmet,QVAPOR,QCLOUD,QRAIN,QICE,QSNOW,QGRAUP,RAINNC,tk,t2m,pressure,u10m,v10m,slp'"
ARWPOST_FREC=21600   # Post processing frequency (seconds)
INPUT_ROOT_NAME='wrfout'
INTERP_METHOD=1

### LETKF setting
OBS=""                                                     # Name of observation folder.
RADAROBS="/OSSE_20140122_DBZ2.5_VR1.0_SO2KM/"              # Name of radar observation folder.
EXP=SCALE_BDY_${DOMAINCONF}_${CONFIGURATION}               # name of experiment

### initial date setting
IDATE=20140124120000     #EXPERIMENT INITIAL DATE
EDATE=20140124120000     #EXPERIMENT END DATE

#### DATA
OBSDIR=${DATADIR}/OBS/$OBS/                                # observations
NRADARS=1
RADAROBSDIR=${DATADIR}/DATA/OBS/$RADAROBS/
TMPDIR=${DATADIR}/TMP/$EXP/                                                             # work directory
OUTPUTDIR=${DATADIR}/input_data/wrf_files/$EXP/                                         # Where results should appear.
GRIBDIR=${HOMEDIR}/DATA/GRIB/FNL/HIRES/ARGENTINA/                # Folder where bdy and inita data gribs are located.
#GRIBDIRSFC=${HOMEDIR}/DATA/GRIB/CFSR/HIRES/ARGENTINASFC/
GRIBTABLE="Vtable.GFS"                                                               # Bdy and init data source Vtable name.
PERTGRIBDIR=${HOMEDIR}/DATA/GRIB/CFSR/HIRES/ARGENTINA/00001/                     # Folder where data for perturbing bdy are located.
PERTGRIBTABLE="Vtable.CFSR2_web"                                                        # Bdy perturbation source vtable name.
GEOG=/home/jruiz/share/LETKF_WRF/wrf/model/GEOG/


#Random dates for boundary perturbations.
ENDPERTDATE=20091231180000                 
INIPERTDATE=20060101000000
PERTREFDATE=$IDATE            #At this date the initial perturbation dates will be taken. This date is used to keep consisntency among the perturbations
                              #used in forecast and analysis experiments. This date must be previous or equal to IDATE.

INPUT_PERT_DATES_FROM_FILE=0          #0 - generate a new set of random dates, 1 - read random dates from a file. 
INI_PERT_DATE_FILE=/home/jruiz/pepe   #List of initial random dates.


#### EXECUTABLES
RUNTIMELIBS=${HOMEDIR}/libs_sparc64/lib/
WRF=${HOMEDIR}/LETKF_WRF/wrf/
LETKF=$WRF/letkf/letkf.exe                     # LETKF module
UPDATEBC=$WRF/model/WRFDA/da_update_bc.exe     
WRFMODEL=$WRF/model/WRFV3.9.1_YAKAIRA/             # WRF model that run in computing nodes.
WRFMODELPPS=$WRF/model/WRFV3.9.1_YAKAIRA/      # WRF model that runs in pps server 
WPS=$WRF/model/WPS3.9.1_YAKAIRA/                # WRF model pre processing utilities (for pps server)
ARWPOST=$WRF/model/ARWpost_YAKAIRA/        # WRF model post processing utilities that run in computing nodes.
SPAWN=$WRF/spawn/
MPIBIN=mpiexec

#### SCRIPTS
UTIL=$WRF/run/util.sh

#### NAMELIST
NAMELISTWRF=$WRF/run/configuration/domain_conf/$DOMAINCONF/namelist.input                      #Namelist for WRF model.
NAMELISTWPS=$WRF/run/configuration/domain_conf/$DOMAINCONF/namelist.wps                        #Namelist for WRF pre processing tools
NAMELISTLETKF=$WRF/run/configuration/letkf_conf/letkf.namelist.$LETKFNAMELIST                  #Namelist for LETKF
NAMELISTARWPOST=$WRF/run/configuration/domain_conf/$DOMAINCONF/namelist.ARWpost                #Namelist for post-processing tools.
NAMELISTOBSOPE=$WRF/run/configuration/letkf_conf/obsope.namelist.$OBSOPENAMELIST               #Namelist for observation operator.
NAMELISTPERTMETEM=$WRF/run/configuration/letkf_conf/pertmetem.namelist.$LETKFNAMELIST          #Namelist for boundary perturbation.

