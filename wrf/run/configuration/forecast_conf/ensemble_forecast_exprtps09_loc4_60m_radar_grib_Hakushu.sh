#KALMAN FILTER CONFIGURATION
DOMAINCONF=CORDOBA_2KBIS               #Define a domain
LETKFNAMELIST=control                   #Define a letkf namelist template

MEMBER=30        #Number of ensemble members.
MAX_DOM=1        #Maximum number of WRF domains.
HOMEDIR=${HOME}/share/   
DATADIR=${HOME}/share/
ANALYSIS=0       #Identify this job as an analysis job.
FORECAST=1       #This is not a forecast job.
INTERPANA=0      #This is used in forecast jobs (but we need to define it here too)
RUN_ONLY_MEAN=0  #This is used in forecast jobs (but we need to define it here too)

USE_ANALYSIS_BC=1 #1 - use analysis as BC , 0 - use forecasts as bc (e.g. global gfs)
                  # if 1 then bc data will be taken from exp_met_em folder in the corresponding INPUT folder.
                  # if 0 then bc data will be taken from for_met_em folder in the corresponding INPUT folder.
                  # default is 1
USE_ANALYSIS_IC=0 #1 - use global analysis as IC, 0 use LETKF analysis as IC
                  #if 0 then profide a LETKF-analysis source (ANALYSIS_SOURC)
                  #default is 0

ANALYSIS_SOURCE=/data10/jruiz/EXPERIMENTS/ANALYSIS_CORDOBA_2KBIS_exprtps09_loc4_60m_radar_grib_Hakushu/

NVERTEXP=27  #Number of vertical levels in initial and boundary conditions input grib data.
NVERTDB=38   #Number of vertical levels in initial and boundary conditions perturbation input grib data.

#AUXILIARY VARIABLE FOR ENSEMBLE SIZE
MM=$MEMBER                      #Variable for iteration limits.
MEANMEMBER=`expr $MEMBER + 1 `  #This is the member ID corresponding to the ensemble mean.

WINDOW=1800       #Forecast initialization frequency (seconds)
GUESFT=7200       #Forecast length (seconds)
WINDOW_START=0    #Window start (seconds from forecast initialization)     
WINDOW_END=$GUESFT    #Window end   (seconds from forecast initialization)
WINDOW_FREC=600   #Output frequency within window (seconds) should be the same as the maximum observation frequency.
ASSIMILATION_FREC=$WINDOW #
NSLOTS=`expr $WINDOW_END \/ $WINDOW_FREC - $WINDOW_START \/ $WINDOW_FREC  + 1 `        #Number of time slots. 
NBSLOT=`expr $ASSIMILATION_FREC \/ $WINDOW_FREC - $WINDOW_START \/ $WINDOW_FREC + 1 `  #Time slot corresponding to the analysis.
if [ $NBSLOT -lt 10 ] ; then
   NBSLOT=0$NBSLOT
fi
SIGMA_OBS="2.0d3"          #NOT USED IN THE FORECAST
SIGMA_OBSV="0.2d0"         #NOT USED IN THE FORECAST
SIGMA_OBSZ="2.0d3"         #NOT USED IN THE FORECAST
SIGMA_OBST="3.0d0"         #NOT USED IN THE FORECAST
GROSS_ERROR="15.0d0"       #NOT USED IN THE FORECAST
COV_INFL_MUL="1.1d0"       #NOT USED IN THE FORECAST
SP_INFL_ADD="0.d0"         #NOT USED IN THE FORECAST
RELAX_ALPHA_SPREAD="0.8d0" #NOT USED IN THE FORECAST
RELAX_ALPHA="0.0d0"        #NOT USED IN THE FORECAST
USE_ADAPTIVE_INFLATION=0   #NOT USED IN THE FORECAST

#DOMAIN AND BOUNDARY DATA

BOUNDARY_DATA_FREQ=21600              #Boundary data frequency. (seconds)
BOUNDARY_DATA_PERTURBATION_FREQ=21600 #Frequency of data used to perturb boundary conditions (seconds)

#POSTPROC CONFIGURATION
OUTLEVS="0.1,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,"      #Level list
OUTVARS="'umet,vmet,pressure,W,QVAPOR,QCLOUD,QRAIN,QICE,QSNOW,QGRAUP,RAINNC,tk,u10m,v10m,slp,mcape,dbz,max_dbz'"  #Variable list.
ARWPOST_FREC=600   # Post processing frequency (seconds)
INPUT_ROOT_NAME='wrfout'
INTERP_METHOD=1
ENABLE_UPP=0

### LETKF setting
OBS=""                                                     # NOT USED IN THE FORECAST
RADAROBS="/OSSE_20140122_DBZ2.5_VR1.0_SO2KM/"              # NOT USED IN THE FORECAST
EXP=FORECAST_${DOMAINCONF}_${CONFIGURATION}                # name of experiment

### initial date setting
IDATE=20140122183000
EDATE=20140122200000

#### DATA
OBSDIR=${HOMEDIR}/DATA/OBS/$OBS/                                                          # NOT USED IN THE FORECAST
NRADARS=1                                                                                 # NOT USED IN THE FORECAST
RADAROBSDIR=${HOMEDIR}/DATA/OBS/$RADAROBS/                                                # NOT USED IN THE FORECAST
TMPDIR=${HOMEDIR}/TMP/$EXP/                                                               # Temporal work directory
OUTPUTDIR=${DATADIR}/EXPERIMENTS/$EXP/                                                    # Where results will be stored.
GRIBDIR=${HOMEDIR}/DATA/GRIB/FNL/HIRES/ARGENTINA/                                         # Folder where bdy and initial grib files are located.
GRIBTABLE="Vtable.GFS"                                                                    # Bdy and init data source Vtable name.
PERTGRIBDIR=${HOMEDIR}/DATA/GRIB/CFSR/HIRES/ARGENTINA/                                    # Folder where data for perturbing bdy are located.
PERTGRIBTABLE="Vtable.CFSR"                                                               # Bdy perturbation source vtable name.
GEOG=${HOMEDIR}/LETKF_WRF/wrf/model/GEOG/                                                 # Folder where WPS GEOG dataset is located.

#INITIAL AND BOUNDARY RANDOM PERTURBATIONS
SCALE_FACTOR="0.05"         #Perturbation scale factor.
RANDOM_SCALE_FACTOR="0.0"   #Random perturbation scale factor.
PERTURB_BOUNDARY=1          #Wheter boundary conditions are going to be perturbed.
PERTURB_BOUNDARY_TYPE=1     #DUMMY
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
PERTREFDATE=20140122000000    #At this date the initial perturbation dates will be taken. This date is used to keep consisntency among the perturbations
                              #used in forecast and analysis experiments. This date must be previous or equal to IDATE.
INPUT_PERT_DATES_FROM_FILE=1  #0 - generate a new set of random dates, 1 - read random dates from a file. 
INI_PERT_DATE_FILE=${HOMEDIR}/DATA/INITIAL_RANDOM_DATES/initial_perturbation_dates_60m  #List of initial random dates.

#### EXECUTABLES
RUNTIMELIBS=/apps/SLES11/opt/netcdf-fortran/4.4.1_intel/lib/:/apps/SLES11/opt/hdf5/1.8.14_intel/lib/:/apps/COMMON/opt/intel/compilers_and_libraries_2017.1.132/linux/lib:/apps/COMMON/opt/intel/compilers_and_libraries_2017.4.196//linux/mpi/intel64/lib
WRF=${HOMEDIR}/LETKF_WRF/wrf/                  # WRF folder (for computing nodes)
LETKF=$WRF/letkf/letkf.exe                     # LETKF module (for computing nodes)
UPDATEBC=$WRF/model/WRFDA/da_update_bc.exe     # Update bc tool (WRFVAR) (for computing nodes)
WRFMODEL=$WRF/model/WRFV3.6.1/                     # WRF model that run in computing nodes.
WRFMODELPPS=$WRF/model/WRFV3.6.1/                  # WRF model that runs in pps server  (usually the same as the one for the computing nodes)
WPS=$WRF/model/WPSV3.6.1/                            # WRF model pre processing utilities (for pps server)
ARWPOST=$WRF/model/ARWpost/                    # WRF model post processing utilities that run in computing nodes.
SPAWN=$WRF/spawn/
MPIBIN=/apps/COMMON/opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/bin64/mpiexec


#### SCRIPTS
UTIL=$WRF/run/util.sh                          # Script containing bash functions that will be used during execution.

#### NAMELIST
NAMELISTWRF=$WRF/run/configuration/domain_conf/$DOMAINCONF/namelist.input           #Namelist for WRF model.
NAMELISTWPS=$WRF/run/configuration/domain_conf/$DOMAINCONF/namelist.wps             #Namelist for WRF pre processing tools
NAMELISTLETKF=$WRF/run/configuration/letkf_conf/letkf.namelist.$LETKFNAMELIST       #Namelist for LETKF
NAMELISTARWPOST=$WRF/run/configuration/domain_conf/$DOMAINCONF/namelist.ARWpost     #Namelist for post-processing tools.
NAMELISTOBSOPE=$WRF/run/configuration/letkf_conf/obsope.namelist.$OBSOPENAMELIST    #Namelist for observation operator.
NAMELISTPERTMETEM=$WRF/run/configuration/letkf_conf/pertmetem.namelist.control      #Namelist for boundary conditions perturbation.

