#Configuration file for the breeding cycle.

NBV=1  #Number of bred vectors to be computed.

#################### System configuration ######################################

LD_LIBRARY_PATH_ADD="/data/opt/netcdf-fortran/4.4.1_intel/lib/:/data/opt/hdf5/1.8.14_intel/lib/:/data/opt/netcdf/4.3.2_intel/lib/"
MY_OMP_NUM_THREADS=4
MY_OMP_STACKSIZE=70M
NODESPERJOB=1

OUTPUTDIR="/home/jruiz/share/exp/${CONFIGURATION}/" #Output root directory where results will be stored.
TMPDIR="${SRCDIR}/../tmp/${CONFIGURATION}"   #TMP directory where computations will take place.

SAVE_HISTORY=1

############## Breeding configuration ##################################


BVFREQ=30  #Rescalling frequency in seconds.

DATABASEROOT="/home/jruiz/share/database_v170124/"

INIT_BV=1                   #1 - initialize bv from random perturbations , 0 do not initialize bv.

INIDATE=20130713051030   #Start date for the breeding cycle.
ENDDATE=20130713053930   #End date for the breeding cycle.

BDYPERT=0             #0 means do not perturb buondary , 1 means perturb boundary.
BDYDATAFREQ=30        #Boundary data frequency in seconds.
BDYDATAROOT="/home/jruiz/share/scale_input_data/wrf_osaka_ensemble_1km_2013_07_13/"

FLENGTH=$BVFREQ        #Forecast length in seconds.
FOUTPUTFREQ=$BVFREQ    #Forecast output freq.
FRESTARTFREQ=$BVFREQ   #Forecast restart output frequency.

# Central trajectory info

BVCENTRALTRAJECTORY="mean"   #Can be one of the ensemble members or the ensemble mean.
BVDATAROOT="/home/jruiz/share/scale_input_data/OsakaPAR_1km_control1000m_smallrandompert_new/"  #Folder for the initial perturbations and for central trayectory.
PRC_NUM_X_CT=2  #Number of subdomains in X for the central trajectory.
PRC_NUM_Y_CT=2  #Number of subdomains in Y for the central trajectory.

NPROCS_CT=`expr $PRC_NUM_X_CT \* $PRC_NUM_Y_CT `

# Breeding in place configuration

BREEDING_IN_PLACE=1                      #1 means enable RUNINPLACE, 0 means disable RUNINPLACE.
BREEDING_IN_PLACE_IT=20                  #How many run in place iterations will be performed when RUNINPLACE is enabled.
BREEDING_IN_PLACE_INIDATE=20130713051030 #First date in which running in place will be applied.
BREEDING_IN_PLACE_ENDDATE=20130713053930 #First date in which running in place will be applied.

if [ $BREEDING_IN_PLACE -eq 0 ] ; then
   BREEDING_IN_PLACE_IT=1  #This is equivalent to not performing breeding_in_place.
fi

######## Domain size section    ##########################
#Domain settings has to be consistent with the settings of the central trayectory.

PRC_NUM_X_IN=4   #Number of subdomains in X to integrate the perturbations.
PRC_NUM_Y_IN=5  #Number of subdomains in Y to integrate the perturbations.

IMAX_IN=45
JMAX_IN=36

NPROCS=`expr $PRC_NUM_X_IN \* $PRC_NUM_Y_IN `

######## Postprocessing section ##########################

Z_LEV_TYPE_IN=plev           
Z_MERGE_OUT_IN=.true.
T_MERGE_OUT_IN=.false.
VNAME_IN='"U","V","T","W","QV","QHYD"'
TARGET_ZLEV_IN="1000,975,950,925,900,875,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150"
ZCOUNT_IN=21
START_STEP_IN=1
END_STEP_IN=2
INC_STEP_IN=1

MAPPROJ_ctl_IN=.true.

DELT_IN='5mn'
STIME_IN='00UTC01Jan2000'

DOMAIN_NUM_IN='01'


