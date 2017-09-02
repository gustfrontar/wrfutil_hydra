
DATABASEROOT="/home/jruiz/share/database_v170124/"

OUTPUTDIR="/home/jruiz/share/exp/${CONFIGURATION}/" #Output root directory where results will be stored.

TMPDIR="${SRCDIR}/../tmp/${CONFIGURATION}"   #TMP directory where computations will take place.

INIDATE=201605200000   #Start date of the run.

BDYDATAFREQ=3600        #Boundary data frequency in seconds.
BDYDATAROOT="/home/jruiz/share/scale_input_data/test_germany_15km_3h_cycle/20160520000000/hist/0001/"

FLENGTH=120       #Forecast length in seconds.
FOUTPUTFREQ=30    #Forecast output freq.
FRESTARTFREQ=30   #Forecast restart output frequency.

######## Domain size section    ##########################

PRC_NUM_X_IN=4
PRC_NUM_Y_IN=2

IMAX_IN=60
JMAX_IN=60

NPROCS=`expr $PRC_NUM_X_IN \* $PRC_NUM_Y_IN `

######## Postprocessing section ##########################

Z_LEV_TYPE_IN=plev            
Z_MERGE_OUT_IN=.true.
T_MERGE_OUT_IN=.false.
VNAME_IN='"U","V","T"'
TARGET_ZLEV_IN="1000,850,700,500,300,200"
ZCOUNT_IN=6
START_STEP_IN=1
END_STEP_IN=2
INC_STEP_IN=1

MAPPROJ_ctl_IN=.true.

DELT_IN='5mn'
STIME_IN='00UTC01Jan2000'

DOMAIN_NUM_IN='01'


