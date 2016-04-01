#!/bin/sh

#Load useful functions
source ../../../run/util.sh

###########################################
#CONFIGURATION
###########################################

INIDATE=20130713050940   #20130712230440  #20130713050940   #20130712230440
ENDDATE=20130713060940   #20130713110440  #20130713051540   #20130713110440
DT=30   #Time interval between two obs (seconds)
#Input time in JST, output time is in UTC
UTC_TO_JST=32400 #Japanese standard time increment.

OBSNAME=osaka_pawr
PROCOBSNAME=${OBSNAME}_super1km_latlon

RAWOBSPATH=$HOME/share/OBS/$OBSNAME
PROCOBSPATH=$HOME/share/OBS/$PROCOBSNAME

TMPDIR=$HOME/data/TMP/RADAR_PREP_$PROCOBSNAME
NAMELIST=$TMPDIR/radarprep.namelist

DX=1000.0d0      #Horizontal resolution
DZ=1000.0d0      #Vertical resolution
MAXRANGE=60d3    #Maximum horizontal range
MAXZ=20d3        #Maximum vertical range 
DBZ_TO_POWER=.true.  #Power transformation flag
DEBUG_OUTPUT=.true.  #Debug ouput flag (binary file to visualize data in grads)
ERROR_REF=5.0d0      #Reflectivity error
ERROR_VR=1.0d0       #Radial velocity error
ID_REF_OBS=4001      #ID for reflectivity observations
ID_VR_OBS=4002       #Id for radial velocity in output file.
#Parameters for osse experiments
OSSE_EXP=.false.     #Identify it as an osse experiment data
SIMULATE_ERROR=.false.  #Wheter we are going to add noise
ERROR_REF_SCLX=5000d0   #SPATIAL SCALE OF ERROR IN METERS
ERROR_VR_SCLX=5000d0    #SPATIAL SCALE OF ERRORS IN METERS
ERROR_REF_SCLZ=5000d0   #VERTICAL SCALE OF ERRORS IN METERS
ERROR_VR_SCLZ=5000d0    #VERTICLA SCALE OF ERRORS IN METERS
LATLON_COORD=.true.     #True output in lat lon, false output in azimuth,range,elev
USE_ATTENUATION=.true.  #True consider attenuation (attenuated pixels are discarded)
ATTENUATION_THRESHOLD=0.01 
USE_QCFLAG=.true.      #True use qc information

#--------------------------------------------------------
#Initialize directories
mkdir -p $TMPDIR
mkdir -p $PROCOBSPATH
cp ../radar_prep/radar_prep.exe $TMPDIR

#Generate namelist
echo "&general                              " > $NAMELIST
echo "DX=$DX                                " >>$NAMELIST
echo "DZ=$DZ                                " >>$NAMELIST
echo "MAXRANGE=$MAXRANGE                    " >>$NAMELIST
echo "MAXZ=$MAXZ                            " >>$NAMELIST
echo "DBZ_TO_POWER=$DBZ_TO_POWER            " >>$NAMELIST
echo "DEBUG_OUTPUT=$DEBUG_OUTPUT            " >>$NAMELIST
echo "ERROR_REF=$ERROR_REF                  " >>$NAMELIST
echo "ERROR_VR=$ERROR_VR                    " >>$NAMELIST
echo "ID_REF_OBS=$ID_REF_OBS                " >>$NAMELIST
echo "ID_VR_OBS=$ID_VR_OBS                  " >>$NAMELIST
echo "LATLON_COORD=$LATLON_COORD            " >>$NAMELIST
echo "USE_ATTENUATION=$USE_ATTENUATION      " >>$NAMELIST
echo "USE_QCFLAG=$USE_QCFLAG                " >>$NAMELIST
echo "ATTENUATION_THRESHOLD=$ATTENUATION_THRESHOLD " >>$NAMELIST 
echo "/                                     " >>$NAMELIST
echo "&osse                                 " >>$NAMELIST
echo "OSSE_EXP=$OSSE_EXP                    " >>$NAMELIST
echo "SIMULATE_ERROR=$SIMULATE_ERROR        " >>$NAMELIST
echo "ERROR_REF_SCLX=$ERROR_REF_SCLX        " >>$NAMELIST
echo "ERROR_VR_SCLX=$ERROR_VR_SCLX          " >>$NAMELIST
echo "ERROR_REF_SCLZ=$ERROR_REF_SCLZ        " >>$NAMELIST
echo "ERROR_VR_SCLZ=$ERROR_VR_SCLZ          " >>$NAMELIST
echo "/                                     " >>$NAMELIST
 
#Data processing
cd $TMPDIR

CDATE=$INIDATE
while [ $CDATE -le  $ENDDATE ]
do

JST_DATE=`date_edit2 $CDATE $UTC_TO_JST `
JST_DATE1=`echo $JST_DATE | cut -c1-8  `
JST_DATE2=`echo $JST_DATE | cut -c9-14 `

#Radar data is sometimes times that are not multiples of 
#DT, model simulations usually works at multiples of DT
#So in this case data is sligthly offset to the nearest multiple
#of DT in order to be used in the assimilation experiments.
DATE_FLOOR=`date_floor $CDATE $DT `


ln -sf $RAWOBSPATH/PAWR_LETKF_INPUTV4_${JST_DATE1}-${JST_DATE2}.dat ./radar_input.grd
ln -sf $PROCOBSPATH/radar_${DATE_FLOOR}.dat ./radarobs.dat
ln -sf $PROCOBSPATH/radar_${DATE_FLOOR}.grd ./super.grd

echo "Input  data: $RAWOBSPATH/PAWR_LETKF_INPUTV4_${JST_DATE1}-${JST_DATE2}.dat "
echo "Output data: $PROCOBSPATH/radar_${DATE_FLOOR}.dat "

./radar_prep.exe 

 
CDATE=`date_edit2 $CDATE $DT `
done #Done over the dates.

#--------------------------------------------------------




