#This script is to interpolate wrf outputs (forecast for example) linearly in time.
#This can be used to increase the frequency of bdy conditions for the scale model.

source ../util.sh
FORECASTINIDATE=20130713000000 #Date of forecast initialization.
INIDATE=20130713053000         #Date of interpolation init  (yyyymmddhhmmss)
ENDDATE=20130713054000         #Date of interpolation end   (yyyymmddhhmmss)
FORECAST_OUTPUT_FREQ=300       #Frequency of forecast output (second)
FORECAST_OUTPUT_FREQ_INT=30    #Desired frequency of forecast output.

PPSSERVER=pps2

DOMAIN=02                      #Domain that will be interpolated.
MEMBER=1024                    #Ensemble size

MAXRUN=10

FORECASTBP=$HOME/data/EXPERIMENTS/FORECAST_OSAKA_1KM_NESTED_osaka_nested_1024m/  #Forecast basepath.

OUTPATH=$HOME/share/FORECAST_OSAKA_1KM_NESTED_1024m_TIMEINT/forecast/

TMPDIR=$HOME/data/TMP/FORECAST_INTERPOLATION/

EXEC=$HOME/share/LETKF_WRF/wrf/interpwrfout/interpwrfout.exe

CDATE=$INIDATE

mkdir -p $OUTPATH/$FORECASTINIDATE/
mkdir -p $TMPDIR
cp $EXEC $TMPDIR

cd $TMPDIR

while [ $CDATE -le $ENDDATE  ] ; do

   #For each interpolation date get the date inmediatly below and above in the forecast output
   LDATE=`date_floor $CDATE $FORECAST_OUTPUT_FREQ `
   UDATE=`date_edit2 $LDATE $FORECAST_OUTPUT_FREQ `  

   LWRFFILE=`wrfout_file_name $LDATE $DOMAIN `
   UWRFFILE=`wrfout_file_name $UDATE $DOMAIN `
  
   CWRFFILE=`wrfout_file_name $CDATE $DOMAIN ` 

   echo "I will interpolate the forecast at date $CDATE "
   echo "Lower date is $LDATE and Upper date is $UDATE "

   M=1
   while [ $M -le $MEMBER  ] ; do

    RUN=1
    while [ $RUN -le $MAXRUN -a $M -le $MEMBER ] ; do
     MEM=` ens_member $M `

     WORKDIR=$TMPDIR/$RUN
     echo "Interpolationg member $MEM "

     echo "mkdir -p $WORKDIR                                                                                             "> tmp$RUN.sh
     echo "cd $WORKDIR                                                                                                   ">> tmp$RUN.sh
     echo "mkdir -p $OUTPATH/$FORECASTINIDATE/$MEM/                                                                      ">> tmp$RUN.sh
     echo "mkdir -p $WORKDIR                                                                                             ">> tmp$RUN.sh
     echo "ln -sf $FORECASTBP/forecast/$FORECASTINIDATE/$MEM/$LWRFFILE $WORKDIR/input_file1.nc                           ">> tmp$RUN.sh
     echo "ln -sf $FORECASTBP/forecast/$FORECASTINIDATE/$MEM/$UWRFFILE $WORKDIR/input_file2.nc                           ">> tmp$RUN.sh 
     echo "cp $FORECASTBP/forecast/$FORECASTINIDATE/$MEM/$LWRFFILE     $OUTPATH/$FORECASTINIDATE/$MEM/$CWRFFILE          ">> tmp$RUN.sh
     echo "ln -sf $OUTPATH/$FORECASTINIDATE/$MEM/$CWRFFILE $WORKDIR/input_fileint.nc                                     ">> tmp$RUN.sh
     echo "ln -sf ../interpwrfout.exe ./                                                                                 ">> tmp$RUN.sh
     echo "./interpwrfout.exe  $CDATE    > interpwrf.log                                                                 ">> tmp$RUN.sh
     chmod 755 tmp$RUN.sh

     ssh $PPSSERVER $TMPDIR/tmp$RUN.sh &

    M=`expr $M + 1 `
    RUN=`expr $RUN + 1 `
    done

    time wait


   done

   CDATE=`date_edit2 $CDATE $FORECAST_OUTPUT_FREQ_INT `

done






