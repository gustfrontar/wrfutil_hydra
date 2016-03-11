#This scripts copies the data from an input experiment with the name convention expected by scale LETKF system.

source ../util.sh

INIDATE=20130713000000         #Date of forecast init.       (yyyymmddhhmmss)
FORECAST_OUTPUT_FREQ=300        #Frequency of forecast output (second)
GUESFT=43200                   #Forecast length (second)
DOMAIN=02                      #Domain that will be used as bdy conditions for the scale model.
MEMBER=10                      #Ensemble size
OPERATION="cp -f "            #Can be copy / move or link

CINIDATE=20130713051000        #Data copy/linking will be performed only within this dates.
CENDDATE=20130713052000

#INBASEPATH=$HOME/data/EXPERIMENTS/FORECAST_OSAKA_1KM_NESTED_osaka_nested_1024m/
#OUTBASEPATH=$HOME/share/LETKF_SCALE/ncepfnl/wrfout_osaka_largens/

INBASEPATH=/home/hp150019/k02016/data/EXPERIMENTS/FORECAST_OSAKA_1KM_NESTED_osaka_nested_1024m/
OUTBASEPATH=$HOME/share/INPUT_SCALE/ncepfnl/wrfout_osaka_largens_5min/

#Forecast end date.


CDATE=$CINIDATE

while [ $CDATE -le $CENDDATE  ] ; do

   M=1
   while [ $M -le $MEMBER  ] ; do

    MEM=` ens_member $M `

    MEMSCALE=`echo $MEM | cut -c2-5 `

    mkdir -p $OUTBASEPATH/$MEMSCALE/ 

    my_file=`wrfout_file_name $CDATE $DOMAIN ` 

    $OPERATION $INBASEPATH/forecast/$INIDATE/$MEM/$my_file $OUTBASEPATH/$MEMSCALE/wrfout_$CDATE 

    M=`expr $M + 1 `
   done

   CDATE=`date_edit2 $CDATE $FORECAST_OUTPUT_FREQ `

done






