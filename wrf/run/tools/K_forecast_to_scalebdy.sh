#This scripts copies the data from an input experiment with the name convention expected by scale LETKF system.

source ../util.sh

INIDATE=20130713000000         #Date of forecast init.       (yyyymmddhhmmss)
FORECAST_OUTPUT_FREQ=30        #Frequency of forecast output (second)
GUESFT=43200                   #Forecast length (second)
DOMAIN=02                      #Domain that will be used as bdy conditions for the scale model.
MEMBER=1024                    #Ensemble size
OPERATION="ln -sf "            #Can be copy / move or link

CINIDATE=20130713051000        #Data copy/linking will be performed only within this dates.
CENDDATE=20130713053000

#INBASEPATH=$HOME/data/EXPERIMENTS/FORECAST_OSAKA_1KM_NESTED_osaka_nested_1024m/
#OUTBASEPATH=$HOME/share/LETKF_SCALE/ncepfnl/wrfout_osaka_largens/

INBASEPATH=$HOME/share/FORECAST_OSAKA_1KM_NESTED_1024m_TIMEINT/
OUTBASEPATH=$HOME/share/INPUT_SCALE/ncepfnl/wrfout_osaka_largens/

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






