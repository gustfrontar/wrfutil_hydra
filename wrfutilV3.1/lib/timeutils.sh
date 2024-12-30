date_floor () {
        if [ $# -lt 2 ]; then
                echo "Usage : date_floor"
                echo "    date_floor date [%Y-%m-%d %hh:%mm:%ss] interval [seconds] "
                echo "    ex) "date_floor "2018-11-20 00:00:00" 3600
                exit 1
        fi

        local DATE=$1      #date format [%Y-%m-%d %hh:%mm:%ss]
        local INTERVAL=$2  #[seconds]

	seconds=$(date -u  +%s -d "$DATE UTC")
	mod=$(($seconds % $INTERVAL))
	echo "$(date -u "+%Y-%m-%d %T" -d "$DATE UTC -$((10#$mod)) seconds")"

}

date_ceil () {
        if [ $# -lt 2 ]; then
                echo "Usage : date_ceil"
                echo "    date_ceil date [%Y-%m-%d %hh:%mm:%ss] interval [seconds] "
                echo "    ex) "date_ceil "2018-11-20 00:00:00" 3600
                exit 1
        fi

        local DATE=$1      #date format [%Y-%m-%d %hh:%mm:%ss]
        local INTERVAL=$2  #[seconds]

        seconds=$(date -u  +%s -d "$DATE UTC")
        mod=$(($seconds % $INTERVAL))
	if [ $mod -eq 0 ] ; then
	   echo "$(date -u "+%Y-%m-%d %T" -d "$DATE UTC")"
	else  
	   mod=$(($INTERVAL-$mod))
	   echo "$(date -u "+%Y-%m-%d %T" -d "$DATE UTC +$((10#$mod)) seconds")"
	fi
}

maximum_common_divisor() (
    if (( $1 % $2 == 0)); then
        echo $2
     else
        maximum_common_divisor $2 $(( $1 % $2 ))
    fi
)


write_step_conf() (
    CALLER=$1  #WPS,GUESS,FORECAST

    if [ $CALLER == "WPS" ] ; then 
       if [ $STEP == 0 ] ; then 
         INI_DATE_FCST=$DA_INI_DATE   
         END_DATE_FCST=$(date -u -d "$INI_DATE_FCST     UTC +$SPIN_UP_LENGTH seconds" +"%Y-%m-%d %T")   
       else 
         INI_DATE_FCST=$(date -u -d "$DA_INI_DATE UTC +$((10#$ANALYSIS_FREQ*(10#$STEP-1)+10#$SPIN_UP_LENGTH)) seconds" +"%Y-%m-%d %T")
         END_DATE_FCST=$(date -u -d "$INI_DATE_FCST UTC +$ANALYSIS_WIN_END seconds" +"%Y-%m-%d %T")
       fi 
    fi
    if [ $CALLER == "GUESS" ] ; then
       if [ $STEP == 0 ] ; then 
         INI_DATE_FCST=$DA_INI_DATE
         END_DATE_FCST=$(date -u -d "$INI_DATE_FCST     UTC +$SPIN_UP_LENGTH seconds" +"%Y-%m-%d %T")
         ANALYSIS_DATE=$(date -u -d "$INI_DATE_FCST     UTC +$SPIN_UP_LENGTH seconds" +"%Y-%m-%d %T")
         FCST_BDY_FREQ=$( maximum_common_divisor $SPIN_UP_LENGTH $BDY_FREQ  ) #Optimization of boundary conditions data.
         FCST_OFREQ=$SPIN_UP_LENGTH
       else
         INI_DATE_FCST=$(date -u -d "$DA_INI_DATE UTC +$((10#$ANALYSIS_FREQ*(10#$STEP-1)+10#$SPIN_UP_LENGTH)) seconds" +"%Y-%m-%d %T")
         END_DATE_FCST=$(date -u -d "$INI_DATE_FCST     UTC +$ANALYSIS_WIN_END seconds" +"%Y-%m-%d %T")
         ANALYSIS_DATE=$(date -u -d "$INI_DATE_FCST     UTC +$ANALYSIS_FREQ    seconds" +"%Y-%m-%d %T")
         FCST_BDY_FREQ=$( maximum_common_divisor $ANALYSIS_FREQ $ANALYSIS_WIN_END ) #Optimization of boundary conditions data.
         FCST_OFREQ=$ANALYSIS_WIN_STEP
       fi
    fi
    if [ $CALLER == "FORECAST" ] ; then
       INI_DATE_FCST=$(date -u -d "$FCST_INI_DATE UTC +$(($FCST_INI_FREQ*$STEP)) seconds" +"%Y-%m-%d %T")
       END_DATE_FCST=$(date -u -d "$INI_DATE_FCST     UTC +$FCST_LEAD_TIME seconds" +"%Y-%m-%d %T")
       FCST_BDY_FREQ=$( maximum_common_divisor $FCST_LEAD_TIME $BDY_FREQ ) #Optimization of boundary conditions data.
    fi

    BDY_INI_DATE=$(date_floor "$INI_DATE_FCST" $BDY_FREQ )
    BDY_END_DATE=$(date_ceil  "$END_DATE_FCST" $BDY_FREQ )
    INI_BDY_DATE=$(date_floor "$INI_DATE_FCST" $BDY_INI_FREQ )
    echo "Generating step.conf"
    echo "# This namelist is authomatically generated  " >  $BASEDIR/conf/step.conf
    echo "STEP=$STEP"                                    >> $BASEDIR/conf/step.conf 
    echo "INI_DATE_FCST=\"$INI_DATE_FCST\""              >> $BASEDIR/conf/step.conf
    echo "END_DATE_FCST=\"$END_DATE_FCST\""              >> $BASEDIR/conf/step.conf  
    echo "BDY_INI_DATE=\"$BDY_INI_DATE\""                >> $BASEDIR/conf/step.conf
    echo "BDY_END_DATE=\"$BDY_END_DATE\""                >> $BASEDIR/conf/step.conf
    echo "INI_BDY_DATE=\"$INI_BDY_DATE\""                >> $BASEDIR/conf/step.conf
    echo "FCST_BDY_FREQ=\"$FCST_BDY_FREQ\""              >> $BASEDIR/conf/step.conf
    echo "FCST_OFREQ=\"$FCST_OFREQ\""                    >> $BASEDIR/conf/step.conf
    echo "ANALYSIS_DATE=\"$ANALYSIS_DATE\""              >> $BASEDIR/conf/step.conf
    echo "EXPTYPE=\"$EXPTYPE\""                          >> $BASEDIR/conf/step.conf


    echo "=========================="
    echo "FORECAST FOR STEP " $STEP
    echo "Starting at " $INI_DATE_FCST
    echo "Ending   at " $END_DATE_FCST
    if [ $CALLER == "GUESS" ] ; then
       echo "Analysis time at " $ANALYSIS_DATE
    fi

)














