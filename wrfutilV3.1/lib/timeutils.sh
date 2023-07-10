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





