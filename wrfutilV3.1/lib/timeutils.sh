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

maximum_common_divisor () {

   local n1=$1
   local n2=$2
   local gcd=0
   if test $n1 -gt $n2 ; then
      local i=1
      while test $i -le $n1 ; do
        local a=`expr $n1 % $i`
        local b=`expr $n2 % $i`
        if test $a -eq 0 -a $b -eq 0 ; then
           if test $gcd -lt $i ; then
             gcd=$i
           fi
        fi
        i=`expr $i + 1`
      done
   fi
   if test $n2 -gt $n1 ; then
      i=1
      while test $i -le $n2 ; do
        a=`expr $n1 % $i`
        b=`expr $n2 % $i`
        if test $a -eq 0 -a $b -eq 0 ; then
           if test $gcd -lt $i ; then
             gcd=$i
           fi
        fi
        i=`expr $i + 1`
      done
   fi

   if test $n2 -eq $n1 ; then
     gcd=$n2
   fi

   echo $gcd
}




