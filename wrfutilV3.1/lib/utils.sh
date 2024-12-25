date_floor () {
        if [ $# -lt 2 ]; then
                echo "Usage : date_floor"
                echo "    date_floor [yyyy][mm][dd][hh][mn] [interval(seconds)]"
                echo "    ex) date_floor 201005051200 3600"
                exit 1
        fi

        local DATE=$1
        local INTERVAL=$2

        cy1=`echo $DATE | cut -c1-4`
        cm1=`echo $DATE | cut -c5-6`
        cd1=`echo $DATE | cut -c7-8`
        ch1=`echo $DATE | cut -c9-10`
        cn1=`echo $DATE | cut -c11-12`
        cs1=`echo $DATE | cut -c13-14`

        seconds1=`date +%s -d"$cy1-$cm1-$cd1 $ch1:$cn1:$cs1 UTC"`
        mod=`expr $seconds1 % $INTERVAL `

        DATEFLOOR=`date_edit2 $DATE -$mod `
        echo $DATEFLOOR
}




