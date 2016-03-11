#!/bin/bash
# =======================================================================
#
#	Utility Shell Finctions for WRF_LETKF
#
#                                                   2010.05.11 M.Kunii
# =======================================================================

# -----------------------------
#    cal_date2min
# -----------------------------
cal_date2sec () {
(

	# calculate total minutes from 1990.01.01.00UTC

	if [ $# -lt 5 ]; then
		echo "Usage : cal_date2sec"
		echo "    cal_date2min [yy] [mm] [dd] [hh] [mn] [ss]"
		exit
	fi

	tyy=$1
	tmm=$2
	tdd=$3
	thh=$4
	tmn=$5
        tss=$6
	
	total_sec=0
        total_min=0

	# --- year ---	
	yy=1990
	while [ ${yy} -lt ${tyy} ]
	do
		tmp1=`expr ${yy} % 4` || test 1 -eq 1
		tmp2=`expr ${yy} % 100`	 || test 1 -eq 1
	
		if [ ${tmp1} -eq 0 -a ${tmp2} -ne 0 ]; then
			days=366
		else
			days=365
		fi
		
		total_min=`expr ${total_min} + ${days} \* 1440` || test 1 -eq 1
		yy=`expr ${yy} + 1`	 || test 1 -eq 1
	done

	# --- month ---
	mm=1
	while [ ${mm} -lt ${tmm} ]
	do
		tmp1=`expr ${yy} % 4` || test 1 -eq 1
		tmp2=`expr ${yy} % 100` || test 1 -eq 1
	
		if [ ${tmp1} -eq 0 -a ${tmp2} -ne 0 ]; then
			flag=1
		else
			flag=0
		fi
			
		case ${mm} in
			1|3|5|7|8|10|12)
				days=31
				;;
			4|6|9|11)
				days=30
				;;
			2)
				days=`expr 28 + ${flag}`
				;;
		esac
	
		total_min=`expr ${total_min} + ${days} \* 1440`
		mm=`expr ${mm} + 1`
	done

	# --- day ---
	dd=1
	while [ ${dd} -lt ${tdd} ]
	do
		total_min=`expr ${total_min} + 1440`
		dd=`expr ${dd} + 1`
	done
	
	# --- hour ---	
	hh=0
	while [ ${hh} -lt ${thh} ]
	do		
		total_min=`expr ${total_min} + 60`
		hh=`expr ${hh} + 1`
	done

	# --- minute ---	
	total_min=`expr ${total_min} + ${tmn}`


        total_sec=`expr ${total_min} \* 60`

        total_sec=`expr ${total_sec} + ${tss}`
	
	echo ${total_sec}
	
)
}


# -----------------------------
#    cal_min2date
# -----------------------------
cal_sec2date () {
(

	# calculate date from total minutes from 1990.01.01.00UTC

	if [ $# -lt 1 ]; then
		echo "Usage : cal_sec2date"
		echo "    cal_sec2date 89411040"
		exit
	fi

	input_sec=$1

	# --- year ---	
	yy=1990
	while :
	do
		tmp1=`expr ${yy} % 4` || test 1 -eq 1
		tmp2=`expr ${yy} % 100` || test 1 -eq 1
	
		if [ ${tmp1} -eq 0 -a ${tmp2} -ne 0 ]; then
			days=366
		else
			days=365
		fi
		
		total_sec=`expr ${total_sec} + ${days} \* 86400`
		if [ ${total_sec} -gt ${input_sec} ]; then
			total_sec=`expr ${total_sec} - ${days} \* 86400` || test 1 -eq 1
			break
		fi
		yy=`expr ${yy} + 1`
	done
	tyy=${yy}

	# --- month ---				
	mm=1
	while :
	do
		tmp1=`expr ${yy} % 4` || test 1 -eq 1
		tmp2=`expr ${yy} % 100` || test 1 -eq 1
	
		if [ ${tmp1} -eq 0 -a ${tmp2} -ne 0 ]; then
			flag=1
		else
			flag=0
		fi
			
		case ${mm} in
			1|3|5|7|8|10|12)
				days=31
				;;
			4|6|9|11)
				days=30
				;;
			2)
				days=`expr 28 + ${flag}`
				;;
		esac			
		
		total_sec=`expr ${total_sec} + ${days} \* 86400`
		if [ ${total_sec} -gt ${input_sec} ]; then
			total_sec=`expr ${total_sec} - ${days} \* 86400` || test 1 -eq 1
			break		
		fi				
		mm=`expr ${mm} + 1`	
	done
	if [ ${mm} -lt 10 ]; then
		tmm="0${mm}"
	else
		tmm=${mm}
	fi
	
	# --- day ---	
	dd=1
	while :
	do
		total_sec=`expr ${total_sec} + 86400`
		if [ ${total_sec} -gt ${input_sec} ]; then
			total_sec=`expr ${total_sec} - 86400` || test 1 -eq 1
			break		
		fi		
		dd=`expr ${dd} + 1`		
	done
	if [ ${dd} -lt 10 ]; then
		tdd="0${dd}"
	else
		tdd=${dd}
	fi	

	# --- hour ---
	hh=0
	while :
	do		
		total_sec=`expr ${total_sec} + 3600`
		if [ ${total_sec} -gt ${input_sec} ]; then
			total_sec=`expr ${total_sec} - 3600` || test 1 -eq 1
			break		
		fi		
		hh=`expr ${hh} + 1`
	done
	if [ ${hh} -lt 10 ]; then
		thh="0${hh}"
	else
		thh=${hh}
	fi

        # --- minute ---
        mn=0
        while :
        do
                total_sec=`expr ${total_sec} + 60`
                if [ ${total_sec} -gt ${input_sec} ]; then
                        total_sec=`expr ${total_sec} - 60` || test 1 -eq 1
                        break
                fi
                mn=`expr ${mn} + 1`
        done
        if [ ${mn} -lt 10 ]; then
                tmn="0${mn}"
        else
                tmn=${mn}
        fi


	# --- second ---		
	ss=`expr ${input_sec} - ${total_sec}` || test 1 -eq 1
	if [ ${ss} -lt 10 ]; then
		tss="0${ss}"
	else
		tss=${ss}
	fi	
	
	echo "${tyy} ${tmm} ${tdd} ${thh} ${tmn} ${tss}"
		
)
}

# -----------------------------
#    date_edit
# -----------------------------
date_edit () {
(

	if [ $# -lt 7 ]; then
		echo "Usage : date_edit"
		echo "    date_edit [yyyy] [mm] [dd] [hh] [mn] [ss] [dt(sec)]"
		echo "    ex) date_edit 201005051200 -180"
		exit
	fi
	
	yy=$1
	mm=$2
	dd=$3
	hh=$4
	mn=$5
        ss=$6
	dt=$7

	cal_date2sec ${yy} ${mm} ${dd} ${hh} ${mn} ${ss} ${dt} > tsec.$$
	read tsec < tsec.$$
	
	tsec2=`expr ${tsec} + ${dt}`
	
	cal_sec2date ${tsec2}
	rm -f tsec.$$	
    
)
}
