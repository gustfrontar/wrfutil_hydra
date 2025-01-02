queue (){
	#In the context of this function:
	#QPROC means how many cores do we need for each ensemble member (or for the job in case a unique member is indicated)
	#TOT_CORES will be computed as the number of available cores  (number of lines in PBS_NODEFILE)
	#MAX_JOBS will be computed as the maximum number of simultaneous jobs ( floor( TOT_CORES / QPROC ) )
	#Groups of up to MAX_JOBS runs will be executed until all the ensemble members are processed. 
        cd $QWORKPATH

	ini_mem=${1}
        end_mem=${2}
        mem_print_format='%0'${#ini_mem}'d'

	#Get the number of available nodes and the number of available procs per nods.
	#Construct the machine files, write the scripts and run them. 

        #1 - Create machine files
	NODES=() 
        while read mynode ; do
	   NODES+=( $mynode )
        done < $PBS_NODEFILE
	TOT_CORES=${#NODES[@]}  #Get the total number of cores
        MAX_JOBS=$(( 10#$TOT_CORES / 10#$QPROC ))  #Floor rounding (bash default)

	echo MAX_JOBS = $MAX_JOBS  

        if [ $QPROC -gt $ICORE ] ; then
           #Round QPROC so that the number of nodes per job is integer.
           QPROC=$(( 10#$ICORE * ( 10#$QPROC / 10#$ICORE ) ))
        fi
        if [ $(( $QPROC % $QSKIP )) -eq 0 ] ; then
           echo "Error: QSKIP=$QSKIP must divide QPROC=$QPROC."
           echo "Check the configuration in machine.con"
           exit 1
        fi

	NPCORE=0        #Counter for the number of cores on current job
	NJOB=1          #Counter for the number of jobs.
	NMEM=$ini_mem   #Counter for the ensemble member
	rm -fr machine.*
        while [ $NMEM -le $end_mem ]; do
	    NCORE=$(( ($NJOB-1)*( 10#$QPROC ) + ( 10#$NPCORE ) ))
	    echo "${NODES[${NCORE}]} " >> machine.$(printf "$mem_print_format" $((10#$NMEM)))
	    NPCORE=$(($NPCORE + $QSKIP ))
	    if [ $NPCORE -ge $QPROC ] ; then
               NMEM=$(($NMEM + 1 ))
	       NJOB=$(($NJOB + 1 ))
	       NPCORE=0
            fi
            if [ $NJOB -gt $MAX_JOBS ] ; then
               NPCORE=0
	       NJOB=1
            fi
        done  

        #2 - Create the scripts
        for IMEM in $(seq -w $ini_mem $end_mem) ; do
		echo "source $BASEDIR/conf/config.env"                                             > ${QPROC_NAME}_${IMEM}.pbs 
                echo "source $BASEDIR/conf/machine.conf"                                          >> ${QPROC_NAME}_${IMEM}.pbs
		echo "ERROR=0                    "                                                >> ${QPROC_NAME}_${IMEM}.pbs
                if [ $QOMP -eq 1 ] ; then
		   echo "export OMP_NUM_THREADS=${QSKIP}"                                         >> ${QPROC_NAME}_${IMEM}.pbs
                else
                   echo "export OMP_NUM_THREADS=1"                                                >> ${QPROC_NAME}_${IMEM}.pbs
                fi  
                echo "MEM=$IMEM "                                                                 >> ${QPROC_NAME}_${IMEM}.pbs
                echo "export MPIEXESERIAL=\"\$MPIEXEC -np 1 -machinefile ../machine.$IMEM \"  "   >> ${QPROC_NAME}_${IMEM}.pbs
         	echo "export MPIEXE=\"\$MPIEXEC -np $(( 10#$QPROC / 10#$QSKIP )) -machinefile ../machine.$IMEM \" "   >> ${QPROC_NAME}_${IMEM}.pbs             
                echo "mkdir ${QWORKPATH}/${IMEM}"                                                 >> ${QPROC_NAME}_${IMEM}.pbs
	        echo "cd ${QWORKPATH}/${IMEM}"                                                    >> ${QPROC_NAME}_${IMEM}.pbs
	        echo "${QSCRIPTCMD}"                                                              >> ${QPROC_NAME}_${IMEM}.pbs
	        echo "if [[ -z \${ERROR} ]] || [[ \${ERROR} -eq 0 ]] ; then"                      >> ${QPROC_NAME}_${IMEM}.pbs
	        echo "touch $PROCSDIR/${QPROC_NAME}_${IMEM}_ENDOK  "                              >> ${QPROC_NAME}_${IMEM}.pbs  
	        echo "fi                                           "                              >> ${QPROC_NAME}_${IMEM}.pbs
        done

        #3 - Run the scripts
        IJOB=1     #Counter for the number of running jobs;
	for IMEM in $(seq -w $ini_mem $end_mem ) ; do
	    bash ${QPROC_NAME}_${IMEM}.pbs &> ${QPROC_NAME}_${IMEM}.out &
            IJOB=$(($IJOB + 1))
            if [ $IJOB -gt $MAX_JOBS ] ; then
	       time wait 	    
               IJOB=1
            fi
        done
        time wait
}

check_proc(){
       
    ini_mem=${1}
    end_mem=${2}
    nmem=$(( $((10#$end_mem))-$((10#$ini_mem))+1))
    check=0 
    for cmiem in $(seq -w $ini_mem $end_mem ) ; do
       if [ -e $PROCSDIR/${QPROC_NAME}_${cmiem}_ENDOK ] ; then
          check=$(($check+1))
       else
          echo "Member ${cmiem} finished with errors"       
       fi
    done
    if [ $check -eq $nmem  ] ; then
       echo "All members successfully processed"
    else 
       echo "Some members failed for ${QPROC_NAME}"
       echo "exiting .... "
       exit 1
    fi
    
}





