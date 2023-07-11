queue (){
	ini_mem=${1}
        end_mem=${2}

	ens_size=$(($end_mem-$ini_mem+1))
	#Get the number of available nodes and the number of available procs per nods.
	#Construct the machine files, write the scripts and run them. 

        #1 - Create machine files
	ICORE=1
	NODES+='null'  
        while read mynode ; do
	   NODES+=( $mynode ) ; ICORE=$(($ICORE+1))
        done < $PBS_NODEFILE
	TOT_CORES=$(($ICORE-1))
        MAX_JOBS=$(( $TOT_CORES / $TPROC ))  #Floor rounding (bash default)
	echo MAX_JOBS = $MAX_JOBS  

	IPCORE=1        #Counter for the number of cores on current job
	IJOB=1          #Counter for the number of jobs.
	QMIEM=$ini_mem  #Counter for the ensemble member
	rm -fr machine.*
	echo $NODES
        while [ $QMIEM -le $end_mem ]; do
	    ICORE=$(( ($IJOB-1)*$QPROC + $IPCORE ))
	    echo "${NODES[${ICORE}]} " >> machine.$(printf "%02d" $QMIEM)
	    echo $ICORE , ${NODES[${ICORE}]} , $IPCORE , $IJOB , $QMIEM
	    IPCORE=$(($IPCORE + 1))
	    if [ $IPCORE -gt $QPROC ] ; then
               QMIEM=$(($QMIEM + 1 ))
	       IJOB=$(($IJOB + 1 ))
	       IPCORE=1
            fi
            if [ $IJOB -gt $MAX_JOBS ] ; then
               QMIEM=$(($QMIEM + 1 ))
               IPCORE=1
	       IJOB=1
            fi
        done  

        #2 - Create the scripts
        for QMIEM in $(seq -w $ini_mem $end_mem) ; do
		echo "source $BASEDIR/conf/config.env"                 > ${QPROC_NAME}_${QMIEM}.pbs 
         	echo "source $BASEDIR/conf/experimento.conf"          >> ${QPROC_NAME}_${QMIEM}.pbs 
		echo "$ENVSET  "                                      >> ${QPROC_NAME}_${QMIEM}.pbs
                echo "MIEM=$QMIEM "                                   >> ${QPROC_NAME}_${QMIEM}.pbs
         	echo "export MPIEXE=\"mpiexec -np ${QPROC} -machinefile ../machine.$QMIEM \" ">> ${QPROC_NAME}_${QMIEM}.pbs                   ## Comando MPIRUN con cantidad de nodos y cores por nodos
	        echo "${QSCRIPTCMD}"                                  >> ${QPROC_NAME}_${QMIEM}.pbs
	        echo "if [[ -z \"${RES}\" ]] || [[ \"${RES}\" -eq "OK" ]] ; then" >> ${QPROC_NAME}_${QMIEM}.pbs
	        echo "touch $PROCSDIR/${QPROC_NAME}_${QMIEM}_ENDOK  " >> ${QPROC_NAME}_${QMIEM}.pbs  #Si existe la variable RES en el script la usamos
	        echo "fi                                            " >> ${QPROC_NAME}_${QMIEM}.pbs
        done

        #3 - Run the scripts
        IJOB=1     #Counter for the number of running jobs;
        QMIEM=$ini_mem
        while [ $QMIEM -le $end_mem ] ; do
	    MEMBER=$(printf "%02d" $QMIEM)
	    bash ${QPROC_NAME}_${MEMBER}.pbs > ${QPROC_NAME}_${MEMBER}.out &
            IJOB=$(($IJOB + 1))
	    QMIEM=$(($QMIEM + 1))
            if [ $IJOB -gt $MAX_JOBS ] ; then
	       time wait 	    
               IJOB=1
            fi
        done
}

check_proc(){
       
    ini_mem=${1}
    end_mem=${2}
    nmem=$(($end_mem-$ini_mem+1))
    check=0
    while [ $check -ne $nmem ] ; do 
       check=0
       for cmiem in $(seq -w $ini_mem $end_mem ) ; do
          if [ -e $PROCSDIR/${QPROC_NAME}_${cmiem}_ENDOK ] ; then
             check=$(($check+1))
          else
             sleep 10 #There is at least one missing member ... wait.
          fi
       done
    done
}





