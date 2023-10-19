queue (){
	#In the context of this function:
	#QPROC means how many cores do we need for each ensemble member (or for the job in case a unique member is indicated)
	#TOT_CORES will be computed as the number of available cores  (number of lines in PBS_NODEFILE)
	#MAX_JOBS will be computed as the maximum number of simultaneous jobs ( floor( TOT_CORES / QPROC ) )
	#Groups of up to MAX_JOBS runs will be executed until all the ensemble members are processed. 
        cd $QWORKPATH

	ini_mem=${1}
        end_mem=${2}
        ens_num_length=${#ini_mem}

	#Get the number of available nodes and the number of available procs per nods.
	#Construct the machine files, write the scripts and run them. 

        #1 - Create machine files
	ICORE=1
	NODES+='null'  
        while read mynode ; do
	   NODES+=( $mynode ) ; ICORE=$(($ICORE+1))
        done < $PBS_NODEFILE
	TOT_CORES=$(($ICORE-1))
        MAX_JOBS=$(( $TOT_CORES / $QPROC ))  #Floor rounding (bash default)
	echo MAX_JOBS = $MAX_JOBS  

	IPCORE=1        #Counter for the number of cores on current job
	IJOB=1          #Counter for the number of jobs.
	IMIEM=$ini_mem  #Counter for the ensemble member
	rm -fr machine.*
	echo $NODES
        while [ $IMIEM -le $end_mem ]; do
	    ICORE=$(( ($IJOB-1)*$QPROC + $IPCORE ))
	    echo "${NODES[${ICORE}]} " >> machine.$(printf '%0${ens_num_length}d' $((10#$IMIEM)))
	    echo $ICORE , ${NODES[${ICORE}]} , $IPCORE , $IJOB , $IMIEM , machine.$(printf '%0${ens_num_length}d' $((10#$IMIEM)))
	    IPCORE=$(($IPCORE + 1))
	    if [ $IPCORE -gt $QPROC ] ; then
               IMIEM=$(($IMIEM + 1 ))
	       IJOB=$(($IJOB + 1 ))
	       IPCORE=1
            fi
            if [ $IJOB -gt $MAX_JOBS ] ; then
               IPCORE=1
	       IJOB=1
            fi
        done  

        #2 - Create the scripts
        for IMIEM in $(seq -w $ini_mem $end_mem) ; do
		echo "source $BASEDIR/conf/config.env"                                             > ${QPROC_NAME}_${IMIEM}.pbs 
		echo "source $BASEDIR/lib/errores.env"                                            >> ${QPROC_NAME}_${IMIEM}.pbs
                echo "source $BASEDIR/conf/machine.conf"                                          >> ${QPROC_NAME}_${IMIEM}.pbs
                echo "source $BASEDIR/conf/$QCONF"                                                >> ${QPROC_NAME}_${IMIEM}.pbs
		test $QTHREAD  && echo "export OMP_NUM_THREADS=${QTHREAD}"                        >> ${QPROC_NAME}_${IMIEM}.pbs
                echo "MIEM=$IMIEM "                                                               >> ${QPROC_NAME}_${IMIEM}.pbs
                echo "export MPIEXESERIAL=\"\$MPIEXEC -np 1 -machinefile ../machine.$IMIEM \"  "  >> ${QPROC_NAME}_${IMIEM}.pbs
         	echo "export MPIEXE=\"mpiexec -np ${QPROC} -machinefile ../machine.$IMIEM \" "    >> ${QPROC_NAME}_${IMIEM}.pbs                   ## Comando MPIRUN con cantidad de nodos y cores por nodos           
	       	test $QWORKPATH &&  echo "cd ${QWORKPATH}/${IMIEM}"                                 >> ${QPROC_NAME}_${IMIEM}.pbs
	        echo "${QSCRIPTCMD}"                                                              >> ${QPROC_NAME}_${IMIEM}.pbs
	        echo "if [[ -z \${res} ]] || [[ \${res} -eq "OK" ]] ; then"                       >> ${QPROC_NAME}_${IMIEM}.pbs
	        echo "touch $PROCSDIR/${QPROC_NAME}_${IMIEM}_ENDOK  "                             >> ${QPROC_NAME}_${IMIEM}.pbs  #Si existe la variable RES en el script la usamos
	        echo "fi                                            "                             >> ${QPROC_NAME}_${IMIEM}.pbs
        done

        #3 - Run the scripts
        IJOB=1     #Counter for the number of running jobs;
        IMIEM=$ini_mem
        while [ $IMIEM -le $end_mem ] ; do
	    MEMBER=$(printf "%02d" $IMIEM)
	    bash ${QPROC_NAME}_${MEMBER}.pbs > ${QPROC_NAME}_${MEMBER}.out  2>&1  &
            IJOB=$(($IJOB + 1))
	    IMIEM=$(($IMIEM + 1))
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





