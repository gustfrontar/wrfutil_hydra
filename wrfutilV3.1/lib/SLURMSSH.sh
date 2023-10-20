queue (){
	#In the context of this function:
	#QPROC means how many cores do we need for each ensemble member (or for the job in case a unique member is indicated)
	#TOT_CORES will be computed as the number of available cores  (number of lines in PBS_NODEFILE)
	#MAX_JOBS will be computed as the maximum number of simultaneous jobs ( floor( TOT_CORES / QPROC ) )
	#Groups of up to MAX_JOBS runs will be executed until all the ensemble members are processed. 
        cd $QWORKPATH

        end_mem=${2}
	ens_size=$(($end_mem-$ini_mem+1))
	mem_print_format='%0'${#ini_mem}'d'
	#Get the number of available nodes and the number of available procs per nods.
	#Construct the machine files, write the scripts and run them. 

        #1 - Create machine files
	list=scontrol show hostname $SLURM_NODELIST  #Get assigned list of nodes.
	NODES=($list)
	TOT_CORES=$(( $ICORE * $INODE ))
        MAX_JOBS=$(( $TOT_CORES / $QPROC ))  #Floor rounding (bash default)
	echo MAX_JOBS = $MAX_JOBS  

	IPCORE=1        #Counter for the number of cores on current job
	IJOB=1          #Counter for the number of jobs.
	QMIEM=$ini_mem  #Counter for the ensemble member
	IPNODE=1        #Counter for the node number in the current job
	IPPROCINNODE=1  #Counter of number of procs used in the current node.
	rm -fr machine.*
	echo $NODES
        while [ $QMIEM -le $end_mem ]; do
	    echo "${NODES[${IPNODE}]} " >> machine.$(printf "$mem_print_format" $QMIEM)
	    echo $IPCORE , ${NODES[${IPNODE}]} , $IPCORE , $IJOB , $QMIEM , machine.$(printf "$mem_print_format" $QMIEM)
	    IPCORE=$(($IPCORE + 1))
	    IPPROCINNODE=$(( $IPPROCINNODE + 1 ))
	    if [ $IPCORE -gt $QPROC ] ; then
               QMIEM=$(($QMIEM + 1 ))
	       IJOB=$(($IJOB + 1 ))
	       IPCORE=1
            fi
            if [ $IJOB -gt $MAX_JOBS ] ; then
               IPCORE=1
	       IJOB=1
	       IPNODE=1
	       IPPROCINNODE=1
            fi
	    if [ $IPPROCINNODE -gt $ICORE ] ; then #We reached the maximum number of cores for this node.
	       IPNODE=$(( $IPNODE + 1 ))
	    fi
        done  

        #2 - Create the scripts
        for QMIEM in $(seq -w $ini_mem $end_mem) ; do
		echo "source $BASEDIR/conf/config.env"                                         > ${QPROC_NAME}_${QMIEM}.pbs 
                echo "source $BASEDIR/conf/machine.conf"                                      >> ${QPROC_NAME}_${QMIEM}.pbs
                echo "source $BASEDIR/conf/$QCONF"                                            >> ${QPROC_NAME}_${QMIEM}.pbs
		test $QTHREAD  && echo "export OMP_NUM_THREADS=${QTHREAD}"                    >> ${QPROC_NAME}_${QMIEM}.pbs
                echo "MIEM=$QMIEM "                                                           >> ${QPROC_NAME}_${QMIEM}.pbs
         	echo "export MPIEXE=\"mpiexec -np ${QPROC} -machinefile ../machine.$QMIEM \" ">> ${QPROC_NAME}_${QMIEM}.pbs                   ## Comando MPIRUN con cantidad de nodos y cores por nodos           
	       	test $QWORKDIR &&  echo "cd ${QWORKDIR}/${QMIEM}"                             >> ${QPROC_NAME}_${QMIEM}.pbs
	        echo "${QSCRIPTCMD}"                                                          >> ${QPROC_NAME}_${QMIEM}.pbs
	        echo "if [[ -z \${res} ]] || [[ \${res} -eq "OK" ]] ; then"                   >> ${QPROC_NAME}_${QMIEM}.pbs
	        echo "touch $PROCSDIR/${QPROC_NAME}_${QMIEM}_ENDOK  "                         >> ${QPROC_NAME}_${QMIEM}.pbs  #Si existe la variable RES en el script la usamos
	        echo "fi                                            "                         >> ${QPROC_NAME}_${QMIEM}.pbs
        done

        #3 - Run the scripts
        IJOB=1     #Counter for the number of running jobs;
        QMIEM=$ini_mem
        while [ $QMIEM -le $end_mem ] ; do
	    MEMBER=$(printf "$mem_print_format" $QMIEM)
	    bash ${QPROC_NAME}_${MEMBER}.pbs > ${QPROC_NAME}_${MEMBER}.out  2>&1  &
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





