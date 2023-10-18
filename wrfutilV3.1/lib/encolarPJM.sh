queue (){
        #IMPORTANT NOTE: this is for FUGAKU. We will assume that each ensemble member use an integer number of nodes to run (at least one).
	#IMPORTANT NOTE 2: There is a limit of 128 simultaneous mpi execution jobs. 

	#In the context of this function:
	#QNODE is the number of nodes associated to each ensemble member (at least 1).
	#TOT_CORES will be computed as the number of available cores  (number of lines in PBS_NODEFILE)
	#MAX_JOBS will be computed as the maximum number of simultaneous jobs ( floor( TOT_CORES / QPROC ) )
	#Groups of up to MAX_JOBS runs will be executed until all the ensemble members are processed. 
        cd $QWORKPATH

	ini_mem=${1}
        end_mem=${2}

	ens_size=$(($end_mem-$ini_mem+1))
        MAX_JOBS=$(( $INODE / $QNODE ))  #Floor rounding (bash default)
	TOT_PROCS=$(( $QPROC * $QNODE ))
        echo MAX_JOBS = $MAX_JOBS 
	if [ $MAX_JOBS -gt 128 ] ; then 
	  echo "WARNING: Maximum number of simultaneous runs is 128, this may produce problems!!!" 
        fi 
 
	#Create the machine files (in MPI VCOORD FILE FUJITSU FORMAT)
	NNODE=0          #Starting node
	IJOB=0           #Starint job
	QMIEM=$ini_mem   #Starting ensemble member
	if [ $INODE -gt 1  ] ; then 
           QPROC=$ICORE #If an application is using more than one node we force it to use all the cores within each node. 
        fi
	rm -fr machine.*
	while [ $QMIEM -le $end_mem ]; do
	   for MY_NODE in $(seq -w $NNODE $(($NNODE + QNODE - 1)) ) ; do
	       for MY_CORE in $(seq -w 1 $QPROC ) ; do
	          echo "(${MY_NODE}) core=1" >>  machine.$(printf "%02d" $QMIEM)
	       done
           done
	   IJOB=$(( $IJOB + 1 ))
	   NNODE=$(( $NNODE + $QNODE )) 
	   QMIEM=$(( $QMIEM + 1 ))
	   if [ $IJOB -ge $MAX_JOBS ] ; then 
	      NNODE=0 
	      IJOB=0 
           fi
        done

        
        #2 - Create the scripts
        for QMIEM in $(seq -w $ini_mem $end_mem) ; do
		echo "export PARALLEL=1              "                                             > ${QPROC_NAME}_${QMIEM}.pbs
		echo "source $BASEDIR/conf/config.env"                                            >> ${QPROC_NAME}_${QMIEM}.pbs 
                echo "source $BASEDIR/conf/machine.conf"                                          >> ${QPROC_NAME}_${QMIEM}.pbs
                echo "source $BASEDIR/lib/errores.env"                                            >> ${QPROC_NAME}_${QMIEM}.pbs
                echo "source $BASEDIR/conf/$QCONF      "                                          >> ${QPROC_NAME}_${QMIEM}.pbs
		test $QTHREAD  && echo "export OMP_NUM_THREADS=${QTHREAD}"                        >> ${QPROC_NAME}_${QMIEM}.pbs
                echo "MIEM=$QMIEM "                                                               >> ${QPROC_NAME}_${QMIEM}.pbs
                echo "export MPIEXESERIAL=\"\$MPIEXEC -np 1 -vcoordfile ../machine.$QMIEM \"  "   >> ${QPROC_NAME}_${QMIEM}.pbs
         	echo "export MPIEXE=\"\$MPIEXEC                -vcoordfile ../machine.$QMIEM \"  ">> ${QPROC_NAME}_${QMIEM}.pbs ## Comando MPIRUN con cantidad de nodos y cores por nodos           
	       	test $QWORKPATH &&  echo "cd ${QWORKPAH}/${QMIEM}"                                 >> ${QPROC_NAME}_${QMIEM}.pbs
	        echo "${QSCRIPTCMD}"                                                              >> ${QPROC_NAME}_${QMIEM}.pbs
	        echo "if [[ -z \${res} ]] || [[ \${res} -eq "OK" ]] ; then"                       >> ${QPROC_NAME}_${QMIEM}.pbs
	        echo "touch $PROCSDIR/${QPROC_NAME}_${QMIEM}_ENDOK  "                             >> ${QPROC_NAME}_${QMIEM}.pbs  #Si existe la variable RES en el script la usamos
	        echo "fi                                            "                             >> ${QPROC_NAME}_${QMIEM}.pbs
        done

        #3 - Run the scripts
        IJOB=1     #Counter for the number of running jobs;
        QMIEM=$ini_mem
        while [ $QMIEM -le $end_mem ] ; do
	    MEMBER=$(printf "%02d" $QMIEM)
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





