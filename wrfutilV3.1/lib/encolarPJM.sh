queue (){
        #IMPORTANT NOTE: this is for FUGAKU. We will assume that each ensemble member use an integer number of nodes to run (at least one).
	#IMPORTANT NOTE 2: There is a limit of 128 simultaneous mpi execution jobs. 

	#In the context of this function:
        cd $QWORKPATH

	ini_mem=${1}
        end_mem=${2}
	mem_print_format='%0'${#ini_mem}'d'
	
	TOT_CORES=$(( 10#$INODE * 10#$ICORE ))                        #Total number of available cores
        MAX_JOBS=$(( 10#$TOT_CORES / 10#$QPROC ))                     #Maximum number of simultaneous jobs

        if [ $QPROC -gt $ICORE ] ; then
           #Round QPROC so that the number of nodes per job is integer.
	   QPROC=$(( 10#$ICORE * ( 10#$QPROC / 10#$ICORE ) ))
        fi
	if [ $(( $QPROC % $QSKIP )) -ne 0 ] ; then 
           echo "Error: QSKIP=$QSKIP must divide QPROC=$QPROC."
	   echo "Check the configuration in machine.con"
	   exit 1
	fi

        echo MAX_JOBS = $MAX_JOBS 
	if [ $MAX_JOBS -gt 128 ] ; then 
	  echo "WARNING: Maximum number of simultaneous runs is 128, this may produce problems!!!" 
        fi 

        NCORE=0
        NNODE=0
        NODES=()
        while [ $NCORE -le $TOT_CORES ] ; do
           NODES+=( $NNODE )
           NCORE=$(($NCORE+1))
           if [ $(( NCORE % ICORE )) -eq 0 ] ; then
              NNODE=$(($NNODE+1))
           fi
        done

        NPCORE=0        #Counter for the number of cores on current job
        NJOB=1          #Counter for the number of jobs.
        NMEM=$ini_mem   #Counter for the ensemble member
        rm -fr machine.*
        while [ $NMEM -le $end_mem ]; do
            NCORE=$(( ($NJOB-1)*( 10#$QPROC ) + ( 10#$NPCORE ) ))
            SMEM=$(printf "$mem_print_format" $((10#$NMEM)))
            echo "(${NODES[$NCORE]}) core=1" >> ./machine.$SMEM
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


        #Set environment and go to workdir for this member.
        for IMEM in $(seq -w $ini_mem $end_mem) ; do	  
           echo "#!/bin/bash  "                                                               > ${QPROC_NAME}_${IMEM}.pbs 
           if [ $QOMP -eq 1 ] ; then 
	      echo "export OMP_NUM_THREADS=${QSKIP}"                                         >> ${QPROC_NAME}_${IMEM}.pbs
	      echo "export PARALLEL=${QSKIP}"                                                >> ${QPROC_NAME}_${IMEM}.pbs
           else 
              echo "export OMP_NUM_THREADS=1"                                                >> ${QPROC_NAME}_${IMEM}.pbs
	      echo "export PARALLEL=1"                                                       >> ${QPROC_NAME}_${IMEM}.pbs
           fi
	   echo "export MEM=$IMEM              "                                             >> ${QPROC_NAME}_${IMEM}.pbs 
           echo "source $BASEDIR/conf/config.env"                                            >> ${QPROC_NAME}_${IMEM}.pbs
           echo "source $BASEDIR/conf/machine.conf"                                          >> ${QPROC_NAME}_${IMEM}.pbs
	   echo "mkdir -p ${QWORKPATH}/\${MEM}"                                              >> ${QPROC_NAME}_${IMEM}.pbs
	   echo "cd ${QWORKPATH}/\${MEM}"                                                    >> ${QPROC_NAME}_${IMEM}.pbs
           #Initialize error count and set MPI execution commands.
           echo "ERROR=0                          "                                          >> ${QPROC_NAME}_${IMEM}.pbs
           echo "export MPIEXESERIAL=\"\$MPIEXEC -np 1 -vcoordfile ${QWORKPATH}/machine.\${MEM} \" "  >> ${QPROC_NAME}_${IMEM}.pbs
           echo "export MPIEXE=\"\$MPIEXEC             -vcoordfile ${QWORKPATH}/machine.\${MEM} \" "  >> ${QPROC_NAME}_${IMEM}.pbs           
	   #Add all job specific commands.
           echo "${QSCRIPTCMD}"                                                              >> ${QPROC_NAME}_${IMEM}.pbs
	   #Indicate if this JOB was successful or not
           echo "if [[ -z \${ERROR} ]] || [[ \${ERROR} -eq 0 ]] ; then"                      >> ${QPROC_NAME}_${IMEM}.pbs
           echo "touch $PROCSDIR/${QPROC_NAME}_\${MEM}_ENDOK  "                              >> ${QPROC_NAME}_${IMEM}.pbs  
           echo "fi 	                                    "                                >> ${QPROC_NAME}_${IMEM}.pbs
        done
      
        #3 - Run the scripts
        IJOB=1     #Counter for the number of running jobs;
	for IMEM in $(seq -w $ini_mem $end_mem ) ; do
	    echo "Launching ${QPROC_NAME}_${IMEM}.pbs"
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





