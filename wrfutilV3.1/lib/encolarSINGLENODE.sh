queue (){
        #IMPORTANT NOTE: This is to run in a single node / server. No machine file will be generated. 

	#In the context of this function:
	#QNODE is the number of nodes should be enforced to be equal to one.
	#TOT_CORES will be computed as the number of available cores  (number of lines in PBS_NODEFILE)
	#MAX_JOBS will be computed as the maximum number of simultaneous jobs ( floor( TOT_CORES / QPROC ) )
	#Groups of up to MAX_JOBS runs will be executed until all the ensemble members are processed. 
        cd $QWORKPATH

	QNODE=1 #Single node version of the function.

	ini_mem=${1}
        end_mem=${2}
	mem_print_format='%0'${#ini_mem}'d'

        MAX_JOBS=$(( $ICORE / $QPROC ))  #Floor rounding (bash default)
	TOT_PROCS=$(( $QPROC ))
        echo MAX_JOBS = $MAX_JOBS 
	if [ $MAX_JOBS -gt 128 ] ; then 
	  echo "WARNING: Maximum number of simultaneous runs is 128, this may produce problems!!!" 
        fi 
 
        #2 - Create the scripts
        for IMIEM in $(seq -w $ini_mem $end_mem) ; do
		echo "export PARALLEL=1              "                                             > ${QPROC_NAME}_${IMIEM}.pbs
		echo "source $BASEDIR/conf/config.env"                                            >> ${QPROC_NAME}_${IMIEM}.pbs 
                echo "source $BASEDIR/conf/machine.conf"                                          >> ${QPROC_NAME}_${IMIEM}.pbs
                echo "source $BASEDIR/lib/errores.env"                                            >> ${QPROC_NAME}_${IMIEM}.pbs
                echo "source $BASEDIR/conf/$QCONF      "                                          >> ${QPROC_NAME}_${IMIEM}.pbs
		test $QTHREAD  && echo "export OMP_NUM_THREADS=${QTHREAD}"                        >> ${QPROC_NAME}_${IMIEM}.pbs
                echo "MIEM=$IMIEM "                                                               >> ${QPROC_NAME}_${IMIEM}.pbs
                echo "export MPIEXESERIAL=\"\$MPIEXEC -np 1                                  \"  ">> ${QPROC_NAME}_${IMIEM}.pbs
         	echo "export MPIEXE=\"\$MPIEXEC -np $QPROC                                   \"  ">> ${QPROC_NAME}_${IMIEM}.pbs            
	       	test $QWORKPATH &&  echo "cd ${QWORKPATH}/${IMIEM}"                               >> ${QPROC_NAME}_${IMIEM}.pbs
	        echo "${QSCRIPTCMD}"                                                              >> ${QPROC_NAME}_${IMIEM}.pbs
	        echo "if [[ -z \${res} ]] || [[ \${res} -eq "OK" ]] ; then"                       >> ${QPROC_NAME}_${IMIEM}.pbs
	        echo "touch $PROCSDIR/${QPROC_NAME}_${IMIEM}_ENDOK  "                             >> ${QPROC_NAME}_${IMIEM}.pbs  
	        echo "fi                                            "                             >> ${QPROC_NAME}_${IMIEM}.pbs
        done

        #3 - Run the scripts
        IJOB=1     #Counter for the number of running jobs;
        IMIEM=$ini_mem
        while [ $IMIEM -le $end_mem ] ; do
	    MEMBER=$(printf "$mem_print_format" $IMIEM)
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





