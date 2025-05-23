queue (){
        #IMPORTANT NOTE: This is to run in a single node / server. No machine file will be generated. 

	#In the context of this function:
	#QNODE/$INODE are assumed to be 1.
	#TOT_CORES will be computed as the number of available cores
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
                echo "export MPIEXESERIAL=\"\$MPIEXEC -np 1 \"  "                                 >> ${QPROC_NAME}_${IMEM}.pbs
         	echo "export MPIEXE=\"\$MPIEXEC -np $(( 10#$QPROC / 10#$QSKIP ))  \" "            >> ${QPROC_NAME}_${IMEM}.pbs             
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





