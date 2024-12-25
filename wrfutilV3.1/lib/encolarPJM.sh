queue (){
        #IMPORTANT NOTE: this is for FUGAKU. We will assume that each ensemble member use an integer number of nodes to run (at least one).
	#IMPORTANT NOTE 2: There is a limit of 128 simultaneous mpi execution jobs. 

	#In the context of this function:
        cd $QWORKPATH

	ini_mem=${1}
        end_mem=${2}
	mem_print_format='%0'${#ini_mem}'d'

        echo MAX_JOBS = $MAX_JOBS 
	if [ $MAX_JOBS -gt 128 ] ; then 
	  echo "WARNING: Maximum number of simultaneous runs is 128, this may produce problems!!!" 
        fi 

        TOT_CORES=$(($INODE*$ICORE))                        #Total number of available cores
        MAX_JOBS=$(( $TOT_CORES / ( $QPROC * $QTHREAD ) ))  #Maximum number of simultaneous jobs

        NCORE=1
        NNODE=0
        NODES+='null'
        while [ $NCORE -le $TOT_CORES ] ; do
           NODES+=( $NNODE ) ; NCORE=$(($NCORE+1))
           if [ $(( $NCORE % $ICORE )) == 0 ] ; then
              NNODE=$(($NNODE+1))
           fi
        done 

        NPCORE=1        #Counter for the number of cores on current job
        NJOB=1          #Counter for the number of jobs.
        NMEM=$ini_mem   #Counter for the ensemble member
        rm -fr machine.*
        while [ $NMEM -le $end_mem ]; do
            NCORE=$(( ($NJOB-1)*( 10#$QPROC * 10#$QTHREAD) + ( 10#$NPCORE * 10#$QTHREAD ) ))
            SMEM=$(printf "$mem_print_format" $((10#$NMEM)))
            echo "(${NODES[$($NCORE)]}) core=1 " >> ${QWORKPATH}/machine.$SMEM
            NPCORE=$(($NPCORE +1 ))
            if [ $NPCORE -gt $QPROC ] ; then
               NMEM=$(($NMEM + 1 ))
               NJOB=$(($NJOB + 1 ))
               NPCORE=1
            fi
            if [ $NJOB -gt $MAX_JOBS ] ; then
               NPCORE=1
               NJOB=1
            fi
        done 

        #Set environment and go to workdir for this member.	
        echo "export PARALLEL=1              "                                             > ${QPROC_NAME}.pbs
	echo "export OMP_NUM_THREADS=${QTHREAD}"                                          >> ${QPROC_NAME}.pbs
	echo "export MEM=\$1                "                                            >> ${QPROC_NAME}.pbs #The ensemble member will be an input to the script.
        echo "source $BASEDIR/conf/config.env"                                            >> ${QPROC_NAME}.pbs
        echo "source $BASEDIR/conf/machine.conf"                                          >> ${QPROC_NAME}.pbs
        echo "source $BASEDIR/lib/errores.env"                                            >> ${QPROC_NAME}.pbs
        echo "source $BASEDIR/conf/$QCONF      "                                          >> ${QPROC_NAME}.pbs
	echo "mkdir -p ${QWORKPATH}/\${MEM}"                                             >> ${QPROC_NAME}.pbs
	echo "cd ${QWORKPATH}/\${MEM}"                                                   >> ${QPROC_NAME}.pbs
        #Initialize error count and set MPI execution commands.
        echo "ERROR=0                          "                                          >> ${QPROC_NAME}.pbs
        echo "export MPIEXESERIAL=\"\$MPIEXEC -np 1 -vcoordfile ${QWORKPATH}/machine.\${MEM} \" "  >> ${QPROC_NAME}.pbs
        echo "export MPIEXE=\"\$MPIEXEC             -vcoordfile ${QWORKPATH}/machine.\${MEM} \" "  >> ${QPROC_NAME}.pbs ## Comando MPIRUN con cantidad de nodos y cores por nodos           
	#Add all job specific commands.
        echo "${QSCRIPTCMD}"                                                              >> ${QPROC_NAME}.pbs
	#Indicate if this JOB was successful or not
        echo "if [[ -z \${ERROR} ]] || [[ \${ERROR} -eq 0 ]] ; then"                      >> ${QPROC_NAME}.pbs
        echo "touch $PROCSDIR/${QPROC_NAME}_\${MEM}_ENDOK  "                             >> ${QPROC_NAME}.pbs  #Si existe la variable RES en el script la usamos
        echo "fi 	                                    "                             >> ${QPROC_NAME}.pbs

        #3 - Run the scripts
        IJOB=1     #Counter for the number of running jobs;
	for IMEM in $(seq -w $ini_mem $end_mem ) ; do
	    bash ${QPROC_NAME}.pbs ${IMEM} &> ${QPROC_NAME}_${IMEM}.out   &
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





