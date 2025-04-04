adjust_format(){

  input_nodelist=${1}

  if [[ "$input_nodelist" =~ "[" ]]; then
    cnode=$(echo $input_nodelist | cut -d "[" -f 1)
    output_nodelist=""
    nodenums=$(echo $input_nodelist | cut -d "[" -f 2 | cut -d "]" -f 1)
    for nodenum in $(echo $nodenums | sed -e "s#\,#\ #g"); do
      echo $nodenum >> check.log
      if [[ "$nodenum" =~ "-" ]]; then
        inodes=$(echo $nodenum | cut -d "-" -f 1)
        inodee=$(echo $nodenum | cut -d "-" -f 2)
      else
        inodes=$nodenum
        inodee=$nodenum
      fi
      cnum=${#inodes}
      for i in $(seq -f %$(printf %02g ${cnum})g $((10#${inodes})) $((10#${inodee})) ) ;do
        if [ "$output_nodelist" == "" ]; then
          output_nodelist=$cnode$inodes
        else
          output_nodelist=$output_nodelist","$cnode$i
        fi
      done
    done
  else
    output_nodelist="$input_nodelist"
  fi

  echo $output_nodelist 
}

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

        TOT_CORES=$(( 10#$INODE * 10#$ICORE ))                        #Total number of available cores
        MAX_JOBS=$(( 10#$TOT_CORES / 10#$QPROC ))                     #Maximum number of simultaneous jobs

        if [ $QPROC -gt $ICORE ] ; then
           #Round QPROC so that the number of nodes per job is integer.
           QPROC=$(( 10#$ICORE * ( 10#$QPROC / 10#$ICORE ) ))
        fi
        if [ $(( $QPROC % $QSKIP )) -ne 0 ] ; then
           echo "Error: QSKIP=$QSKIP must divide QPROC=$QPROC."
           echo "Check the configuration in machine.conf"
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

#        NPCORE=0        #Counter for the number of cores on current job
#        NJOB=1          #Counter for the number of jobs.
#        NMEM=$ini_mem   #Counter for the ensemble member
#        rm -fr machine.*
#        while [ $NMEM -le $end_mem ]; do
#            NCORE=$(( ($NJOB-1)*( 10#$QPROC ) + ( 10#$NPCORE ) ))
#            SMEM=$(printf "$mem_print_format" $((10#$NMEM)))
#            echo "(${NODES[$NCORE]}) core=$QSKIP" >> ./machine.$SMEM
#            NPCORE=$(($NPCORE + $QSKIP ))
#            if [ $NPCORE -ge $QPROC ] ; then
#               NMEM=$(($NMEM + 1 ))
#               NJOB=$(($NJOB + 1 ))
#               NPCORE=0
#            fi
#            if [ $NJOB -gt $MAX_JOBS ] ; then
#               NPCORE=0
#               NJOB=1
#            fi
#        done

	IFS=', ' read -r -a NODES <<< $( adjust_format $SLURM_JOB_NODELIST )
	echo NODES="${NODES[@]}"
	TOT_CORES=$(( $ICORE * $INODE ))
        MAX_JOBS=$(( $TOT_CORES / $QPROC ))  #Floor rounding (bash default)
	echo MAX_JOBS = $MAX_JOBS  

	IPCORE=1        #Counter for the number of cores on current job
	IJOB=1          #Counter for the number of jobs.
	IMEM=$ini_mem  #Counter for the ensemble member
	IPNODE=0        #Counter for the node number in the current job
	IPPROCINNODE=1  #Counter of number of procs used in the current node.
	rm -fr machine.*
        while [ $IMEM -le $end_mem ]; do
	    echo "${NODES[${IPNODE}]} " >> machine.$(printf "$mem_print_format" $((10#$IMEM)))
	    echo $IPCORE , ${NODES[${IPNODE}]} , $IPCORE , $IJOB , $IMEM , machine.$(printf "$mem_print_format" $((10#$IMEM)))
	    IPCORE=$(($IPCORE + 1))
	    IPPROCINNODE=$(( $IPPROCINNODE + 1 ))
	    if [ $IPCORE -gt $QPROC ] ; then
               IMEM=$(($IMEM + 1 ))
	       IJOB=$(($IJOB + 1 ))
	       IPCORE=1
            fi
            if [ $IJOB -gt $MAX_JOBS ] ; then
               IPCORE=1
	       IJOB=1
	       IPNODE=0
	       IPPROCINNODE=1
            fi
	    if [ $IPPROCINNODE -gt $ICORE ] ; then #We reached the maximum number of cores for this node.
	       IPNODE=$(( $IPNODE + 1 ))
	    fi
        done  

        #2 - Create the scripts
        for IMEM in $(seq -w $ini_mem $end_mem) ; do
		echo "source $BASEDIR/conf/config.env"                                                           > ${QPROC_NAME}_${IMEM}.pbs 
                echo "source $BASEDIR/conf/machine.conf"                                                        >> ${QPROC_NAME}_${IMEM}.pbs
		echo "ERROR=0                        "                                                          >> ${QPROC_NAME}_${IMEM}.pbs
                if [ $QOMP -eq 1 ] ; then
                  echo "export OMP_NUM_THREADS=${QSKIP}"                                         >> ${QPROC_NAME}_${IMEM}.pbs
                  echo "export PARALLEL=${QSKIP}"                                                >> ${QPROC_NAME}_${IMEM}.pbs
                else
                  echo "export OMP_NUM_THREADS=1"                                                >> ${QPROC_NAME}_${IMEM}.pbs
                  echo "export PARALLEL=1"                                                       >> ${QPROC_NAME}_${IMEM}.pbs
                fi
                echo "MEM=$IMEM "                                                                             >> ${QPROC_NAME}_${IMEM}.pbs
	        echo "export MPIEXESERIAL=\"\$MPIEXEC -np 1 -machinefile ../machine.$IMEM -bootstrap slurm \" ">> ${QPROC_NAME}_${IMEM}.pbs
         	echo "export MPIEXE=\"mpiexec -np ${QPROC} -machinefile ../machine.$IMEM -bootstrap slurm  \" ">> ${QPROC_NAME}_${IMEM}.pbs                   ## Comando MPIRUN con cantidad de nodos y cores por nodos           
	       	echo "mkdir -p ${QWORKPATH}/${IMEM}"                                                           >> ${QPROC_NAME}_${IMEM}.pbs
	       	echo "cd ${QWORKPATH}/${IMEM}"                                                                 >> ${QPROC_NAME}_${IMEM}.pbs
	        echo "${QSCRIPTCMD}"                                                                            >> ${QPROC_NAME}_${IMEM}.pbs
	        echo "if [[ -z \${ERROR} ]] || [[ \${ERROR} -eq 0 ]] ; then"                                    >> ${QPROC_NAME}_${IMEM}.pbs
	        echo "touch $PROCSDIR/${QPROC_NAME}_${IMEM}_ENDOK  "                                           >> ${QPROC_NAME}_${IMEM}.pbs  #Si existe la variable RES en el script la usamos
	        echo "fi                                            "                                           >> ${QPROC_NAME}_${IMEM}.pbs
        done

        #3 - Run the scripts
        IJOB=1     #Counter for the number of running jobs;
        IMEM=$ini_mem
        while [ $IMEM -le $end_mem ] ; do
	    MEMBER=$(printf "$mem_print_format" $IMEM)
	    bash ${QPROC_NAME}_${MEMBER}.pbs &> ${QPROC_NAME}_${MEMBER}.out  &
            IJOB=$(($IJOB + 1))
	    IMEM=$(($IMEM + 1))
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





