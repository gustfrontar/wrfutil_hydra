queue (){
        ini_mem=${1}
        end_mem=${2}
	ens_num_length=${#ini_mem}
	cd $QWORKPATH
        TPROC=$(( 10#$QNODE * 10#$QPROC ))  #Total number of cores to be used. 	
	for IMIEM in $(seq -w $ini_mem $end_mem ) ; do 
        if [ -z ${QSKIP} ]; then #If QSKIP is not set. assume is one.
           echo "WARNING: QSKIP is not set, assuming es 1 " $QSKIP
           QSKIP=1
        fi
        #Optimally computed QTHREAD
        QTHREAD=$QSKIP
        SPROC=$(( 10#$TPROC / 10#$QSKIP ))  
 



	                        echo "#PBS -j oe "                                    >  ${QPROC_NAME}_${IMIEM}.pbs                   ## Redireccion de flujos de salida y error
	test $QUEUE && 		echo "#PBS -q $QUEUE "                                >> ${QPROC_NAME}_${IMIEM}.pbs                   ## Cola a la que se lo va a encolar
        test $QWALLTIME && 	echo "#PBS -l walltime=${QWALLTIME} "                 >> ${QPROC_NAME}_${IMIEM}.pbs	              ## Tiempo estimado de corrida
                                echo "#PBS -l nodes=${QNODE}:ppn=${QPROC}"            >> ${QPROC_NAME}_${IMIEM}.pbs                   ## Recursos que seran utilizados 
        			echo '#PBS -v BASEDIR'                                >> ${QPROC_NAME}_${IMIEM}.pbs		      ## Indica que debe heredar todo el enviroment
	test $QTHREAD   &&	echo "export OMP_NUM_THREADS=$QTHREAD"                >> ${QPROC_NAME}_${IMIEM}.pbs		
				echo "source $BASEDIR/conf/config.env"                >> ${QPROC_NAME}_${IMIEM}.pbs 
                                echo "source $BASEDIR/conf/machine.conf"              >> ${QPROC_NAME}_${IMIEM}.pbs
	                        echo "source $BASEDIR/lib/errores.env"                >> ${QPROC_NAME}_${IMIEM}.pbs
				echo "source $BASEDIR/conf/$QCONF"                    >> ${QPROC_NAME}_${IMIEM}.pbs                   ## Experiment specific configuration file.
				echo "ERROR=0                    "                    >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "MIEM=$IMIEM "                                   >> ${QPROC_NAME}_${IMIEM}.pbs
				echo "export MPIEXESERIAL=\"$MPIEXEC -np 1\" "        >> ${QPROC_NAME}_${IMIEM}.pbs
                		echo "export MPIEXE=\"$MPIEXEC -machinefile=./machine.${IMIEM} -np ${SPROC}\" " >> ${QPROC_NAME}_${IMIEM}.pbs                   ## Comando MPIRUN con cantidad de nodos y cores por nodos
                      test $QWORKPATH &&  echo "mkdir ${QWORKPATH}/${IMIEM}"          >> ${QPROC_NAME}_${IMIEM}.pbs
                      test $QWORKPATH &&  echo "cd ${QWORKPATH}/${IMIEM}"             >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "NODES+='null'                              "    >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "while read mynode ; do                     "    >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "  NODES+=( \$mynode )                      "    >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "done < \$PBS_NODEFILE                      "    >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "TPROC=$TPROC                               "    >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "SKIP=$QSKIP                                "    >> ${QPROC_NAME}_${IMIEM}.pbs    
                                echo "IPCORE=1                                   "    >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "rm -fr machine.${IMIEM}                    "    >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "while [ \$IPCORE -le $TPROC ]; do          "    >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo " echo "\${NODES[\${IPCORE}]} " >> machine.${IMIEM} " >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo " IPCORE=\$(( 10#\$IPCORE + 10#\$SKIP ))    "    >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "done                                       "    >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "export FORT90L=\"\"                        "    >> ${QPROC_NAME}_${IMIEM}.pbs         
                                echo "export KMP_STACKSIZE=1G                    "    >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "export OMP_STACKSIZE=1G                    "    >> ${QPROC_NAME}_${IMIEM}.pbs
	                        echo "${QSCRIPTCMD}"                                  >> ${QPROC_NAME}_${IMIEM}.pbs
				echo "if   [[ -z \${ERROR} ]] || [[ \${ERROR} -eq 0 ]] ; then " >> ${QPROC_NAME}_${IMIEM}.pbs
	                        echo "touch $PROCSDIR/${QPROC_NAME}_${IMIEM}_ENDOK  " >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "elif [[ -z \${ERROR} ]] || [[ \${ERROR} -gt 0 ]] ; then " >> ${QPROC_NAME}_${IMIEM}.pbs
                                echo "touch $PROCSDIR/${QPROC_NAME}_${IMIEM}_ERROR  " >> ${QPROC_NAME}_${IMIEM}.pbs
				echo "fi 				            " >> ${QPROC_NAME}_${IMIEM}.pbs

        # Encolamiento 
        ${QSUB_ROOT}qsub ${QPROC_NAME}_${IMIEM}.pbs   

        done
}

check_proc(){
       
    ini_mem=${1}
    end_mem=${2}
    nmem=$(( $((10#$end_mem))-$((10#$ini_mem))+1))
    check=0
    while [ $check -ne $nmem ] ; do 
       check=0
       for cmiem in $(seq -w $ini_mem $end_mem ) ; do
          if [ -e $PROCSDIR/${QPROC_NAME}_${cmiem}_ENDOK ] ; then
             check=$(($check+1))
          elif [ -e $PROCSDIR/${QPROC_NAME}_${cmiem}_ERROR ] ; then
             echo "Error: ${QPROC_NAME}_${cmiem} finished with errors"
             echo "Error: Aborting this step"
             exit 1
          else
             sleep 10 #There is at least one missing member ... wait.
          fi
       done
    done
}


name2PID(){
        for file in $(ls  ${BASEDIR}/PROCS/${QDEPEND_NAME} 2>/dev/null)
        do
                jid=$(cat $file 2>/dev/null)
                estado=$(qstat -fx $jid 2> /dev/null |grep job_state |grep -E "R|H|Q" )
                #estado=$(qstat -fx $jid |grep job_state |grep -E "R|H|Q" )
                if ! [[ -z $estado ]]
                then
                         echo $jid
                fi
        done
}

function join_by { local IFS=${LIFS};  echo "$*"; }

