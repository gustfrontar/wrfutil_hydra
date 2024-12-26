queue (){
        ini_mem=${1}
        end_mem=${2}
	ens_num_length=${#ini_mem}
	cd $QWORKPATH

        TOT_CORES=$(( 10#$QNODE * 10#$ICORE ))     #Total number of cores to be used.
        SPROC=$(( 10#$TOT_CORES / 10#$QTHREAD ))   #Number of cores for MPI
         	
	for IMEM in $(seq -w $ini_mem $end_mem ) ; do 


	                        echo "#PBS -j oe "                                    >  ${QPROC_NAME}_${IMEM}.pbs                   ## Redireccion de flujos de salida y error
	test $QUEUE && 		echo "#PBS -q $QUEUE "                                >> ${QPROC_NAME}_${IMEM}.pbs                   ## Cola a la que se lo va a encolar
        test $QWALLTIME && 	echo "#PBS -l walltime=${QWALLTIME} "                 >> ${QPROC_NAME}_${IMEM}.pbs	              ## Tiempo estimado de corrida
                                echo "#PBS -l nodes=${QNODE}:ppn=${ICORE}"            >> ${QPROC_NAME}_${IMEM}.pbs                   ## Recursos que seran utilizados 
        			echo '#PBS -v BASEDIR'                                >> ${QPROC_NAME}_${IMEM}.pbs		      ## Indica que debe heredar todo el enviroment
	test $QTHREAD   &&	echo "export OMP_NUM_THREADS=$QTHREAD"                >> ${QPROC_NAME}_${IMEM}.pbs		
				echo "source $BASEDIR/conf/config.env"                >> ${QPROC_NAME}_${IMEM}.pbs 
                                echo "source $BASEDIR/conf/machine.conf"              >> ${QPROC_NAME}_${IMEM}.pbs
				echo "source $BASEDIR/conf/$QCONF"                    >> ${QPROC_NAME}_${IMEM}.pbs                   ## Experiment specific configuration file.
				echo "ERROR=0                    "                    >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "MEM=$IMEM "                                     >> ${QPROC_NAME}_${IMEM}.pbs
				echo "export MPIEXESERIAL=\"$MPIEXEC -np 1\" "        >> ${QPROC_NAME}_${IMEM}.pbs
                		echo "export MPIEXE=\"$MPIEXEC -machinefile=./machine.${IMEM} -np ${SPROC}\" " >> ${QPROC_NAME}_${IMEM}.pbs                   ## Comando MPIRUN con cantidad de nodos y cores por nodos
                      test $QWORKPATH &&  echo "mkdir ${QWORKPATH}/${IMEM}"           >> ${QPROC_NAME}_${IMEM}.pbs
                      test $QWORKPATH &&  echo "cd ${QWORKPATH}/${IMEM}"              >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "NODES+='null'                              "    >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "while read mynode ; do                     "    >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "  NODES+=( \$mynode )                      "    >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "done < \$PBS_NODEFILE                      "    >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "TPROC=$TOT_CORES                           "    >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "IPCORE=1                                   "    >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "rm -fr machine.${IMEM}                     "    >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "while [ \$IPCORE -le $TOT_CORES ]; do      "    >> ${QPROC_NAME}_${IMEM}.pbs
                                echo " echo "\${NODES[\${IPCORE}]} " >> machine.${IMEM} " >> ${QPROC_NAME}_${IMEM}.pbs
                                echo " IPCORE=\$(( 10#\$IPCORE + 10#\$QTHREAD )) "    >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "done                                       "    >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "export FORT90L=\"\"                        "    >> ${QPROC_NAME}_${IMEM}.pbs         
                                echo "export KMP_STACKSIZE=1G                    "    >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "export OMP_STACKSIZE=1G                    "    >> ${QPROC_NAME}_${IMEM}.pbs
	                        echo "${QSCRIPTCMD}"                                  >> ${QPROC_NAME}_${IMEM}.pbs
				echo "if   [[ -z \${ERROR} ]] || [[ \${ERROR} -eq 0 ]] ; then " >> ${QPROC_NAME}_${IMEM}.pbs
	                        echo "touch $PROCSDIR/${QPROC_NAME}_${IMEM}_ENDOK  " >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "elif [[ -z \${ERROR} ]] || [[ \${ERROR} -gt 0 ]] ; then " >> ${QPROC_NAME}_${IMEM}.pbs
                                echo "touch $PROCSDIR/${QPROC_NAME}_${IMEM}_ERROR  " >> ${QPROC_NAME}_${IMEM}.pbs
				echo "fi 				            " >> ${QPROC_NAME}_${IMEM}.pbs

        # Encolamiento 
        ${QSUB_ROOT}qsub ${QPROC_NAME}_${IMEM}.pbs   

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

