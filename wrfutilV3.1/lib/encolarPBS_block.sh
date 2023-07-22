queue (){
        ini_mem=${1}
        end_mem=${2}
	ens_size=$(($end_mem - $ini_mem + 1 ))

        TPROC=$(( $QNODE * $QPROC ))  #Total number of cores to be used. 	
	for QMIEM in $(seq -w $ini_mem $end_mem ) ; do 

	                        echo "#PBS -j oe "                                    >  ${QPROC_NAME}_${QMIEM}.pbs                   ## Redireccion de flujos de salida y error
	test $QUEUE && 		echo "#PBS -q $QUEUE "                                >> ${QPROC_NAME}_${QMIEM}.pbs                   ## Cola a la que se lo va a encolar
        test $QWALLTIME && 	echo "#PBS -l walltime=${QWALLTIME} "                 >> ${QPROC_NAME}_${QMIEM}.pbs	              ## Tiempo estimado de corrida
                                echo "#PBS -l nodes=${QNODE}:ppn=${QPROC}"            >> ${QPROC_NAME}_${QMIEM}.pbs                   ## Recursos que seran utilizados 
        			echo '#PBS -v BASEDIR'                                >> ${QPROC_NAME}_${QMIEM}.pbs		      ## Indica que debe heredar todo el enviroment
	test $QTHREAD   &&	echo "export OMP_NUM_THREADS=$QTHREAD"                >> ${QPROC_NAME}_${QMIEM}.pbs		
				echo "source $BASEDIR/conf/config.env"                >> ${QPROC_NAME}_${QMIEM}.pbs 
				echo "source $BASEDIR/conf/$EXPCONF"                  >> ${QPROC_NAME}_${QMIEM}.pbs 
                		echo "export MPIEXE=\"$(which mpirun) -np ${TPROC}\" ">> ${QPROC_NAME}_${QMIEM}.pbs                   ## Comando MPIRUN con cantidad de nodos y cores por nodos
	                        echo "$ENVSET  "                                      >> ${QPROC_NAME}_${QMIEM}.pbs			
	                        echo "MIEM=$QMIEM "                                   >> ${QPROC_NAME}_${QMIEM}.pbs 			
	                        echo "${QSCRIPTCMD}"                                  >> ${QPROC_NAME}_${QMIEM}.pbs
				echo "if [[ -z \${res} ]] || [[ \${res} -eq "OK" ]] ; then" >> ${QPROC_NAME}_${QMIEM}.pbs
	                        echo "touch $PROCSDIR/${QPROC_NAME}_${QMIEM}_ENDOK  " >> ${QPROC_NAME}_${QMIEM}.pbs  #Si existe la variable RES en el script la usamos
				echo "fi 				            " >> ${QPROC_NAME}_${QMIEM}.pbs

        # Encolamiento 
        ${QSUB_ROOT}qsub ${QPROC_NAME}_${QMIEM}.pbs   

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

