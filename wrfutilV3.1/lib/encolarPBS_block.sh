queue (){
	                        echo "#PBS -j oe "                                    >  ${QPROC_NAME}_${QMIEM}.pbs                   ## Redireccion de flujos de salida y error
	test $QUEUE && 		echo "#PBS -q $QUEUE "                                >> ${QPROC_NAME}_${QMIEM}.pbs                   ## Cola a la que se lo va a encolar
        test $QWALLTIME && 	echo "#PBS -l walltime=${QWALLTIME} "                 >> ${QPROC_NAME}_${QMIEM}.pbs	              ## Tiempo estimado de corrida

	                        if [ ! -z $QEXCLU ] & [ $QEXCLU -eq 1 ] ; then             
				   echo "#PBS -l nodes=${QNODE}:ppn=${ICORE}"         >> ${QPROC_NAME}_${QMIEM}.pbs                   ## Recursos que seran utilizados (usamos el nodo completo)
                                else 
                                   echo "#PBS -l nodes=1:ppn=${QCORE}"                >> ${QPROC_NAME}_${QMIEM}.pbs                   ## Recursos que seran utilizados (usamos menos de 1 nodo)
				fi
        			echo '#PBS -v BASEDIR'                                >> ${QPROC_NAME}_${QMIEM}.pbs		      ## Indica que debe heredar todo el enviroment
				echo "source $BASEDIR/conf/config.env"                >> ${QPROC_NAME}_${QMIEM}.pbs 
				echo "source $BASEDIR/conf/experimento.conf"          >> ${QPROC_NAME}_${QMIEM}.pbs 
                		echo "export MPIEXE=\"$(which mpirun) -np ${QPROC}\" ">> ${QPROC_NAME}_${QMIEM}.pbs                   ## Comando MPIRUN con cantidad de nodos y cores por nodos
                		echo 'export ARRAYID=$PBS_ARRAYID'                    >> ${QPROC_NAME}_${QMIEM}.pbs      
	                        echo "$ENVSET  "                                      >> ${QPROC_NAME}_${QMIEM}.pbs			
	                        echo "MIEM=$QMIEM "                                   >> ${QPROC_NAME}_${QMIEM}.pbs 			
	                        echo "${QSCRIPTCMD}"                                  >> ${QPROC_NAME}_${QMIEM}.pbs
				echo "if [[ -z \"${RES}\" ]] || [[ \"${RES}\" -eq "OK" ]] ; then" >> ${QPROC_NAME}_${QMIEM}.pbs
	                        echo "touch $PROCSDIR/${QPROC_NAME}_${QMIEM}_ENDOK  " >> ${QPROC_NAME}_${QMIEM}.pbs  #Si existe la variable RES en el script la usamos
				echo "fi                                            " >> ${QPROC_NAME}_${QMIEM}.pbs


##########
## Encolamiento  condicional
##########
        if [[ $EJECUTAR -ne 0 ]]
        then
            ${QSUB_ROOT}qsub ${QPROC_NAME}_${QMIEM}.pbs   
        else
            echo "Recuerde encolar de la siguiente manera:"
            echo "${QSUB_ROOT}qsub ${QPROC_NAME}_${QMIEM}.pbs "
        fi


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

