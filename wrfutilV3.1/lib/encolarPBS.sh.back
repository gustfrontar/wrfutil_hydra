export QSUB_ROOT=/opt/pbs/default/bin/
export AFTEROK="#PBS -W depend=afterok:"
export AFTEROKA="#PBS -W depend=afterokarray:"

queue (){
        QDEPEND=$(getDepend)
        test $QUEUE && 		echo "#PBS -q $QUEUE " >  ${QPROC_NAME}.pbs      						## Cola a la que se lo va a encolar
        test $PCODE && 		echo "#PBS -A $PCODE "  >>  ${QPROC_NAME}.pbs							## Codigo de Usuario/Cuanta si fuera necesario
        test $QPROC_NAME &&	echo "#PBS -N ${QPROC_NAME} "  >>  ${QPROC_NAME}.pbs						## Nombre del Job
        			echo "#PBS -m n " >>  ${QPROC_NAME}.pbs								## Cuando se manda mail (default = nunca)
        test $QWALLTIME && 	echo "#PBS -l walltime=${QWALLTIME} "  >>  ${QPROC_NAME}.pbs					## Tiempo estimado de corrida
        			echo "#PBS -l select=ncpus=${QCORE} "  >>  ${QPROC_NAME}.pbs	## Recursos que seran utilizados
        			#echo "#PBS -l select=${QNODE}:ncpus=${QCORE}:mpiprocs=${QPROC}   "  >>  ${QPROC_NAME}.pbs	## Recursos que seran utilizados
        			echo "#PBS -j oe "  >>  ${QPROC_NAME}.pbs							## Redireccion de flujos de salida y error
        			#echo "#PBS -V "  >>  ${QPROC_NAME}.pbs								## Indica que debe heredar todo el enviroment
        			echo '#PBS -v BASEDIR'  >>  ${QPROC_NAME}.pbs								## Indica que debe heredar todo el enviroment
        			echo $QDEPEND  >>  ${QPROC_NAME}.pbs								## Indica la dependencia del Job
				echo "export WRFUTILDIR=$WRFUTILDIR" >>${QPROC_NAME}.pbs 
				echo "source $WRFUTILDIR/config.env" >>${QPROC_NAME}.pbs 
				echo "source $BASEDIR/experimento.conf" >>${QPROC_NAME}.pbs 
                		echo "export MPIEXE=$(which mpirun) -np ${QPROC}"  >> ${QPROC_NAME}.pbs                         ## Comando MPIRUN con cantidad de nodos y cores por nodos
                		#echo 'export ARRAYID=$PBS_ARRAY_INDEX'  >> ${QPROC_NAME}.pbs                         ## Comando MPIRUN con cantidad de nodos y cores por nodos
                		echo 'export ARRAYID=$PBS_ARRAYID'  >> ${QPROC_NAME}.pbs                         ## Comando MPIRUN con cantidad de nodos y cores por nodos
	#QSCRIPTCMD=${QSCRIPTCMD//__JOBID__/"$PBS_ARRAYID"}
	echo "${QSCRIPTCMD}" >> ${QPROC_NAME}.pbs
#       echo ${QSUB_ROOT}/qsub  ${QPROC_NAME}.pbs
#        PID=$(${QSUB_ROOT}/qsub ${QPROC_NAME}.pbs)
        #[[ $PID =~ (^[0-9]*) ]] && echo ${BASH_REMATCH[1]} > ${BASEDIR}/PROCS/${QPROC_NAME}
        echo $PID > ${BASEDIR}/PROCS/${QPROC_NAME}

}

getDepend(){
        depend=$(join_by $(name2PID))
	test $depend && echo ${QAFTEROK}${depend}
        if ! [[ -z $depend ]]
        then
                test echo ${AFTEROK}${depend}
        else
                echo "#PBS -h "
        fi
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

