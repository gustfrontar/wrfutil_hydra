export QSUB_ROOT=/opt/pbs/default/bin/

### QPROC_NAME QPROC QTHREADS QARRAY QWALLTIME QDEPEND_NAME $QEXCLU $QUEUE  $PCODE 
queue (){
        QDEPEND=$(getDepend)
	nodes=$(((${QPROC}+${ICORE} -1)/${ICORE}))  ## Esto hace division entera parte superior
	cores=$(((${QPROC}+${nodes} -1)/${nodes}))  ## Esto hace division entera parte superior, asi que si no da exacto, reserva mas cores que los necesarios
	mpiproc=${QPROC}
	threads=${QTHREADS}
	qfile=${QPROC_NAME}.pbs
	scrfile=${qfile}
	qAsize=$((10#$MIEMBRO_FIN-10#$MIEMBRO_INI+1))
	concProc="1"
	rm -f $qfile
	touch $qfile
	if [[ $QEXCLU -ne 0 && ( ${QPROC} -lt  ${ICORE}) && $QARRAY ]] ### Hay que compartir nodos exclusivos
	then 
		concProc="$(($ICORE/$QPROC))"
		test  $qAsize -lt $concProc  && concProc=$qAsize
		scrfile=${QPROC_NAME}.scr
		rm -f $scrfile
		touch $scrfile
		chmod 755  $scrfile
	fi
	
###########
## Encabezado PBS
###########
         			echo "#!/bin/bash" >>  $qfile     									## Cola a la que se lo va a encolar
        test $QUEUE && 		echo "#PBS -q $QUEUE " >>  $qfile     									## Cola a la que se lo va a encolar
        test $PCODE && 		echo "#PBS -A $PCODE "  >> $qfile									## Codigo de Usuario/Cuanta si fuera necesario
        test $QPROC_NAME &&	echo "#PBS -N ${QPROC_NAME} "  >> $qfile								## Nombre del Job
        			echo "#PBS -m n " >>  $qfile										## Cuando se manda mail (default = nunca)
        test $QWALLTIME && 	echo "#PBS -l walltime=${QWALLTIME} "  >>  $qfile							## Tiempo estimado de corrida
        test $QDATE && 		echo "#PBS -a ${QDATE} "  >>  $qfile							## Tiempo estimado de corrida
        			echo "#PBS -j oe "  >>  $qfile										## Redireccion de flujos de salida y error
                                echo "#PBS -o $BASEDIR/LOGS/"  >>  $qfile                                                   		## Redireccion de flujos de salida y error
                                echo "#PBS -v WRFUTILDIR"  >>  $qfile                                 			## Indica las variables que se heredan del enviroment
        			#echo "#PBS -V "  >>  $qfile										## Indica que debe heredar todo el enviroment
        			echo "#PBS -l select=${nodes}:ncpus=${cores}:mpiprocs=${mpiproc}:ompthreads=${threads}">>$qfile 	## Recursos que seran utilizados
        test $QARRAY &&         echo "#PBS -J ${QARRAY}:${concProc}"  >> $qfile                                                    			## Genera un array de procesos en PBSPRO
        			echo $QDEPEND  >>  $qfile										## Indica la dependencia del Job
                                echo 'cd $PBS_O_WORKDIR' >> $qfile
	if [[ "$scrfile" != "$qfile" ]]
	then
                                echo 'begin=$PBS_ARRAY_INDEX' >> $qfile
                                echo 'end=$(($PBS_ARRAY_INDEX +'${concProc}' -1))' >> $qfile
                                echo 'if [[ $end > $((10#'$MIEMBRO_FIN')) ]] ; then end=$((10#'$MIEMBRO_FIN')) ; fi' >> $qfile
                                echo 'seq $begin $end | xargs -n1 -P'${concProc}' -I{} $PWD/'$scrfile' {}' >> $qfile
	fi

#########
## Encabezado de Script //Generalmente dependiente del cluster que se use
#########
                                echo 'test ! $WRFUTILCONF &&'" source $BASEDIR/config.env $NOMBRE" >> $scrfile 
                                echo 'test ! $WRFUTILEXP &&'" source $BASEDIR/experimento.conf" >> $scrfile 
                                echo 'test ! $WRFUTILENV &&'" source $BASEDIR/entorno.sh" >> $scrfile
				echo 'export TMPDIR=/glade/scratch/$USER/temp'>> $scrfile 
				echo 'mkdir -p $TMPDIR' >> $scrfile 

#       echo "MPI_SHEPHERD=true" >>  $scrfile							
        if [[ $QUEUE == "share" ]]
        then
                		echo "export MPIEXE="  >>$scrfile 
                		#echo "export MPIEXE=\"mpirun  \`hostname\` -np ${mpiproc}\""  >>$scrfile 
                		#echo "export MPI_USE_ARRAY=false" >>$scrfile 
        else
                		echo "export MPIEXE='/opt/sgi/mpt/mpt-2.15/bin/mpiexec_mpt omplace'" >>$scrfile 
        fi
        #test $(($QARRAY)) -ne 0  &&  echo 'ARRAYID=$PBS_ARRAY_INDEX'  >> $scrfile              #PBSPRO           
        test $(($QARRAY)) -ne 0  &&  echo 'ARRAYID=${1:-$PBS_ARRAY_INDEX}'  >> $scrfile              #PBSPRO           
##########
## Cuerpo de script 
##########
        echo "${QSCRIPTCMD}" >>$scrfile 
        echo 'echo "FIN" >'" ${BASEDIR}/PROCS/${QPROC_NAME}" >> $scrfile 
	chmod 755 $scrfile
	chmod 755 $qfile

##########
## Encolamiento  condicional
##########
        if [[ $EJECUTAR -ne 0 ]]
        then
                ${QSUB_ROOT}/qsub $qfile > ${BASEDIR}/PROCS/${QPROC_NAME}
        else
                echo "Recuerde encolar de la siguiente manera:"
                echo "${QSUB_ROOT}/qsub $qfile > ${BASEDIR}/PROCS/${QPROC_NAME}"
        fi
}

#queue (){
#        QDEPEND=$(getDepend)
#        test $QUEUE && 		echo "#PBS -q $QUEUE " >  ${QPROC_NAME}.pbs      						## Cola a la que se lo va a encolar
#        test $PCODE && 		echo "#PBS -A $PCODE "  >>  ${QPROC_NAME}.pbs							## Codigo de Usuario/Cuanta si fuera necesario
#        test $QPROC_NAME &&	echo "#PBS -N ${QPROC_NAME} "  >>  ${QPROC_NAME}.pbs						## Nombre del Job
#        			echo "#PBS -m n " >>  ${QPROC_NAME}.pbs								## Cuando se manda mail (default = nunca)
#        test $QWALLTIME && 	echo "#PBS -l walltime=${QWALLTIME} "  >>  ${QPROC_NAME}.pbs					## Tiempo estimado de corrida
#        			echo "#PBS -l select=${QNODES}:ncpus=${QPROC}:mpiprocs=$((${QNODES}*${QPROC}/${QTHREADS})):ompthreads=${QTHREADS} "  >>  ${QPROC_NAME}.pbs	## Recursos que seran utilizados
#        			echo "#PBS -j oe "  >>  ${QPROC_NAME}.pbs							## Redireccion de flujos de salida y error
#                                echo "#PBS -o $BASEDIR/LOGS/ "  >>  ${QPROC_NAME}.pbs                                                   ## Redireccion de flujos de salida y error
#                                echo "#PBS -v BASEDIR,NOMBRE,WRFUTILDIR"  >>  ${QPROC_NAME}.pbs                                 ## Indica las variables que se heredan del enviroment
#        			#echo "#PBS -V "  >>  ${QPROC_NAME}.pbs								## Indica que debe heredar todo el enviroment
#        test $QARRAY &&         echo "#PBS -J $QARRAY"  >> ${QPROC_NAME}.pbs                                                    ## Genera un array de procesos en PBSPRO
#        			echo $QDEPEND  >>  ${QPROC_NAME}.pbs								## Indica la dependencia del Job
#                                echo 'test ! $WRFUTILCONF &&'" source $BASEDIR/config.env $NOMBRE" >>  ${QPROC_NAME}.pbs
#                                echo 'test ! $WRFUTILEXP &&'" source $BASEDIR/experimento.conf" >>  ${QPROC_NAME}.pbs
#                                echo 'test ! $WRFUTILENV &&'" source $BASEDIR/entorno.sh" >>  ${QPROC_NAME}.pbs
#				echo 'export TMPDIR=/glade/scratch/$USER/temp'>>  ${QPROC_NAME}.pbs
#				echo 'mkdir -p $TMPDIR' >>  ${QPROC_NAME}.pbs
#
##       echo "MPI_SHEPHERD=true" >>  ${QPROC_NAME}.pbs							
#        if [[ $QUEUE == "share" ]]
#        then
#                		echo "export MPIEXE=\"mpirun  \`hostname\` -np $((${QNODES}*${QPROC}/${QTHREADS}))\""  >> ${QPROC_NAME}.pbs
#                		echo "export MPI_USE_ARRAY=false" >> ${QPROC_NAME}.pbs
#        else
##                		echo "export MPIEXE='/opt/sgi/mpt/mpt-2.15/bin/mpiexec_mpt omplace'" >> ${QPROC_NAME}.pbs
#        fi
#        test $(($QARRAY)) -ne 0  &&  echo 'ARRAYID=$PBS_ARRAY_INDEX'  >> ${QPROC_NAME}.pbs              #PBSPRO           
#####
#        echo "${QSCRIPTCMD}" >> ${QPROC_NAME}.pbs
#        echo 'echo "FIN" >'" ${BASEDIR}/PROCS/${QPROC_NAME}" >> ${QPROC_NAME}.pbs
#        if [[ $EJECUTAR -ne 0 ]]
#        then
##                ${QSUB_ROOT}/qsub ${QPROC_NAME}.pbs > ${BASEDIR}/PROCS/${QPROC_NAME}
#        else
#                echo "Recuerde encolar de la siguiente manera:"
#                echo "${QSUB_ROOT}/qsub ${QPROC_NAME}.pbs > ${BASEDIR}/PROCS/${QPROC_NAME}"
#        fi
#
#
#
#}

getDepend(){
        depend=$(join_by $(name2PID)) # busca el QPID de todos los proceso de los que depende (name2PID) y los concatena (join_by)
        #afterok="#PBS -W depend=afterok$(test -z ${depend##*[]*}  && echo 'array'):" # Averigua si es un Array o no
        afterok="#PBS -W depend=afterok:" # Averigua si es un Array o no
        test $depend && echo ${afterok}${depend} # Si hay dependencia arma la linea del script
}

name2PID(){
	# Nos fijamos en todos los archivos que sean ${QDEPEND_NAME}*
        for file in $(ls  ${BASEDIR}/PROCS/${QDEPEND_NAME} 2>/dev/null)
        do
                jid=$(cat $file 2>/dev/null)
                if ! [[ $jid = "FIN" ]]
                then
			#Ojo con la linea de abajo (qstat) que puede ser muy INEFICIENTE 
                        estado=$(qstat -fx $jid 2> /dev/null |grep job_state |grep -E "R|H|Q|B" )
                        if ! [[ -z $estado ]]
                        then
                                echo $jid
                        fi
                fi
        done

}

function join_by { local IFS=${LIFS};  echo "$*"; }

