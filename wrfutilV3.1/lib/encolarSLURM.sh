export QSUB_ROOT=/usr/bin/
LIFS=":"
### QPROC_NAME QPROC QTHREADS QARRAY QWALLTIME QDEPEND_NAME $QEXCLU $QUEUE  $PCODE  $QDEPENDTYPE
queue (){
        QDEPEND=$(getDepend)
	totproc=$(($QPROC*$QTHREADS))
	nodes=$(((${totproc}+${ICORE} -1)/${ICORE}))  ## Esto hace division entera parte superior
	ppn=$(((${QPROC}+${nodes} -1)/${nodes}))  ## Esto hace division entera parte superior, asi que si no da exacto, reserva mas cores que los necesarios
	mpiproc=${QPROC}
	threads=${QTHREADS}
	qfile=${QPROC_NAME}.${QUEUESYS}
	scrfile=${qfile}
	qAsize=$((10#$MIEMBRO_FIN-10#$MIEMBRO_INI+1))
	concProc="1"
	rm -f $qfile
	touch $qfile
	if [[ $QNEXCLU -ne 0 && ( ${totproc} -lt  ${ICORE}) && $QARRAY && $QEXCLU -eq 0 ]] ### Hay que compartir nodos exclusivos
	then 
		concProc="$(($ICORE/$totproc))"
		test  $qAsize -lt $concProc  && concProc=$qAsize
		scrfile=${QPROC_NAME}.scr
		rm -f $scrfile
		touch $scrfile
		chmod 755  $scrfile
	fi
	
###########
## Encabezado SLURM
###########
         			echo "#!/bin/bash" >>  $qfile     									## Cola a la que se lo va a encolar
        test $QUEUE && 		echo "#SBATCH -p $QUEUE " >>  $qfile     									## Cola a la que se lo va a encolar
        test $PCODE && 		echo "#SBATCH -A $PCODE "  >> $qfile									## Codigo de Usuario/Cuanta si fuera necesario
        test $QPROC_NAME &&	echo "#SBATCH -J ${QPROC_NAME} "  >> $qfile								## Nombre del Job
        			echo "#SBATCH --mail-type=NONE " >>  $qfile										## Cuando se manda mail (default = nunca)
        test $QWALLTIME && 	echo "#SBATCH -t ${QWALLTIME} "  >>  $qfile							## Tiempo estimado de corrida
        test $QDATE && 		echo "#SBATCH --begin=${QDATE} "  >>  $qfile							## Tiempo estimado de corrida
	[[ $QEXCLU != "0" ]] && echo "#SBATCH --exclusive"        >> $qfile                                                     ## PAra que el proceso corra en nodos exclusivos
                                echo "#SBATCH -o $BASEDIR/LOGS/${QPROC_NAME}.out"  >>  $qfile                                                   		## Redireccion de flujos de salida y error
                                echo "#SBATCH -e $BASEDIR/LOGS/${QPROC_NAME}.err"  >>  $qfile                                                   		## Redireccion de flujos de salida y error
#                                echo "#SBATCH --export=WRFUTILDIR"  >>  $qfile                                 			## Indica las variables que se heredan del enviroment
        			#echo "#SBATCH --export=ALL "  >>  $qfile										## Indica que debe heredar todo el enviroment
        			echo "#SBATCH -n ${mpiproc}"   >>$qfile 				## Cantidad de procesos mpi
				echo "#SBATCH -c ${threads}"   >>$qfile 				## Cantidad de threads por proceso
				echo "#SBATCH -N ${nodes}"     >>$qfile 				## Cantidad de nodos requeridos
				echo "#SBATCH --ntasks-per-node=${ppn}">>$qfile 			## Cantidad maxima de procesos por nodos
				
        test $QARRAY &&         echo "#SBATCH -a ${QARRAY}:${concProc}"  >> $qfile                                                    			## Genera un array de procesos en PBSPRO
        			echo $QDEPEND  >>  $qfile										## Indica la dependencia del Job
                                echo 'cd $SLURM_SUBMIT_DIR' >> $qfile
	if [[ "$scrfile" != "$qfile" ]]
	then
                                echo 'begin=$SLURM_ARRAY_TASK_ID' >> $qfile
                                echo 'end=$(($SLURM_ARRAY_TASK_ID +'${concProc}' -1))' >> $qfile
                                echo 'if [[ $end > $((10#'$MIEMBRO_FIN')) ]] ; then end=$((10#'$MIEMBRO_FIN')) ; fi' >> $qfile
                                echo 'seq $begin $end | xargs -n1 -P'${concProc}' -I{} $PWD/'$scrfile' {}' >> $qfile
	fi

#########
## Encabezado de Script //Generalmente dependiente del cluster que se use
#########
				echo "export WRFUTILDIR=$WRFUTILDIR"  >> $scrfile
                                echo "source $WRFUTILDIR/lib/errores.env"  >> $scrfile 
                                #echo 'test ! $WRFUTILCONF &&'" source $BASEDIR/config.env $NOMBRE" >> $scrfile 
                                echo "source $BASEDIR/config.env $NOMBRE" >> $scrfile 
                                #echo 'test ! $WRFUTILEXP &&'" source $BASEDIR/experimento.conf" >> $scrfile 
                                echo " source $BASEDIR/experimento.conf" >> $scrfile 
                                #echo 'test ! $WRFUTILENV &&'" source $BASEDIR/entorno.sh" >> $scrfile

        if [[ $QUEUE == "share" ]]
        then
                		echo "export MPIEXE="  >>$scrfile 
                		#echo "export MPIEXE=\"mpirun  \`hostname\` -np ${mpiproc}\""  >>$scrfile 
                		#echo "export MPI_USE_ARRAY=false" >>$scrfile 
        else
                		echo "export MPIEXE='/home/opt/intel/compilers_and_libraries_2019.1.144/linux/mpi/intel64/bin/mpirun -bootstrap slurm -n ${QPROC}'" >>$scrfile 
        fi
        #test $(($QARRAY)) -ne 0  &&  echo 'ARRAYID=$SLURM_ARRAY_TASK_ID'  >> $scrfile              #PBSPRO           
        test $(($QARRAY)) -ne 0  &&  echo 'ARRAYID=${1:-$SLURM_ARRAY_TASK_ID}'  >> $scrfile              #PBSPRO           
        test $(($QARRAY)) -ne 0  &&  echo 'ARRAYCNT=${1:-$SLURM_ARRAY_TASK_COUNT}'  >> $scrfile              #PBSPRO           
        echo "export OMP_NUM_THREADS=$QTHREADS"  >> $scrfile              #PBSPRO

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
                ${QSUB_ROOT}/sbatch $qfile > ${BASEDIR}/PROCS/${QPROC_NAME}
        else
                echo "Recuerde encolar de la siguiente manera:"
                echo "${QSUB_ROOT}/sbatch $qfile > ${BASEDIR}/PROCS/${QPROC_NAME}"
        fi
}


getDepend(){
        depend=$(join_by $(name2PID)) # busca el QPID de todos los proceso de los que depende (name2PID) y los concatena (join_by)
        afterok="#SBATCH -d ${QDEPENDTYPE}:" # Averigua si es un Array o no
        test $depend && echo ${afterok}${depend} # Si hay dependencia arma la linea del script
}

name2PID(){
	# Nos fijamos en todos los archivos que sean ${QDEPEND_NAME}*
        for file in $(ls  ${BASEDIR}/PROCS/${QDEPEND_NAME} 2>/dev/null)
        do
                jid=$(cat $file 2>/dev/null)
		jid=${jid##* }
                if ! [[ $jid = "FIN" ]]
                then
			#Ojo con la linea de abajo (qstat) que puede ser muy INEFICIENTE 
                        estado=$(squeue -j $jid -o "%t" |grep -v ST |grep -E "R|S|PD|CG|RD" )
                        if ! [[ -z $estado ]]
                        then
                                echo $jid
                        fi
                fi
        done

}

function join_by { local IFS=${LIFS};  echo "$*"; }

