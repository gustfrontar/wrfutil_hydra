spawn () {
	export I_MPI_SPAWN=1
	export FI_MLX_NS_ENABLE=1
	ulimit -s unlimited
	ulimit -l unlimited

        if [ $# -lt 3 ]; then
                echo "Usage : spawn" 
                echo "spawn_serial serial_program base_dir ini_member end_member node core runtime_flags"
                exit 1
        fi
        COMMAND=$1       #The command to be executed
        BASE_DIR=$2      #Base work directory
        INI_MEMBER=$3    #The first ensemble member in the ensemble.
        END_MEMBER=$4    #The last ensemble member in the ensemble.
        PNODE=$5         #Number of nodes allocated to each task.
        PCORE=$6         #Number of cores allocated to each task.
	RUNTIME_FLAGS=$7 #Some programs requires additional arguments. 

        echo "This is spawn_serial function"
        echo "COMMAND=$COMMAND"
        echo "BASE_DIR=$BASE_DIR"
        echo "INI_MEMBER=$INI_MEMBER"
        echo "END_MEMBER=$END_MEMBER"
        echo "PNODE=$PNODE"
        echo "PCORE=$PCORE"
	echo "RUNTIME_FLAGS="$RUNTIME_FLAGS

        #Get the total number of cores to be used for each instance of the process
        #Note that serial processes will be executed in one core only, the rest of them will be idle.
        #This may be useful to avoid overloading the first node and to better distribute the compuational load
        #among nodes for serial processes. 
        PTCORE=$(( $PNODE * $PCORE ))           #Total number of cores to be used by each ensemble member.
        TCORES=$(( $INODE * $ICORE -$PCORE ))   #Total number of cores available to run the ensemble. (from machine.conf)
        MAX_SIM_MEM=$(( $TCORES / $PTCORE ))    #Floor rounding (bash default) Maximumu number of simultaneous members
	                                        #Note: when spawn is run, we need cores for the spawned processes and
						#at least 1 core for the spawn itself. In this script we allocate 
						#PCORE processes for the spawn itself and the rest for the spawned processes.
						#This is why we substract $PCORE in the computation of TCORES
						#Also if MAX_SIM_MEM == 1 we can not spawn it because it means there is room
						#only for the process and not fo spawner itself. Then we run directly the proceess
						#without using spawn. 
	if [ $MAX_SIM_MEM -eq 0 ] ; then 
	   MAX_SIM_MEM=1
	fi
        echo "INODE=$INODE"
        echo "ICORE=$ICORE"
        echo "MAX_SIM_MEM=$MAX_SIM_MEM"
        if [ $MAX_SIM_MEM -gt $END_MEMBER ] ; then
           MAX_SIM_MEM=$END_MEMBER
        fi

        #Run the command
        #This is a serial program so we will need to use the mpi_serial_wrapper.exe 
        ini_mem=$(( $INI_MEMBER ))
        end_mem=$(( $INI_MEMBER + $MAX_SIM_MEM - 1 ))
        exe_group=1
        #Run several instances of metgrid.exe using the spawner.
        #Even for serial instances we may want to launch them with several cores to prevent to many instances
        #running at the same time in the first node and to better distribute them among the available nodes
        while [ $ini_mem -le $END_MEMBER ] ; do
            echo "Executing group number " $exe_group "Ini member = "$ini_mem " End member = "$end_mem
	    if [ $MAX_SIM_MEM -gt 1 ] ; then #Using the spawner we have 2 or more instances to launch.
               echo $MPIEXEC -np 1 ./spawn.exe ./$COMMAND $PTCORE $BASE_DIR $ini_mem $end_mem $RUNTIME_FLAGS
               $MPIEXEC -np 1 ./spawn.exe ./$COMMAND $PTCORE $BASE_DIR $ini_mem $end_mem $RUNTIME_FLAGS
            else #There is only one instance to run. Then call mpi directly without using the spawner.
	       cd $BASE_DIR/$ini_mem
               echo $MPIEXEC -np $PCORE ./$COMMAND 
	    fi
            if [ $? -ne 0 ] ; then
               dispararError 9 $COMMAND
            fi
            #Set ini_mem and end_mem for the next round.
            ini_mem=$(( $ini_mem + $MAX_SIM_MEM ))
            end_mem=$(( $end_mem + $MAX_SIM_MEM ))
            if [ $end_mem -gt $END_MEMBER ] ; then
               end_mem=$END_MEMBER
            fi
            exe_group=$(( $exe_group + 1 )) #This is just to count the number of cycles performed. 
        done
}


