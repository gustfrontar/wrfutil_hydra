spawn () {
        if [ $# -lt 3 ]; then
                echo "Usage : spawn_serial" 
                echo "spawn_serial serial_program base_dir ini_member end_member node core is_serial"
                exit 1
        fi
        COMMAND=$1      #The command to be executed
        BASE_DIR=$2     #Base work directory
        INI_MEMBER=$3   #The first ensemble member in the ensemble.
        END_MEMBER=$4   #The last ensemble member in the ensemble.
        PNODE=$5        #Number of nodes allocated to each task.
        PCORE=$6        #Number of cores allocated to each task.
        IS_SERIAL=$7    # 1- means serial application , 0 - mean mpi application

        echo "This is spawn_serial function"
        echo "COMMAND=$COMMAND"
        echo "BASE_DIR=$BASE_DIR"
        echo "INI_MEMBER=$INI_MEMBER"
        echo "END_MEMBER=$END_MEMBER"
        echo "PNODE=$PNODE"
        echo "PCORE=$PCORE"
        echo "IS_SERIAL=$IS_SERIAL"

        #Get the total number of cores to be used for each instance of the process
        #Note that serial processes will be executed in one core only, the rest of them will be idle.
        #This may be useful to avoid overloading the first node and to better distribute the compuational load
        #among nodes for serial processes. 
        PTCORE=$(( $PNODE * $PCORE )) #Total number of cores to be used by each ensemble member.
        TCORES=$(( $INODE * $ICORE ))   #Total number of cores available to run the ensemble. (from machine.conf)
        MAX_SIM_MEM=$(( $TCORES / $PTCORE ))  #Floor rounding (bash default) Maximumu number of simultaneous members
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
            tot_group_cores=$(( $PTCORE * ( $end_mem - $ini_mem + 1 ) ))
            if [ $IS_SERIAL -eq 1 ] ; then
               mpiexec -np $tot_group_cores ./spawn_serial.exe ./$COMMAND $PTCORE $BASE_DIR $ini_mem $end_mem
               if [ $? -ne 0 ] ; then
                  dispararError 9 $COMMAND
               fi
            else 
               mpiexec -np $tot_group_cores ./spawn.exe ./$COMMAND $PTCORE $BASE_DIR $ini_mem $end_mem
               if [ $? -ne 0 ] ; then
                  dispararError 9 $COMMAND
               fi
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


