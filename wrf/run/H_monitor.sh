
my_path=$HOME/data/EXPERIMENTS/

a=0
while [ $a -eq 0 ] ; do

N_ACTIVE_PRC=0

for i in `ls $my_path ` ; do

    ilog=0
    cont=1

    while [ $cont -eq 1 ] ; do
        ilog=`expr $ilog + 1 `
        if [ ! -e $my_path/$i/log${ilog}.log   ]  ; then
          cont=0 
        fi
    done
    ilog=`expr $ilog - 1 `

    my_log=$my_path/$i/log${ilog}.log


    #echo "Checking $i "

    tmp=`grep "RUNNING IN" $my_log 2>/dev/null `

 if [ $? -eq 0 ] ; then
    my_server=`echo $tmp | cut -c27-28`
    my_pid=`echo $tmp | cut -c42-48`

    my_server="klogin"${my_server}
 
    ssh $my_server "ps -ef " > $HOME/tmp.tmp

    tmp=`grep $my_pid $HOME/tmp.tmp `  2>/dev/null

    if [ $? -eq 0 ] ; then

      echo "Found and active process @$my_server and pid $my_pid :  $i "

      N_ACTIVE_PRC=`expr $N_ACTIVE_PRC + 1 `
    
    fi
 #else
   #echo "I cannot process $my_log , no sever or pid information"

 fi

  
  

done


echo "There are $N_ACTIVE_PRC experiments running in the cluster "

sleep 60

done
