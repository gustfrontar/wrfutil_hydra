#! /bin/bash
PATH=$PATH:$CONDADIR
source activate wrfutil
HORA=$1
export FILEIN=$RUNDIR/wrfout_d01_$(date -d  "$FECHA_INI +$((10#$CICLO)) hours+$((10#$CICLO_MIN)) minutes + $HORA hours" +'%Y-%m-%d_%H_%M')_00
FileToWait=$RUNDIR/wrfout_d01_$(date -d  "$FECHA_INI +$((10#$CICLO)) hours+$((10#$CICLO_MIN)) minutes + $(($HORA+1)) hours" +'%Y-%m-%d_%H_%M')_00
[[ $((10#$PLAZO-$HORA)) -eq 0  ]]  && FileToWait="$RUNDIR/wrfout_termine"
echo "soy $HORA y proceso $FILEIN  pero espero a $FileToWait"
# Espero a que este lista la proxima hora
while [[ ! -f $FileToWait ]]; do sleep 1; done


	# Para el ensamble
python $WRFUTILDIR/bin/python/write_post_smn_ens_oper_sfc.py  
python $WRFUTILDIR/bin/python/write_post_smn_ens_oper_lev.py 

