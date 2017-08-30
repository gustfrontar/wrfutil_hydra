DATAPATH=/tmp  #Carpeta donde estan los wrfout que tienen el problema de los segundos.

execfile=$HOME/round_seconds/round_seconds.exe

cd $DATAPATH

for f in wrfout*; do
  echo "Procesing file -> $f"
 
  new_file_name=${f:0:28}00
  mv $f $new_file_name
  
  ln -sf $new_file_name ./input_file.nc
  $execfile
  
done
