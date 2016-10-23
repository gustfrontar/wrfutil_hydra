DATAPATH=/...../

cd $DATAPATH

for f in wrfout*; do
  echo "Procesing file -> $f"
  
  ln -sf $f ./input_file.nc
  ./round_second.exe 

  wrfout_d03_2014-01-22_16:50:00
  new_file_name=${f:0:28}00

  mv $f $new_file_name
  
done
