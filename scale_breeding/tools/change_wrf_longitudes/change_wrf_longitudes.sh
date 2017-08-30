DATAPATH=/home/jruiz/share/scale_input_data/wrf_europe_50km_2016_05_18/  #Folder containing the unmodified wrf files.

CPATH=`pwd`

execfile=$CPATH/change_wrf_longitudes.exe

cd $DATAPATH

  for f in wrfout*; do

     echo "Changing longitudes of file -> $f"
  
     ln -sf $f ./input_file.nc

     $execfile
  
  done
