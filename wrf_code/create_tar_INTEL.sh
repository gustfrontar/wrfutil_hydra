#This scripts creates the tar files that will be used by the 
#LETKF-WRF assimilation and forecast cycle. 
#The tar files contains the compiled codes as well as all
#the auxiliary files required to run the model. 
BASEDIR=$(pwd)
TARDIR=/home/ra000007/a04037/data/WRF_tar/     #Where the tar files will be created.
WPS=/home/ra000007/a04037/data/WPS4.5_INTEL_NETCDF4      #Path to WPS compilation.
WRF=/home/ra000007/a04037/data/WRF4.5_INTEL_NETCDF4      #Path to WRF compilation.
WRFDA=/home/ra000007/a04037/data/WRFDA4.5_INTEL_NETCDF3  #Path to WRF-DA compilation.

cd ${WPS}
tar  -h --dereference -cvf ./$(basename $WPS).tar *
echo mv $(basename $WPS).tar $TARDIR
mv $(basename $WPS).tar $TARDIR

cd ${WRF}/run
tar  -h --dereference -cvf ./$(basename $WRF).tar *
mv $(basename $WRF).tar  $TARDIR

cd ${WRFDA}/var/da/
tar  -h --dereference -cvf ./$(basename $WRFDA).tar da_update_bc.exe
mv $(basename $WRFDA).tar  $TARDIR


