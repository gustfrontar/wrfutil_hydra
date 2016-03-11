
#This script runs the program model to radar to process the high resolution
#NHM output.

MODELDATAPATH=/home/pmaldonado/WRF_NATURE_RUNV2/
OBSPATH=$MODELDATAPATH/OBS_0.5ms_2db_wt_noatt_lessvert/
RADARFILE=""

ulimit -s unlimited
export KMP_STACKSIZE=1g #For openmp

#Compile the code.

./make_invap_lessvert.sh


export OMP_NUM_THREADS=16

source util.sh

YY=2014
MM=01
DD=22 

DT=300

ihh=15  #15
imm=00  #00
iss=00

fhh=22  #22
fmm=00  #00
fss=00


while [ $ihh$imm$iss -le  $fhh$fmm$fss ]
do

echo $YY $MM $DD $ihh $imm $iss

MODELFILE=$MODELDATAPATH/wrfout_d03_${YY}-${MM}-${DD}_${ihh}:${imm}:${iss}

echo $MODELFILE


ln -sf $MODELFILE ./input_model.nc 

#Run model_to_radar module
./model_to_radar_invap_lessvert.exe

#Copy files in radar format.
mkdir -p $OBSPATH/$YY$MM$DD$ihh$imm/
mv output_radar.grd $OBSPATH/$YY$MM$DD$ihh$imm/rad001.dat


date_edit $YY $MM $DD $ihh $imm $iss $DT > tmp.txt
read YY MM DD ihh imm iss  < tmp.txt
rm -f tmp.txt


done


