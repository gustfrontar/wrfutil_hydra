
#This script runs the program model to radar to process the high resolution
#NHM output.

MODELDATAPATH=/data2/jruiz/NATURE_RUN/NATURE_RUN_1KM_SPECIFIED/BUBBLE_DATA/
RADARFILE=""

#Compile the code.

./make_test.sh

export OMP_NUM_THREADS=5

source util.sh

iYY=2013
iMM=07
iDD=13
ihh=06
imm=00
iss=00

fYY=2013
fMM=07
fDD=13
fhh=06
fmm=00
fss=00

DT=30

mkdir /data2/jruiz/TEST/

while [ $ihh$imm$iss -le  $fhh$fmm$fss ]
do

echo $iYY $iMM $iDD $ihh $imm $iss

MODELFILE=$MODELDATAPATH/wrfout_d01_$iYY-$iMM-${iDD}_$ihh:$imm:$iss

echo $MODELFILE

ln -sf $MODELFILE ./input_model.nc 

#Run model_to_radar module
./test.exe

#Copy files in radar format.
mv output_radar.grd /data2/jruiz/TEST/radar_test.grd


date_edit $iYY $iMM $iDD $ihh $imm $iss $DT > tmp.txt
read iYY iMM iDD ihh imm iss  < tmp.txt
rm -f tmp.txt


done


