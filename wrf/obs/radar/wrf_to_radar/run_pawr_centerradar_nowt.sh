
#This script runs the program model to radar to process the high resolution
#NHM output.

MODELDATAPATH=/data2/jruiz/NATURE_RUN/NATURE_RUN_1KM_SPECIFIED/BUBBLE_DATA/
RADARFILE=""

#Compile the code.

./make_pawr_centerradar_nowt.sh

export OMP_NUM_THREADS=10

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
fhh=07
fmm=00
fss=00

DT=30

mkdir $MODELDATAPATH/OBS_NOWT/

while [ $ihh$imm$iss -le  $fhh$fmm$fss ]
do

echo $iYY $iMM $iDD $ihh $imm $iss

MODELFILE=$MODELDATAPATH/wrfout_d01_$iYY-$iMM-${iDD}_$ihh:$imm:$iss

echo $MODELFILE

ln -sf $MODELFILE ./input_model.nc 

#Run model_to_radar module
./model_to_radar_nowt.exe

#Copy files in radar format.
mkdir -p $MODELDATAPATH/OBS_NOWT/$iYY$iMM$iDD$ihh$imm/
if [ $iss -eq 00 ]
then
mv output_radar.grd $MODELDATAPATH/OBS_NOWT/$iYY$iMM$iDD$ihh$imm/rad01001.dat
else
mv output_radar.grd $MODELDATAPATH/OBS_NOWT/$iYY$iMM$iDD$ihh$imm/rad01002.dat
fi


date_edit $iYY $iMM $iDD $ihh $imm $iss $DT > tmp.txt
read iYY iMM iDD ihh imm iss  < tmp.txt
rm -f tmp.txt


done


