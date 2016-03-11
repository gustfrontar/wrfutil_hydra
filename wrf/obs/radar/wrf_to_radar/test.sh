
#This script runs the program model to radar to process the high resolution
#NHM output.

MODELDATAPATH=/data1/jruiz/NATURE_RUN/NATURE_RUN_1KM_LESSCAPE/BUBBLE_DATA/
RADARFILE=""

#Compile the code.

./make_idealradar.sh

export OMP_NUM_THREADS=10

source util.sh

iYY=2013
iMM=07
iDD=13
ihh=06
imm=30
iss=00

fYY=2013
fMM=07
fDD=13
fhh=06
fmm=30
fss=00

DT=60

#mkdir $MODELDATAPATH/OBS_IDEAL/

while [ $ihh$imm$iss -le  $fhh$fmm$fss ]
do

echo $iYY $iMM $iDD $ihh $imm $iss

MODELFILE=$MODELDATAPATH/wrfout_d01_$iYY-$iMM-${iDD}_$ihh:$imm:$iss

echo $MODELFILE

ln -sf $MODELFILE ./input_model.nc 

#Run model_to_radar module
./model_to_radar_ideal.exe

#Copy files in radar format.
#mkdir -p $MODELDATAPATH/OBS_IDEAL/$iYY$iMM$iDD$ihh$imm$iss/
#if [ $iss -eq 00 ]
#then
#mv radarobs.dat $MODELDATAPATH/OBS_IDEAL/$iYY$iMM$iDD$ihh$imm$iss/rad0101.dat
#else
#mv radarobs.dat $MODELDATAPATH/OBS_IDEAL/$iYY$iMM$iDD$ihh$imm$iss/rad0201.dat
#fi


date_edit $iYY $iMM $iDD $ihh $imm $iss $DT > tmp.txt
read iYY iMM iDD ihh imm iss  < tmp.txt
rm -f tmp.txt


done


