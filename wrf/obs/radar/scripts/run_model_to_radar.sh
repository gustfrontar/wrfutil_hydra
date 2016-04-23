
#This script runs the program model to radar to process the high resolution
#NHM output.

MODELDATAPATH=/data2/jruiz/NATURE_RUN/NATURE_RUN_1KM_OPENBOUNDARY/
RADARFILE=""

#Compile the code.

./make_model_to_radar.sh

export OMP_NUM_THREADS=20

source util.sh

YY=2013
MM=07
DD=13

DT=30

ihh=06
imm=00
iss=00

fhh=06
fmm=00
fss=00

mkdir $MODELDATAPATH/OBS/

while [ $ihh$imm$iss -le  $fhh$fmm$fss ]
do

echo $YY $MM $DD $ihh $imm $iss

MODELFILE=$MODELDATAPATH/wrfout_d01_$YY-$MM-$DD_$ihh:$imm:$iss

ln -sf $MODELFILE ./input_model.nc 

#Run model_to_radar module
./model_to_radar.exe

#Copy files in radar format.
mkdir -p $MODELDATAPATH/OBS/$YY$MM$DD$ihh$imm/
if [ $iss -eq 00 ]
then
mv output_radar.grd $MODELDATAPATH/OBS/$YY$MM$DD$ihh$imm/rad01001.dat
else
mv output_radar.grd $MODELDATAPATH/OBS/$YY$MM$DD$ihh$imm/rad01002.dat
fi


date_edit $YY $MM $DD $ihh $imm $iss $DT > tmp.txt
read YY MM DD ihh imm iss  < tmp.txt
rm -f tmp.txt


done


