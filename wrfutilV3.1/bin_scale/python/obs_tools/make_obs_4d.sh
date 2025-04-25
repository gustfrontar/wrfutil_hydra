#!/bin/bash

BASEDIR=$(pwd)/../../..
source $BASEDIR/conf/config.env
source $BASEDIR/conf/exp.conf
source $BASEDIR/conf/machine.conf
source $BASEDIR/conf/model.conf
source $BASEDIR/conf/letkf.conf
source $BASEDIR/conf/obs.conf


make_obs_4d(){
  filehead=$1
  filetail=$2
  outfile=$3
  radar=$4
#  timets=$(date -u -d "$DA_INI_DATE UTC +$(( (ANALYSIS_FREQ*(STEP-1))+SPIN_UP_LENGTH+ANALYSIS_FREQ )) seconds" +%s)
  timets=$(date -u -d "$ANALYSIS_DATE" +%s)
  winis=$((ANALYSIS_WIN_INI - ANALYSIS_FREQ ))
  wends=$((ANALYSIS_WIN_END - ANALYSIS_FREQ ))

  files="$( /bin/ls -x ${filehead}*${filetail})"
  rm -rf $outfile
  for f in $files ;do
    timef=$(echo $f | sed  -e "s|$filehead||g" |sed -e "s|$filetail||g")  
    timefs=$(date -ud "${timef:0:4}-${timef:4:2}-${timef:6:2} ${timef:8:2}:${timef:10:2}:${timef:12:2}" +%s)
    dif=$((timefs-timets))
    if ((dif > winis)) && ((dif <= wends)) ;then 
      echo ./exe_make_obs_4d $f $outfile $dif $radar
      ./exe_make_obs_4d $f $outfile $dif $radar
#      echo python ./make_obs_4d.py $f $outfile $dif $radar
#      python ./make_obs_4d.py $f $outfile $dif $radar
    fi
  done
}

TMPPATH=$(cd $OBSPATH && pwd | sed  -e "s|_4d||g")
if [ -z  "$TMPPATH" ] ;then
  echo "ERROR : check OBSPATH" $TMPPATH
  exit 1
fi

OBS4DPATH="${TMPPATH}_4d"
mkdir -p $OBS4DPATH

OBS_IN_NUM=${#OBS_LIST[@]}

STEP=1
ANALYSIS_DATE=$(date -u -d "$DA_INI_DATE UTC +$(( SPIN_UP_LENGTH+ANALYSIS_FREQ )) seconds" +"%F %T")
while [ $(date -ud "$ANALYSIS_DATE" +%s) -le $(date -ud  "$DA_END_DATE" +%s) ] ;do 

ANALYSIS_DATE_PFMT=$(date -u -d "$ANALYSIS_DATE" +"%Y%m%d%H%M%S")  #Path/folder format

for j in $(seq $OBS_IN_NUM) ;do
  i=$((j-1))
  OBS_TYPE=${OBS_LIST[$i]}
  mkdir -p $OBS4DPATH/${OBS_TYPE}

  if [  $OBS_TYPE == "RADARC" ] ;then
    for RADARC in ${RADARC_LIST[@]} ;do
      make_obs_4d "$TMPPATH/${OBS_TYPE}/WRF_${OBS_TYPE}_${RADARC}_" ".dat" $OBS4DPATH/${OBS_TYPE}/WRF_${OBS_TYPE}_${RADARC}_${ANALYSIS_DATE_PFMT}.dat 1

    done
  else 
    make_obs_4d "$TMPPATH/${OBS_TYPE}/WRF_${OBS_TYPE}_" ".dat" $OBS4DPATH/${OBS_TYPE}/WRF_${OBS_TYPE}_${ANALYSIS_DATE_PFMT}.dat 0
  fi
 
done
STEP=$((STEP+1))
ANALYSIS_DATE=$(date -u -d "$ANALYSIS_DATE UTC + $ANALYSIS_FREQ seconds" +"%F %T")

done
