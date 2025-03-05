#!/bin/sh 

SRCDIR="/vol0301/data/hp150019/u10335/PREVENIR/result/2km_SSD/noDA"
#SRCDIR="/vol0301/data/hp150019/u10335/PREVENIR/result/2km_SR/config-0"
nmem=20
VARS_FCST='"PW","PREC","SFC_PRES","SFC_TEMP","U10","V10","T2","Q2","MSLP","topo",'

# time 
#timestart="2018-12-14 00:30:00"
#timeend="2018-12-14 02:30:00"
timestart="2019-10-11 06:30:00"
timeend="2019-10-11 10:00:00"
#timeint="1 hour"
timeint="30 minutes"
# fcst
dtype="fcst"
# z or p 
vtype="z"
# xy or ll (latlon)
htype="ll"

# Lat-Lon 
dlat=0.025
dlon=0.025
### SR
lat_south=-34.5
lat_north=-28.5
lon_west=-68.0
lon_east=-60.0
### SSD
#lat_south=-38.5
#lat_north=-32.5
#lon_west=-63.5
#lon_east=-55.5

mkdir -p log
VARS="$VARS_FCST"
mean=""

timenow=$timestart
timef=$(date -ud "$timenow" +%Y%m%d%H%M%S)
timef2=$(date -ud "$timenow" +"%Y%m%d-%H%M%S.000")
icount=0
while [ $(date -ud "$timenow" +%s) -le $(date -ud "$timeend" +%s) ] ; do
for mem in $(seq -f %04g 1 $nmem) $mean ;do
  icount=$((icount+1))
  count=$(printf %04d $icount)

#  if [ "$htype" == "ll" ] ;then
    file_in="history"
   # file_in="histsfc_merge_xy${vtype}"
    file_out="histsfc_merge_${htype}${vtype}"
#  else
#    file_in="history"
#    file_out="histsfc_merge_${htype}${vtype}"
#  fi

  dir="$SRCDIR/$timef/$dtype/$mem"
  cat conf/sno_template.conf | sed -e "s#<--vars-->#$VARS#g" -e "s#<--path_in-->#\"$dir/$file_in\"\,#g"  -e "s#<--path_out-->#\"$dir/$file_out\"\,#g"  -e "s#<--count-->#$count#g"  > conf/sno_${count}.conf 

if [ "$htype" == "ll" ] ;then
cat << EOF >>  conf/sno_${count}.conf 

&PARAM_SNOPLGIN_HGRIDOPE
 SNOPLGIN_hgridope_type="LATLON",
 SNOPLGIN_hgridope_lat_start=${lat_south},
 SNOPLGIN_hgridope_lat_end=${lat_north},
 SNOPLGIN_hgridope_dlat=${dlat},
 SNOPLGIN_hgridope_lon_start=${lon_west},
 SNOPLGIN_hgridope_lon_end=${lon_east},
 SNOPLGIN_hgridope_dlon=${dlon},
/
EOF

fi

done
timenow=$(date -ud "$timeint $timenow" +"%F %T")
timef=$(date -ud "$timenow" +%Y%m%d%H%M%S)
timef2=$(date -ud "$timenow" +"%Y%m%d-%H%M%S.000")
done

echo "total "$icount" jobs submitted."
pjsub --bulk --sparam "1-$icount" exec_sno_bulk.sh
