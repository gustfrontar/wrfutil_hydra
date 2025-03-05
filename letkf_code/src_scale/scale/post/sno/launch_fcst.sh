#!/bin/sh 

#SRCDIR="/vol0301/data/hp150019/u10335/PREVENIR/result/2km_SSD/noDA-rev"
SRCDIR="/vol0301/data/hp150019/u10335/PREVENIR/result/2km_SSD/config-0-rev"
nmem=20
VARS_FCST='"DENS","QV","QC","QR","QI","QS","QG","QHYD","Umet","Vmet","W","T","PRES","RH","topo",'
VARS_FCST_S='"PW","PREC","SFC_PRES","SFC_TEMP","U10","V10","T2","Q2","MSLP","topo",'


# time 
#timestart="2018-12-14 00:30:00"
#timeend="2018-12-14 02:30:00"
timestart="2019-10-11 06:00:00"
timeend="2019-10-11 06:00:00"
#timestart="2019-10-11 06:30:00"
#timeend="2019-10-11 10:00:00"
timeint="30 minutes"
#timeint="1 hour"
# anal/gues/hist/fcst
dtype="fcst"
# z or p 
vtype="z"

# Lat-Lon 
dlat=0.025
dlon=0.025
### SR
#lat_south=-34.5
#lat_north=-28.5
#lon_west=-68.0
#lon_east=-60.0
### SSD
lat_south=-38.5
lat_north=-32.5
lon_west=-63.5
lon_east=-55.5

mkdir -p log
mean=""

dim="3d"
VARS="$VARS_FCST"

# xy or ll (latlon)
for htype in xy ll ;do 


timenow=$timestart
timef=$(date -ud "$timenow" +%Y%m%d%H%M%S)
timef2=$(date -ud "$timenow" +"%Y%m%d-%H%M%S.000")
icount=0
while [ $(date -ud "$timenow" +%s) -le $(date -ud "$timeend" +%s) ] ; do
for mem in $(seq -f %04g 1 $nmem) $mean ;do
  icount=$((icount+1))
  count=$(printf %04d $icount)

  if [ "$htype" == "ll" ] ;then
    file_in="history_merge_xy${vtype}"
    file_out="history_merge_${htype}${vtype}"
  else
    file_in="history"
    file_out="history_merge_${htype}${vtype}"
  fi

  dir="$SRCDIR/$timef/$dtype/$mem"
  cat conf/sno_template.conf | sed -e "s#<--vars-->#$VARS#g" -e "s#<--path_in-->#\"$dir/$file_in\"\,#g"  -e "s#<--path_out-->#\"$dir/$file_out\"\,#g"  -e "s#<--count-->#$count#g"  > conf/sno_3d_${htype}_${count}.conf 


if [ "$htype" == "ll" ] ;then
cat << EOF >>  conf/sno_3d_${htype}_${count}.conf 

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

elif [ "$vtype" == "z" ] ;then
cat << EOF >>  conf/sno_3d_${htype}_${count}.conf 
&PARAM_SNOPLGIN_VGRIDOPE
 SNOPLGIN_vgridope_type = "ZLEV",
 SNOPLGIN_vgridope_lev_num = 30,
 SNOPLGIN_vgridope_lev_data = 500., 1000.,1500.,2000.,2500.,3000.,3500.,4000.,4500.,5000.,5500.,6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,11000.,11500.,12000.,12500.,13000.,13500.,14000.,14500.,15000., 
/
EOF

fi

done
timenow=$(date -ud "$timeint $timenow" +"%F %T")
timef=$(date -ud "$timenow" +%Y%m%d%H%M%S)
timef2=$(date -ud "$timenow" +"%Y%m%d-%H%M%S.000")
done

done ### xy or ll 

timenow=$timestart
timef=$(date -ud "$timenow" +%Y%m%d%H%M%S)
timef2=$(date -ud "$timenow" +"%Y%m%d-%H%M%S.000")
icount=0
while [ $(date -ud "$timenow" +%s) -le $(date -ud "$timeend" +%s) ] ; do
for mem in $(seq -f %04g 1 $nmem) $mean ;do
  icount=$((icount+1))
  count=$(printf %04d $icount)

  dim="2d"
  VARS="$VARS_FCST_S"
  file_in="history"
  file_out="histsfc_merge_${htype}${vtype}"

  dir="$SRCDIR/$timef/$dtype/$mem"
  cat conf/sno_template.conf | sed -e "s#<--vars-->#$VARS#g" -e "s#<--path_in-->#\"$dir/$file_in\"\,#g"  -e "s#<--path_out-->#\"$dir/$file_out\"\,#g"  -e "s#<--count-->#$count#g"  > conf/sno_${dim}_ll_${count}.conf 

cat << EOF >>  conf/sno_${dim}_ll_${count}.conf 
 
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

done
timenow=$(date -ud "$timeint $timenow" +"%F %T")
timef=$(date -ud "$timenow" +%Y%m%d%H%M%S)
timef2=$(date -ud "$timenow" +"%Y%m%d-%H%M%S.000")
done

echo "total "$icount" jobs submitted."
pjsub --bulk --sparam "1-$icount" exec_sno_bulk_fcst.sh
