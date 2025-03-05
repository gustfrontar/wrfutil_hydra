#!/bin/sh 

#SRCDIR="/vol0301/data/hp150019/u10335/PREVENIR/result/2km_SSD/config-0-rev"
SRCDIR="/vol0301/data/hp150019/u10335/PREVENIR/result/2km_SR/config-0-rev"
nmem=40
VARS_RESTART='"MOMX","MOMY","DENS","RHOT","QV","QC","QI","QR","QS","QG","topo",'
#VARS_FCST='"DENS","QV","QC","QR","QI","QS","QG","QHYD","Umet","Vmet","W","T","PRES","RH","PW","PREC","SFC_PRES","SFC_TEMP","U10","V10","T2","Q2","MSLP","topo",'
VARS_FCST='"DENS","QV","QC","QR","QI","QS","QG","QHYD","Umet","Vmet","W","T","PRES","RH","topo",'

# time 
timestart="2018-12-14 00:00:00"
timeend="2018-12-14 03:00:00"
#timestart="2019-10-11 04:00:00"
#timeend="2019-10-11 04:00:00"
#timeint="5 minutes"
timeint="1 hour"
# anal/gues/hist/fcst
dtype="anal"
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
if [ "$dtype" == "anal" ] || [ "$dtype" == "gues" ]; then
    VARS="$VARS_RESTART"
    mean="mean"
elif [ "$dtype" == "fcst" ] || [ "$dtype" == "hist" ]; then
    VARS="$VARS_FCST"
    mean=""
else
  echo "dtype "$dtype" not supported."
  exit 1
fi

timenow=$timestart
timef=$(date -ud "$timenow" +%Y%m%d%H%M%S)
timef2=$(date -ud "$timenow" +"%Y%m%d-%H%M%S.000")
icount=0
while [ $(date -ud "$timenow" +%s) -le $(date -ud "$timeend" +%s) ] ; do
for mem in $(seq -f %04g 1 $nmem) $mean ;do
  icount=$((icount+1))
  count=$(printf %04d $icount)

  if [ "$dtype" == "anal" ] || [ "$dtype" == "gues" ]; then
    if [ "$htype" == "ll" ] ;then
      file_in="init_merge_xy${vtype}_$timef2"
      file_out="init_merge_${htype}${vtype}_$timef2"
    else
      file_in="init_$timef2"
      file_out="init_merge_${htype}${vtype}_$timef2"
    fi
  elif [ "$dtype" == "fcst" ] || [ "$dtype" == "hist" ]; then
    if [ "$htype" == "ll" ] ;then
      file_in="history_merge_xy${vtype}"
      file_out="history_merge_${htype}${vtype}"
    else
      file_in="history"
      file_out="history_merge_${htype}${vtype}"
    fi
  fi

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

elif [ "$vtype" == "z" ]; then
cat << EOF >>  conf/sno_${count}.conf 
&PARAM_SNOPLGIN_VGRIDOPE
 SNOPLGIN_vgridope_type = "ZLEV",
 SNOPLGIN_vgridope_lev_num = 30,
 SNOPLGIN_vgridope_lev_data = 500., 1000.,1500.,2000.,2500.,3000.,3500.,4000.,4500.,5000.,5500.,6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,11000.,11500.,12000.,12500.,13000.,13500.,14000.,14500.,15000., 
/
EOF

elif [ "$vtype" == "p" ]; then
cat << EOF >>  conf/sno_${count}.conf 
&PARAM_SNOPLGIN_VGRIDOPE
 SNOPLGIN_vgridope_type = "PLEV",
 SNOPLGIN_vgridope_lev_num = 19,
 SNOPLGIN_vgridope_lev_data = 100000., 97500., 95000., 92500., 90000., 85000., 80000., 75000., 70000., 65000., 60000., 55000., 50000., 40000., 30000., 25000., 20000., 15000., 10000.,
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
