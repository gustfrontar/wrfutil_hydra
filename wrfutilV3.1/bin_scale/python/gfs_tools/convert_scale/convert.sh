#!/bin/bash

tstart="$1"
fhour="$2"
tend=$(date -ud "$fhour hour $tstart" +"%F %H")
#tend="$2"
gribfile_prefix="$3"
outdir="$4"

#tstart='2016-08-10 06'
#tend='2016-08-10 06'
#gribfile_prefix='/data7/gylien/realtime/ncepgfs/2016081000/gfs'
#outdir='/data7/gylien/realtime/ncepgfs_grads/2016081000'

gribfile_dir="$(dirname ${gribfile_prefix})"
gribfile_base="$(basename ${gribfile_prefix})"
tint=10800
region_levs='1000 925 850 700 500 300 200 100 50 20 10'
region_latlon='280.0:310.5 -60.0:-20.0'

#----

wkdir="$( cd "$( dirname "$0" )" && pwd )"
cd ${wkdir}

mkdir -p ${outdir}


tinit=`echo $gribfile_dir | grep -o '[0-9]\{8\}/[0-9]\{2\}'`
stimefgrads=$(date -ud "${tinit:0:4}-${tinit:4:2}-${tinit:6:2} ${tinit:9:2}" '+%H:%MZ%d%b%Y')

cat ${wkdir}/ctl/land.ctl | sed "s/<STIME>/${stimefgrads}/g" > ${outdir}/land.ctl
cat ${wkdir}/ctl/sfc.ctl | sed "s/<STIME>/${stimefgrads}/g" > ${outdir}/sfc.ctl
cat ${wkdir}/ctl/atm.ctl | sed "s/<STIME>/${stimefgrads}/g" > ${outdir}/atm.ctl
cat ${wkdir}/ctl/reg.ctl | sed "s/<STIME>/${stimefgrads}/g" > ${outdir}/reg.ctl

ifh=0
time="$tstart"
while (($(date -ud "$time" '+%s') <= $(date -ud "$tend" '+%s'))); do

  cfh="f$(printf %03g $ifh)"
  timef=$(date -ud "$time" '+%Y-%m-%d %H')
  timef2=$(date -ud "$time" '+%Y%m%d%H%M%S')
  echo "[$timef]"

 LEVS_IGNORE="0.4 mb|:15 mb:|:40 mb:|:1 mb:|:2 mb:|:3 mb:|:5 mb:|:7 mb:"
#####

  gfile=${gribfile_prefix}.pgrb2.0p50.${cfh}
#  gfilea=${gribfile_prefix}.pgrb2a.0p50.${cfh}
#  gfileb=${gribfile_prefix}.pgrb2b.0p50.${cfh}
  ofile_land=${outdir}/land_${timef2}.grd
  ofile_sfc=${outdir}/sfc_${timef2}.grd
  ofile_atm=${outdir}/atm_${timef2}.grd
  ofile_reg=${outdir}/reg_${timef2}.grd

  tfile="tmp.grb2"
  rm -f ${tfile} ${ofile_land}
  #--Land data
  wgrib2 ${gfile} -match ":LAND:"  -match ":surface:"             -grib ${tfile}
  wgrib2 ${gfile} -match ":TMP:"   -match ":surface:"     -append -grib ${tfile}
#  wgrib2 ${gfile} -match ":TMP:"   -match "below ground:" -append -grib ${tfile}
  wgrib2 ${gfile} -match ":TSOIL:" -match "below ground:" -append -grib ${tfile}
  wgrib2 ${gfile} -match ":SOILW:" -match "below ground:" -append -grib ${tfile}
  wgrib2 ${tfile} -no_header -ieee ${ofile_land} -g2clib 0

  rm -f ${tfile} ${ofile_sfc}
  #--Surface data
#  wgrib2 ${gfile} -match ":PRES:"  -match ":mean sea level:"            -grib ${tfile}
#  wgrib2 ${gfile} -match ":MSLET:"  -match ":mean sea level:"            -grib ${tfile}
  wgrib2 ${gfile} -match ":PRMSL:" -match ":mean sea level:"            -grib ${tfile}
  wgrib2 ${gfile} -match ":PRES:"  -match ":surface:"           -append -grib ${tfile}
  wgrib2 ${gfile} -match ":UGRD:"  -match ":10 m above ground:" -append -grib ${tfile}
  wgrib2 ${gfile} -match ":VGRD:"  -match ":10 m above ground:" -append -grib ${tfile}
  wgrib2 ${gfile} -match ":TMP:"   -match ":2 m above ground:"  -append -grib ${tfile}
  wgrib2 ${gfile} -match ":RH:"    -match ":2 m above ground:"  -append -grib ${tfile}
  wgrib2 ${gfile} -match ":HGT:"   -match ":surface:"           -append -grib ${tfile}
  wgrib2 ${tfile} -no_header -ieee ${ofile_sfc} -g2clib 0

  rm -f ${tfile} ${ofile_atm}
  #--Upper data
  wgrib2 ${gfile} -match ":HGT:"  -match "mb" -not "mb above ground" | grep -vE "${LEVS_IGNORE}" | sort -n -r -t ":" -k5 | wgrib2 -i ${gfile} -grib ${tfile}
  wgrib2 ${gfile} -match ":UGRD:" -match "mb" -not "mb above ground" | grep -vE "${LEVS_IGNORE}" | sort -n -r -t ":" -k5 | wgrib2 -i ${gfile} -append -grib ${tfile}
  wgrib2 ${gfile} -match ":VGRD:" -match "mb" -not "mb above ground" | grep -vE "${LEVS_IGNORE}" | sort -n -r -t ":" -k5 | wgrib2 -i ${gfile} -append -grib ${tfile}
  wgrib2 ${gfile} -match ":TMP:"  -match "mb" -not "mb above ground" | grep -vE "${LEVS_IGNORE}" | sort -n -r -t ":" -k5 | wgrib2 -i ${gfile} -append -grib ${tfile}
  wgrib2 ${gfile} -match ":RH:"   -match "mb" -not "mb above ground" | grep -vE "${LEVS_IGNORE}" | sort -n -r -t ":" -k5 | wgrib2 -i ${gfile} -append -grib ${tfile}
  wgrib2 ${tfile} -no_header -ieee ${ofile_atm} -g2clib 0

  tfile2="tmp.inv"
  rm -f ${tfile} ${tfile2} ${ofile_reg}
  #--Regional commpact data
  for ilev in ${region_levs}; do
    wgrib2 ${gfile} -match ":UGRD:" -match ":$ilev mb:" -not "mb above ground" >> ${tfile2}
  done
  for ilev in ${region_levs}; do
    wgrib2 ${gfile} -match ":VGRD:" -match ":$ilev mb:" -not "mb above ground" >> ${tfile2}
  done
  for ilev in ${region_levs}; do
    wgrib2 ${gfile} -match ":TMP:"  -match ":$ilev mb:" -not "mb above ground" >> ${tfile2}
  done
  for ilev in ${region_levs}; do
    wgrib2 ${gfile} -match ":HGT:"  -match ":$ilev mb:" -not "mb above ground" >> ${tfile2}
  done
  for ilev in ${region_levs}; do
    wgrib2 ${gfile} -match ":RH:"   -match ":$ilev mb:" -not "mb above ground" >> ${tfile2}
  done
  wgrib2 ${gfile} -match ":UGRD:"  -match ":10 m above ground:" >> ${tfile2}
  wgrib2 ${gfile} -match ":VGRD:"  -match ":10 m above ground:" >> ${tfile2}
  wgrib2 ${gfile} -match ":TMP:"   -match ":2 m above ground:"  >> ${tfile2}
  wgrib2 ${gfile} -match ":MSLET:" -match ":mean sea level:"    >> ${tfile2}
  cat ${tfile2} | wgrib2 -i ${gfile} -small_grib ${region_latlon} ${tfile}
  wgrib2 ${tfile} -no_header -little_endian -ieee ${ofile_reg} -g2clib 0

  rm -f ${tfile} ${tfile2}

ifh=$((ifh+(tint/3600)))
time=$(date -ud "$tint second $time" '+%Y-%m-%d %H')
done
