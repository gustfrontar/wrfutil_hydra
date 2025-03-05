#! /bin/bash -x

### Energy & Mass balance ###

echo "+visualize by gnuplot"
rm -f energy.dat mass.dat mass_q.dat

while read -a line
do
   if [ ${line[0]} == "STEP=" ]; then
      echo ${line[1]} ${line[10]} ${line[11]} ${line[12]} ${line[13]} >> energy.dat
      echo ${line[1]} ${line[6]}  ${line[7]} ${line[8]} ${line[9]} >> mass.dat
      echo ${line[1]} ${line[3]}  ${line[4]}  ${line[5]}  ${line[7]} >> mass_q.dat
   fi
done < monitor.peall

gnuplot < ./visualize/energy.plt || exit
gnuplot < ./visualize/mass.plt   || exit
gnuplot < ./visualize/mass_q.plt || exit
rm -f energy.dat mass.dat mass_q.dat



### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

get_range_opt(){
   local range
   if [ $1 == "auto" ]; then
      range=""
   else
      range="--range "$1
   fi
   echo $range
}
get_eddy_opt(){
   local range
   if [ $1 == "off" ]; then
      range=""
   else
      range="--eddy "$1
   fi
   echo $range
}

var__set=(PT W QHYD)
rangeset_slice=(-3:3 auto 0:1.4e-4)
rangeset_column=(auto auto auto)
eddyset=(time off off)
rangeset_snapshot=(auto auto auto)
eddyset_snapshot=(x off off)

time_set=(00000 00900 01800 02700 03600 05400)

i=0
for var in ${var__set[@]}
do
   # time series
   range=$(get_range_opt ${rangeset_slice[${i}]})
   eddy=$(get_eddy_opt ${eddyset[${i}]})
   gpview history.pe\*.nc@${var},y=6000:26000,x=6000:26000,z=0:15000 --nocont --mean x,y ${eddy} --exch ${range} --wsn 2 -sw:ifl=1 || exit
   mv dcl_0001.png slice_${var}.png

   range=$(get_range_opt ${rangeset_column[${i}]})
   eddy=$(get_eddy_opt ${eddyset[${i}]})
   gpview history.pe\*.nc@${var},y=16000,z=0:15000 --nocont --mean z ${eddy} --exch  ${range} --wsn 2 -sw:ifl=1 || exit
   mv dcl_0001.png column_${var}.png

   # snapshot
   for sec in ${time_set[@]}
   do
       range=$(get_range_opt ${rangeset_snapshot[${i}]})
       eddy=$(get_eddy_opt ${eddyset_snapshot[${i}]})
       gpview history.pe\*.nc@${var},y=16000,z=0:15000,time=${sec} ${eddy} ${range} --nocont --wsn 2 -sw:ifl=1 || exit
       mv dcl_0001.png ${var}${sec}sec.png
   done

   let i="${i} + 1"
done

var__set=(RH QC QR)
rangeset_slice=(0:80 0:1.2e-4 0:1.2e-4)
rangeset_column=(auto auto auto)
eddyset_column=(time time time)

i=0
for var in ${var__set[@]}
do
   # time series
   range=$(get_range_opt ${rangeset_slice[${i}]})
   gpview history.pe\*.nc@${var},y=6000:26000,x=6000:26000,z=0:15000 --nocont --mean x,y --exch ${range} --wsn 2 -sw:ifl=1 || exit
   mv dcl_0001.png slice_${var}.png

   range=$(get_range_opt ${rangeset_column[${i}]})
   eddy=$(get_eddy_opt ${eddyset_column[${i}]})   
   gpview history.pe\*.nc@${var},y=16000,z=0:15000 --nocont --mean z ${eddy} --exch ${range} --wsn 2 -sw:ifl=1 || exit
   mv dcl_0001.png column_${var}.png

   let i="${i} + 1"
done
