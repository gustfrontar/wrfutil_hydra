for mem in $(seq -f "%03g" 0 19)
do
  echo $mem
  cp namelist.input.asim.$mem namelist.input.asim.${mem}.ori
  cp namelist.input.ens.$mem namelist.input.asim.${mem}

  sed -i 's/ use_theta_m                         = 1/ use_theta_m                         = 0/' namelist.input.asim.${mem}
done
