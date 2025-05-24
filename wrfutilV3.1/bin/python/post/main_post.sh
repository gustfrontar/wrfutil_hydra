#!/bin/bash
#SBATCH -p mem1      
#SBATCH -n 10      
#SBATCH --mem 300G 
#SBATCH --time 120 

POSTTYPE="GUES"   #Postprocesing type, can be one of: ANAL/GUES/DAFCST

EXPPATH=$HOME/data/data_assimilation_exps/PREVENIR_LOWRESOLUTION_DA_FUGAKU_2024031900/

DATADIR="$EXPPATH/HIST/${POSTTYPE}/"      #Root path of WRFOUT files.
OUTDIR="$EXPPATH/POST/${POSTTYPE}/"       #Root path for postprocessing output.
SRCDIR="$EXPPATH/bin/python/post/src/"    #Postprocessing code dir.

NMEM=60          #Ensemble size
MAXJOBS=30       #Maximum simultaneous jobs
memfmt="%03g"    #Output format for ens members.
anafreq=300      #Analysis frequency (ANAL or GUES only)
forinifreq=10800 #Forecast initialization frequency (DAFCST only)
foroutfreq=1800  #Forecast output frequency (DAFCST only)
forlen=86400     #FOrecast length (DAFCST only)

timestart="2024-03-19 12:00:00" 
timeend="2024-03-20 06:00:00" 

#Create the output directory

mkdir -p $OUTDIR


icount=0
timef=$timestart

#Loops for postprocessing GUES and ANALYSIS
#===============================================================================
if [ $POSTTYPE == "ANAL" ] || [ $POSTTYPE == "GUES"  ] ; then
  for mem in $(seq -f $memfmt $NMEM) ;do 
    ctime=$timestart
    while [ $(date -ud "$ctime" +%s) -le $(date -ud "$timeend" +%s) ] ;do 
      timed=$(date -ud "$ctime" +"%Y%m%d%H%M%S")
      if [ -d "$DATADIR/$timed/" ] ; then

        if [ $POSTTYPE == "ANAL" ] ; then
          file=$DATADIR/$timed/anal$(printf %05g $((10#$mem)))
        elif [ $POSTTYPE == "GUES" ] ; then
          file=$DATADIR/$timed/gues$(printf %05g $((10#$mem)))
        fi 
        echo "Runing posprocessing for file $file"
        icount=$((icount+1))
        out=$OUTDIR/$timed/
        outdate=$(date -ud "$ctime" +"%Y%m%d%_H%M%S")
        mkdir -p $out
        python $SRCDIR/post_ens_wrfout_to_netcdf.py "$file" "$out" "$mem" "$outdate" &>> $out/log_$mem &
        exec_post "$file" "$mem" "$out" "$outdate" &
        if [ $icount == $MAXJOBS ] ;then
          wait
          icount=0
        fi
      else 
        echo "Warning: Could not find analysis at $timed"
        echo "Continue to the next time."
      fi
      ctime=$(date -ud "$ctime UTC + $anafreq seconds" +"%F %T")
    done
  done 
wait
fi

#Loops for postprocessing DAFCST
#===============================================================================
if [ $POSTTYPE == "DAFCST" ] ; then
  inittime=$timestart
  while [ $(date -ud "$inittime" +%s) -le $(date -ud "$timeend" +%s) ] ;do
    timeini=$(date -ud "$inittime" +"%Y%m%d%H%M%S")
    if [ -d "$DATADIR/$timeini/" ] ; then

      for mem in $(seq -f $memfmt $NMEM) ; do
        ctime=$inittime
        forendtime=$(date -ud "$ctime UTC + $forlen seconds" +"%F %T")
        while [ $(date -ud "$ctime" +%s) -le $(date -ud "$forendtime" +%s) ] ; do
          timefor=$(date -ud "$ctime" +"%F_%T")
          file=$DATADIR/$timeini/$mem/wrfout_d01_$timefor
          echo "Runing posprocessing for file $file"
          icount=$((icount+1))
          out=$OUTDIR/$timeini/
          outdate=$(date -ud "$ctime UTC" +"%Y%m%d_%H%M%S")
          mkdir -p $out
          python $SRCDIR/post_ens_wrfout_to_netcdf.py "$file" "$out" "$mem" "$outdate" &>> $out/log_$mem &
          if [ $icount == $MAXJOBS ] ; then
            wait
            icount=0
          fi
          echo $ctime
          ctime=$(date -ud "$ctime UTC + $foroutfreq seconds" +"%F %T")
        done
      done

    else
       echo "Warning: Could not find forecast initialized at $inittime"
       echo "Continue to the next time."
    fi
    inittime=$(date -ud "$inittime UTC + $forinifreq seconds" +"%F %T")
  done

wait
fi


echo "done."
