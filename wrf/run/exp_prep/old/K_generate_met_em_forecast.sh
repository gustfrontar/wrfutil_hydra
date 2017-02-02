#!/bin/bash
#=======================================================================
#   TO GENERATE MET_EM FILES THAT WILL BE USED TO RUN THE DA CYCLE AND
#   TO GENERATE THE PERTURBATIONS. (RUN THIS SCRIPT IN PPS SERVERS)
#=======================================================================
DOMAINCONF=SINLAKU_20K
### directory settings
MAX_DOMAIN=1

WPSDIR=${HOME}/share/WPSINTEL/
GEOGDIR=${HOME}/share/GEOG/
GRIBDIR=${HOME}/share/CFSR/                             # CFSR DATA FOR PERTURBATION GENERATION
OUTPUTDIR=${HOME}/share/INPUT/$DOMAINCONF/              # FINAL DESTINATION OF THE PERTURBATIONS.
TMPDIR=${HOME}/data/TMP/PREPARE_${DOMAINCONF}           # TEMPORARY WORK DIRECTORY.

MAXRUNNING=10 #Maximum simultaneous process running in pps servers
FORECASTMAXLENGTH=72   #Maximum forecast length
FORECASTOUTPUTFREQ=3   #Forecast output frequency.

FORECASTMAXLENGTHSEC=`expr $FORECASTMAXLENGTH \* 3600 `
FORECASTOUTPUTFREQSEC=`expr $FORECASTOUTPUTFREQ \* 3600 `

GRIBSOURCEE=GFS
FORECASTFREQ=6         #Forecast initialization frequency
FORECASTFREQSEC=`expr $FORECASTFREQ \* 3600 `
GRIBDIRE=${HOME}/share/GFS/

IE=20080820000000   #Experiment initial date.
EE=20080930000000   #Experiment end date.

CWD=`pwd`

source ../util.sh
##################################################
#Enviroment variables setting.
MPIBIN=mpiexec
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HOME}/share/libintel/lib/
ulimit -s unlimited
###################################################

###################################################
echo "Setting up work directories.."
echo " >> removing old work directories.."
echo "Removing $TMPDIR"
rm -rf $TMPDIR

echo "Copying files"
mkdir -p $TMPDIR/WPS/
echo "Creating outputdir $OUTPUTDIR"
mkdir -p $OUTPUTDIR


cp -r  $WPSDIR/*             $TMPDIR/WPS/
ln -sf $GEOGDIR              $TMPDIR/WPS/
cp  $CWD/../configuration/$DOMAINCONF/namelist.wps $TMPDIR/WPS/namelist.wps.template

chmod 766 $TMPDIR/WPS/*.sh $TMPDIR/WPS/*.csh
#######################################################################

cd $TMPDIR/WPS

#GENERATE MET_EM FOR THE EXPERIMENT

   CDATE=$IE

   mkdir -p $OUTPUTDIR/for_met_em/

   while [ $CDATE -le $EE ] ; do
   RUNNING=0
   echo "Preparing met_em files corresponding to the forecast initialized at $CDATE"
    while [ $RUNNING -le $MAXRUNNING -a $CDATE -le $EE ] ; do
     if [ ! -e $TMPDIR/$RUNNING/ ] ; then
        mkdir -p $TMPDIR/$RUNNING/
        ln -sf $TMPDIR/WPS/*.exe   $TMPDIR/$RUNNING/
        ln -sf $TMPDIR/WPS/ungrib  $TMPDIR/$RUNNING/
        ln -sf $TMPDIR/WPS/metgrid $TMPDIR/$RUNNING/
        ln -sf $TMPDIR/WPS/geogrid $TMPDIR/$RUNNING/
        ln -sf $TMPDIR/WPS/namelist.wps.template $TMPDIR/$RUNNING/
        ln -sf $TMPDIR/WPS/ungrib/Variable_Tables/Vtable.${GRIBSOURCEE} $TMPDIR/$RUNNING/Vtable
        ln -sf $TMPDIR/WPS/link_grib.csh $TMPDIR/$RUNNING/
        ln -sf $GEOGDIR              $TMPDIR/$RUNNING/
      fi
  
      cd $TMPDIR/$RUNNING/
      rm DATA* *.grib 

      FDATE=`date_edit2 $CDATE $FORECASTMAXLENGTHSEC `

      cp namelist.wps.template namelist.wps
      edit_namelist_wps ./namelist.wps $CDATE $FDATE $FORECASTOUTPUTFREQSEC

      if [ ! -e geo_em.d01.nc ] ; then
         echo "Running geogrid"
         ./geogrid.exe > geogrid.log
      fi


      CDATES=`echo $CDATE | cut -c1-12`
      FDIR=$GRIBDIRE/GFSFORECAST${CDATES}
      TMPDATE=$CDATE
   echo "Generating met_em for forecast $GRIBDIRE/GFSFORECAST${CDATES} "
   if [ -e $FDIR ] ; then #Some forecasts are missing.
      while [ $TMPDATE -le $FDATE ] ; do
         TMPDATE2=`echo $TMPDATE | cut -c1-12`
         ln -sf $FDIR/${TMPDATE2}.grib2 ./${TMPDATE2}.grib
         TMPDATE=`date_edit2 $TMPDATE $FORECASTOUTPUTFREQSEC ` 
      done

      ./link_grib.csh *.grib

      echo "./ungrib.exe > ./ungrib.log              " >  ./tmp.sh
      echo "./metgrid.exe > ./metgrid.log            " >> ./tmp.sh
      echo "mkdir -p $OUTPUTDIR/for_met_em/$CDATE/   " >> ./tmp.sh
      echo "mv met_em* $OUTPUTDIR/for_met_em/$CDATE/ " >> ./tmp.sh
      echo "Running ungrib and metgrid"
      time sh tmp.sh  &
   else
     echo "[Warning]: Forecast for $CDATE is missing."
   fi
      
      CDATE=`date_edit2 $CDATE $FORECASTFREQSEC `
      echo $CDATE
      RUNNING=`expr $RUNNING + 1 `
      echo "CDATE = $CDATE "

    done

    RUNDATES=`expr $RUNNING - 1 `

    time wait
   done
   






















