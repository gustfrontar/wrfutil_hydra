#!/bin/bash
#=======================================================================
#   TO GENERATE MET_EM FILES THAT WILL BE USED TO RUN THE DA CYCLE AND
#   TO GENERATE THE PERTURBATIONS. (RUN THIS SCRIPT IN PPS SERVERS)
#=======================================================================
CONFIGURATION=OSAKA_1KM_NESTED
### directory settings

WPSDIR=${HOME}/share/WPSINTEL/
GEOGDIR=${HOME}/share/GEOG/
GRIBDIR=${HOME}/share/CFSR/                               # CFSR DATA FOR PERTURBATION GENERATION
OUTPUTDIR=${HOME}/share/INPUT/$CONFIGURATION/             # FINAL DESTINATION OF THE PERTURBATIONS.
TMPDIR=${HOME}/data/TMP/PREPARE_${CONFIGURATION}          # TEMPORARY WORK DIRECTORY.

GRIBSOURCEE=GFS
BOUNDARY_DATA_FREQE=21600
BOUNDARY_DATA_FREQO=300         #We can chose to output met_em at a higher frequency using ungrib interp.

GRIBDIRE=$HOME/share/FNL/

IE=20140910180000   #Experiment initial date.
EE=20140911060000   #Experiment end date.

MAX_DOM=2           #How many domains we will have.

CWD=`pwd`

source ../util.sh
##################################################
#Enviroment variables setting.
MPIBIN=mpiexec
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ra000015/a03094/libintel/lib/
ulimit -s unlimited
###################################################

###################################################
echo "Setting up work directories.."
echo " >> removing old work directories.."
echo "Removing $TMPDIR"
rm -rf $TMPDIR

echo "Copying files"
mkdir -p $TMPDIR/WPS/
mkdir -p $OUTPUTDIR
cp -r  $WPSDIR/*             $TMPDIR/WPS/
ln -sf $GEOGDIR              $TMPDIR/WPS/
cp     $CWD/../configuration/$CONFIGURATION/namelist.wps $TMPDIR/WPS/namelist.wps.template

chmod 766 $TMPDIR/WPS/*.sh $TMPDIR/WPS/*.csh
#######################################################################

cd $TMPDIR/WPS

#GENERATE MET_EM FOR THE EXPERIMENT

   mkdir -p $OUTPUTDIR/exp_met_em/

   #while [ $CDATE -le $EE ] ; do
   RUNNING=0
   # while [ $RUNNING -le $MAXRUNNING -a $CDATE -le $EE ] ; do
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
      rm DATA* *.grib2

      cp namelist.wps.template namelist.wps
      edit_namelist_wps ./namelist.wps $IE $EE $BOUNDARY_DATA_FREQO

      if [ ! -e geo_em.d01.nc ] ; then
         echo "Running geogrid"
         ./geogrid.exe > geogrid.log
      fi
      
      CDATE=$IE
      while [ $CDATE -le $EE ] ; do
       CDATES=`echo $CDATE | cut -c1-8`
       CHOUR=`echo $CDATE | cut -c9-10`
       FILE=fnl_${CDATES}_${CHOUR}_00.grib2
       ln -sf $GRIBDIRE/$FILE ./$CDATE.grib2
       CDATE=`date_edit2 $CDATE $BOUNDARY_DATA_FREQE` 
      done
       

      ./link_grib.csh *.grib2

      echo "./ungrib.exe > ./ungrib.log " > ./tmp.sh
      echo "./metgrid.exe ./metgrid.log " >> ./tmp.sh
      sh tmp.sh 
    
    #  CDATE=`date_edit2 $CDATE $BOUNDARY_DATA_FREQE`
    #  echo $CDATE
    #  RUNNING=`expr $RUNNING + 1 `

    #COPY THE DATA TO ITS FINAL DESTINATION
    mv $TMPDIR/$RUNNING/met_em* $OUTPUTDIR/exp_met_em/

   























