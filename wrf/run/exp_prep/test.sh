#!/bin/bash
#=======================================================================
#   TO GENERATE MET_EM FILES THAT WILL BE USED TO RUN THE DA CYCLE AND
#   TO GENERATE THE PERTURBATIONS. (RUN THIS SCRIPT IN PPS SERVERS)
#   This scripts generates met_em to compute the perturbations of the boundary
#   conditions. In a nested domain this is usually performed only at domain 1
#=======================================================================
CONFIGURATION=CORDOBA_2K
### directory settings

WPSDIR=${HOME}/share/LETKF_WRF/wrf/model/WPS/
GEOGDIR=${HOME}/share/LETKF_WRF/wrf/model/GEOG/
OUTPUTDIR=${HOME}/share/INPUT/$CONFIGURATION/             # FINAL DESTINATION OF THE PERTURBATIONS.
TMPDIR=${HOME}/data/TMP/PREPARE_${CONFIGURATION}_DB       # TEMPORARY WORK DIRECTORY.

MAXRUNNING=2 #Maximum simultaneous process running in pps servers

GRIBSOURCE=CFSR
BOUNDARY_DATA_FREQ=21600    #Boundary data perturbations will be output at the oringinal data frequency.
GRIBDIR=$HOME/share/DATA/CFSR_LOWRES/

IE=20080101000000   #Experiment initial date.
EE=20101231180000   #Experiment end date.

MAX_DOM=1


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
rm   -rf $TMPDIR

echo  "Copying files"
mkdir -p $TMPDIR/WPS/
mkdir -p $OUTPUTDIR
cp -r  $WPSDIR/*             $TMPDIR/WPS/
ln -sf $GEOGDIR              $TMPDIR/WPS/
cp     $CWD/../configuration/$CONFIGURATION/namelist.wps $TMPDIR/WPS/namelist.wps.template


chmod 766 $TMPDIR/WPS/*.sh $TMPDIR/WPS/*.csh
#######################################################################
cd $TMPDIR/WPS

#GENERATE MET_EM FOR THE DATABASE
   
   CDATE=$IE
   mkdir -p $OUTPUTDIR/db_met_em/

   while [ $CDATE -le $EE ] ; do
   RUNNING=1
    while [ $RUNNING -le $MAXRUNNING -a $CDATE -le $EE ] ; do
     if [ ! -e $TMPDIR/$RUNNING/ ] ; then
        mkdir -p $TMPDIR/$RUNNING/
        ln -sf $TMPDIR/WPS/*.exe   $TMPDIR/$RUNNING/
        ln -sf $TMPDIR/WPS/ungrib  $TMPDIR/$RUNNING/
        ln -sf $TMPDIR/WPS/metgrid $TMPDIR/$RUNNING/
        ln -sf $TMPDIR/WPS/geogrid $TMPDIR/$RUNNING/
        ln -sf $TMPDIR/WPS/namelist.wps.template $TMPDIR/$RUNNING/
        ln -sf $TMPDIR/WPS/ungrib/Variable_Tables/Vtable.${GRIBSOURCE} $TMPDIR/$RUNNING/Vtable
        ln -sf $TMPDIR/WPS/link_grib.csh $TMPDIR/$RUNNING/
        ln -sf $GEOGDIR              $TMPDIR/$RUNNING/
      fi
  
      cd $TMPDIR/$RUNNING/
      rm DATA* *.grib 

      cp namelist.wps.template namelist.wps
      edit_namelist_wps ./namelist.wps $CDATE $CDATE $BOUNDARY_DATA_FREQ

      if [ ! -e geo_em.d01.nc ] ; then
         echo "Running geogrid"
         ./geogrid.exe > geogrid.log
      fi

      CDATES=`echo $CDATE | cut -c1-10`
      FILE=pgbl00.gdas.${CDATES}.grb2
      if [ ! -e $GRIBDIR/$FILE ] ; then
        #If file does not exist download it from the web.
        download_cfsr $CDATE $CDATE 21600 $GRIBDIR
      fi
    
      echo  ${GRIBDIR}/${FILE} 
      ln -sf ${GRIBDIR}/${FILE} ./$CDATE.grib
      
      ./link_grib.csh *.grib

      echo "./ungrib.exe > ./ungrib.log " > ./tmp.sh
      echo "./metgrid.exe ./metgrid.log " >> ./tmp.sh
      sh tmp.sh &

      CDATE=`date_edit2 $CDATE $BOUNDARY_DATA_FREQ`
      echo $CDATE
      RUNNING=`expr $RUNNING + 1 `

    done

    time wait
    #Copying data
    RUNNING=1
    while [ $RUNNING -le $MAXRUNNING ] ; do
      mv $TMPDIR/$RUNNING/met_em* $OUTPUTDIR/db_met_em/   

      RUNNING=`expr $RUNNING + 1 `  
    done
   done






