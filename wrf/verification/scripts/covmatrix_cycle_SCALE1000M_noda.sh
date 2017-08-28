#!/bin/bash -l
#PBS -l nodes=1:ppn=40
#PBS -S /bin/bash

module load intel

#To compute some components of the ensemble covariance matrix.
BINPATH=/data10/jruiz/LETKF_WRF/wrf/verification/covmatrix/
UTIL=/data10/jruiz/LETKF_WRF/wrf/run/util.sh

DATAPATH=/data10/jruiz/EXPERIMENTS/

EXPERIMENT=OsakaPAR_1km_control1000m_smallrandompert_noda

TYPE="analgp"

INIDATE=20130713052500
ENDDATE=20130713054000

INCREMENT=60   #Time resolution in seconds.

MAX_RUNNING=20        #Number of threads
MEMBER=1000             #Ensemble size
BOOTSTRAP=.false.
BOOTSTRAP_SAMPLES=100
SKIPX=1
SKIPY=1
SKIPZ=1
SMOOTHCOV=.false.
SMOOTHCOVLENGTH=20.0d4
SMOOTHDX=1.0d3
COMPUTEMOMENTS=.false.
MAX_MOMENTS=4
COMPUTEHISTOGRAM=.false.
NBINS=50
NBASE_VARS=4
NCOV_VARS=3
COMPUTEINDEX=.true.
DELTA=40
BASE_VARS="'qhydro','u','v','tk'"
COV_VARS="'qhydro','u','v'"
TMPDIR=/home/jruiz/TMP/$EXPERIMENT/                      #Temporary work directory.

ulimit -s unlimited
export OMP_STACKSIZE=500M
export KMP_STACKSIZE=500M

source $UTIL

#Create temporary folder
mkdir -p $TMPDIR
cp $BINPATH/covariance_matrix.exe $TMPDIR

#Write namelist
my_namelist=$TMPDIR/covariance_matrix.namelist
echo "&general                                                            " >  $my_namelist
echo "nbv=$MEMBER                                                         " >> $my_namelist
echo "npoints=0                                                           " >> $my_namelist
echo "pvarname="'tk','p','u','w','tk','p','u','dbz','w',"                 " >> $my_namelist
echo "bootstrap=$BOOTSTRAP                                                " >> $my_namelist
echo "bootstrap_samples=$BOOTSTRAP_SAMPLES                                " >> $my_namelist
echo "nignore=0                                                           " >> $my_namelist
echo "skipx=$SKIPX                                                        " >> $my_namelist
echo "skipy=$SKIPY                                                        " >> $my_namelist
echo "skipz=$SKIPZ                                                        " >> $my_namelist
echo "smoothcov=$SMOOTHCOV                                                " >> $my_namelist
echo "smoothcovlength=$SMOOTHCOVLENGTH                                    " >> $my_namelist
echo "smoothdx=$SMOOTHDX                                                  " >> $my_namelist
echo "computemoments=$COMPUTEMOMENTS                                      " >> $my_namelist
echo "computehistogram=$COMPUTEHISTOGRAM                                  " >> $my_namelist
echo "computeindex=$COMPUTEINDEX                                          " >> $my_namelist
echo "delta=$DELTA                                                        " >> $my_namelist
echo "nbase_vars=$NBASE_VARS                                              " >> $my_namelist
echo "ncov_vars=$NCOV_VARS                                                " >> $my_namelist
echo "base_vars=$BASE_VARS                                                " >> $my_namelist
echo "cov_vars=$COV_VARS                                                  " >> $my_namelist
echo "/                                                                   " >> $my_namelist


#Link files
mkdir -p $TMPDIR
cd $TMPDIR
export OMP_NUM_THREADS=$MAX_RUNNING
export KMP_NUM_THREADS=$MAX_RUNNING


CDATE=$INIDATE

 while [ $CDATE -le $ENDDATE  ] ; do

   cp  $DATAPATH/$EXPERIMENT/ctl/${TYPE}.ctl  ./input.ctl

   rm $TMPDIR/fcst?????.grd  

   M=1

   while [ $M -le $MEMBER  ] ; do

    MEM=`add_zeros $M 5`
    MEM2=`add_zeros $M 4`
  
    ln -sf $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/${MEM2}.grd ./fcst${MEM}.grd

    M=`expr $M + 1 `

   done


   time ./covariance_matrix.exe #> cov.log 


   mv covindex*.grd       $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/
   mv corrindex*.grd      $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/   
   mv corrdistindex*.grd  $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/

   mv cov_profile*.txt    $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/
   mv corr_profile*.txt   $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/ 
   mv num_profile*.txt    $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/

   mv covmat_x???y???z???.grd $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/
   mv bssprd_x???y???z???.grd $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/
   mv impact_x???y???z???.grd $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/
   mv bsmean_x???y???z???.grd $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/
   mv histgr_x???y???z???.grd $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/

   mv histogram.grd           $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/
   mv minvar.grd              $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/
   mv maxvar.grd              $DATAPATH/$EXPERIMENT/$CDATE/$TYPE/

   CDATE=`date_edit2 $CDATE $INCREMENT `

 done

