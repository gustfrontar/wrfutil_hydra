#PBS -l nodes=1:ppn=40
#PBS -S /bin/bash

#To compute some components of the ensemble covariance matrix.

MAX_RUNNING=40         #Number of threads
CDATE=20130713053400   #The date
MEMBER=100            #Ensemble size
BOOTSTRAP=.false.
BOOTSTRAP_SAMPLES=100
SKIPX=1
SKIPY=1
SKIPZ=1
SMOOTHCOV=.false.
SMOOTHCOVLENGTH=1.0d4
SMOOTHDX=1.0d3
COMPUTEMOMENTS=.true.
MAX_MOMENTS=4
COMPUTEHISTOGRAM=.true.
NBINS=50
ANALYSIS_PATH=/data1/jruiz/EXPERIMENTS/OsakaPAR_1km_control1000m_smallrandompert/            #Ensemble data path
TMPDIR=/data1/jruiz/TMP/test/                     #Temporary work directory.
CTL_PATH=/data1/jruiz/EXPERIMENTS/OsakaPAR_1km_control1000m_smallrandompert/ctl/analgz.ctl   #Ensemble data ctl file.

ulimit -s unlimited
export OMP_STACKSIZE=500M
export KMP_STACKSIZE=500M

source ../../run/util.sh

#Create temporary folder
mkdir -p $TMPDIR
cp ../covmatrix/*.exe $TMPDIR

#Write namelist
my_namelist=$TMPDIR/covariance_matrix.namelist
echo "&general                                                            " >  $my_namelist
echo "nbv=$MEMBER                                                         " >> $my_namelist
echo "npoints=18                                                          " >> $my_namelist
echo "plon=135.48,135.48,135.48,135.602,135.602,135.602,135.602,135.602,   \
          136.29,136.29,136.29,136.29,136.29,135.94,135.94,135.94,135.94,135.94,  " >> $my_namelist                                        
echo "plat=34.7,34.7,34.7,34.71,34.71,34.71,34.71,34.71,                   \
          34.801,34.801,34.801,34.801,34.801,34.88,34.88,34.88,34.88,34.88,       " >> $my_namelist
echo "pvarname=tk,p,u,w,tk,p,u,dbz,w,tk,dbz,p,qv,w,tk,dbz,p,qv,           " >> $my_namelist
echo "plev=5000,5000,5000,1957.5,1957.5,1957.5,1957.5,1957.5, \
            5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,            " >> $my_namelist
echo "dep=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,                            " >> $my_namelist 
echo "error=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,                          " >> $my_namelist
echo "bootstrap=$BOOTSTRAP                                                " >> $my_namelist
echo "bootstrap_samples=$BOOTSTRAP_SAMPLES                                " >> $my_namelist
echo "nignore=9                                                           " >> $my_namelist
echo "ignorevarname='qc','qr','qi','qs','qg','topo','rain','snow','max_dbz' " >> $my_namelist
echo "skipx=$SKIPX                                                        " >> $my_namelist
echo "skipy=$SKIPY                                                        " >> $my_namelist
echo "skipz=$SKIPZ                                                        " >> $my_namelist
echo "smoothcov=$SMOOTHCOV                                                " >> $my_namelist
echo "smoothcovlength=$SMOOTHCOVLENGTH                                    " >> $my_namelist
echo "smoothdx=$SMOOTHDX                                                  " >> $my_namelist
echo "computemoments=$COMPUTEMOMENTS                                      " >> $my_namelist
echo "max_moments=$MAX_MOMENTS                                            " >> $my_namelist
echo "computehistogram=$COMPUTEHISTOGRAM                                  " >> $my_namelist
echo "max_histogram=$NBINS                                                " >> $my_namelist
echo "/                                                                   " >> $my_namelist


#Link files
cd $TMPDIR
export OMP_NUM_THREADS=$MAX_RUNNING
export KMP_NUM_THREADS=$MAX_RUNNING

   cp -sf $CTL_PATH  ./input.ctl
   M=1

   while [ $M -le $MEMBER  ] ; do

    MEM=`add_zeros $M 5`
    MEM2=`add_zeros $M 4`
  
    ln -sf $ANALYSIS_PATH/$CDATE/guesgz/${MEM2}.grd ./fcst${MEM}.grd

    M=`expr $M + 1 `

   done
time ./covariance_matrix.exe 






