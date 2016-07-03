#PBS -l nodes=1:ppn=40
#PBS -S /bin/bash

#To compute some components of the ensemble covariance matrix.

MAX_RUNNING=40         #Number of threads
CDATE=20130713053400   #The date
MEMBER=1000            #Ensemble size
BOOTSTRAP=.false.
BOOTSTRAP_SAMPLES=100
SKIPX=1
SKIPY=1
SKIPZ=1
SMOOTHCOV=.false.
SMOOTHCOVLENGTH=1.0d4
SMOOTHDX=1.0d3
COMPUTEMOMENTS=.false.
COMPUTECOVAR=.false.
COMPUTEINDEX=.true.
DELTA=20
MAX_MOMENTS=4
ANALYSIS_PATH=/data1/jruiz/EXPERIMENTS/OsakaPAR_1km_control1000m_smallrandompert/            #Ensemble data path
TMPDIR=/data1/jruiz/TMP/covariance_matrix_index_1km_1000m/                                   #Temporary work directory.
CTL_PATH=/data1/jruiz/EXPERIMENTS/OsakaPAR_1km_control1000m_smallrandompert/ctl/analgz.ctl   #Ensemble data ctl file.

ulimit -s unlimited
export OMP_STACKSIZE=500M
export KMP_STACKSIZE=500M

source /data1/jruiz/LETKF_WRF/wrf/run/util.sh

#Create temporary folder
mkdir -p $TMPDIR
cp /data1/jruiz/LETKF_WRF/wrf/verification/covmatrix/*.exe $TMPDIR

#Write namelist
my_namelist=$TMPDIR/covariance_matrix.namelist
echo "&general                                                            " >  $my_namelist
echo "nbv=$MEMBER                                                         " >> $my_namelist
echo "npoints=18                                                          " >> $my_namelist
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
echo "computecovar=$COMPUTECOVAR                                          " >> $my_namelist
echo "computeindex=$COMPUTEINDEX                                          " >> $my_namelist
echo "max_moments=$MAX_MOMENTS                                            " >> $my_namelist
echo "delta=$DELTA                                                        " >> $my_namelist
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
  
    ln -sf $ANALYSIS_PATH/$CDATE/analgz/${MEM2}.grd ./fcst${MEM}.grd

    M=`expr $M + 1 `

   done
time ./covariance_matrix.exe #> cov.log 






