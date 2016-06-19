

MAX_RUNNING=5
CDATE=20100806060000
MEMBER=36
ANALYSIS_PATH=$HOME/datos/EXPERIMENTS/ANALYSIS_SA_60KM_control40msa_test/
TMPDIR=$HOME/data/TMP/covariance_matrix/

ulimit -s unlimited
export OMP_STACKSIZE=200M
export KMP_STACKSIZE=200M

source ../../run/util.sh

#Create temporary folder
mkdir -p $TMPDIR
cp ../covmatrix/*.exe $TMPDIR

#Write namelist

my_namelist=$TMPDIR/covariance_matrix.namelist

echo "&general                                     " >  $my_namelist
echo "nbv=$MEMBER                                  " >> $my_namelist
echo "npoints=2                                    " >> $my_namelist
echo "plon=-58,-58,                                " >> $my_namelist                                        
echo "plat= -10, -30,                              " >> $my_namelist
echo "pvarname='umet','tk'                         " >> $my_namelist
echo "plev=500,1000                                " >> $my_namelist
echo "dep=1,1,                                     " >> $my_namelist 
echo "error=1,1,                                   " >> $my_namelist
echo "/                                            " >> $my_namelist


#Link files
cd $TMPDIR

export OMP_NUM_THREADS=$MAX_RUNNING

   cp -sf $ANALYSIS_PATH/anal/$CDATE/plev00001.ctl  ./input.ctl
   M=1
   while [ $M -le $MEMBER  ] ; do

    MEM=`ens_member $M `
  
    ln -sf $ANALYSIS_PATH/anal/$CDATE/plev${MEM}.dat ./fcst${MEM}.grd
   

    M=`expr $M + 1 `
   done


time ./covariance_matrix.exe


for i in `ls covmat*.grd ` ; do

  ctlname=`echo $i | cut -c1-19 `

  cp ./input.ctl ./${ctlname}.ctl

  sed -i 's/plev00001.dat/'${ctlname}.grd'/g'   ./${ctlname}.ctl

done

for i in `ls  impact*.grd ` ; do

  ctlname=`echo $i | cut -c1-19 `

  cp ./input.ctl ${ctlname}.ctl

  sed -i 's/plev00001.dat/'${ctlname}.grd'/g'   ./${ctlname}.ctl

done



