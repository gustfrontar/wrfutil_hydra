

MAX_RUNNING=10
CDATE=20080827180000
MEMBER=40
ANALYSIS_PATH=/home/hp150019/k02016/data/EXPERIMENTS/ANALYSIS_SINLAKU_60K_control40m/
TMPDIR=$HOME/data/TMP/covariance_matrix/



source ../../run/util.sh

mkdir -p $TMPDIR
cp ../*.exe $TMPDIR


#Write namelist

my_namelist=$TMPDIR/covariance_matrix.namelist

echo "&general                                     " >  $my_namelist
echo "nbv=$MEMBER                                  " >> $my_namelist
echo "npoints=5                                    " >> $my_namelist
echo "plon=150,150, 150, 150,150                   " >> $my_namelist                                        
echo "plat= 40, 30,  20,  10,  5                   " >> $my_namelist
echo "plev= 90, 90,  90,  90, 90                   " >> $my_namelist
echo "dep=1,1,1,1                                  " >> $my_namelist 
echo "error=1,1,1,1                                " >> $my_namelist
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


./covariance_matrix.exe


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



