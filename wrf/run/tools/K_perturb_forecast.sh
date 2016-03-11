#This scripts reads the output from a wrf forecast (wrfout format) and perturbs some variables
#with the selected spatial scale. Perturbations are random and umbalanced (suitable only for
#small scale ensemble forecasting and data assimilation).

#Perturbations are smooth in space (horizontal and vertial) as well as in time.

#This perturbed forecast can then be used as initial and boundary conditions for small scale ensembles or DA
#In this case the script is prepared to generate initial and boundary conditions for SCALE-LETKF experiments.

source ../util.sh
FORECASTINIDATE=20130713000000 #Date of forecast initialization.
SOURCEMEMBER=00060             #Perturbations will be generated around a particular member (or may be the ensemble mean or a control forecast)
INIDATE=20130713051000         #Date of interpolation init  (yyyymmddhhmmss)
ENDDATE=20130713054000         #Date of interpolation end   (yyyymmddhhmmss)
FORECAST_OUTPUT_FREQ=30        #Forecast frequency.

PPSSERVER=pps2

DOMAIN=02                      #Domain that will be interpolated.
MEMBER=100                     #How many perturbations are we going to create.

MAXRUN=40

FORECASTBP=$HOME/share/FORECAST_OSAKA_1KM_NESTED_1024m_TIMEINT/forecast/ #Forecast basepath.

OUTPATH=$HOME/share/FORECAST_OSAKA_1KM_NESTED_100m_RANDOMPERT/forecast/  #Where the perturbed forecast will be created.

TMPDIR=$HOME/data/TMP/FORECAST_PERTURBATION/

EXEC=$HOME/share/LETKF_WRF/wrf/add_pert_wrfout/compute_pert_wrfout.exe

TIME_RANGE=`date_diff $ENDDATE $INIDATE`

PERTURB_T=.true.              
PERTURB_Q=.true.       
PERTURB_WIND=.true.         
PERTURB_T_AMP=0.4d0       
PERTURB_Q_AMP=0.4d-3       
PERTURB_WIND_AMP=0.5d0 
PERTURB_T_SCLH=1.0d4  
PERTURB_Q_SCLH=1.0d4      
PERTURB_WIND_SCLH=1.0d4
PERTURB_T_SCLV=0.5d4
PERTURB_Q_SCLV=0.5d4     
PERTURB_WIND_SCLV=0.5d4
DT=$FORECAST_OUTPUT_FREQ            
SCLT=300                         
NTIMES=`expr $TIME_RANGE \/ $FORECAST_OUTPUT_FREQ ` #Ntimes is automathically computed from the total number of times.            
NTIMES=`expr $NTIMES + 1 `

mkdir -p $OUTPATH/$FORECASTINIDATE/
mkdir -p $TMPDIR
cp $EXEC $TMPDIR

#Create namelist
namelist=$TMPDIR/perturb.namelist

echo "&general                            " >  $namelist
echo "PERTURB_T=$PERTURB_T                " >> $namelist
echo "PERTURB_Q=$PERTURB_Q                " >> $namelist
echo "PERTURB_WIND=$PERTURB_WIND          " >> $namelist
echo "PERTURB_T_AMP=$PERTURB_T_AMP        " >> $namelist
echo "PERTURB_Q_AMP=$PERTURB_Q_AMP        " >> $namelist
echo "PERTURB_WIND_AMP=$PERTURB_WIND_AMP  " >> $namelist
echo "PERTURB_T_SCLH=$PERTURB_T_SCLH      " >> $namelist
echo "PERTURB_Q_SCLH=$PERTURB_Q_SCLH      " >> $namelist
echo "PERTURB_WIND_SCLH=$PERTURB_WIND_SCLH" >> $namelist
echo "PERTURB_T_SCLV=$PERTURB_T_SCLV      " >> $namelist
echo "PERTURB_Q_SCLV=$PERTURB_Q_SCLV      " >> $namelist
echo "PERTURB_WIND_SCLV=$PERTURB_WIND_SCLV" >> $namelist
echo "DT=$FORECAST_OUTPUT_FREQ            " >> $namelist
echo "SCLT=$SCLT                          " >> $namelist
echo "NTIMES=$NTIMES                      " >> $namelist
echo "/                                   " >> $namelist

cd $TMPDIR

   M=1
   while [ $M -le $MEMBER  ] ; do

    RUN=1
    while [ $RUN -le $MAXRUN -a $M -le $MEMBER ] ; do
     MEM=` ens_member $M `

     WORKDIR=$TMPDIR/$RUN

     echo "Perturbing member $MEM "

     echo "mkdir -p $WORKDIR                                                                                             ">  tmp$RUN.sh
     echo "cd $WORKDIR                                                                                                   ">> tmp$RUN.sh
     echo "mkdir -p $OUTPATH/$FORECASTINIDATE/$MEM/                                                                      ">> tmp$RUN.sh
     CDATE=$INIDATE
     ITIME=1
      while [ $CDATE -le $ENDDATE ] ; do
       if [ $ITIME -lt 10 ] ; then
         ITIME=0$ITIME
       fi
       if [ $ITIME -lt 100 ] ; then
         ITIME=0$ITIME
       fi

       my_file=`wrfout_file_name $CDATE $DOMAIN `
       echo "cp $FORECASTBP/$FORECASTINIDATE/$SOURCEMEMBER/$my_file  $OUTPATH/$FORECASTINIDATE/$MEM/                     ">> tmp$RUN.sh
       echo "ln -sf $OUTPATH/$FORECASTINIDATE/$MEM/$my_file       $WORKDIR/wrfout${ITIME}.nc                             ">> tmp$RUN.sh
       CDATE=`date_edit2 $CDATE $FORECAST_OUTPUT_FREQ`
       ITIME=`expr $ITIME + 1 `
     done

     echo "cp $namelist $WORKDIR                                                                                         ">> tmp$RUN.sh
     echo "$EXEC  > perturbwrf.log                                                                                       ">> tmp$RUN.sh
     chmod 755 tmp$RUN.sh

     ssh $PPSSERVER $TMPDIR/tmp$RUN.sh       &

    M=`expr $M + 1 `
    RUN=`expr $RUN + 1 `
    done

    time wait


   done








