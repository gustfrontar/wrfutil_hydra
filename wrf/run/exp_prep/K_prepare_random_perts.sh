#!/bin/bash
#. /usr/share/modules/init/sh
#module unload pgi-12.10
#module load intel-2013.1.117

#=======================================================================
#   This script compute the dates that has to be substracted in order 
#   to compute the perturbations for the initial and boundary conditions.
#=======================================================================
MEMBER=5000              #Number of ensemble members
MAX_TIMES=1440           #Maximum number of different states to generate perturbations.
EXP_LENGTH=5             #Number of assimilation cycles in the experiment.
CONFIGURATION=OSAKA_1KM_DOWNSCALLING #SINLAKU_60K
### directory settings
OUTPUTDIR=${HOME}/share/INPUT/$CONFIGURATION/pert_date/       # FINAL DESTINATION OF THE PERTURBATIONS.
BOUNDARY_DATA_FREQ=6
CWD=`pwd`

#PERTURBATION DATA BASE INFO.
PINI=20060101000000
PEND=20091231000000

#EXPERIMENT DATES
EINI=20130713000000
EEND=20130714000000

source ../util.sh
##################################################
#Enviroment variables setting.

###################################################
# Usually do not modify below
#-----------------------------------------------------------------------
echo "Setting up work directories.."
echo " >> removing old work directories.."
#rm -rf $OUTPUTDIR
mkdir -p $OUTPUTDIR
cd $OUTPUTDIR


M=1

while [ $M -lt $MEMBER ]
do

  echo "Computing dates for ensemble member $M"
  MEM=`ens_member $M `

  PROCEED=0
  while [ $PROCEED -eq 0  ]
  do
 
  #Generate random numbers to pick random dates.
  #First date is totally random. Second date is within 5-25 days from the previous date. 
  
  TMPNUM=`expr $EXP_LENGTH \/ 4 `
  TMPTIME=`expr $MAX_TIMES - $TMPNUM - 25 `
  number1=$RANDOM
  let "number1 %=$TMPTIME "
  number2=$RANDOM
  let "number2 %=20"
  number2=`expr $number1 + $number2 + 5`

  PROCEED=1
  #Test to see if number1 and number2 has been chosen before.
  if [ $M -gt 1 ]
  then
    for element in $(seq 1 `expr $M - 1 `  )
    do
     if [ ${hist_number1[$element]} -eq $number1 ] && [ ${hist_number2[$element]} -eq $number2  ]
     then
       PROCEED=0
       echo "Number 1 and number 2 has been chosen before, let's chose another combination"
       echo "If this happens a lot, increase the size of the analysis data base "
     fi
    done
  fi
  if [ $PROCEED -eq 1 ]
  then
  hist_number1[$M]=$number1
  hist_number2[$M]=$number2
  fi

  done

  #Get dates corresponding to these random numbers
  #2009 has been used for the generation of the perturbations.
  INTERVAL=`expr $number1 \* 24 \* 60 \* 60 `
  DATE1=`date_edit2 $PINI $INTERVAL` 

  INTERVAL=`expr $number2 \* 24 \* 60 \* 60 `
  DATE2=`date_edit2 $PINI $INTERVAL` 
 
  echo $DATE1 $DATE2

  CDATE=$EINI

  INC=`expr $BOUNDARY_DATA_FREQ \* 60 \* 60 `
  
  CONT=0
  while [ $CDATE -le $EEND ]
  do
   INC2=`expr $INC \* $CONT `
   CDATE1=`date_edit2 $DATE1 $INC2` 
   CDATE2=`date_edit2 $DATE2 $INC2`
   echo "$CDATE1 $CDATE2 " > $OUTPUTDIR/${CDATE}_${MEM}.dates 
   
   CDATE=`date_edit2 $CDATE $INC `
   CONT=`expr $CONT + 1 `
 done


 M=`expr $M + 1 `
done


