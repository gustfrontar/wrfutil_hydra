#!/bin/bash
# =======================================================================
#
#       Utility Shell Finctions for WRF_LETKF
#
#                                                   2010.05.11 M.Kunii
# =======================================================================

# -----------------------------
#    date_edit
# -----------------------------
date_edit () {
(

        if [ $# -lt 7 ]; then
                echo "Usage : date_edit"
                echo "    date_edit [yyyy] [mm] [dd] [hh] [mn] [dt(min)]"
                echo "    ex) date_edit 201005051200 -180"
                exit
        fi

        yy=$1
        mm=$2
        dd=$3
        hh=$4
        mn=$5
        ss=$6
        dt=$7

        echo $yy-$mm-$dd $hh:$mn:$ss

        seconds=`date +%s -d"$yy-$mm-$dd $hh:$mn:$ss UTC"`

        seconds=`expr $seconds + $dt \* 60 `


        date -u +%Y%m%d%H%M%S -d"@$seconds "


)
}


date_edit2 () {
(

        if [ $# -lt 2 ]; then
                echo "Usage : date_edit"
                echo "    date_edit [yyyy][mm][dd][hh][mn] [dt(sec)]"
                echo "    ex) date_edit 201005051200 -60"
                exit
        fi

        CDATEL=$1
        dt=$2

        cy=`echo $CDATEL | cut -c1-4`
        cm=`echo $CDATEL | cut -c5-6`
        cd=`echo $CDATEL | cut -c7-8`
        ch=`echo $CDATEL | cut -c9-10`
        cn=`echo $CDATEL | cut -c11-12`
        cs=`echo $CDATEL | cut -c13-14`

        seconds=`date +%s -d"$cy-$cm-$cd $ch:$cn:$cs UTC"`

        seconds=`expr $seconds + $dt `


        date -u +%Y%m%d%H%M%S -d"@$seconds "


)
}

date_edit3 () {
(

        if [ $# -lt 2 ]; then
                echo "Usage : date_edit"
                echo "    date_edit [yyyy][mm][dd][hh][mn] [dt(min)]"
                echo "    ex) date_edit 201005051200 -180"
                exit
        fi

        CDATEL=$1
        dt=$2

        cy=`echo $CDATEL | cut -c1-4`
        cm=`echo $CDATEL | cut -c5-6`
        cd=`echo $CDATEL | cut -c7-8`
        ch=`echo $CDATEL | cut -c9-10`
        cn=`echo $CDATEL | cut -c11-12`
        cs=`echo $CDATEL | cut -c13-14`

        seconds=`date +%s -d"$cy-$cm-$cd $ch:$cn:$cs UTC"`

        seconds=`expr $seconds + $dt \* 60 `


        date -u +%Y\ %m\ %d\ %H\ %M\ %S -d"@$seconds "

)
}

#Compute the difference in seconds between two dates.
date_diff () {
(

        if [ $# -lt 2 ]; then
                echo "Usage : date_diff"
                echo "    date_diff [yyyy1][mm1][dd1][hh1][mn1] [yyyy2][mm2][dd2][hh2][mn2]"
                echo "    ex) date_edit 201005051200 201005051230"
                exit 1 
        fi

        local DATE1=$1
        local DATE2=$2

        cy1=`echo $DATE1 | cut -c1-4`
        cm1=`echo $DATE1 | cut -c5-6`
        cd1=`echo $DATE1 | cut -c7-8`
        ch1=`echo $DATE1 | cut -c9-10`
        cn1=`echo $DATE1 | cut -c11-12`
        cs1=`echo $DATE1 | cut -c13-14`

        seconds1=`date +%s -d"$cy1-$cm1-$cd1 $ch1:$cn1:$cs1 UTC"`

        cy2=`echo $DATE2 | cut -c1-4`
        cm2=`echo $DATE2 | cut -c5-6`
        cd2=`echo $DATE2 | cut -c7-8`
        ch2=`echo $DATE2 | cut -c9-10`
        cn2=`echo $DATE2 | cut -c11-12`
        cs2=`echo $DATE2 | cut -c13-14`

        seconds2=`date +%s -d"$cy2-$cm2-$cd2 $ch2:$cn2:$cs2 UTC"`


        echo ` expr $seconds1 - $seconds2 `
)
}

#Given a date and a time interval get the lower closest date which is a multiple of the interval.
date_floor () {
        if [ $# -lt 2 ]; then
                echo "Usage : date_floor"
                echo "    date_floor [yyyy][mm][dd][hh][mn] [interval(seconds)]"
                echo "    ex) date_edit 201005051200 3600"
                exit 1
        fi

        local DATE=$1
        local INTERVAL=$2

        cy1=`echo $DATE | cut -c1-4`
        cm1=`echo $DATE | cut -c5-6`
        cd1=`echo $DATE | cut -c7-8`
        ch1=`echo $DATE | cut -c9-10`
        cn1=`echo $DATE | cut -c11-12`
        cs1=`echo $DATE | cut -c13-14`

        seconds1=`date +%s -d"$cy1-$cm1-$cd1 $ch1:$cn1:$cs1 UTC"`
        mod=`expr $seconds1 % $INTERVAL `

        DATEFLOOR=`date_edit2 $DATE -$mod `
        echo $DATEFLOOR
}


add_zeros() {

local number=$1
local size=$2

local result=`printf "%0${size}d" $number`

echo $result

}



ungrib_file_name () {
local DATE="$1"
local PREFIX="$2"

    cy=`echo $DATE | cut -c1-4`
    cm=`echo $DATE | cut -c5-6`
    cd=`echo $DATE | cut -c7-8`
    ch=`echo $DATE | cut -c9-10`
    cn=`echo $DATE | cut -c11-12`
    cs=`echo $DATE | cut -c13-14`

    echo ${PREFIX}:${cy}-${cm}-${cd}_${ch}:${cn}
}

met_em_file_name () {

local    DATE="$1"
local    DOMAIN="$2"
  
    cy=`echo $DATE | cut -c1-4`
    cm=`echo $DATE | cut -c5-6`
    cd=`echo $DATE | cut -c7-8`
    ch=`echo $DATE | cut -c9-10`
    cn=`echo $DATE | cut -c11-12`
    cs=`echo $DATE | cut -c13-14`

    echo met_em.d${DOMAIN}.${cy}-${cm}-${cd}_${ch}:${cn}:${cs}.nc
}

wrfout_file_name () {

local    DATE="$1"
local    DOMAIN="$2"
  
    cy=`echo $DATE | cut -c1-4`
    cm=`echo $DATE | cut -c5-6`
    cd=`echo $DATE | cut -c7-8`
    ch=`echo $DATE | cut -c9-10`
    cn=`echo $DATE | cut -c11-12`
    cs=`echo $DATE | cut -c13-14`

    echo wrfout_d${DOMAIN}_${cy}-${cm}-${cd}_${ch}:${cn}:${cs}
}

date_in_wrf_format() {


local    DATE="$1"

    cy=`echo $DATE | cut -c1-4`
    cm=`echo $DATE | cut -c5-6`
    cd=`echo $DATE | cut -c7-8`
    ch=`echo $DATE | cut -c9-10`
    cn=`echo $DATE | cut -c11-12`
    cs=`echo $DATE | cut -c13-14`

    echo ${cy}-${cm}-${cd}_${ch}:${cn}:${cs}

}

date_in_scale_format() {

local DATE="$1"

    cy=`echo $DATE | cut -c1-4`
    cm=`echo $DATE | cut -c5-6`
    cd=`echo $DATE | cut -c7-8`
    ch=`echo $DATE | cut -c9-10`
    cn=`echo $DATE | cut -c11-12`
    cs=`echo $DATE | cut -c13-14`

    echo ${cy}${cm}${cd}-${ch}${cn}${cs}.000

}


edit_namelist () {

local    NAMELIST=$1
local    IDATE=$2
local    SCALEDATE=`date_in_scale_format $IDATE `

    iy=`echo $IDATE | cut -c1-4`
    im=`echo $IDATE | cut -c5-6`
    id=`echo $IDATE | cut -c7-8`
    ih=`echo $IDATE | cut -c9-10`
    in=`echo $IDATE | cut -c11-12`
    is=`echo $IDATE | cut -c13-14`

    sed -i 's/SYYYY/'${iy}'/g'                                   $NAMELIST
    sed -i 's/SMM/'${im}'/g'                                     $NAMELIST
    sed -i 's/SDD/'${id}'/g'                                     $NAMELIST 
    sed -i 's/SHH/'${ih}'/g'                                     $NAMELIST
    sed -i 's/SMN/'${in}'/g'                                     $NAMELIST
    sed -i 's/SSS/'${is}'/g'                                     $NAMELIST
    sed -i 's|DATABASEROOT|'${DATABASEROOT}'|g'                  $NAMELIST
    sed -i 's|BDYDATAFREQ|'${BDYDATAFREQ}'|g'                    $NAMELIST
    sed -i 's|BDYNFILES|'${BDYNFILES}'|g'                        $NAMELIST 
    sed -i 's|FOUTPUTFREQ|'${FOUTPUTFREQ}'|g'                    $NAMELIST
    sed -i 's|FLENGTH|'${FLENGTH}'|g'                            $NAMELIST
    sed -i 's|INITDATE|'$SCALEDATE'|g'                           $NAMELIST
    sed -i 's|FRESTARTFREQ|'${FRESTARTFREQ}'|g'                  $NAMELIST
    sed -i 's|PRC_NUM_X_IN|'${PRC_NUM_X_IN}'|g'                  $NAMELIST
    sed -i 's|PRC_NUM_Y_IN|'${PRC_NUM_Y_IN}'|g'                  $NAMELIST
    sed -i 's|IMAX_IN|'${IMAX_IN}'|g'                            $NAMELIST
    sed -i 's|JMAX_IN|'${JMAX_IN}'|g'                            $NAMELIST

}

edit_namelist_postproc () {

local    NAMELIST=$1
local    IDATE=$2

    iy=`echo $IDATE | cut -c1-4`
    im=`echo $IDATE | cut -c5-6`
    id=`echo $IDATE | cut -c7-8`
    ih=`echo $IDATE | cut -c9-10`
    in=`echo $IDATE | cut -c11-12`
    is=`echo $IDATE | cut -c13-14`

    sed -i 's/SYYYY/'${iy}'/g'                      $NAMELIST
    sed -i 's/SMM/'${im}'/g'                        $NAMELIST
    sed -i 's/SDD/'${id}'/g'                        $NAMELIST
    sed -i 's/SHH/'${ih}'/g'                        $NAMELIST
    sed -i 's/SMN/'${in}'/g'                        $NAMELIST
    sed -i 's/SSS/'${is}'/g'                        $NAMELIST
    sed -i 's|Z_LEV_TYPE_IN|'${Z_LEV_TYPE_IN}'|g'   $NAMELIST
    sed -i 's|Z_MERGE_OUT_IN|'${Z_MERGE_OUT_IN}'|g' $NAMELIST
    sed -i 's|T_MERGE_OUT_IN|'${T_MERGE_OUT_IN}'|g' $NAMELIST
    sed -i 's|MAPPROJ_ctl_IN|'${MAPPROJ_ctl_IN}'|g' $NAMELIST
    sed -i 's|VNAME_IN|'${VNAME_IN}'|g'             $NAMELIST
    sed -i 's|TARGET_ZLEV_IN|'${TARGET_ZLEV_IN}'|g' $NAMELIST
    sed -i 's|START_STEP_IN|'${START_STEP_IN}'|g'   $NAMELIST
    sed -i 's|END_STEP_IN|'${END_STEP_IN}'|g'       $NAMELIST
    sed -i 's|INC_STEP_IN|'${INC_STEP_IN}'|g'       $NAMELIST
    sed -i 's|DOMAIN_NUM_IN|'${DOMAIN_NUM_IN}'|g'   $NAMELIST
    sed -i 's|ZCOUNT_IN|'${ZCOUNT_IN}'|g'           $NAMELIST
    sed -i 's|DELT_IN|'${DELT_IN}'|g'               $NAMELIST
    sed -i 's|STIME_IN|'${STIME_IN}'|g'             $NAMELIST

}

edit_namelist_breeding () {

local    NAMELIST=$1
local    INIT_BV=$2

#Note that INIT_BV is passed as an input argumet in function rescale.
#This option should be enambled only at the first breeding cycle.

    sed -i 's|INIT_BV|'${INIT_BV}'|g'                   $NAMELIST
    #We do not perform BREEDING IN PLACE if the flag is set to 0 
    #or if we are initializing the bred vectors.
    if [ $BREEDING_IN_PLACE = 0 -o ${INIT_BV} = ".true." ];then
     sed -i 's|BREEDING_IN_PLACE|.false.|g' $NAMELIST
    else
     sed -i 's|BREEDING_IN_PLACE|.true.|g' $NAMELIST
    fi

}

safe_init_tmpdir () {
#-------------------------------------------------------------------------------
# Safely initialize a temporary directory
#
# Usage: safe_init_tmpdir DIRNAME
#
#   DIRNAME  The temporary directory
#-------------------------------------------------------------------------------

local DIRNAME="$1"


#-------------------------------------------------------------------------------

if [ -z "$DIRNAME" ]; then
  echo "[Warning] $FUNCNAME: '\$DIRNAME' is not set." 
  exit 1
fi

mkdir -p $DIRNAME
res=$? && ((res != 0)) && exit $res

if [ ! -d "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not a directory." 
  exit 1
fi
if [ ! -O "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not owned by you." 
  exit 1
fi

rm -fr $DIRNAME/*
res=$? && ((res != 0)) && exit $res


#-------------------------------------------------------------------------------
}

safe_init_outputdir () {
#-------------------------------------------------------------------------------
# Safely initialize a temporary directory
#
# Usage: safe_init_tmpdir DIRNAME
#
#   DIRNAME  The temporary directory
#-------------------------------------------------------------------------------

local DIRNAME="$1"

#-------------------------------------------------------------------------------

if [ -z "$DIRNAME" ]; then
  echo "[Warning] $FUNCNAME: '\$DIRNAME' is not set." 
  exit 1
fi
if [ -e "$DIRNAME" ]; then
   echo "[Error] $DIRNAME exists: Please remove it manually to avoid data loss"
   exit 1
fi

mkdir -p $DIRNAME
res=$? && ((res != 0)) && exit $res

if [ ! -d "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not a directory." 
  exit 1
fi
if [ ! -O "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not owned by you." 
  exit 1
fi

mkdir -p $DIRNAME/pp_out/
mkdir -p $DIRNAME/init_out/
mkdir -p $DIRNAME/scale_out/
mkdir -p $DIRNAME/logs/



#-------------------------------------------------------------------------------
}

set_tmpdir_breeding () {

local DIRNAME=$1

mkdir -p $DIRNAME/run/

mkdir -p $DIRNAME/rescale/

mkdir -p $DIRNAME/scripts/

mkdir -p $DIRNAME/out/

mkdir -p $DIRNAME/out/latlon/

mkdir -p $DIRNAME/out/logs/

local ibv=1

 while [ $ibv -le $NBV ] ; do

   ibv=`add_zeros $ibv 4`
  
   mkdir -p $DIRNAME/run/${ibv}/pn/

   link_files $DIRNAME/run/${ibv}/pn/

   mkdir -p $DIRNAME/run/${ibv}/pp/

   link_files $DIRNAME/run/${ibv}/pp/

   mkdir -p $DIRNAME/rescale/${ibv}

   ln -sf ${SRCDIR}/../bin/breeding $DIRNAME/rescale/${ibv}

   ln -sf ${SRCDIR}/../run/util.sh  $DIRNAME/scripts/

   cp     ${SRCDIR}/../configuration/${CONFIGURATION}/nml.breeding $DIRNAME/rescale/${ibv}

   ibv=`expr $ibv + 1 `

 done

}

link_files () {

local DIRDEST=$1

#Copy namelist files

cp ${SRCDIR}/../configuration/${CONFIGURATION}/*nml* $DIRDEST

#Link binary files        #######################################

ln -sf ${SRCDIR}/../bin/scale-rm* $DIRDEST

ln -sf ${SRCDIR}/../bin/net2g     $DIRDEST

ln -sf ${DATABASEROOT}/land/param.bucket.conf $DIRDEST

ln -sf ${DATABASEROOT}/rad/PARAG.29 $DIRDEST

ln -sf ${DATABASEROOT}/rad/PARAPC.29 $DIRDEST

ln -sf ${DATABASEROOT}/rad/VARDATA.RM29 $DIRDEST 

ln -sf ${DATABASEROOT}/rad/cira.nc $DIRDEST

ln -sf ${DATABASEROOT}/rad/rad_o3_profs.txt $DIRDEST

ln -sf ${DATABASEROOT}/rad/MIPAS/day.atm $DIRDEST

ln -sf ${DATABASEROOT}/rad/MIPAS/equ.atm $DIRDEST

ln -sf ${DATABASEROOT}/rad/MIPAS/sum.atm $DIRDEST

ln -sf ${DATABASEROOT}/rad/MIPAS/win.atm $DIRDEST


}




set_my_log () {

 local cont=1
 while [ $cont -ge 1  ] ; do
   if [ ! -e $OUTPUTDIR/log${cont}.log ] ; then
     my_log=$OUTPUTDIR/log${cont}.log
     cont=0
   else 
     cont=`expr $cont + 1 `
   fi
 done

}

set_cycle_dates () {

#DEFINE IMORTANT DATES FOR THE CURRENT CYCLE (YYYYMMDDHHNNSS)M
FDATE=`date_edit2 $CDATE $GUESFT `           #Forecast end date
ADATE=`date_edit2 $CDATE $WINDOW `           #Analysis date for the current cycle
WSDATE=`date_edit2 $CDATE $WINDOW_START `    #Assimilation window start date
WEDATE=`date_edit2 $CDATE $WINDOW_END   `    #Assimilation window end   date
#BDYDATE=`date_edit2 $CDATE $DINC   `         #Dummy date for boundary data preparation.


echo ">>> IMPORTANT DATES DEFINED IN THIS CYCLE"

echo ">>> FDATE=   $FDATE "
echo ">>> ADATE=   $ADATE "
echo ">>> WSDATE=  $WSDATE"
echo ">>> WEDATE=  $WEDATE"
#echo ">>> BDYDATE= $BDYDATE"

}


rescale_bv () {

#Prepare initialize bred_vectors_script.

local bvdate=$1

local ibv=1

local maxproc=`expr $NPROCS - 1 `

local maxproc_ct=`expr $NPROCS_CT - 1 `

local init_bv=".false."

local iteration=$2   #In a breeding in place cycle the iteration number that is being performed.

iteration=`add_zeros ${iteration} 4 `

while [ $ibv -le $NBV ] ; do 

  local pp=${ibv}
  local pn=`expr ${ibv} + 1 `

  pp=`add_zeros ${pp}   4`
  pn=`add_zeros ${pn}   4`
  ibv=`add_zeros $ibv   4`
  
  cp ${SRCDIR}/../configuration/${CONFIGURATION}/nml.breeding $TMPDIR/rescale/${ibv}/

  edit_namelist_breeding $TMPDIR/rescale/${ibv}/nml.breeding $init_bv


  local script="$TMPDIR/scripts/rescale_${ibv}.sh"

  echo "#!/bin/bash                                                   " >  $script
  echo "set -x                                                        " >> $script
  echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_ADD:\$LD_LIBRARY_PATH " >> $script
  echo "export OMP_NUM_THREADS=$MY_OMP_NUM_THREADS                    " >> $script
  echo "export OMP_STACKSIZE=$MY_OMP_STACKSIZE                        " >> $script
  echo "cd $TMPDIR/rescale/${ibv}                                     " >> $script

  #Link central trajectory.
  
  local proc=0
  while [ $proc -le $maxproc_ct ] ; do
    proc=`add_zeros $proc 6`

    echo "ln -sf ${BVDATAROOT}/${bvdate}/gues/${BVCENTRALTRAJECTORY}/init.pe${proc}.nc ./ctrl.pe${proc}.nc " >> $script

    proc=`expr $proc + 1 `

  done

  #Link the perturbations.
  local proc=0
  while [ $proc -le $maxproc ] ; do

    proc=`add_zeros $proc 6`

    echo "cp ${TMPDIR}/out/${bvdate}/${ibv}/pp_o.pe${proc}.nc ./pp.pe${proc}.nc " >> $script
    echo "cp ${TMPDIR}/out/${bvdate}/${ibv}/pn_o.pe${proc}.nc ./pn.pe${proc}.nc " >> $script

    if [ $BREEDING_IN_PLACE -eq 1 ] ; then #Copy additional files required by breeding in place.

       bvpdate=`date_edit2 $bvdate -${BVFREQ} ` 
       echo "ln -sf ${TMPDIR}/out/${bvpdate}/${ibv}/pp_r.pe${proc}.nc ./pp_old.pe${proc}.nc " >> $script
       echo "ln -sf ${TMPDIR}/out/${bvpdate}/${ibv}/pn_r.pe${proc}.nc ./pn_old.pe${proc}.nc " >> $script

    fi

    proc=`expr $proc + 1 ` 

  done

  echo "./breeding > $TMPDIR/out/logs/breeding_${bvdate}_${iteration}_${ibv}.log   " >> $script 

  #echo "exit " >> $script #DEBUG

  echo "mv $TMPDIR/rescale/${ibv}/norm.grd ${TMPDIR}/out/${bvdate}/${ibv}/norm_${iteration}.grd " >> $script

  local proc=0
  while [ $proc -le $maxproc ] ; do

    proc=`add_zeros $proc 6 `

    #Copy the rescaled perturbation to output.
    echo "mv ./pp.pe${proc}.nc ${TMPDIR}/out/${bvdate}/${ibv}/pp_r.pe${proc}.nc " >> $script
    echo "mv ./pn.pe${proc}.nc ${TMPDIR}/out/${bvdate}/${ibv}/pn_r.pe${proc}.nc " >> $script

    proc=`expr $proc + 1 `

  done


  echo "#Use: echo qsub  -l nodes=$NODESPERJOB:procs=$NPROCS -Wblock=true " >> $script

  cd $TMPDIR/out/logs/ 
  qsub  -l nodes=$NODESPERJOB:ppn=$NPROCS -Wblock=true $script &
  #bash $script & 

  ibv=`expr $ibv + 1 `

done

  #exit #DEBUG

time wait

}

init_bv () {

#This fuction generates initial perturbations using the boundary data source provided.
#This function uses the breeding module routines to initialize bred vectors. Details about
#random perturbations are specified in nml.breeding namelist.

local bvdate=$1

local ibv=1

local maxproc=`expr $NPROCS - 1 `

local maxproc_ct=`expr $NPROCS_CT - 1 `

#If we are in the first iteration check if the bv will be initialized.
if [ $ITER -eq 1 -a $INIT_BV -eq 1 ] ; then
   local init_bv=".true."
else
   local init_bv=".false."
fi

local  iy=`echo $bvdate | cut -c1-4`
local  im=`echo $bvdate | cut -c5-6`
local  id=`echo $bvdate | cut -c7-8`
local  ih=`echo $bvdate | cut -c9-10`
local  in=`echo $bvdate | cut -c11-12`
local  is=`echo $bvdate | cut -c13-14`


while [ $ibv -le $NBV ] ; do


  local proc=0
  local pn=${ibv}
  local pp=`expr ${ibv} + 1 `

  pp=`add_zeros ${pp}   4`
  pn=`add_zeros ${pn}   4`
  ibv=`add_zeros $ibv   4`

  #Link wrf files.


  if [ $BDYPERT -eq 1 ] ; then
   echo "Perturbed boundary conditions will be used in the initialization of the perturbation"
   echo "PN will use ${pn} member ad bdy, PP will use ${pp} member as bdy"
   link_bdy_data_ens $TMPDIR/run/${ibv}/pp ${pp}
   link_bdy_data_ens $TMPDIR/run/${ibv}/pn ${pn}
  fi
  if [ $BDYPERT -eq 0 ] ; then
   echo "Unperturbed boundary conditions will be used in the initialization of the perturbation"
   link_bdy_data_ens $TMPDIR/run/${ibv}/pp unperturbed
   link_bdy_data_ens $TMPDIR/run/${ibv}/pn unperturbed
  fi

  #Link and edit namelists for scale pp and scale init steps.

  cp ${SRCDIR}/../configuration/${CONFIGURATION}/*nml* $TMPDIR/run/${ibv}/pp
  cp ${SRCDIR}/../configuration/${CONFIGURATION}/*nml* $TMPDIR/run/${ibv}/pn

  edit_namelist $TMPDIR/run/${ibv}/pp/nml.scale_pp  ${bvdate}
  edit_namelist $TMPDIR/run/${ibv}/pn/nml.scale_pp  ${bvdate}

  edit_namelist $TMPDIR/run/${ibv}/pp/nml.scale_init ${bvdate}
  edit_namelist $TMPDIR/run/${ibv}/pn/nml.scale_init ${bvdate}

  cp ${SRCDIR}/../configuration/${CONFIGURATION}/nml.breeding $TMPDIR/rescale/${ibv}/

  edit_namelist_breeding $TMPDIR/rescale/${ibv}/nml.breeding $init_bv

  #Create directory for output
 
  mkdir -p $TMPDIR/out/${bvdate}/${ibv}/

  #Create a simple job script and submit it to the qeue

  local script=$TMPDIR/scripts/init_bv_${ibv}.sh

  echo "#!/bin/bash                                                                                            " >  $script
  echo "set -x                                                                                                 " >> $script
  echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_ADD:\$LD_LIBRARY_PATH                                          " >> $script
  echo "export OMP_NUM_THREADS=$MY_OMP_NUM_THREADS                                                             " >> $script
  echo "export OMP_STACKSIZE=$MY_OMP_STACKSIZE                                                                 " >> $script
  #Run for positive member
  echo "cd $TMPDIR/run/${ibv}/pp                                                                               " >> $script
  echo "mpiexec -np $NPROCS ./scale-rm_pp nml.scale_pp > $TMPDIR/out/logs/scale_pp_${bvdate}_${ibv}pp.log      " >> $script
  echo "mpiexec -np $NPROCS ./scale-rm_init nml.scale_init > $TMPDIR/out/logs/scale_init_${bvdate}_${ibv}pp.log " >> $script
  #Run for negative member
  echo "cd $TMPDIR/run/${ibv}/pn                                                                               " >> $script
  echo "mpiexec -np $NPROCS ./scale-rm_pp nml.scale_pp > $TMPDIR/out/logs/scale_pp_${bvdate}_${ibv}pn.log      " >> $script
  echo "mpiexec -np $NPROCS ./scale-rm_init nml.scale_init > $TMPDIR/out/logs/scale_init_${bvdate}_${ibv}pn.log " >> $script

  echo "cd $TMPDIR/rescale/${ibv}                                                                              " >> $script

  #Move the output

  local proc=0
  while [ $proc -le $maxproc_ct ] ; do
    proc=`add_zeros $proc 6`

    echo "ln -sf ${BVDATAROOT}/${bvdate}/gues/${BVCENTRALTRAJECTORY}/init.pe${proc}.nc $TMPDIR/rescale/${ibv}/ctrl.pe${proc}.nc " >> $script

   proc=`expr $proc + 1 `

  done

  local proc=0
  while [ $proc -le $maxproc ] ; do

    proc=`add_zeros $proc 6`
    local scaledate=`date_in_scale_format $bvdate `
   
    #Copy these perturbations as "_o" which means original or non-rescalled perturbations.

    echo "mv $TMPDIR/run/${ibv}/pp/init_d01_${scaledate}.pe${proc}.nc ./pp.pe${proc}.nc  " >> $script 
    echo "mv $TMPDIR/run/${ibv}/pn/init_d01_${scaledate}.pe${proc}.nc ./pn.pe${proc}.nc  " >> $script 

   proc=`expr $proc + 1 `

  done

  echo "./breeding > $TMPDIR/out/logs/breeding_init_${ibv}.log   " >> $script

  echo "mv $TMPDIR/rescale/${ibv}/norm.grd ${TMPDIR}/out/${bvdate}/${ibv}/norm.grd " >> $script

  local proc=0
  while [ $proc -le $maxproc ] ; do

    proc=`add_zeros $proc 6 `

    #Copy the rescaled perturbation to output.
    echo "mv ./pp.pe${proc}.nc ${TMPDIR}/out/${bvdate}/${ibv}/pp_r.pe${proc}.nc " >> $script
    echo "mv ./pn.pe${proc}.nc ${TMPDIR}/out/${bvdate}/${ibv}/pn_r.pe${proc}.nc " >> $script

    proc=`expr $proc + 1 `

  done

  echo "#Use: echo qsub  -l nodes=$NODESPERJOB:procs=$NPROCS -Wblock=true $script" >> $script

  cd $TMPDIR/out/logs/
  qsub  -l nodes=$NODESPERJOB:ppn=$NPROCS -Wblock=true  $script

  #bash $script &

  ibv=`expr $ibv + 1 `

done

time wait


}


evolve_bv () {


#This fuction generates initial perturbations using the boundary data source provided.
#It also provides land use and topo files.

local bvdate=$1

local iteration=$2  #Breeding in place iteration number.

local ibv=1

local maxproc=`expr $NPROCS - 1 `

local bvdateend=`date_edit2 $bvdate $BVFREQ `

local maxproc=`expr $NPROCS - 1 `

iteration=`add_zeros ${iteration} 4 `


while [ $ibv -le $NBV ] ; do

  local pn=${ibv}
  local pp=`expr ${ibv} + 1 `

  pp=`add_zeros ${pp}   4`
  pn=`add_zeros ${pn}   4`
  ibv=`add_zeros $ibv   4`

  #Link wrf files.

  if [ $BDYPERT -eq 1 ] ; then
   echo "Perturbed boundary conditions will be used in the initialization of the perturbation"
   echo "PN will use ${pn} member ad bdy, PP will use ${pp} member as bdy"
   link_bdy_data_ens $TMPDIR/run/${ibv}/pp ${pp}
   link_bdy_data_ens $TMPDIR/run/${ibv}/pn ${pn}
  fi
  if [ $BDYPERT -eq 0 ] ; then
   echo "Unperturbed boundary conditions will be used in the evolution of the perturbation"
   link_bdy_data_ens $TMPDIR/run/${ibv}/pp unperturbed
   link_bdy_data_ens $TMPDIR/run/${ibv}/pn unperturbed
  fi

  #Link and edit namelists for scale pp and scale init steps.

  cp ${SRCDIR}/../configuration/${CONFIGURATION}/*nml* $TMPDIR/run/${ibv}/pp
  cp ${SRCDIR}/../configuration/${CONFIGURATION}/*nml* $TMPDIR/run/${ibv}/pn

  edit_namelist $TMPDIR/run/${ibv}/pp/nml.scale_pp ${bvdate}
  edit_namelist $TMPDIR/run/${ibv}/pn/nml.scale_pp ${bvdate}

  edit_namelist $TMPDIR/run/${ibv}/pp/nml.scale_init ${bvdate}
  edit_namelist $TMPDIR/run/${ibv}/pn/nml.scale_init ${bvdate}
 
  edit_namelist $TMPDIR/run/${ibv}/pp/nml.scale ${bvdate}
  edit_namelist $TMPDIR/run/${ibv}/pn/nml.scale ${bvdate}

  edit_namelist_postproc $TMPDIR/run/${ibv}/pp/nml.net2g ${bvdate}
  edit_namelist_postproc $TMPDIR/run/${ibv}/pn/nml.net2g ${bvdate}

  edit_namelist_postproc $TMPDIR/run/${ibv}/pp/nml.net2g_latlon ${bvdate}
  edit_namelist_postproc $TMPDIR/run/${ibv}/pn/nml.net2g_latlon ${bvdate}


  #Create directory for output

  mkdir -p $TMPDIR/out/${bvdateend}/${ibv}/grads_pp_o${iteration} #Evolved bred vector grads output.
  mkdir -p $TMPDIR/out/${bvdateend}/${ibv}/grads_pn_o${iteration} #Evolved bred vector grads output.

  mkdir -p $TMPDIR/out/${bvdate}/${ibv}/grads_pp_r${iteration}    #Rescaled bred vector grads output.
  mkdir -p $TMPDIR/out/${bvdate}/${ibv}/grads_pn_r${iteration}    #Rescaled bred vector grads output.

  if [ $SAVE_HISTORY -eq "1" ];then

   mkdir -p $TMPDIR/out/${bvdate}/${ibv}/history_pn_${iteration}   #History files for the negative perturbation.
   mkdir -p $TMPDIR/out/${bvdate}/${ibv}/history_pp_${iteration}   #History files for the positive perturbation.

  fi

  mkdir -p $TMPDIR/out/${bvdateend}/${ibv}/       #For evolved netcdf restart files.

  #Create a simple job script and submit it to the qeue

  local script="$TMPDIR/scripts/evolve_bv_${ibv}.sh"
  local scaledate=`date_in_scale_format $bvdate `
  local scaledateend=`date_in_scale_format $bvdateend `

  echo "#!/bin/bash                                                                                            " >  $script
  echo "set -x                                                                                                 " >> $script
  echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_ADD:\$LD_LIBRARY_PATH                                          " >> $script
  echo "export OMP_NUM_THREADS=$MY_OMP_NUM_THREADS                                                             " >> $script
  echo "export OMP_STACKSIZE=$MY_OMP_STACKSIZE                                                                 " >> $script
  #Run for positive member
  echo "cd $TMPDIR/run/${ibv}/pp                                                                               " >> $script
  echo "mpiexec -np $NPROCS ./scale-rm_pp nml.scale_pp > $TMPDIR/out/logs/scale_pp_${bvdate}_${ibv}pp.log      " >> $script
  echo "mpiexec -np $NPROCS ./scale-rm_init nml.scale_init > $TMPDIR/out/logs/scale_init_${bvdate}_${ibv}pp.log " >> $script
  #Replace the initial conditions from the outer model with the rescaled bred vector from the previous cycle.
  local proc=0
  while [ $proc -le $maxproc ] ; do
    proc=`add_zeros $proc 6`
    echo "rm $TMPDIR/run/${ibv}/pp/init_d01_${scaledate}.pe${proc}.nc                                                 " >> $script
    echo "cp $TMPDIR/out/${bvdate}/${ibv}/pp_r.pe${proc}.nc $TMPDIR/run/${ibv}/pp/init_d01_${scaledate}.pe${proc}.nc  " >> $script
    proc=`expr $proc + 1 `
  done
  echo "mpiexec -np $NPROCS ./scale-rm      nml.scale      >  $TMPDIR/out/logs/scale_${bvdate}_${ibv}pp.log    " >> $script
  echo "mpiexec -np $NPROCS ./net2g         nml.net2g      >  $TMPDIR/out/logs/net2g_${bvdate}_${ibv}pp.log    " >> $script
  if [ $iteration -eq "0001" ] ; then
    echo "mpiexec -np $NPROCS ./net2g  nml.net2g_latlon  >  $TMPDIR/out/logs/net2g_latlon_${bvdate}_${ibv}pn.log " >> $script
  fi
  #Run for negative member
  echo "cd $TMPDIR/run/${ibv}/pn                                                                               " >> $script
  echo "mpiexec -np $NPROCS ./scale-rm_pp nml.scale_pp > $TMPDIR/out/logs/scale_pp_${bvdate}_${ibv}pn.log      " >> $script
  echo "mpiexec -np $NPROCS ./scale-rm_init nml.scale_init > $TMPDIR/out/logs/scale_init_${bvdate}_${ibv}pn.log " >> $script
  #Replace the initial conditions from the outer model with the rescaled bred vector from the previous cycle.
  local proc=0
  while [ $proc -le $maxproc ] ; do
    proc=`add_zeros $proc 6`
    echo "rm $TMPDIR/run/${ibv}/pn/init_d01_${scaledate}.pe${proc}.nc                                                 " >> $script
    echo "cp $TMPDIR/out/${bvdate}/${ibv}/pn_r.pe${proc}.nc $TMPDIR/run/${ibv}/pn/init_d01_${scaledate}.pe${proc}.nc  " >> $script
    proc=`expr $proc + 1 `
  done
  echo "mpiexec -np $NPROCS ./scale-rm      nml.scale      >  $TMPDIR/out/logs/scale_${bvdate}_${ibv}pn.log        " >> $script
  echo "mpiexec -np $NPROCS ./net2g         nml.net2g      >  $TMPDIR/out/logs/net2g_${bvdate}_${ibv}pn.log        " >> $script
  if [ $iteration -eq "0001" -a $ITER -eq 1 ] ; then
    echo "mpiexec -np $NPROCS ./net2g  nml.net2g_latlon  >  $TMPDIR/out/logs/net2g_latlon_${bvdate}_${ibv}pn.log     " >> $script 
    echo "mv $TMPDIR/run/${ibv}/pn/lon*.grd $TMPDIR/out/latlon/                                                      " >> $script
    echo "mv $TMPDIR/run/${ibv}/pn/lat*.grd $TMPDIR/out/latlon/                                                      " >> $script
  fi

  #Move the output
  local proc=0
  while [ $proc -le $maxproc ] ; do

    proc=`add_zeros $proc 6`

    #Copy the evolved perturbations to the output directory.

    echo "mv $TMPDIR/run/${ibv}/pp/restart_d01_${scaledateend}.pe${proc}.nc $TMPDIR/out/${bvdateend}/${ibv}/pp_o.pe${proc}.nc  " >> $script
    echo "mv $TMPDIR/run/${ibv}/pn/restart_d01_${scaledateend}.pe${proc}.nc $TMPDIR/out/${bvdateend}/${ibv}/pn_o.pe${proc}.nc  " >> $script


    if [ $SAVE_HISTORY -eq "1" ] ; then
       #Save full history files for each iteration and also save the configuration files.
       echo "mv $TMPDIR/run/${ibv}/pp/history.pe${proc}.nc $TMPDIR/out/${bvdate}/${ibv}/history_pp_${iteration}/history.pe${proc}.nc " >> $script
       echo "mv $TMPDIR/run/${ibv}/pn/history.pe${proc}.nc $TMPDIR/out/${bvdate}/${ibv}/history_pn_${iteration}/history.pe${proc}.nc " >> $script
       echo "mv $TMPDIR/run/${ibv}/pp/nml.*                $TMPDIR/out/${bvdate}/${ibv}/history_pp_${iteration}/                     " >> $script
       echo "mv $TMPDIR/run/${ibv}/pn/nml.*                $TMPDIR/out/${bvdate}/${ibv}/history_pn_${iteration}/                     " >> $script
    fi

    proc=`expr $proc + 1 `
 
  done

  echo "mv $TMPDIR/run/${ibv}/pp/*${bvdateend}*.grd $TMPDIR/out/${bvdateend}/${ibv}/grads_pp_o${iteration}/      " >> $script
  echo "mv $TMPDIR/run/${ibv}/pn/*${bvdateend}*.grd $TMPDIR/out/${bvdateend}/${ibv}/grads_pn_o${iteration}/      " >> $script

  echo "mv $TMPDIR/run/${ibv}/pp/*${bvdateend}*.ctl $TMPDIR/out/${bvdateend}/${ibv}/grads_pp_o${iteration}/      " >> $script
  echo "mv $TMPDIR/run/${ibv}/pn/*${bvdateend}*.ctl $TMPDIR/out/${bvdateend}/${ibv}/grads_pn_o${iteration}/      " >> $script

  #Copy the rescaled perturbation in grads format.

  echo "mv $TMPDIR/run/${ibv}/pp/*${bvdate}*.grd $TMPDIR/out/${bvdate}/${ibv}/grads_pp_r${iteration}/     " >> $script
  echo "mv $TMPDIR/run/${ibv}/pn/*${bvdate}*.grd $TMPDIR/out/${bvdate}/${ibv}/grads_pn_r${iteration}/     " >> $script
  
  echo "mv $TMPDIR/run/${ibv}/pp/*${bvdate}*.ctl $TMPDIR/out/${bvdate}/${ibv}/grads_pp_r${iteration}/     " >> $script
  echo "mv $TMPDIR/run/${ibv}/pn/*${bvdate}*.ctl $TMPDIR/out/${bvdate}/${ibv}/grads_pn_r${iteration}/     " >> $script
  

  echo "#Use: echo qsub  -l nodes=$NODESPERJOB:procs=$NPROCS -Wblock=true " >> $script

  cd $TMPDIR/out/logs/
  qsub  -l nodes=$NODESPERJOB:ppn=$NPROCS -Wblock=true $script &
  #bash $script & 

  ibv=`expr $ibv + 1 `

done

time wait


}

evolve_rescale_bv () {

local bvdate=$1

echo "Las iteraciones de breeding in place son $BREEDING_IN_PLACE_IT"

local bvdateend=`date_edit2 $bvdate $BVFREQ `

if [ $bvdate -le $BREEDING_IN_PLACE_ENDDATE -a $bvdate -ge $BREEDING_IN_PLACE_INIDATE ] ; then

FINAL_ITER=$BREEDING_IN_PLACE_IT

else

FINAL_ITER=1

fi

echo "We will perform $FINAL_ITER iteratons in this cycle"

echo "Running cycle starting at $bvdate "

local my_iter=1
while [ $my_iter -le $FINAL_ITER ] ; do

 echo "Evolve bred vector for iter = $my_iter "
 evolve_bv $CDATE $my_iter

 echo "Rescale bred vector for iter = $my_iter "
 rescale_bv $FDATE $my_iter

my_iter=`expr $my_iter + 1`
done


}

link_bdy_data ()  {
#Links wrf bdy conditions
local DESTDIR=$1

  local ID=$INIDATE
  local ED=`date_edit2 $ID $FLENGTH`

  local CDATE=${ID}
  local FILECOUNT=0

  while [ $CDATE -le $ED ] ; do
 
   FILECOUNT=`add_zeros $FILECOUNT 5`

   ln -sf ${BDYDATAROOT}/wrfout_${CDATE} $DESTDIR/wrfout_d01_${FILECOUNT}

   FILECOUNT=`expr $FILECOUNT + 1 `
   
   CDATE=`date_edit2 $CDATE $BDYDATAFREQ`
   
 done 

 BDYNFILES=${FILECOUNT}

}

link_bdy_data_ens ()  {
#Links wrf bdy conditions
local DESTDIR=$1
local MEM=$2

  local ID=$INIDATE
  local ED=`date_edit2 $ID $FLENGTH`

  local CDATE=${ID}
  local FILECOUNT=0

  while [ $CDATE -le $ED ] ; do

   FILECOUNT=`add_zeros $FILECOUNT 5`

   ln -sf ${BDYDATAROOT}/${MEM}/wrfout_${CDATE} $DESTDIR/wrfout_d01_${FILECOUNT}

   FILECOUNT=`expr $FILECOUNT + 1 `

   CDATE=`date_edit2 $CDATE $BDYDATAFREQ`

 done

 BDYNFILES=${FILECOUNT}

}










