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
(
        
        if [ $# -lt 2 ]; then
                echo "Usage : date_floor"
                echo "    date_diff [yyyy][mm][dd][hh][mn] [interval(seconds)]"
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
)
}

ens_member () {

local    MEMBER="$1"
local    MEMBER_STR=$MEMBER

    if test $MEMBER -lt 10000
    then
      MEMBER_STR=0$MEMBER_STR
    fi
    if test $MEMBER -lt 1000
    then
      MEMBER_STR=0$MEMBER_STR
    fi
    if test $MEMBER -lt 100
    then
      MEMBER_STR=0$MEMBER_STR
    fi
    if test $MEMBER -lt 10
    then
      MEMBER_STR=0$MEMBER_STR
    fi

    echo $MEMBER_STR
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

edit_namelist_input () {

local    NAMELIST=$1
local    IDATE=$2
local    EDATE=$3
local    OUTPUT_FREQ=$4                      #In seconds
local    BOUNDARY_DATA_FREQ=$5               #In seconds

    iy=`echo $IDATE | cut -c1-4`
    im=`echo $IDATE | cut -c5-6`
    id=`echo $IDATE | cut -c7-8`
    ih=`echo $IDATE | cut -c9-10`
    in=`echo $IDATE | cut -c11-12`
    is=`echo $IDATE | cut -c13-14`

    ey=`echo $EDATE | cut -c1-4`
    em=`echo $EDATE | cut -c5-6`
    ed=`echo $EDATE | cut -c7-8`
    eh=`echo $EDATE | cut -c9-10`
    en=`echo $EDATE | cut -c11-12`
    es=`echo $EDATE | cut -c13-14`

    sed -i 's/SYYYY/'${iy}'/g'                   $NAMELIST
    sed -i 's/SMM/'${im}'/g'                     $NAMELIST
    sed -i 's/SDD/'${id}'/g'                     $NAMELIST
    sed -i 's/SHH/'${ih}'/g'                     $NAMELIST
    sed -i 's/SMN/'${in}'/g'                     $NAMELIST
    sed -i 's/SSS/'${is}'/g'                     $NAMELIST
    sed -i 's/EYYYY/'${ey}'/g'                   $NAMELIST
    sed -i 's/EMM/'${em}'/g'                     $NAMELIST
    sed -i 's/EDD/'${ed}'/g'                     $NAMELIST
    sed -i 's/EHH/'${eh}'/g'                     $NAMELIST
    sed -i 's/EMN/'${en}'/g'                     $NAMELIST
    sed -i 's/ESS/'${es}'/g'                     $NAMELIST
    sed -i 's/OUTPUTFREQ/'${OUTPUT_FREQ}'/g'     $NAMELIST
    sed -i 's/BDYFREQ/'${BOUNDARY_DATA_FREQ}'/g' $NAMELIST

}

edit_namelist_wps () {

local    NAMELIST=$1
local    IDATE=$2
local    EDATE=$3
local    BOUNDARY_DATA_FREQ=$4         #In seconds

    iy=`echo $IDATE | cut -c1-4`
    im=`echo $IDATE | cut -c5-6`
    id=`echo $IDATE | cut -c7-8`
    ih=`echo $IDATE | cut -c9-10`
    in=`echo $IDATE | cut -c11-12`
    is=`echo $IDATE | cut -c13-14`

    ey=`echo $EDATE | cut -c1-4`
    em=`echo $EDATE | cut -c5-6`
    ed=`echo $EDATE | cut -c7-8`
    eh=`echo $EDATE | cut -c9-10`
    en=`echo $EDATE | cut -c11-12`
    es=`echo $EDATE | cut -c13-14`

    sed -i 's/SYYYY/'${iy}'/g'                   $NAMELIST
    sed -i 's/SMM/'${im}'/g'                     $NAMELIST
    sed -i 's/SDD/'${id}'/g'                     $NAMELIST
    sed -i 's/SHH/'${ih}'/g'                     $NAMELIST
    sed -i 's/SMN/'${in}'/g'                     $NAMELIST
    sed -i 's/SSS/'${is}'/g'                     $NAMELIST
    sed -i 's/EYYYY/'${ey}'/g'                   $NAMELIST
    sed -i 's/EMM/'${em}'/g'                     $NAMELIST
    sed -i 's/EDD/'${ed}'/g'                     $NAMELIST
    sed -i 's/EHH/'${eh}'/g'                     $NAMELIST
    sed -i 's/EMN/'${en}'/g'                     $NAMELIST
    sed -i 's/ESS/'${es}'/g'                     $NAMELIST
    sed -i 's/BDYFREQ/'${BOUNDARY_DATA_FREQ}'/g' $NAMELIST

}

edit_namelist_letkf () {

        NAMELIST=$1


    sed -i 's/@NBV@/'${MEMBER}'/g'                                  $NAMELIST
    sed -i 's/@NSLOTS@/'${NSLOTS}'/g'                               $NAMELIST
    sed -i 's/@NBSLOT@/'${NBSLOT}'/g'                               $NAMELIST
    sed -i 's/@SIGMA_OBS@/'${SIGMA_OBS}'/g'                         $NAMELIST
    sed -i 's/@SIGMA_OBSV@/'${SIGMA_OBSV}'/g'                       $NAMELIST
    sed -i 's/@SIGMA_OBSZ@/'${SIGMA_OBSZ}'/g'                       $NAMELIST
    sed -i 's/@SIGMA_OBST@/'${SIGMA_OBST}'/g'                       $NAMELIST
    sed -i 's/@GROSS_ERROR@/'${GROSS_ERROR}'/g'                     $NAMELIST
    sed -i 's/@COV_INFL_MUL@/'${COV_INFL_MUL}'/g'                   $NAMELIST
    sed -i 's/@SP_INFL_ADD@/'${SP_INFL_ADD}'/g'                     $NAMELIST
    sed -i 's/@RELAX_ALPHA_SPREAD@/'${RELAX_ALPHA_SPREAD}'/g'       $NAMELIST
    sed -i 's/@RELAX_ALPHA@/'${RELAX_ALPHA}'/g'                     $NAMELIST

}

edit_namelist_arwpost () {

NAMELIST=$1
local IDATE=$2
local EDATE=$3
local LOCALOUTPUTFREQ=$4 #In seconds


LOCALDATE=$IDATE
while [ $LOCALDATE -le $EDATE ] ; do

    ARWEDATE=$LOCALDATE
    LOCALDATE=` date_edit2 $LOCALDATE $LOCALOUTPUTFREQ `

done

    iy=`echo $IDATE | cut -c1-4`
    im=`echo $IDATE | cut -c5-6`
    id=`echo $IDATE | cut -c7-8`
    ih=`echo $IDATE | cut -c9-10`
    in=`echo $IDATE | cut -c11-12`
    is=`echo $IDATE | cut -c13-14`

    ey=`echo $ARWEDATE | cut -c1-4`
    em=`echo $ARWEDATE | cut -c5-6`
    ed=`echo $ARWEDATE | cut -c7-8`
    eh=`echo $ARWEDATE | cut -c9-10`
    en=`echo $ARWEDATE | cut -c11-12`
    es=`echo $ARWEDATE | cut -c13-14`

    sed -i 's/SYYYY/'${iy}'/g'                   $NAMELIST
    sed -i 's/SMM/'${im}'/g'                     $NAMELIST
    sed -i 's/SDD/'${id}'/g'                     $NAMELIST
    sed -i 's/SHH/'${ih}'/g'                     $NAMELIST
    sed -i 's/SMN/'${in}'/g'                     $NAMELIST
    sed -i 's/SSS/'${is}'/g'                     $NAMELIST
    sed -i 's/EYYYY/'${ey}'/g'                   $NAMELIST
    sed -i 's/EMM/'${em}'/g'                     $NAMELIST
    sed -i 's/EDD/'${ed}'/g'                     $NAMELIST
    sed -i 's/EHH/'${eh}'/g'                     $NAMELIST
    sed -i 's/EMN/'${en}'/g'                     $NAMELIST
    sed -i 's/ESS/'${es}'/g'                     $NAMELIST
    sed -i 's/OUTPUTFREQ/'${LOCALOUTPUTFREQ}'/g'           $NAMELIST
    sed -i 's/INPUT_ROOT_NAME/'${INPUT_ROOT_NAME}'/g'      $NAMELIST
    sed -i 's/OUTPUT_ROOT_NAME/'${OUTPUT_ROOT_NAME}'/g'    $NAMELIST    
    sed -i 's/OUTPUTVARS/'${OUTVARS}'/g'                   $NAMELIST  
    sed -i 's/OUTPUTLEVS/'${OUTLEVS}'/g'                   $NAMELIST 
    sed -i 's/INTERP_METHOD/'${INTERP_METHOD}'/g'          $NAMELIST 
    
}

edit_wrf_wrf () {

local SCRIPT=$1
echo "#!/bin/bash                                           "   > $SCRIPT
echo "set -x                                                "  >> $SCRIPT
echo "WORKDIR=\$1                                           "  >> $SCRIPT
echo "                                                      "  >> $SCRIPT
echo "cd \$WORKDIR                                          "  >> $SCRIPT
if [ $SYSTEM -eq 1 ] ; then 
echo "ulimit -s unlimited                                   "  >> $SCRIPT
fi                     
echo "./wrf.exe                                             "  >> $SCRIPT

chmod 766 $SCRIPT

}

edit_wrf_real () {

local SCRIPT=$1
echo "#!/bin/bash                                           "   > $SCRIPT
echo "set -x                                                "  >> $SCRIPT
echo "WORKDIR=\$1                                           "  >> $SCRIPT
echo "                                                      "  >> $SCRIPT
echo "cd \$WORKDIR                                          "  >> $SCRIPT
if [ $SYSTEM -eq 1 ] ; then 
echo "ulimit -s unlimited                                   "  >> $SCRIPT
fi                     
echo "./real.exe                                            " >> $SCRIPT

chmod 766 $SCRIPT

}

edit_wrf_interpana () {

#Run real using LETKF anal met_em as input to generate wrinput in the forecast grid 
#(only if the LETKF analysis and forecast grids are different)

#Wrfinput and wrfbdy that are previously generated from perturbed gfs outpust are preserved
#as wrfinput_d01.tmp and wrfbdy_d01.tmp and restored at the end.
#The wrfinput_d01 generated from the LETKF analysis file is stored as anal (as in the LETKF cycle).

local SCRIPT=$1
if [ $FORECAST -eq 1 -a $INTERPANA -eq 1 ] ; then  #Forecast and analysis have different grids.
 echo "#!/bin/bash                                           "  > $SCRIPT
 echo "set -x                                                " >> $SCRIPT
 echo "WORKDIR=\$1                                           " >> $SCRIPT
 echo "                                                      " >> $SCRIPT
 echo "cd \$WORKDIR                                          " >> $SCRIPT
 if [ $SYSTEM -eq 1 ] ; then
  echo "ulimit -s unlimited                                  " >> $SCRIPT
 fi
 echo "mv wrfbdy_d01   wrfbdy_d01.tmp                        " >> $SCRIPT
 echo "mv wrfinput_d01 wrfinput_d01.tmp                      " >> $SCRIPT
 echo "cp ../WRF/namelist.input.real2 ./namelist.input       " >> $SCRIPT
 met_em_file=met_em_file_name $CDATE 01
 echo "mv ${met_em_file}.anal ${met_em_file}                 " >> $SCRIPT
 echo "./real.exe                                            " >> $SCRIPT
 echo "cp wrfinput_d01.gfs ./anal                            " >> $SCRIPT
 echo "ln -sf ./wrfinput_d01 ./input1.nc                     " >> $SCRIPT
 echo "ln -sf ./anal         ./input2.nc                     " >> $SCRIPT
 echo "../WRF/merge_wrfinput.exe  > ./merge_wrfinput.log     " >> $SCRIPT
 echo "mv wrfbdy_d01.tmp wrfbdy_d01                          " >> $SCRIPT
 echo "mv wrfinput_d01.tmp wrfinput_d01                      " >> $SCRIPT

 chmod 766 $SCRIPT
fi

}

edit_wrf_pre () {

local SCRIPT=$1
echo "#!/bin/bash                                                         "  > $SCRIPT
echo "set -x                                                              " >> $SCRIPT
echo "WORKDIR=\$1                                                         " >> $SCRIPT
echo "MEM=\$2                                                             " >> $SCRIPT
echo "echo \$WORKDIR                                                      " >> $SCRIPT
echo "cd \$WORKDIR                                                        " >> $SCRIPT
if [ $SYSTEM -eq 0 ] ; then
echo "../WRF/dummy-mpi                                                    " >> $SCRIPT
fi
#MERGE AND UPDATE LATERAL AND LOW BOUNDARY CONDITIONS           
echo "mv wrfinput_d01 wrfinput_d01.gfs                                    " >> $SCRIPT
echo "mv anal wrfinput_d01                                                " >> $SCRIPT
echo "echo \"&control_param                           \" > parame.in      " >> $SCRIPT
echo "echo \"da_file='./wrfinput_d01'                 \" >> parame.in     " >> $SCRIPT
echo "echo \"wrf_bdy_file='./wrfbdy_d01'              \" >> parame.in     " >> $SCRIPT
echo "echo \"wrf_input='./wrfinput_d01.gfs'           \" >> parame.in     " >> $SCRIPT
echo "echo \"debug=.true.                             \" >> parame.in     " >> $SCRIPT
echo "echo \"update_lateral_bdy=.true.                \" >> parame.in     " >> $SCRIPT
echo "echo \"update_low_bdy=.true.                    \" >> parame.in     " >> $SCRIPT
echo "echo \"update_lsm=.true.                        \" >> parame.in     " >> $SCRIPT
echo "echo \"iswater=16                               \" >> parame.in     " >> $SCRIPT
echo "echo \"/                                        \" >> parame.in     " >> $SCRIPT
echo "                                                                    " >> $SCRIPT
echo "./da_update_bc.exe > daupdatebc\${MEM}.log                          " >> $SCRIPT

chmod 766 $SCRIPT

}

edit_wrf_post () {

SCRIPT=$1

echo "#!/bin/bash                                                  " > $SCRIPT 
echo "WORKDIR=\$1                                                  " >> $SCRIPT
echo "MEM=\$2                                                      " >> $SCRIPT
echo "cd \$WORKDIR                                                 " >> $SCRIPT

#Insert DUMMY MPI CALL
if [ $SYSTEM -eq 0 ] ; then
echo "../WRF/dummy-mpi                                             " >> $SCRIPT
fi
if [ $SYSTEM -eq 1 ]; then
echo "ulimit -s unlimited                                          " >> $SCRIPT
fi

echo "cat ./rsl.error.* > ./wrf\${MEM}.log                         " >> $SCRIPT
# --- RENAME OUTPUT FOR ANALYSIS
if [ $ANALYSIS -eq 1 ] ; then

  local CDATEL=$WSDATE

  local LOCAL_OUTFREC=$WINDOW_FREC

  local it=1

  while [ ${CDATEL} -le ${WEDATE} ] ; do

  local itm=$it
  if [ ${it} -lt 10 ]
  then
  itm="0${itm}"
  fi


  local_file=` wrfout_file_name $CDATEL $DOMAIN`
  if [ $it -eq $NBSLOT ] ; then
   
  echo "cp $local_file ../LETKF/gues\${MEM}                          " >> $SCRIPT
  fi

  echo "mv $local_file ../LETKF/gs${itm}\${MEM}                      " >> $SCRIPT

  CDATEL=`date_edit2 $CDATEL $LOCAL_OUTFREC `
  it=`expr ${it} + 1`
  done

  local MMS=` ens_member $MEANMEMBER `

  echo "if [  \$MEM -eq $MM ] ; then " >> $SCRIPT
  echo "cp ../LETKF/gs${NBSLOT}\${MEM} ../LETKF/gues${MMS}             " >> $SCRIPT
  echo "cp ../LETKF/gs${NBSLOT}\${MEM} ../LETKF/gs${NBSLOT}${MMS}      " >> $SCRIPT
  echo "fi" >> $SCRIPT

fi

# --- RENAME OUTPUT FOR FORECAST
if [ $FORECAST -eq 1 ] ; then

  local CDATEL=$DATE

  local LOCAL_OUTFREC=$WINDOW_FREC

  local it=1

  while [ ${CDATEL} -le ${FDATE} ] ; do

  local itm=$it
  if [ ${it} -lt 10 ]
  then
  itm="0${itm}"
  fi


  local_file=` wrfout_file_name $CDATEL $DOMAIN`

  echo "mv $local_file ../LETKF/gs${itm}\${MEM}                      " >> $SCRIPT

  CDATEL=`date_edit2 $CDATEL $LOCAL_OUTFREC `
  it=`expr ${it} + 1`
  done

fi


chmod 766 $SCRIPT


}

safe_rm_tmpdir () {
#-------------------------------------------------------------------------------
# Safely remove a temporary directory
#
# Usage: safe_rm_tmpdir DIRNAME
#
#   DIRNAME  The temporary directory
#-------------------------------------------------------------------------------

local DIRNAME="$1"

#-------------------------------------------------------------------------------

if [ -z "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '\$DIRNAME' is not set." >&2
  exit 1
fi
if [ ! -e "$DIRNAME" ]; then
  return 0
fi
if [ ! -d "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not a directory." >&2
  exit 1
fi
if [ ! -O "$DIRNAME" ]; then
  echo "[Error] $FUNCNAME: '$DIRNAME' is not owned by you." >&2
  exit 1
fi

rm -fr $DIRNAME
res=$? && ((res != 0)) && exit $res

#-------------------------------------------------------------------------------
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

mkdir -p $DIRNAME/LETKF
mkdir -p $DIRNAME/SCRIPTS
if [ $SYSTEM -eq 0 ] ; then
mkdir -p $DIRNAME/SPAWN
fi
mkdir -p $DIRNAME/WRF
mkdir -p $DIRNAME/add_pert
mkdir -p $DIRNAME/OBS
#Aditional executables for forecast jobs.
if [ $FORECAST -eq 1 ] ; then
mkdir -p $DIRNAME/WPS
mkdir -p $DIRNAME/wrf_to_wps
fi
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
if [ -e "$DIRNAME" -a $RESTART -eq 0 ]; then
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

#rm -fr $DIRNAME/*
#res=$? && ((res != 0)) && exit $res

# SET OUTPUT DIRECTORY
if [ $RESTART -eq 0  ] ; then

 echo "This is a new experiment -> Building output directory from scracth"

 mkdir -p $DIRNAME/gues
 mkdir -p $DIRNAME/anal
 mkdir -p $DIRNAME/configuration 

 if [ $FORECAST -eq 1 ] ; then
   mkdir -p $DIRNAME/forecast/
 fi

else

 if [ ! -e $DIRNAME/gues/ -o ! -e $DIRNAME/anal/ ] ; then
  echo "[Error] This is a restart run but OUTPUTDIR=$DIRNAME does not exist "
  exit 1
 fi
 echo "[Warning] This is a restart experiment -> Using the previous output directory (data can be partially overwritten) "

fi


#-------------------------------------------------------------------------------
}

generate_vcode () {
#Generate vcode files.

local TMPDIRL="$1"               #Temporary directory

NODE=0
JOB=1

if [ $RUN_ONLY_MEAN -ne 1 ] ; then

 while [ $JOB -le $MM -a $JOB -le $MAX_BACKGROUND_JOBS ] ; do

 MNODE=1
  vcoord_file=$TMPDIRL/vcoord_${JOB}
  rm -fr $vcoord_file
  while [ $MNODE -le $NODES_PER_MEMBER ]; do
   echo "( $NODE ) " >> $vcoord_file
   MNODE=`expr $MNODE + 1 `
   NODE=`expr $NODE + 1 `
  done
  JOB=`expr $JOB + 1 `

 done

else

  JOB=1
  MNODE=1
  vcoord_file=$TMPDIRL/vcoord_${JOB}
  rm -fr $vcoord_file
  while [ $MNODE -le $NODES_PER_MEMBER ]; do
   echo "( $NODE ) " >> $vcoord_file
   MNODE=`expr $MNODE + 1 `
   NODE=`expr $NODE + 1 `
  done

fi

}

copy_data () {

#COPY LETKF

cp $LETKF $TMPDIR/LETKF/letkf.exe
cp $NAMELISTLETKF $TMPDIR/LETKF/letkf.namelist.template

#COPY WRF
cp $WRFMODEL/run/*         $TMPDIR/WRF/
rm -f $TMPDIR/WRF/*.exe $TMPDIR/WRF/namelist.input
cp $WRFMODEL/main/wrf.exe  $TMPDIR/WRF/
cp $WRFMODEL/main/real.exe $TMPDIR/WRF/
cp $ARWPOST/ARWpost.exe    $TMPDIR/WRF/
cp -r $ARWPOST/src         $TMPDIR/WRF/
cp $UPDATEBC               $TMPDIR/WRF/da_update_bc.exe
cp $NAMELISTWRF            $TMPDIR/WRF/namelist.input.template

if [ $SYSTEM -eq 0 ];then
#COPY SPAWN
cp $SPAWN/dummy-mpi $TMPDIR/SPAWN
cp $SPAWN/spawn     $TMPDIR/SPAWN
fi

#COPY AND COMPILE PERTURBATION COMPUTATION CODE
cp $WRF/add_pert/*  $TMPDIR/add_pert
#ssh $PPSSERVER " cd $TMPDIR/add_pert && ./make_compute_pert_metem.sh >  $TMPDIR/add_pert/compile.log "

#COPY BASH SCRIPTS
cp $UTIL        $TMPDIR/SCRIPTS/util.sh
chmod 766 $TMPDIR/SCRIPTS/*.sh

if [ $FORECAST -eq 1 -a $INTERPANA -eq 1 ] ; then
#COPY WPS
cp $NAMELISTWPS   $TMPDIR/WPS/namelist.wps.template
cp $WPS/* $TMPDIR/WPS/
cp $WRF_TO_WPS/wrf_to_wps.exe     $TMPDIR/wrf_to_wps/
cp $WRF_TO_WPS/merge_wrfinput.exe $TMPDIR/WRF/
fi

}

set_cycle_dates () {

#DEFINE IMORTANT DATES FOR THE CURRENT CYCLE (YYYYMMDDHHNNSS)
FDATE=`date_edit2 $CDATE $GUESFT `           #Forecast end date
ADATE=`date_edit2 $CDATE $WINDOW `           #Analysis date for the current cycle
WSDATE=`date_edit2 $CDATE $WINDOW_START `    #Assimilation window start date
WEDATE=`date_edit2 $CDATE $WINDOW_END   `    #Assimilation window end   date
BDYDATE=`date_edit2 $CDATE $DINC   `         #Dummy date for boundary data preparation.


echo ">>> IMPORTANT DATES DEFINED IN THIS CYCLE"

echo ">>> FDATE=   $FDATE "
echo ">>> ADATE=   $ADATE "
echo ">>> WSDATE=  $WSDATE"
echo ">>> WEDATE=  $WEDATE"
echo ">>> BDYDATE= $BDYDATE"

}

set_total_nodes () {

if [ $RUN_ONLY_MEAN -ne 1  ] ; then

 if [ $MM -le $MAX_BACKGROUND_JOBS ] ; then
  TOTAL_NODES=`expr $NODES_PER_MEMBER \* $MM `                     #Total number of nodes to be used.
 else
  TOTAL_NODES=`expr $NODES_PER_MEMBER \* $MAX_BACKGROUND_JOBS `
 fi

else

 TOTAL_NODES=$NODES_PER_MEMBER

fi

if [ $SYSTEM -eq  0 ] ; then

   if [ $TOTAL_NODES -ge 3841 ] ; then
      QEUE="large"
   fi
   if [ $TOTAL_NODES -lt 3841 ] ; then
      QEUE="small" 
   fi

fi



PROC_PER_MEMBER=`expr $NODES_PER_MEMBER \* $PROC_PER_NODE `         #Total number of procs per ensemble member.
PROC_FOR_LETKF=$TOTAL_NODES                                         #Total number of nodes for LETKF
TOTAL_PROC=`expr $TOTAL_NODES \* $PROC_PER_NODE `                   #Total number of procs to be used.

echo ">>> NODES THAT WILL BE USED IN THIS EXPERIMENT"
echo ">>> TOTAL_NODES=         $TOTAL_NODES "
echo ">>> NODES_PER_MEMBER=    $NODES_PER_MEMBER"
echo ">>> PROC_FOR_LETKF=      $PROC_FOR_LETKF"

}

save_configuration () {

local MYSCRPT=$1
local DESTDIR=$OUTPUTDIR/configuration/

#Save the current configuration files in the output directory.
if [ $FORECAST -eq 1 ] ; then

 cp -r $CDIR/configuration/forecast_conf/${CONFIGURATION}.sh  $DESTDIR  #Save experiment conf.

else

 cp -r $CDIR/configuration/analysis_conf/${CONFIGURATION}.sh  $DESTDIR  #Save experiment conf.

fi

cp -r $CDIR/configuration/machine_conf/${MCONFIGURATION}.sh  $DESTDIR  #Save machine conf.
cp -r $CDIR/configuration/$DOMAINCONF                        $DESTDIR  #Save domain conf.
cp -r $MYSCRIPT                                              $DESTDIR  #Save main script.

}


generate_run_forecast_script_k_background () {
#CREATE THE SCRIPT TO BE SUBMITED TO K COMPUTER

      local_script=$1
      if [ $RUN_ONLY_MEAN -eq 1 ] ; then
         local INIMEMBER=$MM  #Run only the last member.
      else
         local INIMEMBER=1
      fi

      #Prepare the script for bulk job.
      echo "#!/bin/bash    "                                     > $local_script
      echo "#PJM --rsc-list \"node=${TOTAL_NODES}\"              ">> $local_script
      echo "#PJM --rsc-list \"rscgrp=${QEUE}\"                   ">> $local_script
      echo "#PJM --rsc-list \"elapse=${ELAPSE}\"                 ">> $local_script
      echo "#PJM --mpi \"proc=${TOTAL_NODES}           \"        ">> $local_script
      #echo "#PJM --mpi \"shape=${TOTAL_NODES}       \"           ">> $local_script
      echo "#PJM --mpi assign-online-node                        ">> $local_script
      echo "#PJM --stg-transfiles all                            ">> $local_script
      echo "#PJM -S                                              ">> $local_script
      #PROGRAMS AND SCRIPTS
      echo "#PJM --stgin \"$RUNTIMELIBS/*            ./lib/      \" ">> $local_script
    if [ $ANALYSIS -eq 1  ] ; then
      echo "#PJM --stgin \"${TMPDIR}/SCRIPTS/*       ./SCRIPTS/  \" ">> $local_script
    fi
      #Generate staging list.
      #UPLOAD SCRIPTS
      M=$INIMEMBER
      while [ $M -le $MM ] ; do
       MEM=`ens_member $M`
       echo "#PJM --stgin \"$TMPDIR/ENSINPUT/${MEM}/*  ./WRF$MEM/ \"   ">> $local_script
       M=`expr $M + 1 `
      done
      #UPLOAD ANALYSIS FILES (ONLY IF ITERATION IS GREATHER THAN ONE)
      if [ $ITER -gt 1 ] ; then
       M=$INIMEMBER
       while [ $M -le $MM ] ; do
        MEM=`ens_member $M`
        echo "#PJM --stgin \"$OUTPUTDIR/anal/${CDATE}/anal$MEM  ./WRF$MEM/anal \" ">> $local_script
        M=`expr $M + 1 `
       done
      fi
     #COPY WRF MODEL AND DUMMY_MPI
      echo "#PJM --stgin \"${TMPDIR}/WRF/*                ./WRF/ \" ">> $local_script
      echo "#PJM --stgin \"${TMPDIR}/SPAWN/dummy-mpi      ./WRF/ \" ">> $local_script
     #STAGEOUT ANALYSIS AND GUES
     if [ $ANALYSIS -eq 1 ] ; then
      M=$INIMEMBER
      while [ $M -le $MEANMEMBER ] ; do
        MEM=`ens_member $M`
        #echo "#PJM --stgout   \"./LETKF/gues${MEM}          ${RESULTDIRG}/gues${MEM} \" ">> $local_script
        echo "#PJM --stgout   \"./WRF$MEM/*.log             ${RESULTDIRG}/           \" ">> $local_script
        M=`expr $M + 1 `
      done
       echo "#PJM --stgout   \"./LETKF/NOUT*                     ${RESULTDIRA}/ \"     ">> $local_script
       echo "#PJM --stgout   \"./LETKF/*                         ${TMPDIR}/CURRENT_LETKF/ \"     ">> $local_script
     fi
     #STAGEOUT FORECASTS
     if [ $FORECAST -eq 1 ] ; then
      M=$INIMEMBER
      while [ $M -le $MEANMEMBER ] ; do
        MEM=`ens_member $M`
        echo "#PJM --stgout   \"./LETKF/gs??${MEM}   ${RESULTDIRG}/        \" ">> $local_script
        echo "#PJM --stgout   \"./WRF$MEM/*.log      ${RESULTDIRG}/        \" ">> $local_script
        M=`expr $M + 1 `
      done

     fi
       
      echo "BASEDIR=\`pwd\`                           ">> $local_script
      echo ". /work/system/Env_base                   ">> $local_script
      echo "if [ -d "./lib" ] ; then                  ">> $local_script
      echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\${BASEDIR}/lib:/opt/local/globus/lib:/opt/FJSVpxtof/sparc64fx/lib64:/opt/FJSVtclang/GM-1.2.0-15/lib64">> $local_script
      echo "fi                                      ">> $local_script


     M=$INIMEMBER
     while [  $M -le $MM ];do
       MEM=`ens_member $M `
       echo "ln -sf \${BASEDIR}/WRF/* \${BASEDIR}/WRF${MEM}/ " >> $local_script
       echo "cp \${BASEDIR}/WRF/namelist.input.real  \${BASEDIR}/WRF${MEM}/namelist.input " >> $local_script
       M=`expr $M + 1 `
     done

     
     #SCRIPT
     echo "export PARALLEL=${PROC_PER_NODE}        " >> $local_script
     echo "export OMP_NUM_THREADS=${PROC_PER_NODE} " >> $local_script
  
     M=$INIMEMBER
     while [  $M -le $MM ];do
      JOB=1
      while [  $JOB -le $MAX_BACKGROUND_JOBS -a $JOB -le $MM ];do
      MEM=`ens_member $M `
         echo "mpiexec -np ${NODES_PER_MEMBER} --vcoordfile ./WRF/vcoord_${JOB} ./SCRIPTS/WRF_REAL.sh \${BASEDIR}/WRF${MEM}/ &  " >> $local_script
         JOB=`expr $JOB + 1 `
         M=`expr $M + 1 `
      done
      echo "time wait " >> $local_script
     done

    if [ $FORECAST -eq 1 -a $INTERPANA -eq 1  ] ; then
    #Only for ensemble forecast and for the case in which the LETKF analysis grid an the forecast grids are different.
     M=$INIMEMBER
     while [  $M -le $MM ];do
      JOB=1
      while [  $JOB -le $MAX_BACKGROUND_JOBS -a $JOB -le $MM ];do
      MEM=`ens_member $M `
         echo "mpiexec -np ${NODES_PER_MEMBER} --vcoordfile ./WRF/vcoord_${JOB} ./SCRIPTS/WRF_INTERPANA.sh \${BASEDIR}/WRF${MEM}/ &  " >> $local_script
         JOB=`expr $JOB + 1 `
         M=`expr $M + 1 `
      done
      echo "time wait " >> $local_script
     done
    fi
  
     if [ $ITER -gt 1 ];then
     M=$INIMEMBER
     while [  $M -le $MM ];do
      JOB=1
      while [  $JOB -le $MAX_BACKGROUND_JOBS -a $JOB -le $MM ];do
      MEM=`ens_member $M `
         echo "mpiexec -np 1 --vcoordfile ./WRF/vcoord_${JOB} ./SCRIPTS/WRF_PRE.sh \${BASEDIR}/WRF${MEM} ${MEM} &  " >> $local_script
         JOB=`expr $JOB + 1 `
         M=`expr $M + 1 `
      done
      echo "time wait " >> $local_script
     done      
     fi
        
     M=$INIMEMBER
     while [  $M -le $MM ];do
       MEM=`ens_member $M `
       echo "cp \${BASEDIR}/WRF/namelist.input.wrf  \${BASEDIR}/WRF${MEM}/namelist.input " >> $local_script
       M=`expr $M + 1 `
     done

 
     M=$INIMEMBER
     while [  $M -le $MM ];do
      JOB=1
      while [  $JOB -le $MAX_BACKGROUND_JOBS -a $JOB -le $MM ];do
      MEM=`ens_member $M `
         echo "mpiexec -np ${NODES_PER_MEMBER} --vcoordfile ./WRF/vcoord_${JOB} ./SCRIPTS/WRF_WRF.sh \${BASEDIR}/WRF${MEM} &  " >> $local_script
         JOB=`expr $JOB + 1 `
         M=`expr $M + 1 `
      done
      echo "time wait " >> $local_script
     done 

     M=$INIMEMBER
     while [  $M -le $MM ];do
      JOB=1
      while [  $JOB -le $MAX_BACKGROUND_JOBS -a $JOB -le $MM ];do
      MEM=`ens_member $M `

         echo "mpiexec -np 1 --vcoordfile ./WRF/vcoord_${JOB} ./SCRIPTS/WRF_POST.sh \${BASEDIR}/WRF${MEM} $MEM &  " >> $local_script
         JOB=`expr $JOB + 1 `
         M=`expr $M + 1 `
      done
      echo "time wait " >> $local_script
     done 

}

get_conventional_observations () {

local ADATES=`echo $ADATE | cut -c1-10`  #Short version of analysis date (YYYYMMDDHH)

  if [ -e $OBSDIR/obs$ADATES/t-3.dat  ] ; then
   cp $OBSDIR/obs$ADATES/t-3.dat $TMPDIR/OBS/obs01.dat
  else
   touch $TMPDIR/OBS/obs01.dat
  fi
  if [ -e $OBSDIR/obs$ADATES/t-2.dat  ] ; then
   cp $OBSDIR/obs$ADATES/t-2.dat $TMPDIR/OBS/obs02.dat
  else
   touch $TMPDIR/OBS/obs02.dat
  fi
  if [ -e $OBSDIR/obs$ADATES/t-1.dat  ] ; then
   cp $OBSDIR/obs$ADATES/t-1.dat $TMPDIR/OBS/obs03.dat
  else
   touch $TMPDIR/OBS/obs03.dat
  fi
  if [ -e $OBSDIR/obs$ADATES/t.dat  ] ; then
   cp $OBSDIR/obs$ADATES/t.dat $TMPDIR/OBS/obs04.dat
  else
   touch $TMPDIR/OBS/obs04.dat
  fi
  if [ -e $OBSDIR/obs$ADATES/t+1.dat  ] ; then
   cp $OBSDIR/obs$ADATES/t+1.dat $TMPDIR/OBS/obs05.dat
  else
   touch $TMPDIR/OBS/obs05.dat
  fi
  if [ -e $OBSDIR/obs$ADATES/t+2.dat  ] ; then
   cp $OBSDIR/obs$ADATES/t+2.dat $TMPDIR/OBS/obs06.dat
  else
   touch $TMPDIR/OBS/obs06.dat
  fi
  if [ -e $OBSDIR/obs$ADATES/t+3.dat  ] ; then
   cp $OBSDIR/obs$ADATES/t+3.dat $TMPDIR/OBS/obs07.dat
  else
   touch $TMPDIR/OBS/obs07.dat
  fi
}

get_radar_observations () {

local CDATEL=$WSDATE

local itradar=1

while [ $itradar -le $NRADARS  ] ; do
if [ $itradar -lt 10 ] ; then
  itradar=0$itradar
fi

  local it=1
  while [ ${CDATEL} -le ${WEDATE}  ] ; do
   if [ $it -lt 10 ] ; then
    it=0$it
   fi

  OBSFILE=$OBSDIR/RADAR${itrad}_${CDATEL}.dat
  echo $OBSFILE
  if [ -e $OBSFILE ] ; then
   cp -f $OBSFILE ./rad${CSLOT}01.dat
  fi

  it=`expr ${it} + 1 `
done

itradar=`expr ${itradar} + 1`
done

}

perturb_met_em () {

local local_script=$1
local EXEC=$TMPDIR/add_pert/compute_pert_metem.exe
if [ $RUN_ONLY_MEAN -eq 1 ] ; then
   local INIMEMBER=$MM  #Run only the last member.
else
   local INIMEMBER=1
fi

echo "#!/bin/bash                                                               "  > $local_script            
echo "set -x                                                                    " >> $local_script
#This script perturbs met_em files.
#To avoid perturbing again the same met_em file for the same ensemble member
#First the script checks wether the perturbed met_em file exists and if it do not
#exist the script creates the file.
echo "MEM=\$1                       #Ensemble member                            " >> $local_script
echo "SCALE_FACTOR=\$2               #Scale factor of perturbation amplitude.   " >> $local_script
echo "RANDOM_SCALE_FACTOR=\$3        #Amplitude for random perturbations.       " >> $local_script
echo "INIDATE=$CDATE                   #INITIAL DATE                            " >> $local_script
echo "ENDDATE=$BDYDATE                 #END DATE                                " >> $local_script
echo "BOUNDARY_DATA_FREQ=$BOUNDARY_DATA_FREC #Boundary data frequency (seconds) " >> $local_script
echo "BOUNDARY_DATA_PERTURBATION_FREQ=$BOUNDARY_DATA_PERTURBATION_FREQ          " >> $local_script
echo "WORKDIR=$TMPDIR/ENSINPUT/\$MEM/  #Temporary work directory                " >> $local_script
echo "PERTMETEMDIR=$PERTMETEMDIR    #Met em data base for perturbations         " >> $local_script
echo "PERTDATEDIR=$PERTDATEDIR      #Dates database for perturbations           " >> $local_script
echo "METEMDIR=$METEMDIR            #Met em for curren experiment dates         " >> $local_script
echo "EXEC=$EXEC                     #Executable for perturbation computation.  " >> $local_script
echo "source $TMPDIR/SCRIPTS/util.sh                                            " >> $local_script
echo "ulimit -s unlimited                                                       " >> $local_script
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_ADD:\$LD_LIBRARY_PATH             " >> $local_script
echo "export PATH=$LD_PATH_ADD:$PATH                                            " >> $local_script
echo "mkdir -p \$WORKDIR                                                       " >> $local_script
echo "cd \$WORKDIR                                                             " >> $local_script
echo "rm -fr \$WORKDIR/met_em*                                                 " >> $local_script

#################################################
#   CYCLE TO CREATE THE PERTURBATIONS
#################################################
echo "CDATE=\$INIDATE                                                          " >> $local_script
echo "while [  \$CDATE -le \$ENDDATE ] ; do                                    " >> $local_script
echo "CFILE=\`met_em_file_name \$CDATE 01\`                                    " >> $local_script
   #Get the dates in the perturbation data that are closer to the current date.
echo "LDATE=\`date_floor \$CDATE \$BOUNDARY_DATA_PERTURBATION_FREQ \`          " >> $local_script 
echo "UDATE=\`date_edit2 \$LDATE \$BOUNDARY_DATA_PERTURBATION_FREQ \`          " >> $local_script
   #Get the dates that we will use to create the perturbation and link the corresponding met_em files
echo "   read CDATE1 CDATE2 < \$PERTDATEDIR/\${LDATE}_\$MEM.dates              " >> $local_script
echo "   TMPFILE=\`met_em_file_name \$CDATE1 01 \`                             " >> $local_script
echo "   ln -sf \$PERTMETEMDIR/\$TMPFILE ./input_filea1.nc                     " >> $local_script
echo "   TMPFILE=\`met_em_file_name \$CDATE2 01 \`                             " >> $local_script
echo "   ln -sf \$PERTMETEMDIR/\$TMPFILE ./input_filea2.nc                     " >> $local_script
echo "   read CDATE1 CDATE2 < \$PERTDATEDIR/\${UDATE}_\$MEM.dates              " >> $local_script
echo "   TMPFILE=\`met_em_file_name \$CDATE1 01 \`                             " >> $local_script
echo "   ln -sf \$PERTMETEMDIR/\$TMPFILE ./input_fileb1.nc                     " >> $local_script
echo "   TMPFILE=\`met_em_file_name \$CDATE2 01 \`                             " >> $local_script
echo "   ln -sf \$PERTMETEMDIR/\$TMPFILE ./input_fileb2.nc                     " >> $local_script
   #Copy the unperturbed met_em file (this one will be modified).
echo "   cp \$METEMDIR/\$CFILE .                                               " >> $local_script
echo "   ln -sf ./\$CFILE ./ctrl_met_em.nc                                     " >> $local_script
echo "   chmod 766 ./ctrl_met_em.nc                                            " >> $local_script
   #Get the time dinstance in seconds between the current time and LDATE
echo "   TIMEDISTANCE=\`date_diff \$CDATE \$LDATE \`                           " >> $local_script
   #Run the program 
   # ctrl_met_em.nc = ctrl_met_em.nc + scale_factor *[ ( input_file1.nc - input_file2.nc ) ]
   # the [] indicates a time interpolation between LDATE and UDATE.
echo "   $EXEC \$SCALE_FACTOR \$RANDOM_SCALE_FACTOR \$TIMEDISTANCE \$BOUNDARY_DATA_PERTURBATION_FREQ " >> $local_script
echo "CDATE=\`date_edit2 \$CDATE \$BOUNDARY_DATA_FREQ \`                       " >> $local_script
echo "done                                                                     " >> $local_script
#We are done!
chmod 766 $local_script


M=$INIMEMBER
while [ $M -le $MM ] ; do
  
  RUNNING=0
  while [ $RUNNING -le $MAX_RUNNING -a $M -le $MM ] ; do
    MEM=`ens_member $M `
    #Do not perturb the ensemble mean run.
    if [ $M -eq $MEANMEMBER ] ; then
       TMP_SCALE_FACTOR="0.0"
       TMP_RANDOM_SCALE_FACTOR="0.0"
    else
       TMP_SCALE_FACTOR=$SCALE_FACTOR
       TMP_RANDOM_SCALE_FACTOR=$RANDOM_SCALE_FACTOR
    fi
    ssh $PPSSERVER " $local_script $MEM $TMP_SCALE_FACTOR $TMP_RANDOM_SCALE_FACTOR > $TMPDIR/add_pert/perturb_met_em${MEM}.log  2>&1 " &
    RUNNING=`expr $RUNNING + 1 `
    M=`expr $M + 1 `
  done
  time wait
done

}


wrf_to_met_em () {
#This function transforms WRF output into met_em files
#This function should be run after running perturb_met_em
local local_script=$1

if [  $INTERPANA -eq 1 ] ; then

if [ $RUN_ONLY_MEAN -eq 1 ] ; then
   local INIMEMBER=$MM  #Run only the last member.
else
   local INIMEMBER=1
fi


#We only perform this operation if interpana is 1.
#This means that the input analysis grid is different from the forecast grid 
#Usually the forecast grid has a larger resolution.

echo ">>> [Warning]: This function works only if analysis files are in mercator or lat-lon grids"
local EXEC=$TMPDIR/wrf_to_wps/wrf_to_wps.exe
local APATH=$ANALYSIS_SOURCE/anal/$CDATE/


echo "#!/bin/bash                                                               "  > $local_script
echo "set -x                                                                    " >> $local_script
#This script perturbs met_em files.
#To avoid perturbing again the same met_em file for the same ensemble member
#First the script checks wether the perturbed met_em file exists and if it do not
#exist the script creates the file.
echo "MEM=\$1                      #Ensemble member                             " >> $local_script
echo "WORKDIR=$TMPDIR/ENSINPUT/\$MEM/wrf_to_wps #Temporary work directory       " >> $local_script
echo "EXEC=$EXEC                   #Executable                                  " >> $local_script
echo "source $TMPDIR/SCRIPTS/util.sh                                            " >> $local_script
echo "ulimit -s unlimited                                                       " >> $local_script
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_ADD:\$LD_LIBRARY_PATH             " >> $local_script
echo "export PATH=$LD_PATH_ADD:$PATH                                            " >> $local_script
echo "mkdir -p \$WORKDIR                                                        " >> $local_script
echo "cd \$WORKDIR                                                              " >> $local_script

#################################################
#  FROM LETKF_ANALYSIS TO UNGRIB BINARY FORMAT
#################################################
echo "AFILE=\$APATH/anal${MEM}                                                 " >> $local_script
echo "ln -fs \$AFILE ./input.nc                                                " >> $local_script
echo "$EXEC                                                                    " >> $local_script
echo "ungrib_file=ungrib_file_name $CDATE DATA                                 " >> $local_script
echo "ln -fs ./output.grd \$ungrib_file                                        " >> $local_script
#################################################
#  FROM UNGRIB BINARY FORMAT TO MET_EM
#################################################
echo "ln -fs \$TMPDIR/WPS/* ./                                                 " >> $local_script
echo "cp namelist.wps.template namelist.wps                                    " >> $local_script
echo "edit_namelist_wps ./namelist.wps $CDATE $CDATE $BOUNDARY_DATA_FREC       " >> $local_script
echo "$MPIBIN -np 1 ./metgrid.exe                                              " >> $local_script
echo "metgrid_file=met_em_file_name $CDATE 01                                  " >> $local_script
echo "mv \$metgrid_file ../\${metgrid_file}.anal                               " >> $local_script

chmod 766 $local_script

M=$INIMEMBER
while [ $M -le $MM ] ; do

  RUNNING=0
  while [ $RUNNING -le $MAX_RUNNING -a $M -le $MM ] ; do
    MEM=`ens_member $M`
    #Do not perturb the ensemble mean run.
    ssh $PPSSERVER " $local_script $MEM  > $TMPDIR/wrf_to_wps/wrf_to_wps{MEM}.log  2>&1 " &
    RUNNING=`expr $RUNNING + 1 `
    M=`expr $M + 1 `
  done
  time wait
done

fi


}

arw_postproc () {

#This script perturbs met_em files.
#To avoid perturbing again the same met_em file for the same ensemble member
#First the script checks wether the perturbed met_em file exists and if it do not
#exist the script creates the file.

# FOR THE ANALYSIS CASE

if [ $ANALYSIS -eq 1 ] ; then

  local CDATE=$ADATE                     #INITIAL DATE
  local WORKDIR=$TMPDIR/ENSINPUT/        #Temporary work directory
  local ANALDIR=$OUTPUTDIR/anal/$CDATE/
  local GUESDIR=$OUTPUTDIR/gues/$CDATE/
                                         
  ARWPOST_FREC=$WINDOW_FREC

  INPUT_ROOT_NAME=tmpin
  OUTPUT_ROOT_NAME=tmpout
  cp $NAMELISTARWPOST $WORKDIR/namelist.ARWpost
  edit_namelist_arwpost $WORKDIR/namelist.ARWpost $CDATE $CDATE $ARWPOST_FREC

  echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_ADD:\$LD_LIBRARY_PATH " >  ${WORKDIR}/tmp.sh
  echo "export PATH=$LD_PATH_ADD:$PATH                                " >> ${WORKDIR}/tmp.sh
  if [ $SYSTEM -eq  1 ] ; then
     echo " ulimit -s unlimited                                       " >> ${WORKDIR}/tmp.sh
  fi
  echo "MEM=\$1                                                       " >> ${WORKDIR}/tmp.sh
  echo "mkdir -p ${WORKDIR}/\${MEM}                                   " >> ${WORKDIR}/tmp.sh
  echo "cd ${WORKDIR}/\${MEM}                                         " >> ${WORKDIR}/tmp.sh
  echo "DATADIR=${ANALDIR}                                            " >> ${WORKDIR}/tmp.sh
  echo "PREFIX=anal                                                   " >> ${WORKDIR}/tmp.sh
  echo "ln -sf \${DATADIR}/\${PREFIX}\${MEM} ./tmpin                  " >> ${WORKDIR}/tmp.sh
  echo "ln -sf \${DATADIR}/plev\${MEM}.dat ./tmpout.dat               " >> ${WORKDIR}/tmp.sh
  echo "ln -sf \${DATADIR}/plev\${MEM}.ctl ./tmpout.ctl               " >> ${WORKDIR}/tmp.sh
  echo "ln -sf $ARWPOST/src .                                         " >> ${WORKDIR}/tmp.sh
  echo "ln -sf $WORKDIR/namelist.ARWpost ./namelist.ARWpost           " >> ${WORKDIR}/tmp.sh
  echo "$ARWPOST/ARWpost.exe > \${DATADIR}/arwpost\${MEM}.log         " >> ${WORKDIR}/tmp.sh
  echo "sed -i 's/tmpout/plev'\${MEM}'/g' \${DATADIR}/plev\${MEM}.ctl " >> ${WORKDIR}/tmp.sh

  echo "DATADIR=${GUESDIR}                                            " >> ${WORKDIR}/tmp.sh
  echo "PREFIX=gues                                                   " >> ${WORKDIR}/tmp.sh
  echo "ln -sf \${DATADIR}/\${PREFIX}\${MEM} ./tmpin                  " >> ${WORKDIR}/tmp.sh
  echo "ln -sf \${DATADIR}/plev\${MEM}.dat ./tmpout.dat               " >> ${WORKDIR}/tmp.sh
  echo "ln -sf \${DATADIR}/plev\${MEM}.ctl ./tmpout.ctl               " >> ${WORKDIR}/tmp.sh
  echo "ln -sf $ARWPOST/src .                                         " >> ${WORKDIR}/tmp.sh
  echo "ln -sf $WORKDIR/namelist.ARWpost ./namelist.ARWpost           " >> ${WORKDIR}/tmp.sh
  echo "$ARWPOST/ARWpost.exe > \${DATADIR}/arwpost\${MEM}.log         " >> ${WORKDIR}/tmp.sh
  echo "sed -i 's/tmpout/plev'\${MEM}'/g' \${DATADIR}/plev\${MEM}.ctl " >> ${WORKDIR}/tmp.sh

  chmod 766 ${WORKDIR}/tmp.sh

  M=1
  while [ $M -le $MEANMEMBER ] ; do
    RUNNING=0
    while [ $RUNNING -le $MAX_RUNNING -a $M -le $MEANMEMBER ] ; do
      MEM=`ens_member $M`
      ssh $PPSSERVER " ${WORKDIR}/tmp.sh $MEM " & 
      RUNNING=`expr $RUNNING + 1 `
      M=`expr $M + 1 `
    done
    time wait
  done

fi

if [ $FORECAST -eq  1 ] ; then 
#This script perturbs met_em files.
#To avoid perturbing again the same met_em file for the same ensemble member
#First the script checks wether the perturbed met_em file exists and if it do not
#exist the script creates the file.

  local IDATE=$CDATE                   #INITIAL DATE
  local EDATE=$FDATE                   #FINAL DATE
  local WORKDIR=$TMPDIR/ENSINPUT/      #Temporary work directory
  local GUESDIR=$RESULTDIRG

  if [ $RUN_ONLY_MEAN -eq 1 ] ; then
     local INIMEMBER=$MM  #Run only the last member.
  else
     local INIMEMBER=1
  fi

  ARWPOST_FREC=$WINDOW_FREC

  INPUT_ROOT_NAME=tmpin
  OUTPUT_ROOT_NAME=tmpout
  cp $NAMELISTARWPOST $WORKDIR/namelist.ARWpost.template
  edit_namelist_arwpost $WORKDIR/namelist.ARWpost $CDATE $CDATE $ARWPOST_FREC

  echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_ADD:\$LD_LIBRARY_PATH " >  ${WORKDIR}/tmp.sh
  echo "export PATH=$LD_PATH_ADD:$PATH                                " >> ${WORKDIR}/tmp.sh
  if [ $SYSTEM -eq  1 ] ; then
     echo " ulimit -s unlimited                                       " >> ${WORKDIR}/tmp.sh
  fi
  echo "MEM=\$1                                                       " >> ${WORKDIR}/tmp.sh
  echo "cd ${WORKDIR}/\${MEM}                                         " >> ${WORKDIR}/tmp.sh
  echo "ln -sf $TMPDIR/ARWPOST/src .                                  " >> ${WORKDIR}/tmp.sh

  local TMPDATE=$IDATE
  it=1
  while [ $TMPDATE -le $EDATE  ] ; do
   if [ $it -lt 10 ] ; then
     it=0$it
   fi
   echo "cp namelist.ARWpost.template namelist.ARWpost                                       " >> ${WORKDIR}/tmp.sh
   echo "edit_namelist_arwpost $WORKDIR/namelist.ARWpost $TMPDATE $TMPDATE $ARWPOST_FREC     " >> ${WORKDIR}/tmp.sh
   echo "ln -sf ${GUESDIR}/gs${it}\${MEM} ./tmpin                                            " >> ${WORKDIR}/tmp.sh
   echo "ln -sf ${GUESDIR}/\$localfile    ./tmpin                                            " >> ${WORKDIR}/tmp.sh
   echo "ln -sf ${GUESDIR}/plev${TMPDATE}_\${MEM}.dat ./tmpout.dat                           " >> ${WORKDIR}/tmp.sh
   echo "ln -sf ${GUESDIR}/plev${TMPDATE}_\${MEM}.ctl ./tmpout.ctl                           " >> ${WORKDIR}/tmp.sh
   echo "ln -sf $ARWPOST/src .                                                               " >> ${WORKDIR}/tmp.sh
   echo "$ARWPOST/ARWpost.exe > ${GUESDIR}/arwpost\${MEM}.log                                " >> ${WORKDIR}/tmp.sh
   echo "sed -i 's/tmpout/plev${TMPDATE}_'\${MEM}'/g' ${GUESDIR}/plev${TMPDATE}_\${MEM}.ctl  " >> ${WORKDIR}/tmp.sh

   TMPDATE=date_edit2 $TMPDATE $WINDOW_FREC
   it=`expr $it + 1 `
  done

  chmod 766 ${WORKDIR}/tmp.sh

  M=$INIMEMBER
  while [ $M -le $MM ] ; do
    RUNNING=0
    while [ $RUNNING -le $MAX_RUNNING -a $M -le $MM ] ; do
      MEM=`ens_member $M`
      ssh $PPSSERVER " ${WORKDIR}/tmp.sh $MEM  " &
      RUNNING=`expr $RUNNING + 1 `
      M=`expr $M + 1 `
    done
    time wait
  done

fi

}

# ------------------------------------
#   FUNCTION : pjsub_end_check
# ------------------------------------
sub_and_wait() {
input_shell=$1
n=$2
if [ $# -eq 1 ]; then
        n=$RANDOM
fi

        request_log=$1
        id=`pjsub -z jid ${input_shell}`
        echo "${input_shell} SUMBITTED, ID= $id"

#        len=`echo ${#id}`
#        lenm=`expr ${len} - 1`
#        id=`echo ${id} | cut -c2-${lenm}`

        echo "WAITING FOR THE JOB..."
        qsub_end=0
        while [ ${qsub_end} -eq 0 ]
        do
                pjstat > qstat.log${n}
                grep ${id} qstat.log${n} > grep.log${n}

                if [ ! -s grep.log${n} ]; then
                        qsub_end=1
                        echo "JOB FINISHED"
                else
                        sleep 30
                fi
        done

        rm -f qstat.log${n} grep.log${n}
}

check_analysis () {
#TO DO, add checks over forecasts.

local M=1
local cycle_error=0
while [ $M -le $MEMBER ] ; do
 local MEM=`ens_member $M `
 grep "SUCCESS COMPLETE WRF" ${RESULTDIRG}/wrf${MEM}.log > null
 if [ $? -ne 0 ] ; then
   echo "[Error]: WRF for ensemble member $MEM" 
   echo "====================================="
   echo "SHOWING LAST PART OF wrf${MEM}.log     "
   tail ${RESULTDIRG}/wrf${MEM}.log 
   cycle_error=1
 fi

 if [ $ITER -ne 1 ] ; then
  grep  "Update_bc completed successfully" ${RESULTDIRG}/daupdatebc${MEM}.log > null
  if [ $? -ne 0 ] ; then
   echo "[Error]: WRF da update bc for ensemble member $MEM"
   echo "====================================="
   echo "SHOWING LAST PART OF dapudatebc${MEM}.log     "
   tail ${RESULTDIRG}/daupdatebc${MEM}.log
   cycle_error=1
  fi
 fi

 grep  "Successful completion of ARWpost" ${RESULTDIRG}/arwpost${MEM}.log > null
 if [ $? -ne 0 ] ; then
   echo "[Error]: ARWPOST for gues ensemble member $MEM"
   echo "====================================="
   echo "SHOWING LAST PART OF arwpost${MEM}.log     "
   tail ${RESULTDIRG}/arwpost${MEM}.log
   cycle_error=1
 fi

 grep  "Successful completion of ARWpost" ${RESULTDIRA}/arwpost${MEM}.log > null
 if [ $? -ne 0 ] ; then
   echo "[Error]: ARWPOST for anal ensemble member $MEM"
   echo "====================================="
   echo "SHOWING LAST PART OF arwpost${MEM}.log     "
   tail ${RESULTDIRA}/arwpost${MEM}.log
   cycle_error=1
 fi

 M=`expr $M + 1 `
done

local MMS=`ens_member $MEANMEMBER`
if [ ! -e ${RESULTDIRA}/anal$MMS ] ; then
  echo "[Error]: Cannot find analysis ensemble mean."
  cycle_error=1
fi

grep  "PARTIAL OBSERVATIONAL DEPARTURE" ${RESULTDIRA}/NOUT-000 > null
if [ $? -ne 0  ] ; then
 echo "[Error]: LETKF do not finish properly."
 tail ${RESULTDIRA}/NOUT-000
 cycle_error=1
fi

if [ $cycle_error -eq 1 ] ; then
   echo "CICLE ABNORMAL END -> STOP EXECUTION "
   exit 1
else
   echo "CICLE NORMAL END -> CONTINUE EXECUTION"
fi

}


