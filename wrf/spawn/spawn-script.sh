#!/bin/bash
# ===================================================================
# ensfcst.sh
#   This script runs the WRF model
# ===================================================================
set -e

. spawn-script.input ### INPUT PARAMETERS

num_mpi=${MPI_PER_MEMBER}
window_length=6
ft_lag=3
flag_init=0

ft=`expr ${window_length} + ${ft_lag}`

source ../util.sh

# ==============================
#    COPY INPUT DATA
# ==============================

#Create the output directory.

min=`expr ${window_length} \* 60`
read AY AM AD AH AMN <<EOF
`date_edit ${syy} ${smm} ${sdd} ${shh} ${smn} ${min}`
EOF

if test $FLAG -eq 0 ; then
    #INITIAL STEP USE INPUT DATA AS IT IS.
    flag_init=1
    
elif test $FLAG -eq 1 ; then
    #FURTHER STEPS USE PREVIOUS ANALYSIS.
    echo "Overwrite Initial Data" ${MEM}
    ls -lh ../init_merge wrfinput_d01 anal.grd
    for file in "../init_merge" "wrfinput_d01" "anal.grd" ; do
	if [ ! -e $file ] ; then
	    echo $file "not found for" ${MEM} 1>&2
	    exit 1
	fi
    done
    ../init_merge wrfinput_d01 anal.grd > log.init_merge${MEM} #2> /dev/null
    
    # UPDATE BOUNDARY CONDITIONS.
    echo "Update lower boundary condition" ${MEM}
    cat << EOF > parame.in
&control_param
da_file='./wrfinput_d01'
wrf_input='./wrfinput_d01.org'
debug=.true.
update_lateral_bdy=.false.
update_low_bdy=.true.
iswater=16
/
EOF
    ../da_update_bc.exe > log.da_update_low${MEM} #2> /dev/null
    #  update_bc_lateral
    echo "Update lateral boundary condition" ${MEM}
    cat << EOF > parame.in
&control_param
da_file='./wrfinput_d01'
wrf_bdy_file='./wrfbdy_d01'
debug=.true.
update_lateral_bdy=.true.
update_low_bdy=.false.
iswater=16
/
EOF
    ../da_update_bc.exe > log.da_update_lat${MEM} #2> /dev/null
    flag_init=1
    
elif test $FLAG -eq 2 ; then
    #FURTHER STEPS USE PREVIOUS ANALYSIS.
    echo "Overwrite Initial Data" ${MEM}
    ls -lh ../init_merge wrfinput_d01 anal.grd
    for file in "../init_merge" "wrfinput_d01" "anal.grd" ; do
	if [ ! -e $file] ; then
            echo $file "not found" 1>&2
            exit 1
        fi
    done
    ../init_merge wrfinput_d01 anal.grd > log.init_merge${MEM} #2> /dev/null
    
    # UPDATE BOUNDARY CONDITIONS (only lateral)
    echo "Update lateral boundary condition" ${MEM}
    cat << EOF > parame.in
&control_param
da_file='./wrfinput_d01'
wrf_bdy_file='./wrfbdy_d01'
debug=.true.
update_lateral_bdy=.true.
update_low_bdy=.false.
iswater=16
/
EOF
    ../da_update_bc.exe > log.da_update_lat${MEM} #2> /dev/null
    flag_init=1
fi
# --- check initial fields ---
if test ${flag_init} -eq 0 ; then
    echo "xxx failed to prepare initial field : ${MEM} xxx" 1>&2
    exit 1
fi


# ------------------------------
#    STEP 320: WRF
# ------------------------------ 
# -- remove logs ---

rm -f rsl_error.* rsl_out.*
rm -f wrfout_d01*
# --- run the pgm ---
#ulimit -s unlimited
export OMP_NUM_THREADS=${WRF_NUM_THREADS}
export PARALLEL=${WRF_NUM_THREADS}
./wrf.exe ### SPAWNED

# --- end transaction ---
rm -f log_error* log_out*
for i in `seq 0 \`expr ${num_mpi} - 1\`` ; do
    ii=`printf %04d ${i}`
    cat rsl.error.${ii} >> log_error.${MEM}
    cat rsl.out.${ii}   >> log_out.${MEM}
done

# ------------------------------
#     CONVERT
# ------------------------------ 

for it in `seq 3 ${ft}` ; do
    diff_hour=${it}
    min=`expr ${diff_hour} \* 60`
    read vyy vmm vdd vhh vmn <<EOF
`date_edit ${syy} ${smm} ${sdd} ${shh} ${smn} ${min}`
EOF

    ip=`expr ${it} - 2`
    ip=`printf "%02d" ${ip}`
    # --- input file ---
    if [ ! -e wrfout_d01_${vyy}-${vmm}-${vdd}_${vhh}:${vmn}:00 ] ; then
	echo wrfout_d01_${vyy}-${vmm}-${vdd}_${vhh}:${vmn}:00 "not found" 1>&2
	exit 1
    fi
    # --- run the pgm ---
    ../nc2grd wrfout_d01_${vyy}-${vmm}-${vdd}_${vhh}:${vmn}:00 const.grd gs${ip}${MEM}.grd > log.convert_${ip} 2>&1 ### nc2grd_argv version
done

#DONE!
