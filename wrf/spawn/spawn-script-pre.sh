#!/bin/bash
# ===================================================================
# Originally taken from ensfcst.sh
# Merge LETKF analysis with input data of WRF
# Shigenori Otsuka
# ===================================================================
set -e

### DUMMY MPI PROGRAM TO INVOKE MPI_Init
../dummy-mpi

. spawn-script.input ### INPUT PARAMETERS

flag_init=0
if test $FLAG -eq 0 ; then
    #INITIAL STEP USE INPUT DATA AS IT IS.
    flag_init=1
    
elif test $FLAG -ge 1 ; then
    #FURTHER STEPS USE PREVIOUS ANALYSIS.
    echo "Overwrite Initial Data" ${MEM}
    for file in "../init_merge" "wrfinput_d01" "anal.grd" ; do
	if [ ! -e $file ] ; then
	    echo $file "not found for" ${MEM} 1>&2
	    exit 1
	fi
    done
    ../init_merge wrfinput_d01 anal.grd > log.init_merge${MEM}

    if test $FLAG -eq 1 ; then
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
	../da_update_bc.exe > log.da_update_low${MEM}
    fi

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
fi
# --- check initial fields ---
if test ${flag_init} -eq 0 ; then
    echo "xxx failed to prepare initial field : ${MEM} xxx" 1>&2
    exit 1
fi



