#!/bin/bash
# ===================================================================
# Originally taken from ensfcst.sh
# Convert WRF output to GrADS format
# Shigenori Otsuka
# ===================================================================
set -e

### DUMMY MPI PROGRAM TO INVOKE MPI_Init
../dummy-mpi

. spawn-script.input ### INPUT PARAMETERS
. ../util.sh

# --- end transaction ---
rm -f log_error* log_out*
for i in `seq 0 \`expr ${MPI_PER_MEMBER} - 1\`` ; do
    ii=`printf %04d ${i}`
    cat rsl.error.${ii} >> log_error.${MEM}
    cat rsl.out.${ii}   >> log_out.${MEM}
done

# ------------------------------
#     CONVERT
# ------------------------------ 

SWINDOW=`expr ${GUESFT} - ${WINDOW}`
counter=1
for diffmin in `seq ${SWINDOW} ${GUESINTV} ${GUESFT}` ; do
    read vyy vmm vdd vhh vmn <<EOF
`date_edit ${syy} ${smm} ${sdd} ${shh} ${smn} ${diffmin}`
EOF

    # --- input file ---
    if [ ! -e wrfout_d01_${vyy}-${vmm}-${vdd}_${vhh}:${vmn}:00 ] ; then
	echo wrfout_d01_${vyy}-${vmm}-${vdd}_${vhh}:${vmn}:00 "not found" 1>&2
	exit 1
    fi
    # --- run the pgm ---
    ../nc2grd wrfout_d01_${vyy}-${vmm}-${vdd}_${vhh}:${vmn}:00 const.grd gs`printf %02d ${counter}`${MEM}.grd > log.convert_${counter} 2>&1 ### nc2grd_argv version
    counter=`expr ${counter} + 1`
done

#DONE!
