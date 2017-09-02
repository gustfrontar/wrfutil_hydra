#!/bin/bash

#Pourpose: Run the scale model in Hakushu

CONFIGURATION="germany_3km"

SRCDIR=`pwd`

source $SRCDIR/../configuration/${CONFIGURATION}/config.sh

source $SRCDIR/util.sh 

#Load modules in Hakushu   #######################################

. /usr/share/modules/init/sh  #This enables the module alias

module load common_intel/hdf5
module load common_intel/netcdf
module load common_intel/netcdf-fortran
module load intel-2013.1.117
module load mpi/intel-4.0

#Create tmpdir            #######################################

safe_init_tmpdir $TMPDIR

#Create outputdir         #######################################

#safe_init_outputdir $OUTPUTDIR


#Link configuration files #######################################

link_files $TMPDIR

cd $TMPDIR

ln -sf $DATABASEROOT ./database

ln -sf $BDYDATAROOT/* .

#Edit namelist_pp

edit_namelist $TMPDIR/nml.scale_pp

#Run scale+_pp             #######################################

cd $TMPDIR

mpiexec -np $NPROCS ./scale-rm_pp nml.scale_pp > scale_pp.log

#Link boundary data         #######################################

#link_bdy_data $TMPDIR

#Edit namelist_init

#edit_namelist $TMPDIR/nml.scale_init

#Run scale_init             #######################################

#cd $TMPDIR

#mpiexec -np $NPROCS ./scale-rm_init nml.scale_init > scale_init.log

#Edit namelist_scale

#edit_namelist $TMPDIR/nml.scale

#Run scale model            #######################################

#cd $TMPDIR

#mpiexec -np $NPROCS ./scale-rm nml.scale > scale_rm.log


#Edit namelist post         #######################################

#edit_namelist_postproc $TMPDIR/nml.net2g


#Run postprocessing         #######################################

#mpiexec  -np  $NPROCS  ./net2g  nml.net2g  > net2g.log





