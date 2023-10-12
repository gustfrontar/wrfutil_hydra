Code compilation:

Step 0
======

The follwing code should be installed before hand
WRF, WPS, WRFDA

WRF, WRFDA and WPS code is available from the following git repositories.
https://github.com/wrf-model/WRF
https://github.com/wrf-model/WPS

This guide provide information on how to compile WRF, WPS and WRFDA from scratch
https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php

Also the geographycal statical data has to be downloaded from this location
https://www2.mmm.ucar.edu/wrf/users/download/get_sources_wps_geog.html

Additionally to WRF code, the following programs needs to be compiled.
There are a couple of compilation scripts available for some compilers and 
computer environments. If these do not correspond to your environment you will
have to edit the compilation scripts before runing a successfull compilation.

LETKF -> /letkf_code/src/wrf/letkf/make_letkf_$COMPILER.sh 
ENSEMBLE SPAWNER -> /letkf_code/src/ensemble_spawn/make_spawner_$COMPILER.sh [ Only required to run the scripts that use spawn to distribute the ensemble members ]
WRF_TO_WPS -> /letkf_code/src/wrf/wrf_to_wps/make_wrf_to_wps_$COMPILER.sh 


Step 1
======

Edit and run the script>
/wrf_code/create_tar_$COMPILER

This will create tar files with the executables and required files to run
WRF, WPS and WRFDA. This tar files will be used to create experiments using the
scripts under the wrfutil folder.


Step 2
======

Go to the /wrfutilv$VERSION  folder.

Under the "conf" folder you will find subfolders containing configuration files templates. 

Copy or modify one of these templates according to your environment.

config.env -> contains the experiment name and the location of the "tar" files containing the compiled code. 

machine.conf -> contains the machine characteristics and the resources to be required by the experiment.
                It also contain some machine dependent variables that need to be set to run the experiments.
                If your system is not listed you will have to find out which additional settings are required
                by your system in order to run the different steps. 

model.conf   -> contains the basic settings of the WRF model domain. Size, location, resolution, etc. 
                All other configurations will be taken from the namelist templates located in the "namelists" folder.
          

assimilation.conf -> in an assimilation experiment (this is defined in the config.env file) this file controls the assimilation frequency
                     the start date, end date, etc. 

forecast.conf     -> in a forecast experiment (this is defined in the config.env file) this file controls the forecast initialization frequency
                     , forecast lead time, etc.

letkf.conf        -> in an assimilation experiment this file controls some aspects of the Local Ensemble Transform Based assimilation as
                     localization scales, type of observation assimilated, inflation factors, etc.


Some detail about the meaning of different flags defined in the configuration files are provided in the configuration files.


Step 3
======

To create a experiment use the scripts located in  /wrfutilv$VERSION/bin/

createExpASIM.sh $CONF_TEMPLATE $EXP_PATH -> creates a data assimilation experiment
createExpFCST.sh $CONF_TEMPLATE $EXP_PATH -> creates a forecast experiment

$EXP_PATH is the absolute PATH where your experiment will be created (e.g. /home/user/my_experiments/my_experiment_name/). 
This PATH will be "stored" in the file /conf/config.env within your experiment path.

$CONF_TEMPLATE should be the name of one of the folders under /wrfutil$VERSION/conf/ you can create as many configuration
templates as you whish, or you can use the ones provided in this repository.

Once th folder $EXP_PATH is created, a copy of the bash scripts and the tar files will be done so all the required information will be 
copyed to the new folder. The only data that is not copied to the experiment folder is the GFS/WRF boundary conditions data
and the WRF static geographycal data (mainly because of their size). 


Step 4
======

Prepare boundary conditions data for your experiments. For the moment the outhermost boundary conditions for the experiment
came from the Global Ensemble Forecasting System (NCEP). Data is freely available through the Amazon Web Services repository.
The script located in wrfutil$VERSION/bin/python/get_archived_gefs_data.py can be used to downlad this data for a particular
time. 
More information about this dataset can be found here:
https://aws.amazon.com/marketplace/pp/prodview-qumzmkzc2acri#resources

Step 5
======

a) Run a simple forecast experiment. 

Runing an ensemble forecast experiment initialized from the GEFS using the WRF model. 

You can use the main_FCST.sh script. This script will run WPS using the GEFS data to create the initial and boundary conditions.
If required it will expand the ensemble of boundary conditions using one of the available methods. 
Then it will run real.exe and wrf.exe for each ensemble members. Ensemble members will be distributed into different nodes
according to the settings in machine.conf 

NOTE: In FUGAKU the expansion of the ensemble to add boundary and initial perturbations can not be performed in the computation nodes.
So you need to run main_WPS.sh first to generate the original met_em... files from the GEFS ensemble in the computation nodes.
Then you need to run main_Pert.sh in the login or FUGAKU pps nodes, and the main_FCST.sh script to run real and wrf (with RUN_WPS and RUN_PERT = 0 in the forecast.conf 
configuration file). 

b) Run a simple data assimilation experiment. 





















