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

Additional notes and configuration files will be included in the wrf_code folder.
A folder for each machine should be there containing configure files specific for that
machine as well as some additional considerations for the compilation of WRF on that machine.

Also the geographycal static data has to be downloaded from this location
https://www2.mmm.ucar.edu/wrf/users/download/get_sources_wps_geog.html

Additionally to WRF code, the following programs needs to be compiled.
There are a couple of compilation scripts available for some compilers and 
computer environments. If these do not correspond to your environment you will
have to edit the compilation scripts before runing a successfull compilation.

These additional codes (LETKF and associated tools, WRF_TO_WPS and associated tools and,
PERT_MET_EM and associated tools) can be compiled using the make_all.sh script under
the directory letkf_code. This script include compilation options for some machines.
You can use these compilation options as templates for your own system. 

Binaries are packed into tar files in the letkf_code/bin directory.


Step 1
======

Compile the LETKF code and assimilation cycle tools
using the script /letkf_code/make_all.sh 

You may need to edit the make_all.sh script to create a set of enviroment variables
compatible with your compiler and system. A few instances of compiler options and flags
are available as a reference. 

Also, the first time you run this script check that variables MAKE_LETKF, MAKE_PERT_MET_EM and MAKE_WRF_TO_WPS
are all set to true. Binary files will be created and added to tar files located in  /letkf_code/bin/
The tar file name will contain the ENVIRONMENT name (e.g. INTEL_FUGAKU, FUJITSU_FUGAKU or the name
chosen for your custom environment)

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

obs.conf          -> in an assimilation experiment this file contains the path where the observation files are stored. It also contains two bash 
                     arrays, one containing the list of observation types to be assimilated. The other one contains the list of radars to be assimilated
                     (in case of radar data are being assimilated)

More details about the meaning of different flags defined in the configuration files are provided in the configuration files.


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
The script located in wrfutil$VERSION/bin/python/gfs_tools/get_archived_gefs_data.py can be used to downlad this data for a particular
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

There are different ways to run the scripts:
-Submitting a single job that performs all the taskts in the computation nodes ( QUEUESYS in machine.conf = PJM [Fugaku] , PBSSSH [ PBS cluster] , SLURMSSH[ slurm cluster ]
 In this case, select the corresponding option in machine.conf. To submit the job (e.g. the main_FCAST.sh script to perform all the forecasts) use the /wrfutil$VERSION/bin/queue_sub.sh script.
 E.g. ./queue_sub main_FCST.sh  . This will submitt a single job the to queue that will attempt to run all the requested forecasts. 
 If the job is interupted due to a time limit issue, you can resubmitt the job again using the same command. The script will authomatically resume the runs from the last successfull run.
 (Note: if the job could not complete a single run, for example a single forecast initialization, then longer wall times should be allocated for the job before submitting it again to the queue)
 Once the script is runing additional scripts to run individual members in the computation nodes will be created in the $BASEDIR/WRF/ in this folder also numerated folders will be created where the 
 individual ensemble member runs will be performed (refer to the run_WRF.sh script for more details).
 A similar procedure is used for the preparation of initial and boundary conditions using the WPS tools (refer to the run_WPS.sh script for more details). 

-Submitting individual jobs to the computation nodes ( QUEUESYS in machine.conf = PBS_block ):
 TODO

-Running the scripts on a stand alone server (whitout a queue system, QUEUESYS in machine.conf = SSH ):
 TODO


b) Run a simple data assimilation experiment. 





















