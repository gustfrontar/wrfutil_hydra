#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 16:58:28 2023

@author: jruiz

"""

import module_pert_met_em as mpm
import numpy as np
import os
from multiprocessing import Pool

np.random.seed(10) #Fix the random seed (generated perturbation weigths are fixed)

conf = dict()
#Get value of environment variables. 
conf['BasePathOri'] = __BASE_PATH_ORI__                  #The base path for the original ensemble member met_em files.
conf['BasePathOut'] = __BASE_PATH_OUT__
conf['TargetEnsembleSize']= __TARGET_ENSEMBLE_SIZE__     #The number of ensemble member in the target ensemble.
conf['ActualEnsembleSize']= __ACTUAL_ENSEMBLE_SIZE__     #The number of ensemble members in the original ensemble.
conf['PertAmp']= __PERT_AMP__                            #The perturbation randomization amplitud (for RejuGauss and RejuUniform )
conf['PertType'] = __PERT_TYPE__                         #The type of perturbation transformation (RejuGauss , RejuUniform , Specular 
conf['VarList'] = __VAR_LIST__                           #List containing the variable names to be perturbed. 
conf['NetCDF4'] = __NETCDF4__                            #If TRUE, NETCDF4 output will be forced independently of the input format.
                                                         #if FALSE, then the output format will be the same as the input format.
conf['UseCp'] = True                                     #Use CP to create the new met_em files before computing the perturbations. 
conf['ReduceMemUse'] = True                              #Do not store the expaneded ensemble in memory. 
conf['FileSearchString'] = 'met_em*.nc*'

#==============================================================================
#Get the file lists.
#EnsFileList is a list of list. There are as many list as files for each ens member (times x procs).
#each list contains as many elements as ensemble members.
#TargetEnsFileList is a list of list with the same structure as the EnsFileList but corresponding to the 
#target ensemble. This files will be created later.  
EnsFileList , TargetEnsFileList = mpm.get_file_lists( conf )

conf = mpm.get_netcdf_type( conf , EnsFileList[0][0] )
print( 'The base data model is ',conf['DataModel'] )


#Get the perturbation trnasform matrix
Tm = mpm.get_pert_matrix( conf )
    
#For each met_em file (this includes times and eventually procs) read, transform and expand, and write
#the ensemble perturbations.
def pert_file_type( my_ens_file_list , my_target_ens_file_list , my_Tm , conf ) :
     #Create ensemble files for the output ensemble
     print('Processing file type: ', my_ens_file_list[0] )
     template_file = my_ens_file_list[0]
     mpm.create_ens_files( my_target_ens_file_list , template_file , conf )
     for my_var in conf['VarList']:
         print('Processing variable :',my_var)
         print('Reading ...')
         my_ens = mpm.read_ens( my_ens_file_list , my_var , conf )                #Read
         if conf['ReduceMemUse'] :
            print('Transforming and writing ...')
            mpm.trans_and_write_ens( my_ens , my_target_ens_file_list , my_var , Tm , conf , new_mean = None )       #Transform and write.
         else   :
            print('Transforming ...')
            my_ens = mpm.trans_ens( my_ens , my_Tm , conf , new_mean = None )        #Transform
            print('Checking ...')
            my_ens = mpm.check_ens( my_ens , my_var )                                #Check
            print('Writing ...')
            mpm.write_ens( my_target_ens_file_list , my_var , my_ens , conf )        #Write
         

     #for my_var in conf['VarList']:
     #    my_ens = mpm.read_ens( my_ens_file_list , my_var , conf )                #Read
     #    my_ens = mpm.trans_ens( my_ens , my_Tm , conf , new_mean = None )        #Transform
     #    my_ens = mpm.check_ens( my_ens , my_var )                                #Check
     #    mpm.write_ens( my_target_ens_file_list , my_var , my_ens , conf )        #Write

THREADS_NUM=int(os.environ['OMP_NUM_THREADS'])
pool = Pool(processes=THREADS_NUM)
pool.starmap(pert_file_type,[(EnsFileList[p],TargetEnsFileList[p],Tm,conf) for p in range(len(EnsFileList))])

