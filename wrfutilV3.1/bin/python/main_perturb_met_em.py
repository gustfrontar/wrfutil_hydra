#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 16:58:28 2023

@author: jruiz

"""

import module_pert_met_em as mpm
import numpy as np
import os

np.random.seed(10) #Fix the random seed (generated perturbation weigths are fixed)

conf = dict()
#Get value of environment variables. 
#The base path for the original ensemble member met_em files.
conf['BasePathOri'] = __BASE_PATH_ORI__

conf['BasePathOut'] = __BASE_PATH_OUT__

#The number of ensemble member in the target ensemble.
conf['TargetEnsembleSize']= __TARGET_ENSEMBLE_SIZE__
#The number of ensemble members in the original ensemble.
conf['ActualEnsembleSize']= __ACTUAL_ENSEMBLE_SIZE__
#The perturbation randomization amplitud (for RejuGauss and RejuUniform )
conf['PertAmp']= __PERT_AMP__
#The type of perturbation transformation (RejuGauss , RejuUniform , Specular )
conf['PertType'] = __PERT_TYPE__
#TODO check variable list.
#The list of variables that whose perturbations will be transformed.
conf['VarList'] = __VAR_LIST__

#["PRES", "SM", "ST", "GHT", "SKINTEMP",
#           "ST100200", "ST040100", "ST010040", "ST000010",
#           "SM100200", "SM040100", "SM010040", "SM000010",
#           "PSFC", "RH", "UU", "VV", "TT", "PMSL" ]

conf['FileSearchString'] = 'met_em*.nc*'

#==============================================================================
#Get the file lists.
#EnsFileList is a list of list. There are as many list as files for each ens member (times x procs).
#each list contains as many elements as ensemble members.
#TargetEnsFileList is a list of list with the same structure as the EnsFileList but corresponding to the 
#target ensemble. This files will be created later.  
EnsFileList , TargetEnsFileList = mpm.get_file_lists( conf )


#Get the perturbation trnasform matrix
Tm = mpm.get_pert_matrix( conf )
    
#For each met_em file (this includes times and eventually procs) read, transform and expand, and write
#the ensemble perturbations.
#TODO this loop can be paralelized over EnsFileList (times x procs)



for  ilist , file_list_in in enumerate( EnsFileList ) :

     file_list_out = TargetEnsFileList[ilist]    
     #Create ensemble files for the output ensemble
     template_file = EnsFileList[ilist][0]
     mpm.create_ens_files( file_list_out , template_file , conf )
    
     for my_var in conf['VarList']:
         my_ens = mpm.read_ens( file_list_in , my_var , conf )          #Read
         my_ens = mpm.trans_ens( my_ens , Tm , conf , new_mean = None ) #Transform
         my_ens = mpm.check_ens( my_ens , my_var )                      #Check
         mpm.write_ens( file_list_out , my_var , my_ens , conf )        #Write





