# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import netCDF4 as nc
import glob
import os
import sys
import gc

def get_file_lists ( conf ) :
    #Generate the list with the original ensemble files and the target ensemble files.
    #Returns a lists of lists, one list for each time containing all ensemble members. 
    BasePathOri = conf['BasePathOri']
    BasePathOut = conf['BasePathOut']
    Nt = conf['TargetEnsembleSize']
    Ne = conf['ActualEnsembleSize']
    FileSearchString = conf['FileSearchString']

    #Get the list of met_em files.
    FileList = glob.glob( BasePathOri + '/01/' + FileSearchString )
    FileList.sort()

    #Generate the list of original met_ems
    EnsFileList = []
    for my_file in FileList :
        my_list = []
        for iens in range(Ne) :
           file_path = BasePathOri + '/'+ "{:02d}".format(iens+1)  +'/' + os.path.basename(my_file)
           if os.path.isfile( file_path ) :
              my_list.append( file_path )
           else : 
              print('Error: File not found ' + file_path )
              sys.exit(1)
        
        EnsFileList.append( my_list )
            
    #Generate the list of output met_ems
    TargetEnsFileList = []
    for my_file in FileList :
        my_list = []
        for iens in range(Nt) :
           my_list.append( BasePathOut + '/'+ "{:02d}".format(iens+1)  +'/' + os.path.basename(my_file) )
        TargetEnsFileList.append( my_list )
        
    return EnsFileList , TargetEnsFileList 

def get_netcdf_type( conf , TestFile ) :

    my_ds = nc.Dataset( TestFile , 'r' )

    conf['DataModel'] = my_ds.data_model

    return conf

def read_ens( file_list , variable , conf ) :
    #file_list a list containing the full path to the files belonging to the ensemble for a particular time.
    #variable is the variable to be accessed
    #conf the general configuration dictionary. 
    
    ens_size = len( file_list )
    
    for ifile , my_file  in enumerate( file_list ) :
        my_ds = nc.Dataset( my_file , 'r' )
        my_var = my_ds[variable][:]
        if ifile == 0 :            
           my_var_shape = my_var.shape
           ens = np.zeros(np.append(my_var_shape,ens_size))
        ens[...,ifile] = my_var
        
        my_ds.close()
    gc.collect()    
    return ens

def write_ens( file_list , variable , ens , conf ) :
    
    for ifile , my_file in enumerate( file_list ) :

        my_ds = nc.Dataset( my_file , 'r+' )
        my_ds[variable][:] = ens[...,ifile]
        my_ds.close()
    

def create_ens_files( file_list , template_file , conf ) :

    if conf['UseCp'] : 
       #Use system copy to generate the new files. Note that netcdf type will be inherited from the original netcdf files. 
       for my_file in file_list :
          #Create the file directory if it does not exist.
          os.makedirs( os.path.dirname( my_file ) , exist_ok=True) 
          copycmd = 'cp ' + template_file + ' ' + my_file 
          os.system(copycmd)
    else   :
       #This function creates all the files specified in file_list using template_file as a template.
       template_ds = nc.Dataset( template_file )
       for my_file in file_list :
          #Create the file directory if it does not exist.
          os.makedirs( os.path.dirname( my_file ) , exist_ok=True)
          #Create the file and transfer template information into the file.
          if conf['NetCDF4'] :
             my_ds = nc.Dataset( my_file , 'w' , format='NETCDF4' )            #Force NETCDF4 Output (independent of input data format)
          else               :
             my_ds = nc.Dataset( my_file , 'w' , format=conf['DataModel'] )    #Use the input data format as output data format.

          my_ds.setncatts( template_ds.__dict__)
          for dname , the_dim in template_ds.dimensions.items():
             my_ds.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
          for v_name, varin in template_ds.variables.items():
             if conf['NetCDF4'] or conf['DataModel'] == 'NETCDF4' :
                my_var = my_ds.createVariable(v_name, varin.datatype, varin.dimensions , compression='zlib' )
             else                                                 :
                my_var = my_ds.createVariable(v_name, varin.datatype, varin.dimensions )
             my_var.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
             #if v_name not in conf['VarList'] : 
             my_var[:] = varin[:]
          my_ds.close()
       template_ds.close()
    gc.collect()   

    
def get_pert_matrix( conf ) :
    
    #This routine generates a transformation matrix that produce the desired perturbation.
    Nt = conf['TargetEnsembleSize']
    Ne = conf['ActualEnsembleSize']
    
    
    if conf['PertType'] == 'RejuGauss' :
        Tm = np.random.randn(Ne,Nt)
        Tm = np.matmul( Tm , (np.eye(Nt) - (1.0/Nt)*np.ones((Nt,Nt)) ) )
        Tm = np.sqrt( conf['PertAmp'] / (Nt-1) ) * Tm
        iens = 0
        for ii in range( Nt ) :
             Tm[iens,ii] = Tm[iens,ii] + 1.0
             iens = iens + 1
             if iens == Ne :
                 iens = 0  
    
    elif conf['PertType'] == 'RejuUniform' :
        #All elements of the transformation matrix are positive. 
        Tm = np.random.rand(Ne,Nt)
        Tm = np.matmul( Tm , (np.eye(Nt) - (1.0/Nt)*np.ones((Nt,Nt)) ) )
        Tm = np.sqrt( conf['PertAmp'] / (Nt-1) ) * Tm
        #We force all elements of the transformation matrix to be possitive. 
        #The mean perturbation should still be equal to 0. 
        Tm = Tm - np.min(Tm)
        iens = 0
        for ii in range( Nt ) :
             Tm[iens,ii] = Tm[iens,ii] + 1.0
             iens = iens + 1
             if iens == Ne :
                 iens = 0  
    
    elif conf['PertType'] == 'Specular' :
         if Nt != 2 * Ne :
             print('Warning: This method is optimal when the target ensemble size is equal to the actual ensemble size')
         Tm = np.zeros((Ne,Nt))
         iens = 0
         amp = 1.0
         for ii in range( Nt ) :
             Tm[iens,ii] = amp
             iens = iens + 1
             if iens == Ne :
                 iens = 0
                 amp = -1.0 * amp

    return Tm


# Transforms the ensemble perturbation according to the selected mode.
def trans_ens( ens , Tm , conf , new_mean = None ) :
    #ens a numpy array containing data for each ensemble member. 
    #the ensemble dimension is always the last dimension.
    #new_mean is for the case we would like to recenter the ensemble mean.
    
    ens_shape = ens.shape
    ens_size =  ens_shape[-1]              #The last dimension of ens defines the original ensemble size.
    ens_size_target = np.shape( Tm )[1]    #The columns of Tm defines the size of the target ensemble.
    ens_shape_target = np.array(ens_shape)
    ens_shape_target[-1] = ens_size_target    
    #Reshape the ensemble (files are different variables, columns are ensemble members)
    ens = ens.reshape( np.prod( ens_shape[0:-1] ) , ens_size )
    #Get the ensemble mean and the perturbations.
    ens_mean = np.mean( ens , axis = 1 )
    ens = ens - np.tile( ens_mean , ( ens_size , 1 ) ).T 
    #Transform the perturbations given the transformation matrix
    ens = np.matmul( ens , Tm )
    if new_mean is not None :
       ens = ens + np.tile( new_mean , ( ens_size_target , 1 ) ).T
    else :
       ens = ens + np.tile( ens_mean , ( ens_size_target , 1 ) ).T
    #Reshape the ensemble to recover its original strucutre.
    ens = ens.reshape( ens_shape_target )
    return ens



# Transform the ensemble computing one new memeber at at time to be more memory efficient.
def trans_and_write_ens( ori_ens , my_target_file_list , my_var , Tm , conf , new_mean = None ) :
    #ens a numpy array containing data for each ensemble member. 
    #the ensemble dimension is always the last dimension.
    #new_mean is for the case we would like to recenter the ensemble mean.

    ori_ens_shape = ori_ens.shape
    var_shape  = ori_ens.shape[0:-1]
    ori_ens_size = ori_ens_shape[-1]  #The number of original ensemble members.

    print('Reshaping the ensemble')
    ori_ens = ori_ens.reshape( np.prod( var_shape ) , ori_ens_size )             #Reshape ensemble members (1 member per column of a matrix)
    ori_ens_mean = np.mean( ori_ens , axis = 1 )                                 #Compute the ensemble mean
    ori_ens = ori_ens - np.tile( ori_ens_mean , ( ori_ens_size , 1 ) ).T         #Substract the mean to obtain the ensemble perturbations.

    if new_mean is not None :
        new_mean = new_mean.reshape( np.prod( var_shape , 1 ) )                  #Reshape the recentering mean (if any)

    for imember , file_member in enumerate( my_target_file_list ) :
        print('Processing member ',imember+1,' to be written at ',file_member )
        if imember == 0 :
           print('New member weigths ',Tm[:,imember] ) 
        #Generate the new member
        my_member = np.matmul( ori_ens , Tm[:,imember] )                  #Compute the new perturbations according to the selected method.
        if new_mean is not None :
           my_member = my_member + new_mean                               #Recenter the new perturbations around a new mean
        else :
           my_member = my_member + ori_ens_mean                           #Recenter the new perturbations around the original mean
        #Check the new member
        my_member = check_ens( my_member , my_var )                       #Do some sanity checks.

        #Write the new member to the corresponding file.
        print('Writing the file')
        my_member = my_member.reshape( var_shape )
        my_ds = nc.Dataset( file_member , 'r+' )
        my_ds[my_var][:] = np.copy( my_member )
        my_ds.close()
        gc.collect()

    return 

def check_ens( ens , variable ) :
    #Variable specific checks (enforce physically meaningful ranges)
    if variable == 'HR' :
        ens[ens < 0.0] = 0.0
        ens[ens > 1.0] = 1.0
    if variable == 'SPECHUMD' :
        ens[ens < 0.0] = 0.0
    if variable == 'QR' :
        ens[ens < 0.0] = 0.0
    if variable == 'QG' :
        ens[ens < 0.0] = 0.0
    if variable == 'QS' :
        ens[ens < 0.0] = 0.0
    if variable == 'QC' :
        ens[ens < 0.0] = 0.0
    if variable == 'QI' :
        ens[ens < 0.0] = 0.0
    if variable == 'QH' :
        ens[ens < 0.0] = 0.0
        
    return ens


