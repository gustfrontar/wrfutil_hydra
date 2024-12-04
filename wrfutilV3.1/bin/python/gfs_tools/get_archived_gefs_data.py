#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
import wget

RemoteServer = "https://noaa-gefs-pds.s3.amazonaws.com"
LocalDataDir = '/home/jruiz/salidas/GFSDATA/'
DataType = 'ENS'  #ENS o DET

Members = ['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30']
Cycles  = ['06','12','18']
LeadTimes=['000','003','006','009','012','015','018','021','024','027','030','033','036']
Dates = ['20240319']


#Build download and destination list.

if DataType == 'ENS'  :
    
    file_origa_list = []
    file_origb_list = []
    file_dest_list  = []

    for my_date in Dates :
        
        for my_cycle in Cycles :
            
            for my_member in Members :
                if my_member == '00' :
                    prefix = 'gec'
                else                 :
                    prefix = 'gep'
                
                for my_lead in LeadTimes :
                   
                    file_origa_list.append( RemoteServer + '/gefs.' + my_date + '/' + my_cycle + '/atmos/pgrb2ap5/' + prefix + my_member + '.t' + my_cycle + 'z.pgrb2a.0p50.f' + my_lead  )
                    file_origb_list.append( RemoteServer + '/gefs.' + my_date + '/' + my_cycle + '/atmos/pgrb2bp5/' + prefix + my_member + '.t' + my_cycle + 'z.pgrb2b.0p50.f' + my_lead  )
                    file_dest_list.append( LocalDataDir + '/gefs.' + my_date + '/' + my_cycle + '/pgrb2b/0' + my_member + '/' )
                    

print('A total number of ',len(file_origa_list),' files will be downloaded')


#We download the data.
for ifile , my_filea in enumerate( file_origa_list ) :
    print('Downloading file ',my_filea )
    file_desta = file_dest_list[ifile] + '/' + os.path.basename( my_filea )  
    #Check wether the file is already there.
    if not os.path.isfile( file_desta ) :
       os.makedirs( file_dest_list[ifile] , exist_ok=True)
       filename = wget.download( my_filea , out= file_dest_list[ifile] )
    my_fileb = file_origb_list[ifile] 
    file_destb = file_dest_list[ifile] + '/' + os.path.basename( my_fileb )
    if not os.path.isfile( file_destb ) :
       os.makedirs( file_dest_list[ifile] , exist_ok=True)
       filename = wget.download( my_fileb , out= file_dest_list[ifile] )
    #Merge the two files   
    merged_file = file_dest_list[ifile] + '/' + os.path.basename( my_fileb ).replace('pgrb2b','pgrb2')
    if not os.path.isfile( merged_file ) : 
       os.system('cat ' + file_desta + ' ' + file_destb + ' > ' + merged_file ) 
    
    
    
    
