#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
import wget

RemoteServer = "https://noaa-gefs-pds.s3.amazonaws.com"
LocalDataDir = '/vol0004/ra000007/data/u10335/PREVENIR/wrfutil_hydra/GFSDATA/'
DataType = 'ENS'  #ENS o DET

Members = ['00','01','02']
Cycles  = ['00']
LeadTimes=['000','003','006']
Dates = ['20240320']

opt_scale = True

#Build download and destination list.

#Format convecsion to GrADS for SCALE (with wgrib2)   
if opt_scale == True : 
   for ymd in Dates :
       for hh in Cycles :    
           datef=ymd[0:4]+"-"+ymd[4:6]+"-"+ymd[6:8]+" "+ hh
           for mem in Members :
               if mem == "00" :
                  ename= mem+'/gec'+mem
               else :
                  ename= mem+'/gep'+mem
              
               tstart=datef
               fhour=LeadTimes[-1]
               gribfile=LocalDataDir + '/gefs.' + ymd + '/' + hh + '/pgrb2b/' + ename + '.t' + hh + 'z'
               outdir= LocalDataDir + '/gefs.' + ymd + '/' + hh + '/grads_scale/' + mem
               os.system('./convert_scale/convert.sh ' +  '"' + tstart + '" ' +  fhour + ' ' + gribfile + ' ' + outdir)
