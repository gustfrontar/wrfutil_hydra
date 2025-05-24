# -*- coding: utf-8 -*-
#################################################
# Author: DMSR, Servicio Meteorologico Nacional #
# Create: 02/2022 - P. Maldonado                #
#################################################
import os, re
from datetime import datetime, timedelta

#################
# GET FUNCTIONS #
#################
def get_fcst_dates(start_date, length):
   ''' start_date: str, length: int in minutes'''

   fini = datetime.strptime(start_date, '%Y-%m-%d_%H:%M:%S')
   fend = fini + timedelta(minutes=int(length))

   ini = fini.strftime('%Y-%m-%d_%H:%M:%S')
   end = fend.strftime('%Y-%m-%d_%H:%M:%S')
   lead = length//60

   return ini, end, lead

def get_exp_names(filename, type_):

   type_ = type_.upper()

   if 'anal' in filename:
      sname = 'ANAL'
      lname = 'Analysis'
      if 'MEMBER' == type_:
         type_ = filename[-4:]
      flength = 0
   elif any(x in filename for x in ('gues', 'GUES')):
      sname = 'GUES'
      lname = 'First Guess'
      if 'MEMBER' == type_:
         type_ = (re.search('/(\d{2})/', filename).group(1)).zfill(4)
      #### flength = int(os.environ['ANALISIS'])
         flength=3600
   else:
      print('Error: filename not set')
      sname, lname, type_, flength = None, None, None, None

   # DET
   #   sname = ''
   #   lname = 'Deterministic'
   #   type_ = (re.search('/(\d{2})/', filename).group(1)).zfill(4)
   #   flength = int(os.environ['PLAZO'])

   # ENS
   #   sname = ''
   #   lname = 'Ensemble'
   #   type_ = (re.search('/(\d{2})/', filename).group(1)).zfill(4)
   #   flength = int(os.environ['PLAZO'])

   return sname, lname, type_, flength

