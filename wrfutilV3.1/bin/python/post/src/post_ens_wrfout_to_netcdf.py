# -*- coding: utf-8 -*-
#################################################
# Author: DMSR, Servicio Meteorologico Nacional #
# Create: 02/2022 - P. Maldonado                #
#################################################
import sys, os, time, re, glob
import numpy as np
import xarray as xr
import models_io.util_post as utpost
import util_common as utcommon
import models_io.catalog_variables as catalog_vars


def main():
   tini = time.time()
   print('--- Hello from {} ---'.format(os.path.basename(__file__)))

   # Parse input parameters 
   FILELIST = sys.argv[1].split(',')

   PATHOUT = sys.argv[2]
 
   MEMBER=sys.argv[3]

   OUTDATE=sys.argv[4]

   # Set variables
   '''
   DIMS = {'3D': {'vars': {'Q': 'PRESSURE', 
                           'GEOPT': 'PRESSURE', 
                           'T': 'PRESSURE', 
                           'Umet': 'PRESSURE', 
                           'Vmet': 'PRESSURE', 
                           'W': 'PRESSURE', 
                           'QRAIN': 'PRESSURE', 
                           'QGRAUP': 'PRESSURE', 
                           'QSNOW': 'PRESSURE', 
                  'levs': {'PRESSURE': np.array([1000., 975., 950., 925., 900., 850., 800., 750., 700., 650., 600., 550., 500., 400., 300., 250., 200., 150.])}
   		},
           '2D': {'vars': {'PSFC': None, 
                           'SLP': None, 
                           'T2': None, 
                           'TD2': None, 
                           'Q2': None, 
                           'PP': None, 
                           'DBZ': 'Z_AGL', 
                           'MDBZ': None, 
                           'MCAPE': None, 
                           'CIN': None, 
                           'Umet10': None, 
                           'Vmet10': None, 
                           'GRAUPELNC': None, 
                           'LCL': None, 
                           'LFC': None, 
                           'SRH_1000': None, 
                           'SRH_3000': None, 
                           'TPE': 'PRESSURE'},
                #   'levs': {'Z_AGL': np.arange(1000., 11000., 1000.)}
                   'levs': {'Z_AGL': np.arange(1000., 11000., 1000.), 'PRESSURE': np.array([850.])}
          		}
   	}

   '''

   DIMS = {'3D': {'vars': {'Q': 'Z_AGL',
                          'PRESSURE': 'Z_AGL',
                          'T': 'Z_AGL',
                          'Umet': 'Z_AGL',
                          'Vmet': 'Z_AGL',
                          'W': 'Z_AGL',
                          'QCLOUD': 'Z_AGL',
                          'QRAIN': 'Z_AGL',
                          'QGRAUP': 'Z_AGL',
                          'QSNOW': 'Z_AGL',
                          'DBZ': 'Z_AGL'},
                 'levs': {'Z_AGL': np.array([250.,500.,1000.,2000.,3000.,4000.,5000.,6000.,8000.,10000.,12000.,14000.])} 
                },
           '2D': {'vars': {'MCAPE': None,
                          'PSFC': None,
                          'T2': None,
                          'CIN': None,
                          'Q2': None,
                          'PP': None,
                          'MDBZ': None,
                          'UH':None,
                          'Umet10': None,
                          'Vmet10': None},
                 'levs': {'Z_AGL': None}}
         }

#   '''

   # Get experiment type properties
   #sname, lname, member, flength = utcommon.get_exp_names(FILENAME, 'member')
   #exp = {'title': None, 'member': member, 'flength': flength}
   #print(lname, exp['member'])
   print(f"desde el python {FILELIST}")

   #filelist = sorted(glob.glob(FILENAME))
   #print(f"lista de archivos encontrados {filelist}")

   for dim in DIMS: 
      #exp['title'] = ('{} {}').format(lname, dim)

      # Create SMN dataset
      ds = utpost.smn_netcdf_from_list(FILELIST, DIMS[dim]['vars'], catalog_vars, DIMS[dim]['levs'], 'WRF')

      # Compress data
      ds = utpost.lossy_compression(ds, 12)

      # Save to netcdf 
      fileout = '{}/WRF.D{}.T{}.M{}.nc'.format(PATHOUT, OUTDATE , dim, MEMBER)
      ds.to_netcdf(fileout)

   tiempo = float('{:.4f}'.format(time.time()-tini))
   print('Execution Time: {} seconds'.format(tiempo))


if __name__ == "__main__":
   main()
