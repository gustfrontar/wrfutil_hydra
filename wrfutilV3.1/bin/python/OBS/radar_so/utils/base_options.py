import argparse
import os
from configparser import ConfigParser
import shlex
import numpy as np

class BaseOptions():
    """
       Configuración del control de calidad de radar

       A través de esta clase se puede acceder a la lista de archivos a procesar, se
       especifican los radares que se van a controlar, y las fechas que se van a tener
       en cuenta.
       También se especifica el nombre de los campos y la configuración de los filtros.

    """

    def __init__(self):
        """Reset the class; indicates the class hasn't been initialized"""
        self.initialized = True
        self.dataroot = os.environ['SO_DATAROOT']
        self.instrument_list = os.environ['SO_INSTRUMENT_LIST'].split(',')
        self.output_path = os.environ['SO_OUTPUTROOT']

        self.so_grid = [float(x) for x in os.environ['SO_GRID'].split(',')]
        self.so_vars = [os.environ['SO_VARS'].split(',')]
        self.var_opt = [float(x) for x in os.environ['SO_CREF_OPT'].split(',')]
        self.reference_date = os.environ['SO_REFERENCE_DATE']

        self.window_length = int(os.environ['SO_WINDOW_LENGTH'])
        self.window_type = 'centered'
        self.min_nobs = int(os.environ['SO_MIN_NOBS'])


    def print_options(self,opt):
        """ Print and save options """
        message = ''
        message += '----------------- Options ---------------\n'
        for k, v in sorted(vars(opt).items()):
            comment = ''
            default = self.parser.get_default(k)
            if v != default:
                comment = '\t[default: %s]' % str(default)
            message += '{:>25}: {:<30}{}\n'.format(str(k), str(v), comment)
        message += '----------------- End -------------------'
        print(message)
