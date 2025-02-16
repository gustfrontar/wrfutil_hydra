import numpy as np
import mod_plot as mp

conf=dict()
conf['data_type']  = ['ANAL','GUES']
conf['members']    = ['00001','emean']
conf['data_path']  = '../../../HIST/'
conf['plot_path']  = '../../../PLOT/'
conf['plot_type'] = ['ctt','maxrefw']  #['ctt','maxrefw']
conf['nworkers']  = 20

###########################################################


mp.plot_loop( conf )


