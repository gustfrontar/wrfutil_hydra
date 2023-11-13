import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors
import matplotlib.patches as mpatches
#from cartopy import crs
#from cartopy.feature import NaturalEarthFeature, COLORS
from netCDF4 import Dataset
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords , vinterp , destagger , extract_times )
import glob


def ploter_loop( conf , plot_types ) :

 #Plot forecasts.   
 if 'FCST' in conf['exp_type']  :
   base_data_dir = conf['histdir'] + '/FCST/' 
   attribs=dict()
   attribs['exp_type']='FCST'

   #Get the list of date paths
   date_list = glob.glob( base_data_dir + '/*' )
   for my_date in date_list  :
       attribs['date'] = os.path.basename( my_date ) 
       print( 'Ploting date ', attribs['date']  )

       #Get the list of ensemble member paths
       mem_list = glob.glob( my_date + '/*' )
       for my_mem in mem_list :
          attribs['mem'] = os.path.basename( my_mem ) 
          print( 'Ploting member ', attribs['mem'] )
          file_list = glob.glob( my_mem + '/wrfout*' )
          for my_file in file_list :
            attribs['file'] = os.path.basename( my_file ) 
            attribs['file_path'] = my_file
            print('Ploting file ', attribs['file'] )
         
            for my_type in conf['plot_type_list'] :
                get_figure( conf , plot_types[ my_type ] , attribs )
 
 #Plot analysis 
 if 'ANAL' in conf['exp_type']  :
   base_data_dir = conf['histdir'] + '/ANAL/'
   attribs=dict()
   attribs['exp_type']='ANAL'

   #Get the list of date paths
   date_list = glob.glob( base_data_dir + '/*' )
   for my_date in date_list  :
      attribs['date'] = os.path.basename( my_date )
      print( 'Ploting date ', attribs['date']  )

      #Get the list of ensemble member files
      file_list = glob.glob( my_date + '/anal*' )
      file_list.append( my_date + '/guesemean' )  #Add the gues ensemble mean to the list. 
      for my_file in file_list :
         attribs['mem'] = os.path.basename( my_file )[4:]
         print( 'Ploting member ', attribs['mem'] )
         attribs['file'] = os.path.basename( my_file )
         attribs['file_path'] = my_file
         print('Ploting file ', attribs['file'] )

         for my_type in conf['plot_type_list'] :
             get_figure( conf , plot_types[ my_type ] , attribs )

                


def get_figure( conf , plot_type , attribs ) :

   if plot_type['pvar1'] is not None :  #Lets plot variable 1 as shaded.

      ncfile = Dataset(attribs['file_path'])
      pvar1 = getvar(ncfile,plot_type['pvar1'])
      if plot_type['pvar1vlev'] is not None :
         if plot_type['pvar2'] == 'W' :
            pvar2=destagger(pvar2,0)
         pvar1 = vinterp( ncfile , field=pvar1 , vert_coord='ght_msl',interp_levels=plot_type['pvar1vlev'],extrapolate=True,log_p=True)

      lat , lon = latlon_coords( pvar1 ) 
      pvar1 = to_np( pvar1 )

   if plot_type['pvar2'] is not None :  #Lets plot variable 2 as countour 
      if ~ ncfile._isopen  :
        ncfile = Dataset(attribs['file_path'])
      pvar2 = getvar(ncfile,plot_type['pvar2'])
      lat , lon = latlon_coords( pvar2 )
      
      if plot_type['pvar2vlev'] is not None :
         if plot_type['pvar2'] == 'W' : 
            pvar2=destagger(pvar2,0)
         pvar2 = vinterp( ncfile , field=pvar2 , vert_coord='ght_msl',interp_levels=plot_type['pvar2vlev'],extrapolate=True,log_p=True) 
      pvar2 = np.squeeze( to_np( pvar2 ) )

   lon = to_np( lon ) ; lat = to_np( lat )
   plot_date = np.datetime_as_string( extract_times( ncfile , timeidx=0 , meta=False ) , unit='s')    #os.path.basename( attribs['file_path'] )[11:]  
   print( plot_date )
   plot_title=''
   file_name = str( plot_date ) + '.png'
   fig = plt.figure(figsize=(8,6))
   #ax = plt.axes(projection=cart_proj) #create axes
   # Download and add the states and coastlines
   #states = NaturalEarthFeature(category="cultural", scale="10m",
   #                                      facecolor="none",name="admin_1_states_provinces_shp")
   #ax.add_feature(states, linewidth=.5, edgecolor="black")
   #ax.coastlines('50m', linewidth=0.8)
   if plot_type['pvar2'] is not None :  #Lets plot variable 1 as shade
      print( plot_type['pvar1'] + ':' , pvar1.min() , pvar1.max() )
      plt.contourf( lon , lat , pvar1 ,  cmap=plot_type['pvar1cmap'] , levels=plot_type['pvar1clev'] , extend='both')
      # Add a color bar
      plt.colorbar(shrink=.98)
      if plot_type['pvar1vlev'] is not None :
         file_name =  str( plot_type['pvar1vlev'][0] ) + '_' + file_name
         plot_title = '@' + str( plot_type['pvar1vlev'][0] ) + plot_title 
      file_name = plot_type['pvar1'] + '_' + file_name
      plot_title = plot_type['pvar1'] + ' ' + plot_title 
   if plot_type['pvar2'] is not None :  #Lets plot variable 2 as countour
      
      print( plot_type['pvar2'] + ':' , pvar2.min() , pvar2.max() )
      plt.contour( lon , lat , pvar2 , levels=plot_type['pvar2clev'] )
      if plot_type['pvar2vlev'] is not None :
         file_name =  str( plot_type['pvar2vlev'][0] ) + '_' + file_name
         plot_title = '@' + str( plot_type['pvar2vlev'][0] )+ ' ' + plot_title 
      file_name = plot_type['pvar2'] + '_' + file_name
      plot_title = plot_type['pvar2'] + ' ' + plot_title  
   # Set the map bounds
   plt.xlim((np.min(lon),np.max(lon)))
   plt.ylim((np.min(lat),np.max(lat)))
   # Add the gridlines
   plt.grid()
   file_name=file_name.replace(':','_')
   file_name=file_name.replace('T','_')
   file_path=conf['plotdir'] + '/' + attribs['exp_type'] + '/' + attribs['date'] + '/' + attribs['mem'] + '/' 
   os.makedirs(file_path , exist_ok=True)
   plt.title( plot_title )
   plt.savefig( file_path + file_name , dpi=None, facecolor='w',edgecolor='w' )
   plt.close()
