import numpy as np
import os
import copy
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.colors as mcolors
from matplotlib.colors import from_levels_and_colors
import matplotlib.patches as mpatches
import cartopy
from cartopy import crs
from cartopy.feature import NaturalEarthFeature, COLORS
from netCDF4 import Dataset
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords , vinterp , destagger , extract_times )
import glob
from multiprocessing import Pool

states_provinces = NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none',edgecolor='white')


def plot_loop( conf ) :

#This function standarizes the loop over the plot types, data types,
#dates and ensemble members. 
#This function paralelizes figure generation. We generate a configuration 
#dictionary for each plot and the call the ploting function for each configuration.



   conf_list = []
   #Generate a set of configuration dictionaries (one for each plot)

   for my_type in conf['data_type'] :
      out_path = conf['plot_path'] + my_type + '/'
      os.makedirs( out_path , exist_ok=True )

      #Generamos el archivo con la lat y lon de la reticula
      conf['grid_file'] = out_path + '/lat_lon.npz'
      if ~os.path.isfile( conf['grid_file'] ) :
         #Read some wrf variable for a wrf file to get the grid
         file_list = glob.glob( conf['data_path'] + '/' + my_type + '/*/*' + conf['members'][0]  ) 
         wrf_ds = Dataset( file_list[0]  ) 
         data    = getvar( wrf_ds , "slp" , timeidx=None )
         lats , lons = latlon_coords( data )
         lats = to_np( lats ) ; lons = to_np( lons )
         cart_proj = get_cartopy( data )
         wrf_ds.close()
         np.savez_compressed( conf['grid_file'] , lons = lons , lats = lats , cart_proj = cart_proj , allow_pickle=True )


      for my_mem in conf['members'] :
         file_list = glob.glob( conf['data_path'] + '/' + my_type + '/*/*' + my_mem  )

         for my_file in file_list :

             for my_plot in conf['plot_type'] :

                print(' Ploting data in file : '+ my_file )

                conf['date'] = my_file.split('/')[-2]
                conf['type'] = my_plot
                conf['wrf_file']  = my_file
                conf['fig_file']  = out_path + '/' + conf['type'] + conf['date'] + '_' + my_mem + '.png'
                conf['data_file'] = out_path + '/' + conf['type'] + conf['date'] + '_' + my_mem + '.npz'

                conf_list.append( copy.deepcopy( conf ) )

   #Launch the paralel processing of the figures.
   p = Pool(conf['nworkers'])
   p.map( plot_fig , conf_list )
          


def plot_fig( conf ) :
#Make different plot based on the plot type.

   data = np.load( conf['grid_file'] , allow_pickle=True )
   lats = data['lats']
   lons = data['lons']
   cart_proj = data['cart_proj'].item()

   fig = plt.figure(figsize=(15,10))
   ax1 = fig.add_subplot(1, 1, 1, projection=crs.PlateCarree()  )

   # Agrego las latitudes y longitudes. Con esto de x_inline e y_inline = False hace que 
   # las etiquetas de los ejes no se encuentren en el medio de la figura.
   gl = ax1.gridlines(draw_labels=True, color='black', linestyle='dotted', x_inline = False, y_inline = False)
   gl.top_labels = False
   gl.bottom_labels = True
   gl.right_labels = False  # Oculta etiquetas de latitud en la parte derecha
   gl.left_labels = True  # Muestra etiquetas de latitud en la parte izquierda

   #PLOT TYPE DEPENDENT SECTION

   # ---- CTT ( Cloud Top Temperature ) ---- 
   if conf['type'] == 'ctt' :
   #Get the required data
      if os.path.isfile( conf['data_file'] ) :
         data_dict = np.load( conf['data_file'] , allow_pickle=True )
         data = data_dict['ctt']
      else :
         wrf_ds = Dataset( conf['wrf_file']  )
         data    = getvar( wrf_ds , 'ctt' , units = 'DegC' , timeidx=None )
         data = to_np( data )
         np.savez_compressed( conf['data_file'] , ctt = data )
         wrf_ds.close()

      colormap = get_colormap( 'ctt' )
      contourf = plt.contourf(lons, lats, data , cmap = colormap 
            , transform=crs.PlateCarree()
            , levels=np.arange(-97, 25), extend='both')
      plt.colorbar(contourf , label="T (°C)", shrink = .8)
      ax1.set_title(f'Cloud top temperature: ' + conf['date'] )
      mapedgecolor = 'white'

   # ---- Maximum reflectivity and maximum vertical velocity ----
   elif conf['type'] == 'maxrefw' :

      if os.path.isfile( conf['data_file'] ) :
         data_dict = np.load( conf['data_file'] )
         data  = data_dict['mdbz']
         data2 = data_dict['w'] 
      else :
         wrf_ds = Dataset( conf['wrf_file']  )
         data    = getvar( wrf_ds , 'mdbz' , timeidx=None )
         data2   = np.max( getvar( wrf_ds , 'wa' , timeidx=None ) , 0 )
         data = to_np( data ) ; data2 = to_np( data2 )
         np.savez_compressed( conf['data_file'] , mdbz = data , w = data2 )
         wrf_ds.close()

      colormap = get_colormap( 'dbz' )
      contourf=plt.contourf(lons,lats,data,
               transform=crs.PlateCarree(),
               cmap=colormap , levels=np.arange(5., 75., 5.), extend='both')
      contour=plt.contour(lons,lats, data2 ,
               transform=crs.PlateCarree(),levels=[5,10,20,30])
      plt.colorbar(contourf , label="DbZ", shrink = .8)
      ax1.set_title(f'Ref. + max. W: ' + conf['date'] )
      mapedgecolor = 'black' 

   # Para graficar las costas
   ax1.add_feature(cartopy.feature.COASTLINE, edgecolor = mapedgecolor , linewidth=0.4)
   # Para graficar las fronteras
   ax1.add_feature(cartopy.feature.BORDERS, edgecolor = mapedgecolor , linewidth=0.4)
   ax1.add_feature(states_provinces, edgecolor = mapedgecolor , linewidth=0.4)

   print('Saving figure in : ' + conf['fig_file'] )
   plt.savefig(conf['fig_file'])
   plt.close()



def get_colormap( plottype ) :
#Return the colormap corresponding to the plot type

    if plottype == 'ctt' :
       #Genero una barra de colores personalizada
       # Definir los colores de la paleta "Enhanced IR"
       colors = [
       (0.0, 'darkmagenta'),
       (0.1, 'white'),
       (0.2, 'black'),
       (0.27, 'red'),
       (0.33, 'yellow'),
       (0.37, 'lime'),
       (0.42, 'darkblue'),
       (0.50, 'cyan'),
       (0.55, 'white'),
       (1.0, 'black')]

       # Crear el colormap personalizado
       return mcolors.LinearSegmentedColormap.from_list('Enhanced_IR', colors)

    elif plottype == 'dbz' :
       dbz_levels = np.arange(5., 75., 5.)
       dbz_rgb = np.array([[4,233,231],[1,159,244], [3,0,244],
                [2,253,2], [1,197,1],
                [0,142,0], [253,248,2],
                [229,188,0], [253,149,0],
                [253,0,0], [212,0,0],
                [188,0,0],[248,0,253],
                [152,84,198]], np.float32) / 255.0

       colormap , dbz_norm = from_levels_and_colors(dbz_levels, dbz_rgb, extend='max')
       return colormap

    else :
       print('Error: Not recognized colormap')


    










