#print __doc__

# Author: Juan Ruiz (jruiz@cima.fcen.uba.ar)
# License: BSD 3 clause

# History: Juan Ruiz created 2017.

def main_qc( filename , options ) :

   import numpy as np
   import matplotlib.pyplot as plt
   import pyart
   import netCDF4

#  import qctools

   order_ref=False
   order_v  =False

   radar = pyart.io.read(filename)

   #Get the nyquist velocity
   nyquistv=radar.get_nyquist_vel(0,check_uniform=True)

   #===================================================
   # DEALIASING 
   #===================================================

   if options['ifdealias'] : 

     #Uso una de las funciones de dealiasing de pyart con los parametros por defecto
     winddealias=pyart.correct.region_dealias.dealias_region_based(radar, interval_splits=options['interval_splits'], interval_limits=None, skip_between_rays=options['skip_between_rays'], skip_along_ray=options['skip_along_ray'], centered=True, nyquist_vel=None, check_nyquist_uniform=True, gatefilter=False, rays_wrap_around=True, keep_original=False, set_limits=True, vel_field='V', corr_vel_field=None)

     #Add wind dealias to the radar objetc.
     radar.fields['Vda'] = winddealias
     radar.fields['Vda']['coordinates']=radar.fields['V']['coordinates']
     radar.fields['Vda']['units']=radar.fields['V']['units']
     radar.fields['Vda']['long_name']=radar.fields['V']['long_name']
     radar.fields['Vda']['standard_name']=radar.fields['V']['standard_name']


   #===================================================
   # RHO HV FILTER
   #===================================================

   if options['rhofilter']   :



   #===================================================
   # ECHO TOP FILTER
   #===================================================


   if options['ifetfilter']  :
     if ( ~ order_ref ) : 
        [ ref_array , ref_az , ref_level , ref_time , ref_az_exact ]=order_variable( radar , options['reflectivity_var_name'] )
        order_ref=True

     [ echo_top .... ]=compute_echo_top( ref_array , options ) 

   return radar

   

def order_variable ( radar , var_name )  :  

   import numpy as np

   #order_var es la variable ordenada con los azimuths entre 0 y 360 (si hay rayos repetidos se promedian).
   #order_azimuth es el azimuth "aproximado" utilizando 0 como azimuth inicial y avanzando en intervalos regulares e iguales a la resolucion
   #del azimuth en grados.
   #levels son los angulos de elevacion.
   #azimuth_exact es un array que contiene para cada nivel el valor exacto del azimuth que corresponde a cada rayo. 

   ray_angle_res = np.unique( radar.ray_angle_res['data'] )
   if( np.size( ray_angle_res ) >= 2 )  :
      print('Warning: La resolucion en azimuth no es uniforme en los diferentes angulos de elevacion ')
      print('Warning: El codigo no esta preparado para considerar este caso y puede producir efectos indesaedos ')
   ray_angle_res=np.nanmean( ray_angle_res )
   print 'The resolution in azimuth is: %5.3f' % ( ray_angle_res )


   levels=np.unique(radar.elevation['data'])
   azimuth=radar.azimuth['data']
   time=radar.time['data']

   order_azimuth=np.arange(0.0,360.0-ray_angle_res,ray_angle_res) #Asuming quasi regular azimuth location
   naz=np.size(order_azimuth)
   nel=np.size(levels)

   var=radar.fields[var_name]['data']

   nr=var.shape[1]

   #Allocate arrays
   order_var=np.zeros((naz,nr,nel))
   order_time=np.zeros((naz,nel)) 
   azimuth_exact=np.zeros((naz,nel))

   order_var[:]=np.nan
   order_time[:]=np.nan
   azimuth_exact[:]=np.nan

   

   for ilev in range(0, nel) :

      levmask= radar.elevation['data'] == levels[ilev] 
      #Find the azimuths corresponding to the current elevation.
      azlev=azimuth[ levmask ]
      timelev=time[ levmask ]
      #Get variabile values corresponding to the current elevation
      varlev=var[ levmask , : ]

      #For the first azimuth which is a special case because it contains zero.
      az_index=np.logical_or( azlev <= ray_angle_res/2.0 , azlev >= 360 - ray_angle_res/2.0 )
     
      if ( np.sum(az_index) > 0 ) : 
         order_var[0,:,ilev] = np.nanmean( varlev[az_index,:] , 0 )
         order_time[0,ilev] = np.nanmean( timelev[ az_index ] )
         azimuth_exact[0,ilev] = np.nanmean( azlev[ az_index ] )

      #Para los que vienen despues.
      for iaz in range(1,naz) :
         #Search for all the rays that are close to the current azimuth and level.
         az_index=np.logical_and( azlev <= order_azimuth[iaz] + ray_angle_res/2.0 , azlev >= order_azimuth[iaz] - ray_angle_res/2.0 )
         if( np.sum( az_index ) > 0 ) :
            order_var[iaz,:,ilev] = np.nanmean( varlev[az_index,:] , 0 )
            order_time[iaz,ilev] = np.nanmean( timelev[ az_index ] )
            azimuth_exact[iaz,ilev] = np.nanmean( azlev[ az_index ] )
            azimuth_exact[iaz,ilev] = np.nanmean( azlev[ az_index ] )

   return order_var , order_azimuth , levels , order_time , azimuth_exact

