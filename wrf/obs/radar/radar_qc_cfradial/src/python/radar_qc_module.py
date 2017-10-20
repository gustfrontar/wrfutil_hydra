#print __doc__

# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause


def main_qc( filename , options ) :

   sys.path.append('../fortran')

   import numpy as np
   import matplotlib.pyplot as plt
   import pyart
   import common_qc_tools  as qct  #Fortran code routines.
   import netCDF4

#  Constant parameters

   #Codigos que permitan identificar el control que actuo en cada pixel.

   #Para la reflectividad
   QCCODE_ATTENUATION = 10
   QCCODE_SPECKLE     = 11
   QCCODE_TEXTURE     = 12
   QCCODE_RHOFILTER   = 13
   QCCODE_SIGN        = 14
   QCCODE_BLOCKING    = 15

   #Para la velocidad radial
   QCCODE_DEALIAS     = 30

   #El codigo de los datos buenos para reflectividad y velocidad radial.
   QCCODE_GOOD        = 0


   name_v=options['dv_var_name']
   name_ref=options['reflectivity_var_name']

   radar = pyart.io.read(filename)

   #Get the nyquist velocity
   nyquistv=radar.get_nyquist_vel(0,check_uniform=True)

   #===================================================
   # RESHAPE VARIABLES
   #===================================================
      #From time,range -> azimuth,range,elevation

      [ ref_array , az , level , time , index , az_exact ]=order_variable( radar , name_ref )
      [ v_array  ]=order_variable( radar , name_v )

      na=ref_array.shape[0]
      nr=ref_array.shape[1]
      ne=ref_array.shape[2]

      cref_array=ref_array   #Define the corrected reflectivity array.
      cv_array  =v_array     #Define the corrected radial velocity array.


   #===================================================
   # INITIALIZE QC FLAGS
   #===================================================

      #TODO if present V, initialize for V
      #TODO if present DBZ, initialize for DBZ
      #Estos seran arrays que contengan para cada pixel
      #un codigo que indique si es un dato bueno o si 
      #no paso alguno de los contrles que indique cual fue
      qcref_array=np.zeros(cref_array.shape)
      qcv_array  =np.zeros(v_array.shape)

   #===================================================
   # GEOREFERENCE RADAR DATA
   #===================================================

      #Use pyart rutines to compute lat,lon and z at each grid point
      #Reorder data
      [ altitude_array ]=order_variable( radar , 'altitude' )
      [ longitude_array ]=order_variable( radar , 'longitude' )
      [ latitude_array ]=order_variable( radar , 'latitude' ) 

   #===================================================
   # DEALIASING 
   #===================================================

   if options['ifdealias'] : 

     #Uso una de las funciones de dealiasing de pyart con los parametros por defecto
     winddealias=pyart.correct.region_dealias.dealias_region_based(radar, interval_splits=options['interval_splits'], interval_limits=None, skip_between_rays=options['skip_between_rays'], skip_along_ray=options['skip_along_ray'], centered=True, nyquist_vel=None, check_nyquist_uniform=True, gatefilter=False, rays_wrap_around=True, keep_original=False, set_limits=True, vel_field=name_v, corr_vel_field=None)

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

     #qct.echo_top(reflectivity=cref_array,heigth,rrange,na,nr,ne,nx,ny,nz,output_data_3d,output_data_2d)



   
   #===================================================
   # SPECKLE FILTER
   #===================================================

   if options['ifspfilter']  :

     #Compute the number pixels with reflectivities over spfiltertr sourrounding each pixels in the box defined by nx,ny,nz.
     qct.speckle_filter(var=cref_array,na=na,nr=nr,ne=ne,nx=options['spfilternx'],ny=options['spfilterny'],nz=options['spfilternz'],threshold=options['spfiltertr'],speckle=speckle)

     #Set the pixels with values below the threshold as undef. 
     cref_array[ speckle < options['spfiletertr']  ]=options['undef'] 

   
   #===================================================
   # ATTENUATION FILTER
   #===================================================



   #===================================================
   # TOPOGRAPHY BLOCKING FILTER
   #===================================================



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

   if ( var_name == 'altitude' ) :
      var=radar.gate_altitude['data']
   else if( var_name == 'longitude' ) :
      var=radar.gate_longitude['data'] 
   else if( var_name == 'latitude'  ) :
      var=radar.gate_latitude['data'] 
   else
      var=radar.fields[var_name]['data']

   nr=var.shape[1]

   #Allocate arrays
   order_var=np.zeros((naz,nr,nel))
   order_time=np.zeros((naz,nel)) 
   order_index=np.zeros((naz,nel))   #This variable can be used to convert back to the azimuth - range array
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
            order_index[iaz,ilev] = az_index[0]  #If multiple levels corresponds to a single azimuth / elevation chose the first one.

   return order_var , order_azimuth , levels , order_time , order_index , azimuth_exact


   def order_variable_inv (  radar , var , order_index )  :

       import numpy as np
   
       #Esta funcion es la inversa de la funcion order variable. Es decir que toma un array ordenado como azimuth , range y elevation y lo vuelve
       #a ordenar como azimuth-elevation y range. Es decir el orden original que se encuentra en los archivos con formato cfradial y que heredan los objetos radar de pyart.

       #var es la variable ordenada como var(azimuth,range,elevation)
       #order_index (azimuth,elevation) contiene la posicion original de los haces que fueron asignados a cada azimuth y elevacion por la funcion order_variable.
       #nr numero de puntos en la direccion del rango.
       #nb numero de beams. 

       na=var.shape[0]
       nr=var.shape[1]
       ne=var.shape[2]

       nb=radar.azimuth['data'].shape[0]

       output_var=np.ones((na,nb)) * options['undef']
       

       for ia in range(0,na)  :

          for ie in range(0,ne)  :
       
             output_var[ia,order_index[ia,ie]]=var[ia,:,ie] 

        #TODO ver como convertir la salida en un masked array.


   return output_var


