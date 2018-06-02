#print __doc__
# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause


def main_qc( filename , options ) :

   import sys
   import time
   sys.path.append('../fortran')

   import numpy as np
   import numpy.ma as ma
   import matplotlib.pyplot as plt
   import pyart
   from common_qc_tools  import qc  #Fortran code routines.
   from common_qc_tools  import qc_const
   import netCDF4

   import matplotlib.pyplot as plt

   qc_const.undef = options['undef']   #This set the undef value for the fortran routines.
   undef          = options['undef']   #This set the undef value for the rest of the script.

   output=dict() #Initialize output dictionary.

   computed_etfilter=False   #Flag to indicate if echo top has been computed already.

#  Constant parameters

   #Codigos que permitan identificar el control que actuo en cada pixel.

   #Para la reflectividad
   QCCODE_ATTENUATION = 10
   QCCODE_SPECKLE     = 11
   QCCODE_TEXTURE     = 12
   QCCODE_RHOFILTER   = 13
   QCCODE_SIGN        = 14
   QCCODE_BLOCKING    = 15
   QCCODE_ECHOTOP     = 16
   QCCODE_ECHODEPTH   = 17

   #Para la velocidad radial
   QCCODE_DEALIAS     = 30

   #El codigo de los datos buenos para reflectividad y velocidad radial.
   QCCODE_GOOD        = 0

   #Shortcut to variable names
   name_v=options['v_name']
   name_ref=options['ref_name']
   name_rho=options['rho_name']

   radar = pyart.io.read(filename)

   #Get the nyquist velocity
   nyquistv=radar.get_nyquist_vel(0,check_uniform=True)

   #===================================================
   # RESHAPE VARIABLES
   #===================================================
   #From time,range -> azimuth,range,elevation

   startt=time.time()


   if name_ref in radar.fields :
        start=time.time()

        [ output['ref'] , output['az'] , output['level'] , output['time'] , output['index'] , output['az_exact'] ]=order_variable( radar , name_ref , undef )
        na=output['ref'].shape[0]
        nr=output['ref'].shape[1]
        ne=output['ref'].shape[2]
        output['na']=na
        output['nr']=nr
        output['ne']=ne

        output['cref'] = np.zeros(output['ref'].shape) 

        output['cref'][:] = output['ref']                 #Initialize the corrected reflectivity array.

        output['cref'][ output['cref'] == undef ]=options['norainrefval']
        output['ref'] [ output['ref']  == undef ]=options['norainrefval']

        output['qcref'] = np.zeros(output['cref'].shape)  #Set the qc flag array to 0.

        end=time.time()

        print("The elapsed time in {:s} is {:2f}".format("ref -> az,r,el",end-start) )

 
   if name_v in radar.fields  : 

        start=time.time()

        [ output['v'] , output['az'] , output['level'] , output['time'] , output['index'] , output['az_exact']  ]=order_variable( radar , name_v , undef  )
 
        na=output['ref'].shape[0]
        nr=output['ref'].shape[1]
        ne=output['ref'].shape[2]
        output['na']=na
        output['nr']=nr
        output['ne']=ne
 
        output['cv'] = np.zeros(output['v'].shape) 

        output['cv'][:] = output['v']                     #Initialize the corrected doppler velocity array

        output['qcv'] = np.zeros(output['v'].shape)       #Set the qc flag array to 0.

        end=time.time()

        print("The elapsed time in {:s} is {:2f}".format("v -> az,r,el",end-start) )

   #===================================================
   # GEOREFERENCE RADAR DATA
   #===================================================

   #Use pyart rutines to compute x, y and z at each grid point
   #Reorder data
   start=time.time()

   #dm is a dummy variable

   [ output['altitude'] , dm , dm , dm , dm , dm ] = order_variable( radar , 'altitude' , undef )
   [ output['x']        , dm , dm , dm , dm , dm ] = order_variable( radar , 'x' , undef ) 
   [ output['y']        , dm , dm , dm , dm , dm ] = order_variable( radar , 'y' , undef )

   #Compute distance to radar
   output['distance']=np.power( np.power(output['x'],2)+np.power(output['y'],2) , 0.5 )

   end=time.time()

   print("The elapsed time in {:s} is {:2f}".format("x,y,z -> az,r,el",end-start) )

   #===================================================
   # DEALIASING 
   #===================================================

   if ( options['ifdealias'] and (name_v in radar.fields) ) : 
     
     start=time.time()

     #Uso una de las funciones de dealiasing de pyart con los parametros por defecto
     winddealias=pyart.correct.region_dealias.dealias_region_based(radar, interval_splits=options['interval_splits'],interval_limits=None, 
                 skip_between_rays=options['skip_between_rays'], skip_along_ray=options['skip_along_ray'],centered=True,nyquist_vel=None,
                 check_nyquist_uniform=True,gatefilter=False,rays_wrap_around=True,keep_original=False,set_limits=True,
                 vel_field=name_v,corr_vel_field=None)

     #Add wind dealias to the radar objetc.
     radar.fields['Vda']                  = winddealias
     radar.fields['Vda']['coordinates']   = radar.fields[name_v]['coordinates']
     radar.fields['Vda']['units']         = radar.fields[name_v]['units']
     radar.fields['Vda']['long_name']     = radar.fields[name_v]['long_name']
     radar.fields['Vda']['standard_name'] = radar.fields[name_v]['standard_name']

     #Re-order dealiased wind data.
     [ output['cv'] , output['az'] , output['level'] , output['time'] , output['index'] , output['az_exact']  ]=order_variable( radar , 'Vda' , undef  )

     end=time.time()

     print("The elapsed time in {:s} is {:2f}".format("dealiasing",end-start) )

   #===================================================
   # DEALIASING BORDER FILTER
   #===================================================

   #TODO

   #===================================================
   # RHO HV FILTER
   #===================================================

   if options['ifrhofilter']   :

      start=time.time()

      nx=options['rhofilternx']
      ny=options['rhofilterny']
      nz=options['rhofilternz']

      if rho_name in radar.fields  :

         [ output['rho'] , dm , dm , dm , dm , dm  ]=order_variable( radar , rho_name , undef  )

         output['rho_smooth']=box_functions_2d(var=output['rho'],na=na,nr=nr,ne=ne,nx=nx,ny=nz,nz=nz,operation='MEAN',threshold=0.0)

         output['cref'][ output['rho_smooth'] < options['rhofiltertr'] ] = undef
         output['qcref'][output['rho_smooth'] < options['rhofiltertr'] ] = QCCODE_RHOFILTER
 
      else   :
         display('Warning: could not perform RHO-HV filter because rho was not found on this file')

      if [ not options['rhofilter_save'] ] :
          output['rho_smooth']=0
          output['rho']=0

      end=time.time()
 
      print("The elapsed time in {:s} is {:2f}".format("rho filter",end-start) )


   #===================================================
   # ECHO TOP FILTER  
   #===================================================


   if ( options['ifetfilter'] ) :
  
     start=time.time()

     if ( not computed_etfilter )     :

        tmp_z=np.zeros((nr,ne))
        tmp_d=np.zeros((nr,ne))

        tmp_z[:]=output['altitude'][0,:,:]    #Store only one RHI section of the altitutde (we will assume that the altitude of a pixel is independent
                                              #of the azimuth).
        tmp_d[:]=output['distance'][0,:,:]    #Store only one RHI section of the distance to radar (we will assume that the altitude of a pixel is independent
                                              #of the azimuth). 
        tmp_max_z=np.zeros((na,nr,ne))

        for ii in range(0,output['ne'])     :       #Estimate the maximum radar data height assoicated with each gate.
           tmp_max_z[:,:,ii]=output['altitude'][:,:,output['ne']-1]
 
        [tmp_data_3d,tmp_data_2d]=qc.echo_top(reflectivity=output['cref'],heigth=tmp_z,rrange=tmp_d,na=na,nr=nr,ne=ne
                                                ,nx=options['etfilternx'],ny=options['etfilterny'],nz=options['etfilternz'])
        computed_etfilter = True  #In case we need any of the other variables computed in this routine.

        tmp_ref=0 #Unset tmp_ref to free some memory.
        

     #Mask the pixels with echo top values under the threshold.
     #Do not mask tose pixels where the echo top threshold in below the maximum radar height (i.e. the heigth of the maximum elevation)
     output['cref'][ np.logical_and( tmp_max_z > options['etfiltertr'] , tmp_data_3d[:,:,:,0] < options['etfiltertr'] ) ]   = undef
     output['qcref'][ np.logical_and( tmp_max_z > options['etfiltertr'] , tmp_data_3d[:,:,:,0] < options['etfiltertr'] )  ] = QCCODE_ECHOTOP

     #If requested store the auxiliary fields and data in the output dictionary.
     if  ( options['etfilter_save'] )     : 
        output['echo_top'] = tmp_data_3d[:,:,:,0] 


     end=time.time()

     print("The elapsed time in {:s} is {:2f}".format("echo-top filter",end-start) )

   #===================================================
   # ECHO DEPTH FILTER 
   #===================================================


   if ( options['ifedfilter']  ) :

     start=time.time()

     if ( not computed_etfilter )     :

        tmp_z=np.zeros((nr,ne))
        tmp_d=np.zeros((nr,ne))

        tmp_z[:]=output['altitude'][0,:,:]    #Store only one RHI section of the altitutde (we will assume that the altitude of a pixel is independent
                                              #of the azimuth).
        tmp_d[:]=output['distance'][0,:,:] #Store only one RHI section of the distance to radar (we will assume that the altitude of a pixel is independent
                                              #of the azimuth). 
        tmp_max_z=np.zeros((na,nr,ne))

        for ii in range(0,output['ne'])     :       #Estimate the maximum radar data height assoicated with each gate.
           tmp_max_z[:,:,ii]=output['altitude'][:,:,output['ne']-1]

        [tmp_data_3d,tmp_data_2d]=qc.echo_top(reflectivity=output['cref'],heigth=tmp_z,rrange=tmp_d,na=na,nr=nr,ne=ne
                                                ,nx=options['edfilternx'],ny=options['edfilterny'],nz=options['edfilternz'])
        computed_etfilter = True  #In case we need any of the other variables computed in this routine.

     #Mask the pixels with echo depth values under the threshold.
     #Do not mask tose pixels where the echo depth threshold in below the maximum radar height (i.e. the heigth of the maximum elevation)
     output['cref'] [ np.logical_and( tmp_max_z > options['edfiltertr'] , tmp_data_3d[:,:,:,2] < options['edfiltertr'] ) ] = undef
     output['qcref'][ np.logical_and( tmp_max_z > options['edfiltertr'] , tmp_data_3d[:,:,:,2] < options['edfiltertr'] ) ] = QCCODE_ECHODEPTH

     #If requested store the auxiliary fields and data in the output dictionary.
     if  ( options['edfilter_save'] )     :
        output['echo_depth'] = tmp_data_3d[:,:,:,2] 

     end=time.time()

     print("The elapsed time in {:s} is {:2f}".format("echo-depth filter",end-start) )
 
   #===================================================
   # SPECKLE FILTER
   #===================================================

   if options['ifspfilter']  :

       start=time.time()

       #TODO poner todo en terminos de la estructura output
       #Compute the number pixels with reflectivities over spfiltertr sourrounding each pixels in the box defined by nx,ny,nz.
       nx=options['spfilternx']
       ny=options['spfilterny']
       nz=options['spfilternz']
       tr=options['spfilterreftr']
       output['speckle']=qc.box_functions_2d(datain=output['cref'].data,na=na,nr=nr,ne=ne,boxx=nx,boxy=ny,boxz=nz,operation='COUN',threshold=tr) 

       #Set the pixels with values below the threshold as undef. 
       output['cref'][ output['speckle'] < options['spfiltertr']  ] = undef
       
       output['qcref'][ output['speckle'] < options['spfiltertr'] ] = QCCODE_SPECKLE

       #If the field is not included in the output then set it to 0.
       if [ not options['spfilter_save'] ] : 
          output['speckle']=0   

       end=time.time()

       print("The elapsed time in {:s} is {:2f}".format("speckle filter",end-start) )

   #===================================================
   # ATTENUATION FILTER
   #===================================================

   if options['ifattfilter']   :
      
      start=time.time()

      beaml=radar.range['data'][1]-radar.range['data'][0] #Get beam length

      output['attenuation']=qc.get_attenuation(var=output['cref'],na=na,nr=nr,ne=ne,beaml=beaml,cal_error=options['attcalerror'])

      #Set the pixels with values below the threshold as undef. 
      output['cref'][ output['attenuation'] < options['attfiltertr']  ] = undef

      output['qcref'][ output['attenuation'] < options['attfiltertr'] ] = QCCODE_ATTENUATION

      if [ not options['attfilter_save'] ] :
          output['attenuation']=0


      end=time.time()

      print("The elapsed time in {:s} is {:2f}".format("attenuation filter",end-start) )

   #===================================================
   # TOPOGRAPHY BLOCKING FILTER
   #===================================================

   #TODO

   #===================================================
   # DOPPLER NOISE FILTER
   #===================================================

   #TODO

   #===================================================
   # INTERFERENCE FILTER
   #===================================================

   #TODO

   #===================================================
   # LOW DOPPLER VOLOCITY FILTER
   #===================================================

   #TODO

   #===================================================
   # ADD CORRECTED DATA TO RADAR OBJECT
   #===================================================

   #TODO

   #===================================================
   # END
   #===================================================


   endt=time.time()

   print("The elapsed time in {:s} is {:2f}".format("the entire QC",endt-startt) )

   return radar , output

#===========================================================================================================
# OTRAS FUNCIONES CONTENIDAS EN ESTE MODULO
#===========================================================================================================   

def order_variable ( radar , var_name , undef )  :  

   import numpy as np
   import numpy.ma as ma

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
   #print ('The resolution in azimuth is: %5.3f' % ( ray_angle_res ) )


   levels=np.unique(radar.elevation['data'])
   azimuth=radar.azimuth['data']
   time=radar.time['data']

   order_azimuth=np.arange(0.0,360.0-ray_angle_res,ray_angle_res) #Asuming quasi regular azimuth location
   naz=np.size(order_azimuth)
   nel=np.size(levels)

   if ( var_name == 'altitude' ) :
      var=radar.gate_altitude['data']
   elif( var_name == 'longitude' ) :
      var=radar.gate_longitude['data'] 
   elif( var_name == 'latitude'  ) :
      var=radar.gate_latitude['data'] 
   elif( var_name == 'x' )         :
      var=radar.gate_x['data']
   elif( var_name == 'y' )         : 
      var=radar.gate_y['data']
   else  :
      var=radar.fields[var_name]['data'].data


      var[ var == undef ] = np.nan

   nr=var.shape[1]

   #Allocate arrays
   order_var    = undef + np.zeros((naz,nr,nel))
   order_time   =np.zeros((naz,nel)) 
   order_index  =undef + np.zeros((naz,nel))   #This variable can be used to convert back to the azimuth - range array
   azimuth_exact=undef + np.zeros((naz,nel))

   order_var[:]     = undef 
   order_time[:]    = undef 
   azimuth_exact[:] = undef

   

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

   order_var[ np.isnan( order_var ) ]= undef
   order_index[ np.isnan( order_index ) ]=undef

   return order_var , order_azimuth , levels , order_time , order_index , azimuth_exact


   def order_variable_inv (  radar , var , order_index  )  :

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

             if ( not order_index[ia,ie] == undef )  :
       
                output_var[ia,order_index[ia,ie]]=var[ia,:,ie] 

   return output_var


