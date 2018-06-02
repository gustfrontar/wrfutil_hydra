import numpy as np
import datetime as dt
import matplotlib 
import matplotlib.pyplot as plt
import os
#import smooth_2d as s2d
from scipy import signal
from scipy import misc
import warnings        #To supress certain warnings.

from mpl_toolkits.basemap import Basemap

#Some module constant section.

default_figure_size=(10,7.5)

def data_diff(data_pp,data_pn):

    data_diff=dict()

    for key in data_pp:
    
       data_diff[key]=data_pp[key]-data_pn[key]
    
    return data_diff

def norm_bv( data_p , data_n , norm_type='None' , smooth='None' , sigma='None',cutoff='None' , xi='None' , xe='None' , yi='None' , ye='None', zi='None', ze='None'):
#Computes the bv in a similar way as the breeding fortran module does.
#Sigma is the standard deviation of the gaussian filter in units of grid points.
#Smooth is a filter option. Currently gaussian filter is supported.

    if norm_type == 'None' : #UVT, UV or T 
       norm_type = 'UVT'

    if smooth == 'None'    : #None, Gaussian or Lanczos
       filter_norm=False
    else                   :
       filter_norm=True 

    if sigma == 'None' or sigma==0  : 
       sigma = 1
       print('Warning : Sigma = 1 ')

    if cutoff == 'None'    : #Convolution function size for the Gaussian filter 
       cutoff=10

    tmp=np.array( data_p['U'] )
    tmp_shape = np.shape( tmp )
    nx=tmp_shape[0]
    ny=tmp_shape[1]
    nz=tmp_shape[2]
   
    norm=np.zeros((nx,ny,nz))

    if norm_type == 'UVT' :
       
       norm= norm + np.power( data_p['U'] - data_n['U'] , 2. )

       norm= norm + np.power( data_p['V'] - data_n['V'] , 2. )

       norm= norm + np.power( data_p['T'] - data_n['T'] , 2. )

    if norm_type == 'UV' :

       norm= norm + np.power( data_p['U'] - data_n['U'] , 2. )

       norm= norm + np.power( data_p['V'] - data_n['V'] , 2. )

    if norm_type == 'W' :

       norm= norm + np.power( data_p['W'] - data_n['W'] , 2. )

    if norm_type == 'T'  :

       norm= norm + np.power( data_p['T'] - data_n['T'] , 2. )

    if norm_type == 'QV'  :

       norm= norm + np.power( data_p['QV'] - data_n['QV'] , 2. )

    if norm_type == 'QHYD'  :

       norm= norm + np.power( data_p['QHYD'] - data_n['QHYD'] , 2. )

    if filter_norm :
      if smooth == 'Gaussian'  :
         filter_size=np.around( 2*cutoff*sigma , 0 ).astype(int)

         gaussian_conv_func=np.zeros([ 2*filter_size+1 , 2*filter_size+1 ])
         for ii in range(0,2*filter_size+1) : 
            for jj in range(0,2*filter_size+1) :
               gaussian_conv_func[ii,jj] = np.exp( -0.5*(np.power( ii-filter_size,2 ) + np.power( jj-filter_size, 2) ) /np.power(sigma,2)   )
       
         #Normalize the convolving function.
         gaussian_conv_func=gaussian_conv_func/np.sum( gaussian_conv_func)

         #a=plt.figure
         #plt.pcolor( gaussian_conv_func )
         #plt.show(a)


         for iz in range(0,nz) :

             norm[:,:,iz]=signal.convolve2d(norm[:,:,iz],gaussian_conv_func, boundary='symm', mode='same')

      #elif smooth == 'Lanczos'  :
      #   for iz in range(0,nz) :
      #       mask=np.ones((nx,ny))
      #       norm[:,:,iz]=s2d.filter_2d(inputvar2d=norm[:,:,iz],mask=mask,lambdaf=sigma,ctyp=smooth,nx=nx,ny=ny)
      else                      :   #TODO add lanczos 2D filter option.

         print('Smooth option not recognized, we will not smooth the norm')

    norm=np.power(norm,0.5)

    norm_mean , norm_max , norm_min = get_regional_average(norm,xi=xi,xe=xe,yi=yi,ye=ye,zi=zi,ze=ze)

    return norm_mean , norm_max , norm_min  , norm  #Generate a tuple as output.

def norm_bv_2( data_p , data_n , norm_type='None' , smooth='None' , sigma='None',cutoff='None' ):
#Computes the bv in a similar way as the breeding fortran module does.
#Sigma is the standard deviation of the gaussian filter in units of grid points.
#Smooth is a filter option. Currently gaussian filter is supported.

    if norm_type == 'None' : #UVT, UV or T 
       norm_type = 'UVT'

    tmp=np.array( data_p['U'] )
    tmp_shape = np.shape( tmp )
    nx=tmp_shape[0]
    ny=tmp_shape[1]
    nz=tmp_shape[2]

    norm=np.zeros((nx,ny,nz))

    if norm_type == 'UVT' :

       norm= norm + np.power( data_p['U'] - data_n['U'] , 2. )

       norm= norm + np.power( data_p['V'] - data_n['V'] , 2. )

       norm= norm + np.power( data_p['T'] - data_n['T'] , 2. )

    if norm_type == 'UV' :

       norm= norm + np.power( data_p['U'] - data_n['U'] , 2. )

       norm= norm + np.power( data_p['V'] - data_n['V'] , 2. )

    if norm_type == 'W' :

       norm= norm + np.power( data_p['W'] - data_n['W'] , 2. )

    if norm_type == 'T'  :

       norm= norm + np.power( data_p['T'] - data_n['T'] , 2. )

    if norm_type == 'QV'  :

       norm= norm + np.power( data_p['QV'] - data_n['QV'] , 2. )

    if norm_type == 'QHYD'  :

       norm= norm + np.power( data_p['QHYD'] - data_n['QHYD'] , 2. )

    norm=np.power( norm , 0.5 )

    return norm  #Generate a tuple as output.

def filter_field( data , smooth='None' , sigma='None', cutoff='None' ):

    if smooth == 'None'    : #None, Gaussian or Lanczos
       filter_norm=False
    else                   :
       filter_norm=True

    if sigma == 'None' or sigma==0  :
       sigma = 1
       print('Warning : Sigma = 1 ')

    if cutoff == 'None'    : #Convolution function size for the Gaussian filter 
       cutoff=10


    if( np.size( np.shape( data ) ) >= 3 ) :
      nz=data.shape[2]
    else :
      nz=1

    data_s=np.zeros( np.shape(data) )

    if smooth == 'Gaussian'  :
       filter_size=np.round( 2*cutoff*sigma , 0 )
  
       gaussian_conv_func=np.zeros(( 2*filter_size.astype(int) +1 , 2*filter_size.astype(int) +1 ))
       for ii in range(0,2*filter_size.astype(int) + 1 ) :
          for jj in range(0,2*filter_size.astype(int)+1) :
             gaussian_conv_func[ii,jj] = np.exp( -0.5*(np.power( ii-filter_size.astype(int),2 ) + np.power( jj-filter_size.astype(int), 2) ) /np.power(sigma,2)   )

       #Normalize the convolving function.
       gaussian_conv_func=gaussian_conv_func/np.sum( gaussian_conv_func )

       if ( nz == 1 ) :
          data_s[:,:]=signal.convolve2d(data[:,:],gaussian_conv_func, boundary='symm', mode='same')
       else :
          for iz in range(0,nz) :
             data_s[:,:,iz]=signal.convolve2d(data[:,:,iz],gaussian_conv_func, boundary='symm', mode='same')

    #elif smooth == 'Lanczos'  :
    #   for iz in range(0,nz) :
    #       mask=np.ones((nx,ny))
    #       norm[:,:,iz]=s2d.filter_2d(inputvar2d=norm[:,:,iz],mask=mask,lambdaf=sigma,ctyp=smooth,nx=nx,ny=ny)
    else                      :   #TODO add lanczos 2D filter option.

       print('Smooth option not recognized, we will not smooth the norm')

    return data_s  #Generate a tuple as output.


                                                                           

def lat_lon_to_i_j(lonfield,latfield,lonlist,latlist) :
#Gets the i,j which is closer to a particular lat and lon give a latfield, lonfield.
   npoints=latlist.size

   i=np.zeros(latlist.shape)
   j=np.zeros(latlist.shape)
 
   for ipoint in range(0,npoints) :
 
     dist=np.power( latfield - latlist[ipoint] , 2.0 ) + np.power( lonfield - lonlist[ipoint] , 2.0 )

     #Get the indices of the minimum
     i[ipoint] , j[ipoint] = np.unravel_index(dist.argmin(), dist.shape)

   return i , j 

def growth_rate_bv( norm_o , norm_r , xi='None' , xe='None' , yi='None' , ye='None', zi='None', ze='None') :

   nx = np.shape( norm_o )[0]
   ny = np.shape( norm_o )[1]
   nz = np.shape( norm_o )[2]

   nregs=xi.size  #Cuantas regiones tenemos definidas.

   growth_rate=norm_o[:,:,:] - norm_r[:,:,:]  #3d growing rate field.

   growth_rate_mean=np.zeros(nregs)
   growth_rate_max=np.zeros(nregs)
   growth_rate_min=np.zeros(nregs)

   growth_rate_mean , growth_rate_max , growth_rate_min = get_regional_average(growth_rate,xi=xi,xe=xe,yi=yi,ye=ye,zi=zi,ze=ze)

   return growth_rate_mean , growth_rate_max , growth_rate_min , growth_rate

def get_regional_average( var , xi='None' , xe='None' , yi='None' , ye='None', zi='None', ze='None') :
   #Given a 3D variable and a set of regions, get the mean, max and min of the variable over the grid points
   #contained in these regions.
   #Region is defined as a 3D rectangular section of the grid (in grid space). 

   nx = np.shape( var )[0]
   ny = np.shape( var )[1]
   if ( np.size ( np.shape ( var ) ) >= 3 ) :
      nz = np.shape( var )[2]
   else :
      nz = 1

   if xi == 'None'  :
      xi=np.ndarray(1)
      xi[:]=0
   if xe == 'None'  :
      xe=np.ndarray(1)
      xe[:]=(nx - 1)
   if yi == 'None'  :
      yi=np.ndarray(1)
      yi[:]=0
   if ye == 'None'  :
      ye=np.ndarray(1)
      ye[:]=ny-1
   if zi == 'None'  :
      zi=np.ndarray(xi.shape)
      zi[:]=0
   if ze == 'None'  :
      ze=np.ndarray(xi.shape)
      ze[:]=(nz-1)

   nregs=xi.size  #Get number of regions

  
   var_mean=np.zeros(nregs)
   var_max=np.zeros(nregs)
   var_min=np.zeros(nregs)

   for ireg in range(0,nregs) :
 
      if ( nz > 1 ) :
       
        tmp=var[xi[ireg]:xe[ireg]+1,yi[ireg]:ye[ireg]+1,zi[ireg]:ze[ireg]+1] #Get subregion for var.

      else : 
 
        tmp=var[xi[ireg]:xe[ireg]+1,yi[ireg]:ye[ireg]+1]

      #Np.nan.. functions produce warnings with the input consist only of nan  values. This warning is ignored. 
      with warnings.catch_warnings()    :
        warnings.simplefilter("ignore", category=RuntimeWarning)

        var_mean[ireg]=np.nanmean( tmp )

        var_max[ireg]=np.nanmax( tmp )
  
        var_min[ireg]=np.nanmin( tmp )

   return var_mean , var_max , var_min 


def plot_bv(bv,data,lon,lat,plotvars,plotlevels,plotbasedir,mycolorbar='seismic',range='centered',figextension='None',fontsize='None',offset='None',figsize='None',debug=False,ndatalevels='None',date='None'):

#Example of optional inputs with default values.
    if offset == 'None' :
      offset=0

#End of example

    #Create output directory
    if not os.path.exists(plotbasedir):
      os.makedirs(plotbasedir)

    tmp=np.shape(lon)

    nx=tmp[0]
    ny=tmp[1]


    for key in bv     :

      for var in plotvars : 
 
         if var == key     :  #We will plot a figure.
        
             plotvar=bv[key][offset:nx-offset,offset:ny-offset,:]
             plotlon=lon[offset:nx-offset,offset:ny-offset]
             plotlat=lat[offset:nx-offset,offset:ny-offset]

             plotvarname='bv_' + key 
             plot_var_levels(plotvar,plotlon,plotlat,plotlevels,plotbasedir,plotvarname,mycolorbar=mycolorbar,range=range,figextension='None',fontsize='None',offset='None',figsize='None',debug=False,date=date) 



def plot_state(state,lon,lat,plotvars,plotlevels,plotbasedir,mycolorbar='coolwarm',range='maxmin',figextension='None',fontsize='None',offset='None',figsize='None',debug=False,ndatalevels='None',date='None'):

#Example of optional inputs with default values.
    if offset == 'None' :
      offset=0

#End of example

    #Create output directory
    if not os.path.exists(plotbasedir):
      os.makedirs(plotbasedir)

    tmp=np.shape(lon)

    nx=tmp[0]
    ny=tmp[1]


    for key in state     :

      for var in plotvars :

         if var == key     :  #We will plot a figure.

             plotvar=state[key][offset:nx-offset,offset:ny-offset,:]
             plotlon=lon[offset:nx-offset,offset:ny-offset]
             plotlat=lat[offset:nx-offset,offset:ny-offset]

             plotvarname='state_' + key
             plot_var_levels(plotvar,plotlon,plotlat,plotlevels,plotbasedir,plotvarname,varcontour='None',mycolorbar=mycolorbar,range=range,figextension='None',fontsize='None',offset='None',figsize='None',debug=False,date=date)


def plot_var_levels(var,lon,lat,plotlevels,plotbasedir,varname,varcontour='None',clevels='None',ndatalevels='None',range='centered',scale_max='None',scale_min='None',figextension='None',fontsize='None',offset='None',figsize='None',debug=False,date='None',mycolorbar='seismic') :

#Example of optional inputs with default values.
    if date=='None'           :
       got_date=False
       date=''
    else                      :
       got_date=True

    if ndatalevels=='None' :
       ndatalevels=10

    if figextension == 'None' :
      figextension='png'

    if fontsize == 'None' :
      fontsize=20

    if offset == 'None' :
      offset=0

    if figsize == 'None' :
       figsize=default_figure_size

    if varcontour== 'None' :
       plot_contour=False
    else :
       plot_contour=True
       #If we have a contour var, then define the clevels for plotting it.
       if clevels == 'None' :
          if ( np.size( varcontour.shape ) >= 3 ) :
            data_range_max=np.max(abs(np.mean(varcontour,2)))
          else :
            data_range_max=np.max(abs( varcontour ) )

          data_range_min=-data_range_max
          delta_data=( data_range_max - data_range_min ) / 5
		  #clevels=np.zeros( ndatalevels )
          clevels=np.arange( data_range_min , data_range_max , delta_data )

    #Create output directory
    if not os.path.exists(plotbasedir):
      os.makedirs(plotbasedir)

    #Set font size
    matplotlib.rcParams.update({'font.size': fontsize})

    tmp=np.shape(lon)

    nx=tmp[0]
    ny=tmp[1]

    true_lon=np.mean(lon) #For projection.
    true_lat=np.mean(lat) #For projection.
    #Get upper right and lower left lats and lons for prejection.
    ll_lat=lat[0,0]
    ur_lat=lat[nx-1,ny-1]
    ll_lon=lon[0,0]
    ur_lon=lon[nx-1,ny-1]

    for level in plotlevels :

        fig=plt.figure(1,figsize=figsize)

        
        #Basemap plotting section (maybe this can be sent to an independent function)
        #Set the map projection
        m = Basemap(projection='stere',lon_0=true_lon,lat_0=90.0,lat_ts=true_lat,\
          llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,   \
          rsphere=6371200.,resolution='h',area_thresh=10000)
        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()
        # draw parallels.
        parallels = np.arange(-90.0,90,0.5)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)

        # draw meridians
        meridians = np.arange(0.,360.,0.5)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

        if ( np.size( np.shape(var) ) >= 3 ) :
           tmpplot=var[offset:nx-offset,offset:ny-offset,level]
        else :
           tmpplot=var[offset:nx-offset,offset:ny-offset]

        tmplon=lon[offset:nx-offset,offset:ny-offset]

        tmplat=lat[offset:nx-offset,offset:ny-offset]

        lonplot , latplot = m(tmplon,tmplat)
 
        if range == 'centered' or range == 'None'  :

           var_range_max=np.nanmax(abs(tmpplot))
           var_range_min=-var_range_max

        elif range == 'maxmin'                     :
 
           var_range_max=np.nanmax(tmpplot)
           var_range_min=np.nanmin(tmpplot) 
        elif range == 'fixed'                      :
           if scale_max == 'None' :
              var_range_max=np.nanmax(tmpplot)
              var_range_min=np.nanmin(tmpplot)
           else                   :
              var_range_max=scale_max
              var_range_min=scale_min
        elif range == 'positive'                   :
              var_range_min=0
              var_range_max=np.nanmax(tmpplot)  


        tmpplot = np.ma.masked_invalid(tmpplot)

        m.pcolor(lonplot,latplot,tmpplot,cmap=plt.cm.get_cmap(mycolorbar),vmax=var_range_max,vmin=var_range_min)

        m.colorbar()

        if plot_contour   :
 
          if( np.size( np.shape( varcontour ) ) >= 3 ) :
             tmpplot=np.nanmean(varcontour[offset:nx-offset,offset:ny-offset,:],2)
          else :
             tmpplot=varcontour[offset:nx-offset,offset:ny-offset]
          matplotlib.rcParams['contour.negative_linestyle'] = 'solid'  #Forces negative lines to be solid too.
          m.contour(lonplot,latplot,tmpplot,levels=clevels,colors='k',linewidths=2,inline=1,fmt='%1.1f',fontsize=12)

        my_title=varname + ' at level ' + str(level)
        plt.ylabel('Lat')
        plt.xlabel('Lon')

        
        if got_date   :
           my_title=my_title + ' (' + date + ')'
           plt.title(my_title)
           print( 'Generationg the following figure : ' + 'Figure_' + varname + '_' + date + '_' + str(level) + '.' + figextension )
           plt.savefig( plotbasedir + '/Figure_' + varname + '_' + date + '_' + str(level) + '.' + figextension )
        else          :
           plt.title(my_title)
           print( 'Generationg the following figure : ' + 'Figure_' + varname + '_' + str(level) + '.' + figextension )
           plt.savefig( plotbasedir + '/Figure_' + varname + '_' + str(level) + '.' + figextension )


        if debug   :
           plt.show(fig)
           plt.close(fig)
        else       :
           plt.close(fig)





def plot_var_ave(var,lon,lat,plotbasedir,varname,range='centered',varcontour='None',clevels='None',ndatalevels='None',levels='None',figextension='None',fontsize='None',offset='None',figsize='None',debug=False,date='None',mycolorbar='seismic') :

#Example of optional inputs with default values.

    tmp=np.shape(var)

    nx=tmp[0]
    ny=tmp[1]
    nz=tmp[2]

    if date=='None'           :
       got_date=False
       date=''
    else                      :
       got_date=True

    if ndatalevels=='None' :
       ndatalevels=10

    if figextension == 'None' :
      figextension='png'

    if fontsize == 'None' :
      fontsize=20

    if offset == 'None' :
      offset=0

    if figsize == 'None' :
       figsize=default_figure_size
 
    if levels == 'None'  :
       levels=np.arange(0,nz)

    if varcontour == 'None' :
       plot_contour=False
    else                    :
       plot_contour=True
       
       if clevels == 'None' :
          data_range_max=np.nanmax(abs(np.mean(varcontour,2)))
          data_range_min=-data_range_max
          delta_data=( data_range_max - data_range_min )
          clevels=np.zeros( ndatalevels )
          clevels=np.arange( data_range_min , data_range_max , delta_data )

#End of example

    #Create output directory
    if not os.path.exists(plotbasedir):
      os.makedirs(plotbasedir)

    #Set font size
    matplotlib.rcParams.update({'font.size': fontsize})

    true_lon=np.mean(lon) #For projection.
    true_lat=np.mean(lat) #For projection.
    #Get upper right and lower left lats and lons for prejection.
    ll_lat=lat[0,0]
    ur_lat=lat[nx-1,ny-1]
    ll_lon=lon[0,0]
    ur_lon=lon[nx-1,ny-1]

    #Plot averaged norm.
    fig=plt.figure(1,figsize=figsize)


    #Basemap plotting section (maybe this can be sent to an independent function)
    #Set the map projection
    m = Basemap(projection='stere',lon_0=true_lon,lat_0=90.0,lat_ts=true_lat,\
      llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,   \
      rsphere=6371200.,resolution='h',area_thresh=10000)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    # draw parallels.
    parallels = np.arange(-90.0,90,0.5)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)

    # draw meridians
    meridians = np.arange(0.,360.,0.5)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

    #Plot vertical average of the variable.
    tmpplot=np.nanmean(var[offset:nx-offset,offset:ny-offset,:],2) 

    tmplon=lon[offset:nx-offset,offset:ny-offset]

    tmplat=lat[offset:nx-offset,offset:ny-offset]

    lonplot , latplot = m(tmplon,tmplat)

    if range == 'centered' or range == 'None'  :

       var_range_max=np.nanmax(abs(tmpplot))
       var_range_min=-var_range_max

    elif range == 'maxmin'                     :

       var_range_max=np.nanmax(tmpplot)
       var_range_min=np.nanmin(tmpplot)

    tmpplot = np.ma.masked_invalid(tmpplot)

    m.pcolor(lonplot,latplot,tmpplot,cmap=plt.cm.get_cmap(mycolorbar),vmax=var_range_max,vmin=var_range_min)

    m.colorbar()

    if plot_contour   :
      tmpplot=np.nansum(varcontour[offset:nx-offset,offset:ny-offset,:],2)
      matplotlib.rcParams['contour.negative_linestyle'] = 'solid'  #Forces negative lines to be solid too.
      m.contour(lonplot,latplot,tmpplot,levels=clevels,colors='k',linewidths=2,inline=1,fmt='%1.1f',fontsize=12)

    my_title= varname + ' vertical average  ' 
    if got_date   :
      my_title=my_title + ' (' + date + ')'


    plt.title( my_title )

    #plt.grid(True)

    plt.ylabel('Lat')
    plt.xlabel('Lon')

    print( 'Generationg the following figure : ' + 'Figure_' + varname + '_' + date + '_avez.' + figextension )
    plt.savefig( plotbasedir + 'Figure_' + varname + '_' + date + '_avez.' + figextension )

    if debug   :
       plt.show(fig)
    else       :
       plt.close(fig)

    #X-Z -average plot

    fig=plt.figure(1,figsize=figsize)
    tmpplot=np.nanmean(var[offset:nx-offset,offset:ny-offset,:],0) 

    lonplot=np.nanmean(lon[offset:nx-offset,offset:ny-offset],0)


    if range == 'centered' or range == 'None'  :

       var_range_max=np.nanmax(abs(tmpplot))
       var_range_min=-var_range_max

    elif range == 'maxmin'                     :

       var_range_max=np.nanmax(tmpplot)
       var_range_min=np.nanmin(tmpplot)

    var_range_max=np.nanmax(abs(tmpplot))
    var_range_min=-np.nanmax(abs(tmpplot))

    levplot , lonplot = np.meshgrid( levels , lonplot )

    tmpplot = np.ma.masked_invalid(tmpplot)

    plt.pcolor(lonplot,levplot,tmpplot,cmap=plt.cm.get_cmap(mycolorbar),vmax=var_range_max,vmin=var_range_min)
 
    plt.colorbar()

    my_title= varname + ' Y average  ' 
    if got_date   :
       my_title=my_title + ' (' + date + ')'

    
    plt.title( my_title )

    plt.ylim(np.min(levplot),np.max(levplot))
    plt.xlim(np.min(lonplot),np.max(lonplot))

    plt.grid(True)

    plt.ylabel('Height')
    plt.xlabel('Lon')

    if plot_contour   :
      tmpplot=np.nanmean(varcontour[offset:nx-offset,offset:ny-offset,:],0)
      matplotlib.rcParams['contour.negative_linestyle'] = 'solid'  #Forces negative lines to be solid too.
      plt.contour(lonplot,levplot,tmpplot,levels=clevels,colors='k',linewidths=2,inline=1,fmt='%1.1f',fontsize=12)


    print( 'Generationg the following figure : ' + 'Figure_' + varname + '_' + date + '_avey.' + figextension )
    plt.savefig( plotbasedir + 'Figure_' + varname + '_' + date + '_avey.' + figextension )

    if debug   :
       plt.show(fig)
    else       :
       plt.close(fig)

    #Y-Z -average plot
    fig=plt.figure(1,figsize=figsize)
    tmpplot=np.nanmean(var[offset:nx-offset,offset:ny-offset,:],1)

    latplot=np.nanmean(lat[offset:nx-offset,offset:ny-offset],1)
    levplot , latplot = np.meshgrid( levels , latplot )

    if range == 'centered' or range == 'None'  :

       var_range_max=np.nanmax(abs(tmpplot))
       var_range_min=-var_range_max

    elif range == 'maxmin'                     :

       var_range_max=np.nanmax(tmpplot)
       var_range_min=np.nanmin(tmpplot)

    tmpplot = np.ma.masked_invalid(tmpplot)

    plt.pcolor(latplot,levplot,tmpplot,cmap=plt.cm.get_cmap(mycolorbar),vmax=var_range_max,vmin=var_range_min)

    plt.colorbar()

    my_title= varname + ' X average at '
    if got_date   :
      my_title=my_title + ' (' + date + ')'

    plt.title( my_title )

    plt.ylim(np.min(levplot),np.max(levplot))
    plt.xlim(np.min(latplot),np.max(latplot))

    plt.grid(True)

    plt.ylabel('Height')
    plt.xlabel('Lat')


    if plot_contour   :
      tmpplot=np.nanmean(varcontour[offset:nx-offset,offset:ny-offset,:],1)
      matplotlib.rcParams['contour.negative_linestyle'] = 'solid'  #Forces negative lines to be solid too.
      plt.contour(latplot,levplot,tmpplot,levels=clevels,colors='k',linewidths=2,inline=1,fmt='%1.1f',fontsize=12)


    print( 'Generationg the following figure : ' + 'Figure_' + varname + '_' + date + '_avex.' + figextension )
    plt.savefig( plotbasedir + 'Figure_' + varname + '_' + date + '_avex.' + figextension )

    if debug   :
       plt.show(fig)
    else       :
       plt.close(fig)



def plot_norm_timeseries(norm_o,norm_i,norm_r,var_name,reg_name,plotbasedir,ibv,figprefix='None',figsize='None',figextension='None',debug=False)  :

 #Norm_o es la perturbacion evolucionada en t 
 #Norm_i es de donde partio la perturbacion evolucioanda en t-1
 #Norm_r es la perturbacion o rescalada en t.
 #Debido al bip, la serie de BV no se conecta en el tiempo.

 if reg_name == 'None' :
    plot_reg_name = False
 if figsize == 'None' :
    figsize=default_figure_size
 if figextension == 'None' :
    figextension='.png'

 #Create output directory
 if not os.path.exists(plotbasedir):
   os.makedirs(plotbasedir)


#Creates plots of bv norm 
 tmp_shape=norm_o[var_name[0]].shape
 ntimes=tmp_shape[0]
 nbv=tmp_shape[1]
 nregs=tmp_shape[2]

 #Create time series.
 time_serie=np.zeros([ntimes*4,nregs])
 times=np.zeros(ntimes*4)

 for myvar in var_name     :
  for ii in range(1,ntimes)  :
   #Creamos un grafico tipo serrucho para representar la evolucion de la norma.
   time_serie[4*ii,:]=norm_i[myvar][ii-1,ibv,:]
   time_serie[4*ii+1,:]=norm_o[myvar][ii,ibv,:]
   time_serie[4*ii+2,:]=norm_r[myvar][ii,ibv,:]
   time_serie[4*ii+3,:]=np.nan

   times[4*ii]=ii-1
   times[4*ii+1]=ii
   times[4*ii+2]=ii
   times[4*ii+3]=ii

  for ireg in range(0,nregs) : 
   iregstr="%04d" % ( ireg + 1 )

   fig=plt.figure(1,figsize=figsize)

   #for it in range(0,niter)     :
   plt.plot(times,time_serie[:,ireg],'-')

   plt.ylabel('Norm')
   plt.xlabel('Time')

   if debug == True   :
      plt.show()

   print( 'Generationg the following figure : ' + 'Figure_' + figprefix + myvar + 'reg' + iregstr + figextension )
   plt.savefig( plotbasedir + 'Figure_' + figprefix + myvar + 'reg' + iregstr + figextension )

   plt.close(fig)
    
#plt.plot(times,meandisttrack[:,0],'r-',linewidth=3,label='CTRL')
#ax.fill_between(times, meandisttrack[:,0]+errorbar[:,0],meandisttrack[:,0]-errorbar[:,0], facecolor='red', alpha=0.1)

#plt.plot(times,meandisttrack[:,1],'b-',linewidth=3,label='PE')
#ax.fill_between(times, meandisttrack[:,1]+errorbar[:,1],meandisttrack[:,1]-errorbar[:,1], facecolor='blue', alpha=0.1)

#plt.legend(fontsize=14,loc='upper right')

















    

