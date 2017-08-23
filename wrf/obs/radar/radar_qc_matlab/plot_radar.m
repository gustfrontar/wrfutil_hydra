
clear all
close all

path='/home/jruiz/Dropbox/DATA/DATOS_RADAR/ANGUIL/20100111/120/';

files=dir([path '*.nc']);

for ifile=1:size(files,1)
    
    filename=[ path files(ifile).name ];
    
    %Variables
    ref = read_netcdf_var(filename,'dBZ',true);
    doppler = read_netcdf_var(filename,'V',true);
    phidp = read_netcdf_var(filename,'PhiDP',true);
    %kdp = read_netcdf_var(filename,'KDP',true);
    rhohv = read_netcdf_var(filename,'RhoHV',true);
    
    %Dimensions
    elevation  = read_netcdf_var(filename,'elevation',false);
    azimuth    = read_netcdf_var(filename,'azimuth',false);
    range      = read_netcdf_var(filename,'range',false);
    time       = read_netcdf_var(filename,'time',false);
    
    %Radar location
    radar_lon = read_netcdf_var(filename,'longitude',false);
    radar_lat = read_netcdf_var(filename,'latitude',false);
    radar_altitude = read_netcdf_var(filename,'altitude',false);
    
    %Scan time 
    time_coverage_start = read_netcdf_var(filename,'time_coverage_start',false)';
    time_coverage_start([14 17])='_';
    time_coverage_start=time_coverage_start(1:20)
    
   
    [order_ref , az , el]= order_variable(ref,azimuth,elevation) ;
    [order_doppler , az , el] = order_variable(doppler,azimuth,elevation) ;
    [order_phidp , az , el] = order_variable(phidp,azimuth,elevation) ;
    %[order_kdp , az , el]= order_variable(kdp,azimuth,elevation) ;
    [order_rhohv , az , el] = order_variable(rhohv,azimuth,elevation) ;
    
    
   
    [latitude , longitude , z]=georeference_radar_data(el,az,range,radar_altitude,radar_lon,radar_lat);
    
    order_doppler( order_doppler < -40 )=NaN;
    order_ref    ( order_ref    == -32 )=NaN;
    order_rhohv  ( abs(order_rhohv)  <  0.01 )=NaN;
    %order_kdp    ( order_kdp    <  -12.0  )=NaN;
    order_phidp  ( order_phidp    < 0.0 )=NaN;
    
    for iel=1:length(el)
       mkdir([ path '/elev_' num2str(iel) ]);
       figure
       subplot(2,2,1)
          pcolor(longitude(:,:,iel),latitude(:,:,iel),order_ref(:,:,iel));shading flat
          
          caxis([0 80]);
          colorbar
          title('Reflectivity')
          
       subplot(2,2,2)
          pcolor(longitude(:,:,iel),latitude(:,:,iel),order_doppler(:,:,iel));shading flat
          
          caxis([-30 30]);
          colorbar 
          title('Doppler')
          
       subplot(2,2,3)
          pcolor(longitude(:,:,iel),latitude(:,:,iel),order_rhohv(:,:,iel));shading flat
          
          caxis([-1 1]);
          colorbar
          title('RhoHV')
          
       subplot(2,2,4)
          pcolor(longitude(:,:,iel),latitude(:,:,iel),order_phidp(:,:,iel));shading flat
          
          caxis([0 360]);
          colorbar
          title('PhiDP')
          
       figurefile=[path '/elev_' num2str(iel) '/plot_' time_coverage_start '.png'];
       
       print('-dpng',figurefile);
       
       close all
    end
    
end