

function [ radar , ref , dv , dvc , phidp , rhohv , kdp , present_ref , present_dv , present_dvc , present_phidp , present_rhohv , present_kdp ] = read_radar_netcdf(filename , minref )


    [present_ref]=test_netcdf_var(filename,'dBZ');
    [present_dv]=test_netcdf_var(filename,'V');
    [present_dvc]=test_netcdf_var(filename,'Vda');
    [present_phidp]=test_netcdf_var(filename,'PhiDP');
    [present_rhohv]=test_netcdf_var(filename,'RhoHV');
    [present_kdp]=test_netcdf_var(filename,'KDP');

%Get dimensions
    elevation  = double( read_netcdf_var(filename,'elevation',false) ) ;
    azimuth    = double( read_netcdf_var(filename,'azimuth',false) ) ;
    range      = double( read_netcdf_var(filename,'range',false) ) ;
    time       = double( read_netcdf_var(filename,'time',false) ) ;

%Get Radar location
    radar.lon = read_netcdf_var(filename,'longitude',false);
    radar.lat = read_netcdf_var(filename,'latitude',false);
    radar.altitude = read_netcdf_var(filename,'altitude',false);

%Get Scan time
    radar.time_coverage_start = read_netcdf_var(filename,'time_coverage_start',false)';
    radar.time_coverage_end = read_netcdf_var(filename,'time_coverage_end',false)';
    radar.time_coverage_start([14 17])='_';
    radar.time_coverage_start=radar.time_coverage_start(1:20);
    radar.time_coverage_end([14 17])='_';
    radar.time_coverage_end=radar.time_coverage_end(1:20);



%Get Variable Data
    if( present_ref )
       input_ref = read_netcdf_var(filename,'dBZ',true);
       input_ref ( isnan( input_ref ) )= minref ;
    else
       input_ref = NaN( size(range,1), size(elevation,1) );
    end
    if( present_dv )
       input_dv = read_netcdf_var(filename,'V',true);
    else
       input_dv = NaN( size(range,1), size(elevation,1) );
    end
    if( present_phidp )
       input_phidp = read_netcdf_var(filename,'PhiDP',true);
    else
       input_phidp = NaN( size(range,1), size(elevation,1) );
    end
    if( present_kdp )
       input_kdp = read_netcdf_var(filename,'KDP',true);
    else
       input_kdp = NaN( size(range,1), size(elevation,1) );
    end
    if( present_rhohv )
       input_rhohv = read_netcdf_var(filename,'RhoHV',true);
    else
       input_rhohv = NaN( size(range,1), size(elevation,1) );
    end
    if( present_dvc )
       input_dvc = read_netcdf_var(filename,'Vda',false);
       input_dvc( input_dvc == max(max(max( input_dvc) )) )=NaN;
       input_dvc( input_dvc == min(min(min( input_dvc) )) )=NaN;
    else
       input_dvc = NaN( size(range,1), size(elevation,1) );
    end

    ray_angle_res = read_netcdf_var(filename,'ray_angle_res',false);
    if( length( unique( ray_angle_res ) ) >= 2 )
      display(['Warning: La resolucion en azimuth no es uniforme en los diferentes angulos de elevacion '])
      display(['Warning: El codigo no esta preparado para considerar este caso y puede producir efectos indesaedos '])
    end
    radar.ray_angle_res=nanmean(ray_angle_res);
    display(['La resolucion en azimuth es: ' num2str(radar.ray_angle_res)])
  
    %Llamamos a las rutinas que ordenan los datos.
    %Las rutinas que ordenan los datos usaran una reticula regular en rango, azimuth y elevation.
    %esta reticula regular nominal esta definida en base a la resolucion en grados que indica el radar.
    %la ubicacion exacta de cada rayo se guarda en el arreglo az_exact que puede usarse en toda aquella 
    %rutian que requiera geolocalizar los datos (por ejemplo para la preparacion de los archivos para asimilacion).
   
    [ref , az , el , ti , az_exact ]= order_variable(input_ref,azimuth,elevation,time,radar.ray_angle_res);
    [dv , az , el , ti , az_exact ] = order_variable(input_dv,azimuth,elevation,time,radar.ray_angle_res);
    [phidp , az , el , ti , az_exact ] = order_variable(input_phidp,azimuth,elevation,time,radar.ray_angle_res);
    [rhohv , az , el , ti , az_exact ] = order_variable(input_rhohv,azimuth,elevation,time,radar.ray_angle_res);
    [dvc , az , el , ti , az_exact ] = order_variable(input_dvc,azimuth,elevation,time,radar.ray_angle_res);
    [kdp , az , el , ti , az_exact ] = order_variable(input_kdp,azimuth,elevation,time,radar.ray_angle_res);

    radar.ne=length(el);
    radar.na=length(az);
    radar.nr=length(range);

    radar.elevation=el;
    radar.azimuth=az;
    radar.azimuth_exact=az_exact;
    radar.range=range;
    radar.time=ti;


    radar.beam_wid_v=read_netcdf_var(filename,'radar_beam_width_v',false);
    radar.beam_wid_h=read_netcdf_var(filename,'radar_beam_width_h',false);

end
