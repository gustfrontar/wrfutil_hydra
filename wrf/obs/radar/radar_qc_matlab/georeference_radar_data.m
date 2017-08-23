function [radar]=georeference_radar_data(radar)

%Radar is a structure. Function outputs are added to this same structure.
%Sounding is a structure with fields, temperature, pressure, moisture and
%refractivity index that we will need to compute the height of the radar
%beam. At input only sounding.filename has to be defined.

d2r=pi/180; %From deg to radian.
Re=6370e3;  %Earth radius in (m).
Ns=1.21;     
ke=(4/3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEIGTH COMPUTATION FOR EACH RANGE AND ELEVATION
% If no temperature profile is available standard WS-88 RADAR formula is
% used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We do not have sounding information. Will assume standard atmospheric conditions.    
   
   radar.Z=NaN(radar.na,radar.nr,radar.ne);
   for ii=1:radar.nr;
       for kk=1:radar.ne;
              %Two possibilities the second probably more accurate than the
              %first one, but both assume a standard constant value for the
              %refractivity index.
              %radar.height(ii,kk)=radar.altitude + radar.radius(ii)*sin(radar.elev(kk)*d2r) +(radar.radius(ii)^2)/(2*Ns*Re); %
              radar.Z(:,ii,kk)=radar.altitude + sqrt(radar.range(ii)^2+(ke*Re)^2+2*radar.range(ii)*ke*Re*sin(radar.elevation(kk)*d2r))-ke*Re;
       end
   end

radar.height=squeeze(radar.Z(1,:,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LAT LON computation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radar.latitude=NaN(radar.na,radar.nr,radar.ne);
radar.longitude=NaN(radar.na,radar.nr,radar.ne);


%We are assuming round earth with radius Re
%We are assuming approximately sraight beam propagation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radar.Rh=NaN(radar.nr,radar.ne);


%Esta parte del script esta modificada de forma tal de usar el azimuth exacto para calcular
%la longitud y la latitud. 
for kk=1:radar.ne
  for ii=1:radar.na
    %The curvature of the radar beam is not taken into account.
    distance_on_earth=ke*Re*asin(radar.range.*cos(radar.elevation(kk)*d2r)/(ke*Re));
    radar.Rh(:,kk)=squeeze(distance_on_earth); %Radar horizontal range.

    distance_on_earth_km=distance_on_earth/1000;
    az_vector=radar.azimuth_exact(ii,kk)*ones(size(distance_on_earth_km));
    
    %[distance_matrix , az_matrix]=meshgrid(distance_on_earth/1000,radar.azimuth);
    
    [radar.longitude(ii,:,kk) , radar.latitude(ii,:,kk)]=ll_arc_distance(radar.lon,radar.lat,distance_on_earth_km,az_vector);

  end
end










