

function [u , v , z_vad , u_error , v_error , uo , vo , data_coverage , fit]=compute_vad(radar , wind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes a vad wind profile that will be used for wind qc.
% the computation follows the algorithm proposed by Gao and Drogemeier 2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% radar, an structure with the information about the radar grid.
% wind, radial wind data in the original radar grid ( azimuth , range ,
% elevation );
% Outuput
% u , v a vertical profile of u and v that results from the combination of
% data coming from different ranges and elevation angles.
% 

tic

dvrdtitathreshold= 1e2;   %Derivadas mayores que este umbral no van a ser usadas para el VAD.
vertical_resolution=100;  %Vertical resolution of the resulting vad profile in meters.
data_coverage_threshold=0.0; %Ranges with less that 10% of valid data points won't be used for VAD.
%Note: Some mesoscale circulations can produce such gradients, this does
%not mean that they will be removed, the won't be considered for the 
%computation of the large scale wind component at this point.
max_vad_h=22000;             %Maximum height at wich VAD will be computed.
compute_from_vr=true;      %Wether vad is computed from VR or from DVR/DTITA
max_elevation_angle=80;    %Do not consider levels over this elevation angle.


%Find the angles where VAD algorithm will be applied.
index_angles=radar.elevation < max_elevation_angle;
wind=wind(:,:,index_angles);
elevation_angles=radar.elevation(index_angles);

deg2rad=pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute radial veloctity derivatives with respect to azimuth and remove
% outliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( compute_from_vr )

f=wind;
tita=deg2rad*radar.azimuth;

else
dtita=radar.azimuth(2)-radar.azimuth(1);

f=(wind-circshift(wind,[1 0 0]))/(dtita*deg2rad);

%tita=deg2rad*(radar.azimuth+circshift(radar.azimuth,[0 1]))/2;
%tita(1)=deg2rad*(2*radar.azimuth(1)-dtita)/2;
tita=deg2rad*radar.azimuth;

%Remove outliers (due to aliasing or due to small scale features);
%f( abs(f) > dvrdtitathreshold )=NaN;

plot( f(:,200,20) )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit a trigonometric function to find u and v for each range and elevation
% angle.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ alfa , beta , data_coverage , fit ] = compute_trigonometric_fit( f , tita );

if ( compute_from_vr )
    uo=alfa;
    vo=beta;
else  
  uo=beta;
  vo=-alfa;  
end

for ii=1:length(elevation_angles)
   uo(:,ii)=uo(:,ii)/cos(deg2rad*elevation_angles(ii));
   vo(:,ii)=vo(:,ii)/cos(deg2rad*elevation_angles(ii)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine the information provided by each elevation and range to get a
% single profile.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Remove values that have been obtained with only a few pixels.
tmp_index= data_coverage < data_coverage_threshold ;

uo( tmp_index ) = NaN;
vo( tmp_index ) = NaN;
data_coverage(tmp_index)=0;

%Compute a weigth for each range and elevation angle tacking into account
%the number of valid points used in each case and the elevation angle.
weight=data_coverage;
%As the elevation angle increases the proyection of the horizontal wind on
%radial velocity is reduced and errors can be greather.
%for ie=1:length(elevation_angles);
%    weight(:,ie)=weight(:,ie)*cos(elevation_angles(ie)*deg2rad);
%end

%Generate a common vertical grid.

z_vad=0:vertical_resolution:max_vad_h;

uo_vad=NaN(length(z_vad),length(elevation_angles));
vo_vad=NaN(length(z_vad),length(elevation_angles));
w_vad=NaN(length(z_vad),length(elevation_angles));

%Interpolate uo, vo and weights to the common vertical coordinate.
for ie=1:length(elevation_angles)
   uo_vad(:,ie)=interp1(squeeze(radar.height(:,ie)),squeeze(uo(:,ie)),z_vad);
   vo_vad(:,ie)=interp1(squeeze(radar.height(:,ie)),squeeze(vo(:,ie)),z_vad);
   w_vad(:,ie) =interp1(squeeze(radar.height(:,ie)),squeeze(weight(:,ie)),z_vad);   
end

%Now we have several uo and vo estimates for each level. We can compute the
%mean and the standard deviation to have a unique VAD profile for each
%heigth.

tmp_nan=isnan( uo_vad ) | isnan( vo_vad ) | isnan(w_vad);
uo_vad(tmp_nan)=NaN;
vo_vad(tmp_nan)=NaN;
w_vad(tmp_nan)=NaN;

u=nanmean( uo_vad .* w_vad , 2 )./nanmean( w_vad , 2 );
v=nanmean( vo_vad .* w_vad , 2 )./nanmean( w_vad , 2 );

u_error = ( squeeze( nanmean(  (uo_vad.^2) .* w_vad ,2 )./nanmean( w_vad ,2) ) - u.^2 ).^0.5;
v_error = ( squeeze( nanmean(  (vo_vad.^2) .* w_vad ,2 )./nanmean( w_vad ,2) ) - v.^2 ).^0.5;

time=toc;

display(['VAD was computed in ' num2str(time) 'seconds.'])

end