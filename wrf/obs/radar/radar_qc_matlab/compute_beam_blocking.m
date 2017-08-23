

function [max_blocking , blocking, correction]=compute_beam_blocking(radar,terrain)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the blocking of the terrain.
% max_blocking is the along radius maximum blocking encountered by the
% radar beam in its propagation away from the radar antena.
% blocking is the estimated local blocking at each radar gate. 
% The blocking is estimated as the percentage of the main beam blocked by
% the terrain (the Gaussian distribution of the power within the main beam
% is not taken into account).
% Note: As radar height might strongly depend on the meteorological
% conditions it is advaisable to compute them using information of the
% vertical profile of temperature and moisture in order to get accurate
% results of the blocking percentage.
% A simple correction for the blocking is comptued following Fulton et al
% 1998 (WAF) that has been developed for C band radars. It is not clear if
% this correction can be applied as it is for X band radars.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the max_vertical_beam_extent as a function of range and elevation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[max_vertical_beam_extent]=compute_max_vertical_beam_extent(radar);

max_vertical_beam_extent=repmat(max_vertical_beam_extent,[1 1 size(terrain,1)]);
max_vertical_beam_extent=permute(max_vertical_beam_extent,[3 1 2]);

% Compute the height over the terrain and normalize it by the vertical
% extent of the beam (so the beam is circula and has radius equal to one).

terrain=(radar.Z-terrain)./max_vertical_beam_extent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Precompute the blocking of a circular beam of unity radius.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z_ref=[-1:0.1:1];
b_ref=( -z_ref.*sqrt( 1-z_ref.^2) - asin(z_ref) + pi/2 )/pi;

z_ref=[-9e9 z_ref 9e9];
b_ref=[1 b_ref 0];

blocking=interp1(z_ref,b_ref,terrain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loop over radius to find cumulative blocking (i.e. max blocking reached
%along a ray.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_blocking=zeros(radar.na,radar.nr,radar.ne);
max_blocking(:,1,:)=blocking(:,1,:);
for ir=2:radar.nr
    
    max_blocking(:,ir,:)=max_blocking(:,ir-1,:);
    tmp_data=max_blocking(:,ir,:);
    tmp_data2=blocking(:,ir,:);
    blocking_index= tmp_data2 >= tmp_data;
    tmp_data( blocking_index )=tmp_data2( blocking_index );
    max_blocking(:,ir,:)=tmp_data;
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Precompute radar correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b_ref=[0 0.1 0.3 0.4 0.5 0.6 1];
c_ref=[0 0   1   2   3   0   0]; %Reflectivity values that has to be added after the blocking.

%time=toc;
%display(['Terrain blocking was computed in ' num2str(time) ' seconds'])

correction=interp1(b_ref,c_ref,max_blocking);

end