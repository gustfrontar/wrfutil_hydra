
clear all
close all

%This script reads in terrain data in geotiff format provided by 
%http://gdem.ersdac.jspacesystems.or.jp/outline.jsp

%ASTER GDEM Ver. 2 was produced in cooperation with the Japan-US ASTER Science Team using ASTER data which had been acquired from the beginning of the observation until the end of August, 2010. 
%Since the validation of ASTER GDEM and the preparation of distribution site were finished, ASTER GDEM is released now and its distribution started. 
%For the validation result and quality of ASTER GDEM, refer to the ASTER GDEM Project page. 

% Juan Ruiz 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE GRID OF INTEREST... 1 DEGREE SURROUNDING OSAKA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

osaka_lat=34.823;
osaka_lon=135.523;

lat_radar=osaka_lat+(-1:1/3600:1);
lon_radar=osaka_lon+(-1:1/3600:1);

max_osaka_lon=max(osaka_lon);
min_osaka_lon=min(osaka_lon);

max_osaka_lat=max(osaka_lat);
min_osaka_lat=min(osaka_lat);

[lon_radar_mat lat_radar_mat]=meshgrid(lon_radar,lat_radar);
osaka_topography=zeros(size(lat_radar_mat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Interpolate digital terrain information into the desired area.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We will read different data tiles and interpolate them to the desired
%area.

geotiff_data_path='../../terrain_data/Tiles_201402031107/';

geofile=[geotiff_data_path 'ASTGTM2_N33E133/ASTGTM2_N33E133_num.tif'];

infofile = geotiffinfo(geofile);
max_lat=max(infofile.CornerCoords.Lat);
min_lat=min(infofile.CornerCoords.Lat);
max_lon=max(infofile.CornerCoords.Lon);
min_lon=min(infofile.CornerCoords.Lon);

[topo] = geotiffread(geofile);
topo=double(topo);
[nx ny]=size(topo);
lat=min_lat+((1:ny)-1)*(max_lat-min_lat)/(ny-1);
lon=min_lon+((1:nx)-1)*(max_lon-min_lon)/(nx-1);
[lon_mat lat_mat]=meshgrid(lon,lat);

tmp_topo=interp2(lon_mat,lat_mat,topo,lon_radar_mat,lat_radar_mat,'bilinear');







