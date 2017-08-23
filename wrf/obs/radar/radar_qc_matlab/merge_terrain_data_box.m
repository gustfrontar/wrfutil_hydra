clear all
close all

%This script reads in terrain data in geotiff format provided by 
%http://gdem.ersdac.jspacesystems.or.jp/outline.jsp

%ASTER GDEM Ver. 2 was produced in cooperation with the Japan-US ASTER Science Team using ASTER data which had been acquired from the beginning of the observation until the end of August, 2010. 
%Since the validation of ASTER GDEM and the preparation of distribution site were finished, ASTER GDEM is released now and its distribution started. 
%For the validation result and quality of ASTER GDEM, refer to the ASTER GDEM Project page. 

%In this version the data is interpolated to a lower resolution (100m x 100m) closer to the model grid and also to the radar pulse length.

% Juan Ruiz 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE GRID OF INTEREST... 1 DEGREE SURROUNDING OSAKA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

osaka_lat=34.823;
osaka_lon=135.523;
Re=6370e3;
degdx=Re*pi/180;
radar_r_res=100; %Radar pixel resolution in r direction.

lat_radar=osaka_lat+(-1:radar_r_res/degdx:1);
lon_radar=osaka_lon+(-1:radar_r_res/(degdx*cos(osaka_lat*pi/180)):1);

max_osaka_lon=max(osaka_lon);
min_osaka_lon=min(osaka_lon);

max_osaka_lat=max(osaka_lat);
min_osaka_lat=min(osaka_lat);

[lon_radar_mat lat_radar_mat]=meshgrid(lon_radar,lat_radar);
osaka_meantopography=NaN(size(lat_radar_mat));
osaka_maxtopography=NaN(size(lat_radar_mat));
osaka_mintopography=NaN(size(lat_radar_mat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Interpolate digital terrain information into the desired area.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We will read different data tiles and interpolate them to the desired
%area.

for it=33:35 %Tile loop
 for jt=33:36 %Tile loop

   geotiff_data_path='../../terrain_data/Tiles_201402031107/';
   geofile=[geotiff_data_path 'ASTGTM2_N' num2str(it) 'E1' num2str(jt) '/ASTGTM2_N' num2str(it) 'E1' num2str(jt) '_dem.mat'];

   load(geofile,'lon','lat','topo');
   topo=flipdim(topo,1);
   lon_tile=lon;
   lat_tile=lat;
   [ny nx]=size(lon_radar_mat);

   tmp_meantopo=NaN(size(osaka_meantopography));
   tmp_maxtopo=NaN(size(osaka_meantopography));
   tmp_mintopo=NaN(size(osaka_meantopography));

     %Get tile corners.
      max_lon_tile=max(lon_tile);
      min_lon_tile=min(lon_tile);
      max_lat_tile=max(lat_tile);
      min_lat_tile=min(lat_tile);
      [nx_tile ny_tile]=size(topo);
      x_res_tile=(lon_tile(2)-lon_tile(1));
      y_res_tile=(lat_tile(2)-lat_tile(1));
      x_res_radar=(lon_radar(2)-lon_radar(1));
      y_res_radar=(lat_radar(2)-lat_radar(1));

     for ilat=1:ny
          for ilon=1:nx
         %Define box limits assuming lat-lon regular grid.
         lat_s=lat_radar_mat(ilat,ilon)-0.5*y_res_radar;
         lat_n=lat_radar_mat(ilat,ilon)+0.5*y_res_radar;
         lon_w=lon_radar_mat(ilat,ilon)-0.5*x_res_radar;
         lon_e=lon_radar_mat(ilat,ilon)+0.5*x_res_radar;

         %Compute points that fall within the desidered box
          j_s=(lat_s-min_lat_tile)/y_res_tile+1;
          j_e=(lat_n-min_lat_tile)/y_res_tile+1;
          i_s=(lon_w-min_lon_tile)/x_res_tile+1;
          i_e=(lon_e-min_lon_tile)/x_res_tile+1;
        
          j_e=floor(j_e);
          i_e=floor(i_e);
          i_s= ceil(i_s);
          j_s= ceil(j_s); 
         %Check if there are tile points for this box:
          compute_topo=true;
          if ( j_e < 1 | j_s > ny_tile | i_e < 1 | i_s > nx_tile )
             compute_topo=false;   %The tile does not cover this grid point.
          end   
         if ( compute_topo ) 
         %Obtain topography statistics for this box
           if( j_e > ny_tile);j_e=ny_tile;end
           if( j_s < 1      );j_s=1      ;end
           if( i_e > nx_tile);i_e=nx_tile;end
           if( i_s < 1      );i_s=1      ;end
         %*****************************************************************
         tmp_meantopo(ilat,ilon)=nan_safe_ave(nan_safe_ave(topo(j_s:j_e,i_s:i_e)));
         tmp_maxtopo(ilat,ilon) =max         (         max(topo(j_s:j_e,i_s:i_e)));
         tmp_mintopo(ilat,ilon) =min         (         min(topo(j_s:j_e,i_s:i_e)));

         end
         %*****************************************************************
         
        
                  
          end
      end

%tmp_topo=interp2(lon_mat,lat_mat,topo,lon_radar_mat,lat_radar_mat,'bilinear');

osaka_meantopography(~isnan(tmp_meantopo))=0;
tmp_meantopo(isnan(tmp_meantopo))=0;
osaka_meantopography=osaka_meantopography+tmp_meantopo;

osaka_maxtopography(~isnan(tmp_maxtopo))=0;
tmp_maxtopo(isnan(tmp_maxtopo))=0;
osaka_maxtopography=osaka_maxtopography+tmp_maxtopo;

osaka_mintopography(~isnan(tmp_mintopo))=0;
tmp_mintopo(isnan(tmp_mintopo))=0;
osaka_mintopography=osaka_mintopography+tmp_mintopo;

display('The current coverage of the desired area is:...')
100-100*sum(sum(isnan(osaka_meantopography)))/numel(osaka_meantopography)

 end
end

terrain.meanheight=osaka_meantopography;
clear osaka_meantopography
terrain.maxheight=osaka_maxtopography;
clear osaka_maxtopography
terrain.minheight=osaka_mintopography;
clear osaka_mintopography
terrain.lat=lat_radar;
clear lat_radar;
terrain.lon=lon_radar;
clear lon_radar;

save('osaka_terrain.mat','terrain');


