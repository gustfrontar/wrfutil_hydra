function [echo_top , echo_base , echo_depth , max_dbz , max_dbz_z , tmp_elevation , tmp_range ]=compute_multiple_echotop(radar,data,zthresh,tmp_elevation,tmp_range)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Compute a local echo top and echo base for each radar gate (taking into
% account multiple tops as in the case of low level echo overrun by cirrus
% clouds.
% Juan Ruiz 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Necho_top=4;    %Maximum number of echo tops that can be detected.


%Allocate space

echo_top=NaN(size(data));
echo_base=NaN(size(data));
echo_depth=NaN(size(data));
max_dbz=NaN(size(data));
max_dbz_z=NaN(size(data));
%vertical_ref_grad=NaN(size(data));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERP EACH LEVEL REFLECTIVITY TO THE FIRST LEVEL GRID.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_z=20000;
dz=100;

tmp_radar_r=radar.Rh(:,1);
tmp_radar_z=0:dz:max_z;
tmp_radar_data=NaN(radar.na,radar.nr,length(tmp_radar_z));

%Replace NaN values with vertical interpolation (where possible)
for ie=2:radar.ne-1
    tmp_data = data(:,:,ie);
    tmp_data_po= data(:,:,ie+1);
    tmp_data_mo= data(:,:,ie+1);
    nan_index=isnan( tmp_data );
    
    tmp_data( nan_index )= 0.5 * ( tmp_data_po( nan_index) + tmp_data_mo( nan_index ) );
    
    data(:,:,ie) = tmp_data ;
end

%vertical_ref_grad=NaN(radar.na,radar.nr,length(tmp_radar_z));
[tmp_radar_r tmp_radar_z]=meshgrid(squeeze(tmp_radar_r),squeeze(tmp_radar_z));

%data(isnan(data))=radar.replacerefmissing;   %Replace missing values prior
%tmp_radar_data2=NaN(radar.na,radar.nr,length(new_z));

[range elevation]=meshgrid(radar.range,radar.elevation);

    %Compute elevation and range of the new grid points so I can
    %interpolate them using regular grid functions.
    if( isempty(tmp_elevation) || isempty(tmp_range) )

    tmp_elevation=griddata(radar.Rh',radar.height',elevation,tmp_radar_r,tmp_radar_z,'linear');
    tmp_range    =griddata(radar.Rh',radar.height',range,tmp_radar_r,tmp_radar_z,'linear');
    end
%figure;pcolor(tmp_elevation);shading flat
%figure;pcolor(tmp_range);shading flat
                                                                                        
warning off %Supress warning about NaN in interpolated fields.
%Loop over elevations

not_nan=~isnan(tmp_elevation) & ~isnan(tmp_range) ;

tmpe=tmp_elevation(not_nan);
tmpr=tmp_range(not_nan);

% figure;pcolor(range);shading flat
% figure;pcolor(elevation);shading flat

for ia=1:radar.na
 tmp=squeeze(tmp_radar_data(ia,:,:))';    
 tmp2=interp2(range,elevation,squeeze(data(ia,:,:))',tmpr,tmpe,'bilinear');
 tmp(not_nan)=tmp2;
 tmp_radar_data(ia,:,:)=tmp';
end

% 
% figure
% pcolor(squeeze(tmp_radar_data(ia,:,:))');shading flat

%figure;pcolor(range);shading flat
%figure;pcolor(elevation);shading flat
%figure;pcolor(range,elevation,squeeze(data(ia,:,:))');shading flat
%figure;pcolor(squeeze(tmp_radar_data(ia,:,:))');shading flat

%sazimuth=280;
%pcolor(squeeze(tmp_radar_data(sazimuth,:,:)));shading flat
%Loop to analize the vertical profile of reflectivity.
tmp_echo_top_3d=NaN(size(tmp_radar_z));
tmp_echo_base_3d=NaN(size(tmp_radar_z));
tmp_max_dbz_3d=NaN(size(tmp_radar_z));
tmp_max_dbz_z_3d=NaN(size(tmp_radar_z));
%tmp_vertical_ref_grad_3d=NaN(size(tmp_radar_z));

for ia=1:radar.na
    for ir=1:radar.nr

    [tmp_echo_top,tmp_echo_base,tmp_max_dbz,tmp_max_dbz_z] ...
          =compute_echo_top_singlecolumn(tmp_radar_z(:,1),tmp_radar_data(ia,ir,:),ones(size(tmp_radar_data,3),1),Necho_top,zthresh);
          for itop=1:Necho_top
             if( ~isnan(tmp_echo_top(itop)) && ~isnan(tmp_echo_base(itop)) )
              echo_top_index=ceil(tmp_echo_top(itop)/dz);
              echo_base_index=floor(tmp_echo_base(itop)/dz);
      
              tmp_echo_top_3d(echo_base_index:echo_top_index,ir)=tmp_echo_top(itop);
              tmp_echo_base_3d(echo_base_index:echo_top_index,ir)=tmp_echo_base(itop);
              tmp_max_dbz_3d(echo_base_index:echo_top_index,ir)=tmp_max_dbz(itop);
              tmp_max_dbz_z_3d(echo_base_index:echo_top_index,ir)=tmp_max_dbz_z(itop);
              %tmp_vertical_ref_grad_3d(ir,echo_base_index:echo_top_index)=tmp_vertical_ref_grad(itop);
             end
          end 
    end
   %Interpolate back to the radar grid using nearest neighbor interpolation.
   echo_top(ia,:,:)         =interp2(tmp_radar_r,tmp_radar_z,tmp_echo_top_3d,radar.Rh',radar.height','nearest')';
   echo_base(ia,:,:)        =interp2(tmp_radar_r,tmp_radar_z,tmp_echo_base_3d,radar.Rh',radar.height','nearest')';
   max_dbz(ia,:,:)          =interp2(tmp_radar_r,tmp_radar_z,tmp_max_dbz_3d,radar.Rh',radar.height','nearest')';
   max_dbz_z(ia,:,:)        =interp2(tmp_radar_r,tmp_radar_z,tmp_max_dbz_z_3d,radar.Rh',radar.height','nearest')';
   %vertical_ref_grad(ia,:,:)=interp2(tmp_radar_z,tmp_radar_r,tmp_vertical_ref_grad_3d,radar.Rh,radar.height,'nearest'); 

end
% 
% figure
% pcolor(squeeze(tmp_echo_top_3d(:,:)));shading flat
% 
% figure
% pcolor(squeeze(echo_top(ia,:,:))');shading flat


echo_depth=echo_top-echo_base;
end
    
    












