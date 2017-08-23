function [echo_top_3d echo_depth_3d echo_depth_size_3d total_max_dbz_3d total_max_dbz_level_3d max_dbz_3d max_dbz_level_3d]=compute_echotop3d(radar,reflectivity,thresholds,snx,sny)

% This function computes the echo top and interpolate it to the grid
% corresponding to each elevation in order to perform QC operation.
% The inputs are: radar, the structure containing radar info
% (geolocalization, etc), reflectivity a 3D array with reflectivity data in
% the radar original grid (azimuth, r, elevation) and thresholds is a
% vector of thresholds, echo top will be computed and interpolated for all
% the selected threshholds.
% snx and sny control the amount of box based horizontal smooth. Smoothing
% will be performed using a movil box average of size (2*nx+1)*(2*ny+1);
% Output: The output is the echo_top_3d which is the echo top interpolated
% to the grid of each radar elevation for QC opeartions. It also provides
% the echo_depth which is the difference between echo_top and echo_base
% (take into account that depending on the way that echo_base and echo_top
% are computed this could mean different things).
% 
% By Juan Ruiz 2014

nt=length(thresholds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   COMPUTE ECHO TOP AND ECHO BASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo_top_3d=NaN(radar.na,radar.nr,radar.ne,nt);
echo_depth_3d=NaN(radar.na,radar.nr,radar.ne,nt);
max_dbz_3d=NaN(radar.na,radar.nr,radar.ne,nt);
max_dbz_level_3d=NaN(radar.na,radar.nr,radar.ne,nt);
total_max_dbz_3d=NaN(radar.na,radar.nr,radar.ne,nt);
total_max_dbz_level_3d=NaN(radar.na,radar.nr,radar.ne,nt);
echo_depth_size_3d=NaN(radar.na,radar.nr,radar.ne,nt);

[echo_top echo_depth total_max_dbz total_max_dbz_level max_dbz max_dbz_level]=compute_echotop(radar,reflectivity,thresholds);
echo_depth=echo_top-echo_depth;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PERFORM BOX AVERAGE SMOOTHING (IF REQUESTED)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[echo_depth_size]=compute_boxcount(echo_depth,snx,sny,0,1000);
[echo_top]=compute_boxmean(echo_top,snx,sny,0);
[echo_depth]=compute_boxmean(echo_depth,snx,sny,0);
[max_dbz]=compute_boxmean(max_dbz,snx,sny,0);
[max_dbz_level]=compute_boxmean(max_dbz_level,snx,sny,0);
[total_max_dbz]=compute_boxmean(total_max_dbz,snx,sny,0);
[total_max_dbz_level]=compute_boxmean(total_max_dbz_level,snx,sny,0);

echo_top(isnan(echo_top))=0;
echo_depth(isnan(echo_depth))=0;
max_dbz(isnan(max_dbz))=0;
max_dbz_level(isnan(max_dbz_level))=0;
total_max_dbz(isnan(total_max_dbz))=0;
total_max_dbz_level(isnan(total_max_dbz_level))=0;
echo_depth(isnan(echo_depth))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INTERPOLATE FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo_depth_3d(:,:,1,:)=echo_depth;
echo_top_3d(:,:,1,:)=echo_top;
max_dbz_3d(:,:,1,:)=max_dbz;
max_dbz_level_3d(:,:,1,:)=max_dbz_level;
total_max_dbz_3d(:,:,1,:)=total_max_dbz;
total_max_dbz_level_3d(:,:,1,:)=total_max_dbz_level;
echo_depth_size_3d(:,:,1,:)=echo_depth_size;

%Interpolate echo_top

   %Artificially extend data to radar location
   extended_field=NaN(size(echo_top_3d,1),size(echo_top_3d,2)+1,size(echo_top_3d,4));
   extended_range=NaN(size(radar.Rh,1)+1,size(radar.Rh,2));
   extended_field(:,2:end,:)=squeeze(echo_top_3d(:,:,1,:));
   extended_range(2:end,:)=radar.Rh;
   extended_field(:,1,:)=repmat(squeeze(nanmean(extended_field(:,2,:),1)),[radar.na 1]);
   extended_range(1,:)=0;
   
   for kk=1:radar.ne
      tmp_field=permute(extended_field(:,:,:),[2 1 3]);
      echo_top_3d(:,:,kk,:)=permute(interp1(squeeze(extended_range(:,1)),tmp_field,squeeze(radar.Rh(:,kk)),'linear'),[2 1 3]);
   end
   
%Interpolate echo_base


   %Artificially extend data to radar location
   extended_field=NaN(size(echo_depth_3d,1),size(echo_depth_3d,2)+1,size(echo_depth_3d,4));
   extended_range=NaN(size(radar.Rh,1)+1,size(radar.Rh,2));
   extended_field(:,2:end,:)=squeeze(echo_depth_3d(:,:,1,:));
   extended_range(2:end,:)=radar.Rh;
   extended_field(:,1,:)=repmat(squeeze(nanmean(extended_field(:,2,:),1)),[radar.na 1]);
   extended_range(1,:)=0;
   
   for kk=1:radar.ne
      tmp_field=permute(extended_field(:,:,:),[2 1 3]);
      echo_depth_3d(:,:,kk,:)=permute(interp1(squeeze(extended_range(:,1)),tmp_field,squeeze(radar.Rh(:,kk)),'linear'),[2 1 3]);
   end
   
 %Interpolate max_dbz


   %Artificially extend data to radar location
   extended_field=NaN(size(max_dbz_3d,1),size(max_dbz_3d,2)+1,size(max_dbz_3d,4));
   extended_range=NaN(size(radar.Rh,1)+1,size(radar.Rh,2));
   extended_field(:,2:end,:)=squeeze(max_dbz_3d(:,:,1,:));
   extended_range(2:end,:)=radar.Rh;
   extended_field(:,1,:)=repmat(squeeze(nanmean(extended_field(:,2,:),1)),[radar.na 1]);
   extended_range(1,:)=0;
   
   for kk=1:radar.ne
      tmp_field=permute(extended_field(:,:,:),[2 1 3]);
      max_dbz_3d(:,:,kk,:)=permute(interp1(squeeze(extended_range(:,1)),tmp_field,squeeze(radar.Rh(:,kk)),'linear'),[2 1 3]);
   end  
   
      
 %Interpolate max_dbz_level


   %Artificially extend data to radar location
   extended_field=NaN(size(max_dbz_level_3d,1),size(max_dbz_level_3d,2)+1,size(max_dbz_level_3d,4));
   extended_range=NaN(size(radar.Rh,1)+1,size(radar.Rh,2));
   extended_field(:,2:end,:)=squeeze(max_dbz_level_3d(:,:,1,:));
   extended_range(2:end,:)=radar.Rh;
   extended_field(:,1,:)=repmat(squeeze(nanmean(extended_field(:,2,:),1)),[radar.na 1]);
   extended_range(1,:)=0;
   
   for kk=1:radar.ne
      tmp_field=permute(extended_field(:,:,:),[2 1 3]);
      max_dbz_level_3d(:,:,kk,:)=permute(interp1(squeeze(extended_range(:,1)),tmp_field,squeeze(radar.Rh(:,kk)),'linear'),[2 1 3]);
   end  
   
   %Interpolate total_max_ref


   %Artificially extend data to radar location
   extended_field=NaN(size(total_max_dbz_3d,1),size(total_max_dbz_3d,2)+1,size(total_max_dbz_3d,4));
   extended_range=NaN(size(radar.Rh,1)+1,size(radar.Rh,2));
   extended_field(:,2:end,:)=squeeze(total_max_dbz_3d(:,:,1,:));
   extended_range(2:end,:)=radar.Rh;
   extended_field(:,1,:)=repmat(squeeze(nanmean(extended_field(:,2,:),1)),[radar.na 1]);
   extended_range(1,:)=0;
   
   for kk=1:radar.ne
      tmp_field=permute(extended_field(:,:,:),[2 1 3]);
      total_max_dbz_3d(:,:,kk,:)=permute(interp1(squeeze(extended_range(:,1)),tmp_field,squeeze(radar.Rh(:,kk)),'linear'),[2 1 3]);
   end  
   
      
 %Interpolate total_max_ref_level


   %Artificially extend data to radar location
   extended_field=NaN(size(total_max_dbz_level_3d,1),size(total_max_dbz_level_3d,2)+1,size(total_max_dbz_level_3d,4));
   extended_range=NaN(size(radar.Rh,1)+1,size(radar.Rh,2));
   extended_field(:,2:end,:)=squeeze(total_max_dbz_level_3d(:,:,1,:));
   extended_range(2:end,:)=radar.Rh;
   extended_field(:,1,:)=repmat(squeeze(nanmean(extended_field(:,2,:),1)),[radar.na 1]);
   extended_range(1,:)=0;
   
   for kk=1:radar.ne
      tmp_field=permute(extended_field(:,:,:),[2 1 3]);
      total_max_dbz_level_3d(:,:,kk,:)=permute(interp1(squeeze(extended_range(:,1)),tmp_field,squeeze(radar.Rh(:,kk)),'linear'),[2 1 3]);
   end  
   
   %Interpolate echo_top_size


   %Artificially extend data to radar location
   extended_field=NaN(size(echo_depth_size_3d,1),size(echo_depth_size_3d,2)+1,size(echo_depth_size_3d,4));
   extended_range=NaN(size(radar.Rh,1)+1,size(radar.Rh,2));
   extended_field(:,2:end,:)=squeeze(echo_depth_size_3d(:,:,1,:));
   extended_range(2:end,:)=radar.Rh;
   extended_field(:,1,:)=repmat(squeeze(nanmean(extended_field(:,2,:),1)),[radar.na 1]);
   extended_range(1,:)=0;
   
   for kk=1:radar.ne
      tmp_field=permute(extended_field(:,:,:),[2 1 3]);
      echo_depth_size_3d(:,:,kk,:)=permute(interp1(squeeze(extended_range(:,1)),tmp_field,squeeze(radar.Rh(:,kk)),'linear'),[2 1 3]);
   end  
   
   
   
 time=toc;
 
 display(['Echo top and depth data is computed and interpolated in ' num2str(time) ' seconds'])


