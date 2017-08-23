function [wind , u_vad , v_vad , uo , vo , tmp_wind , vr_background , sample_coverage]=compute_wind_qc(radar,wind,nx,ny,nz,nx2,ny2,nz2,vad_threshold,nanomaly_threshold);
% This function performs the wind quality control. Currently the algorithm
% is done in 3 steps:
% 1) First compute vad from wind fiel at each range and elevation angle.
% 2) Temporarily remove the points that does not fit the VAD.
% 3) Compute the local mean with the points that does not fit the VAD
% removed.
% 4) Remove the poinst that are far from the local mean.

%VAD estimate provide a global constrain for the QC (i.e. is this point
%similar to the large scale circulation?) If not do not remove it yet (it
%might be a small scale feature).
%The local mean provides a local quality constrain based on the globa one.
%This algorithm should be furtherly tested to see how it works.


tic

deg2rad=pi/180;

vr_background=NaN(radar.na,radar.nr,radar.ne);
min_data_coverage=0.1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First obtain vad (and local elev and range vad) for 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[u_vad , v_vad , z_vad , u_error , v_error , uo , vo , sample_coverage , fit]=compute_vad( radar , wind );

if( size(vo,2) < radar.ne)
vo(:,end+1:radar.ne)=NaN;
uo(:,end+1:radar.ne)=NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Smooth uo and vo and remove isolated levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uo( sample_coverage < min_data_coverage )=NaN;
vo( sample_coverage < min_data_coverage )=NaN;

%count_wind=compute_boxcount(uo,1,1,0,-1e9);
%tmp_index=count_wind < 0.50; % Ranges and azimuths that are isolated.
%uo(tmp_index)=NaN;
%vo(tmp_index)=NaN;

%figure
%pcolor(uo);shading flat

for iter=1:2
for ie=1:radar.ne
   tmp=uo(:,ie);
   tmp=tmp(~isnan(tmp));
   uomean=mean(tmp);
   uostd =std(tmp);
   tmp_index=abs( uo(:,ie)-uomean ) > 3.0*uostd;
   uo(tmp_index,ie)=NaN;
   vo(tmp_index,ie)=NaN;
   
   tmp=vo(:,ie);
   tmp=tmp(~isnan(tmp));
   vomean=mean(tmp);
   vostd =std(tmp);
   tmp_index=abs(vo(:,ie) - vomean ) > 3.0*vostd;
   uo(tmp_index,ie)=NaN;
   vo(tmp_index,ie)=NaN;
end
end

%figure
%pcolor(uo);shading flat

%Compute mean
tmp_uo=compute_boxmeannc(uo,4,2,0);
tmp_vo=compute_boxmeannc(vo,4,2,0);

%Remove values that are different from the local mean
tmp_index=( abs(tmp_uo - uo) > 5 | abs(tmp_vo-vo) > 5 );
uo(tmp_index)=tmp_uo(tmp_index); %Replace the value with the local mean.
vo(tmp_index)=tmp_vo(tmp_index);

%Recompute mean
tmp_uo=compute_boxmeannc(uo,4,2,0);
tmp_vo=compute_boxmeannc(vo,4,2,0);

count_wind=compute_boxcount(uo,1,1,0,-1e9);
tmp_index=count_wind < 0.50; % Ranges and azimuths that are isolated.
uo(tmp_index)=NaN;
vo(tmp_index)=NaN;

% 
% %Replace missing values with the local mean if its not missing too. This
% %extends the ranges and elevation angles with uo and vo values and
% %maximizes the number of good data retained by the QC algorithm.
uo(isnan(uo))=tmp_uo(isnan(uo));
vo(isnan(vo))=tmp_vo(isnan(vo));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate a background field using the VAD information.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 elev=repmat(radar.elevation',[radar.nr 1]);
 cosaz=cos(radar.azimuth*deg2rad);
 sinaz=sin(radar.azimuth*deg2rad);
 cosel=cosd(elev);

for ii=1:radar.na
    vr_background(ii,:,:)=vo.*cosel*cosaz(ii)+uo.*cosel*sinaz(ii);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRST COMPUTE LOCAL MEAN REMOVING POINTS THAT DOES NOT FIT TO VAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp_index= abs(vr_background-wind ) > vad_threshold | isnan(vr_background) ;

tmp_wind=wind;
tmp_wind(tmp_index)=NaN;

%figure
%pcolor(tmp_wind(:,:,6));shading flat


%This mean function computes the local mean around a point 
%using the sourrounding points but not
%the point itself.
[localmean localstd]=compute_boxmean_V2(tmp_wind,nx,ny,nz,nx2,ny2,nz2);
%figure
%pcolor(localmean(:,:,6));shading flat

%tmp_index= abs ( wind - localmean ) > mean_threshold | isnan(localmean);
localstd ( localstd==0 )= NaN;
standard_dvc=abs(tmp_wind-localmean)./localstd;
tmp_index= standard_dvc > nanomaly_threshold | isnan(standard_dvc);

%Remove the points where the difference between the local mean and the
%value is greather than the indicated threshold.
wind( tmp_index ) = NaN;


time=toc;
display(['Time to complete the wind quality control: ' num2str(time) ' seconds'])
end


