function [echo_top echo_base total_max_ref total_max_ref_level max_ref max_ref_level max_radar_heigth]=compute_echotop(radar,data,zthresh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% data, a 3D array (usually azimuth,r,elev) containing a radar field
% (reflectivity). 
% radar estructure contining geolocation information for data.
% zthresh the reflectivity thresholds that will be used for echo top
% computation.
% Output:
% echo_top, top of echoes between zmin and zmax and that pass certain
% checks.
% echo_base, base of echoes between zmin and zmax and that pass certain
% checks.
% max_ref maximum reflectivity in the column (between zmin and zmax) and
% corresponding to echos that pass the continuity checks.
% max_ref_level height corresponding to the max_ref detected.
% max_radar_height: at a given radius (corresonding to the lowest radar
% elevation) the maximum heigth of the volume scan.
% total_max_ref is the maximum reflectivity in the column without tacking
% into account the requirements to detect an echo (vertical extension
% height, etc).
% total_max_ref_level is the level where the total_max_ref is found.

% Al 2D outputs are given in the grid corresponding to the lowest
% elevation. This is because this elevation has the largest range.

%The interpolation approach used in this function migth produce spurius
%results close to the radar (associated with high elevation angles) in this
%part the resulting vertical resolution will be poor although in the
%original data is very high. This can be solved performing also a vertical
%interpolation to a regular height grid. tmp_radar_data and tmp_radar_z can
%be maped to a regular R-Z grid. 

% Juan Ruiz 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Some function parameters.
zmin=0; %Only echoes above this level will be considered in echo top computation.
zmax=20000; %Echoes above this heigth won't be considered for echo top computation.

continuity_threshold_base=1; %how many consecutive points over the threshold will be needed to initiate a cloud.
continuity_threshold_top=1;  %how many points below the reflectivity threshold do we need to be sure that the cloud is over.
mem_size=max([continuity_threshold_base continuity_threshold_top]);

%For example an echo is detected if we have at least mem_size consecutive
%gates with reflectivity over the selected threshold. An echo topo is
%detected if we have at least mem_size consecutive gates with reflectivity
%below the threshold (and we have detected an echo below that level);

Necho_top=3;    %Maximum number of echo tops that can be detected.
Necho_base=3;   %Maximum number of echo bases that can be detected.

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERP EACH LEVEL REFLECTIVITY TO THE FIRST LEVEL GRID.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_radar_data=NaN(size(data));
tmp_radar_z=NaN(size(data));

data(isnan(data))=radar.replacerefmissing;   %Replace missing values prior
                                             %to data interpolation.

warning off %Supress warning about NaN in interpolated fields.
%Loop over elevations
for ie=1:radar.ne
tmp_data=permute(squeeze(data(:,:,ie)),[2 1]);
tmp_z=permute(squeeze(radar.Z(:,:,ie)),[2 1]);

tmp_radar_data(:,:,ie)=permute(interp1(radar.Rh(:,ie),tmp_data,radar.Rh(:,1)),[2 1]);
tmp_radar_z(:,:,ie)=permute(interp1(radar.Rh(:,ie),tmp_z,radar.Rh(:,1)),[2 1]);

end
warning on

max_radar_heigth=max(tmp_radar_z,2);

%Echo top and echo base has to be within zmin and zmax.
tmp_radar_data( tmp_radar_z < zmin | tmp_radar_z > zmax )=0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  COMPUTE THE MAXIMUM HEIGHT OF ECHOES ACCORDING TO THE THRESHOLDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nt=length(zthresh);

%Allocate memory for output variables.
echo_top=NaN(radar.na,radar.nr,nt);
echo_base=NaN(radar.na,radar.nr,nt);

max_ref=NaN(radar.na,radar.nr,nt);
max_ref_level=NaN(radar.na,radar.nr,nt);
total_max_ref=zeros(radar.na,radar.nr,nt);
total_max_ref_level=zeros(radar.na,radar.nr,nt);

%We will use TMP_RADAR_DATA and TMP_RADAR_Z the reflectivity and height of
%each level interpolated to the lowest level grid.


%Loop over thresholds to compute echo top according to each threshold.

for it=1:nt
    
    over_zmin=false(radar.na,radar.nr);
    over_zmax=false(radar.na,radar.nr);
    %Loop over the elevations.
    
    %Allocate temporary arrays.
    tmp_echo_top=NaN(radar.na,radar.nr);
    tmp_echo_base=NaN(radar.na,radar.nr);
    tmp_max_ref=NaN(radar.na,radar.nr);
    tmp_max_ref_level=NaN(radar.na,radar.nr);
    base_detected=false(radar.na,radar.nr);
    top_detected=false(radar.na,radar.nr);
    previus_level_ref=zeros(radar.na,radar.nr,mem_size);
    previus_level_z  =zeros(radar.na,radar.nr,mem_size);
    tmp_total_max_ref=zeros(radar.na,radar.nr);
    tmp_total_max_ref_level=zeros(radar.na,radar.nr);
    
    for ie=1:radar.ne
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get the total_max_ref and total_max_ref_level
        tmp_data=tmp_radar_data(:,:,ie);
        tmp_height=tmp_radar_z(:,:,ie);
        tmp_index= tmp_data > tmp_total_max_ref;
        tmp_total_max_ref(tmp_index)=tmp_data(tmp_index);
        tmp_total_max_ref_level(tmp_index)=tmp_height(tmp_index);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        previus_level_ref(:,:,1:end-1)=previus_level_ref(:,:,2:end);
        previus_level_ref(:,:,end)=tmp_radar_data(:,:,ie);
        
        previus_level_z(:,:,1:end-1)=previus_level_z(:,:,2:end);
        previus_level_z(:,:,end)=tmp_radar_z(:,:,ie);
        
        
        %Can we initiate a cloud? we require at least
        %continuity_threshold_base consecutive values above the threshold
        %and over zmin.
        init_index= sum(previus_level_ref(:,:,end-continuity_threshold_base+1:end) > zthresh(it),3)== continuity_threshold_base & (~base_detected | top_detected);
        init_index2= sum(previus_level_ref(:,:,end-continuity_threshold_base+1:end) > zthresh(it),3)== continuity_threshold_base & ~base_detected;
        base_detected( init_index2 )= true;
        top_detected( init_index ) = false; %Redefine echo top.
        tmp_echo_top( init_index ) = NaN;
        
        %Set the level of cloud initiation as the echo_base;
         tmp_height=squeeze(previus_level_z(:,:,end-continuity_threshold_base+1));
        %Note that if echo top is redefined, then echo_base is also
        %redefined.
         tmp_echo_base( init_index )=tmp_height( init_index );
        

            
        %If we have a base then start the computation of max_reflectivity
        %within the cloud (not the total cloud reflectivity that does not
        %perform any continuity check).

         tmp_data=squeeze(max(previus_level_ref(:,:,end-continuity_threshold_base+1:end),[],3)); %Tmp_data is redefined;
         tmp_max_ref( init_index2 )=tmp_data( init_index2 );
         %Local loop to find the level of maximum reflectivity.
            for zz=1:mem_size
               tmp_data=squeeze(previus_level_ref(:,:,end-zz+1));
               tmp_height=squeeze(previus_level_z(:,:,end-zz+1));
               tmp_index3=tmp_data==tmp_max_ref & init_index2 ;
               tmp_max_ref_level(tmp_index3)=tmp_height(tmp_index3);
            end
        
            
         %Continue with the computation of max ref and max ref level within
         %the cloud until echo_top is reached.
         tmp_data=squeeze(previus_level_ref(:,:,end));
         tmp_height=squeeze(previus_level_z(:,:,end));
         tmp_index=( base_detected & tmp_data > tmp_max_ref & ~top_detected );
         tmp_max_ref( tmp_index ) = tmp_data( tmp_index );
         tmp_max_ref_level( tmp_index ) = tmp_height( tmp_index );
        
         %Detect echo top
         %We require at least 3 consecutive levels where the reflectivity is below the
         %threshold.

         tmp_height=squeeze(previus_level_z(:,:,end-continuity_threshold_top+1));
         top_index= base_detected & sum(previus_level_ref(:,:,end-continuity_threshold_top+1:end) > zthresh(it),3) == 0 & ~ top_detected ;
         top_detected(top_index)=true;
         tmp_echo_top( top_index )=tmp_height( top_index );

      %  end
        
        if ( ie == radar.ne )
        %What happens if we reach the last level and the condition for echo
        %top is not detected.
           tmp_index= base_detected & ~ top_detected ;
           %Consider the gates where a base has been detected but the top
           %cannot be detected.
           tmp_height=squeeze(previus_level_z(:,:,end));
           %In this case assign the highest z possible for this particular
           %gate.
           tmp_echo_top( tmp_index )=tmp_height( tmp_index );
        
        end
    end
    
    
    echo_top(:,:,it)=tmp_echo_top;
    echo_base(:,:,it)=tmp_echo_base;
    max_ref(:,:,it)=tmp_max_ref;
    total_max_ref(:,:,it)=tmp_total_max_ref;
    max_ref_level(:,:,it)=tmp_max_ref_level;
    total_max_ref_level(:,:,it)=tmp_total_max_ref_level;
    
end

   
time=toc;
    
display(['Echo top was computed in ' num2str(time) ' seconds'])    
    












