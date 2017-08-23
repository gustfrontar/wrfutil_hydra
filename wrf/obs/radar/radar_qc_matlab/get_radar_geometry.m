
clear all
close all

undef=-9999;

path='/home/jruiz/Dropbox/DATA/DATOS_RADAR/ANGUIL/20100111/120/';

files=dir([path '*.nc_dealiased']);

zthresh=5; %Threshold in dBZ for echo top computation.
snx=5;
sny=5;

minref=0;


error_ref=5; %Reflectivity error in DBZ.
error_dv=2;  %Doppler velocity error in m/s.

id_ref_obs=4001;
id_dv_obs=4002;

%Increments for cartesian grid (m).
dx=2000;
dz=2000;
maxz=25e3;  %Maximum height for the observations.

%QC
rhohvt=0.8; %Rhohv threshold for QC
echo_top_t=3000;%Echo top threshold (m).


radar_type=1;

radar.undef=undef;

ifile=3;


    
    filename=[ path files(ifile).name ];
    
    [radar , ref , dv , dvc , phidp , rhohv ]=read_radar_netcdf( filename );
    
    rhohv(rhohv < 0)=NaN;
    [rhohvsmooth]=compute_boxmean(rhohv,1,1,1);
    
    
    ref(ref < minref)=minref;
    
    %Georeference will be performed for each time due to possible changes
    %in observation strategy. 
    
    [ radar ]=georeference_radar_data( radar );
    
    radar.replacerefmissing=minref;
    
    radar.error_ref=error_ref;
    radar.error_dv =error_dv ;
    radar.id_ref_obs= id_ref_obs;
    radar.id_dv_obs = id_dv_obs;
    radar.radar_type=radar_type;
    
    %Topography blocking.
    %Interpolar una topografia del WRF o de alguna otra fuente...
   
    
    
 