

function [radar , ref , wind , qcflag , attenuation]=read_radar_data_seq(filename,endian)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function writes a radar volume in an unformat binary.
% The output format is based on sequential fortran binary IO to speed up
% LETKF IO.
% This version takes into account only reflectivity and wind data. Other
% radars may include polarimetric variables.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
undef=-327.68;

fid=fopen(filename,'r',endian);

%Write radar characteristics header.
radar.beam_wid_r=100.0;  %Radial resolution in meters.
radar.lambda=3;
radar.nvars=2;  %Currently wind and reflectivity.

%TIME HEADER
fread(fid,1,'int');%Fortran sequential header
radar.year=fread(fid,1,'float');
radar.month=fread(fid,1,'float');
radar.day=fread(fid,1,'float');
radar.hour=fread(fid,1,'float');
radar.minute=fread(fid,1,'float');
radar.second=fread(fid,1,'float');  
fread(fid,1,'int');%Fortran sequential header

%RADAR LOCATION AND CHARACTERISTICS HEADER
fread(fid,1,'int');%Fortran sequential header
radar.longitude=fread(fid,1,'float');
radar.latitude=fread(fid,1,'float');
radar.altitude=fread(fid,1,'float');
radar.beam_wid_h=fread(fid,1,'float');
radar.beam_wid_v=fread(fid,1,'float');
radar.beam_wid_r=fread(fid,1,'float');  
radar.lambda=fread(fid,1,'float');     %Wave length in cm
radar.undef=fread(fid,1,'float');
fread(fid,1,'int');%Fortran sequential header

%DATA SIZE HEADER
fread(fid,1,'int');%Fortran sequential header
radar.na=double(fread(fid,1,'int'));
radar.nr=double(fread(fid,1,'int'));
radar.ne=double(fread(fid,1,'int'));
radar.nvars=double(fread(fid,1,'int'));
fread(fid,1,'int');%Fortran sequential header

%AZIMUTH
fread(fid,1,'int');%Fortran sequential header
radar.azimuth=fread(fid,radar.na,'float');
fread(fid,1,'int');%Fortran sequential header

%RANGE
fread(fid,1,'int');%Fortran sequential header
radar.radius=fread(fid,radar.nr,'float');
fread(fid,1,'int');%Fortran sequential header

%ELEVATION
fread(fid,1,'int');%Fortran sequential header
radar.elev=fread(fid,radar.ne,'float');
fread(fid,1,'int');%Fortran sequential header

%GLOBAL ATTENUATION FACTOR
fread(fid,1,'int');%Fortran sequential header
radar.attenuation_factor=fread(fid,1,'float');
fread(fid,1,'int');%Fortran sequential header

%Prealocate data to be read.

ref=NaN(radar.na,radar.nr,radar.ne);
wind=NaN(radar.na,radar.nr,radar.ne);
qcflag=NaN(radar.na,radar.nr,radar.ne);
attenuation=NaN(radar.na,radar.nr,radar.ne);

%Read main block of reflectivity data

for ie=1:radar.ne
   fread(fid,1,'int');%Fortran sequential header
   ref(:,:,ie)=fread(fid,[radar.na radar.nr],'float');
   fread(fid,1,'int');%Fortran sequential header
end

%Read main block of wind data

for ie=1:radar.ne
   fread(fid,1,'int');%Fortran sequential header
   wind(:,:,ie)=fread(fid,[radar.na radar.nr],'float');
   fread(fid,1,'int');%Fortran sequential header
end

%Read main block of attenuation

% for ie=1:radar.ne
%    fread(fid,1,'int');%Fortran sequential header
%    attenuation(:,:,ie)=fread(fid,[radar.na radar.nr],'float');
%    fread(fid,1,'int');%Fortran sequential header
% end
% 
% %Read main block of qcflag
% 
% for ie=1:radar.ne
%    fread(fid,1,'int');%Fortran sequential header
%    qcflag(:,:,ie)=fread(fid,[radar.na radar.nr],'float');
%    fread(fid,1,'int');%Fortran sequential header
% end



fclose(fid);

time=toc;
display(['Data was read in ' num2str(time) 'seconds'])

%Replace undef with NaN;

wind(wind == radar.undef)=NaN;
ref(ref==radar.undef)=NaN;
qcflag(qcflag==radar.undef)=NaN;
attenuation(attenuation==radar.undef)=NaN;

end



