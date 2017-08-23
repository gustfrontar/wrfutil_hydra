

function []=write_radar_data_seq(radar,ref,wind,attenuation,qcflag,filename,endian)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function writes a radar volume in an unformat binary.
% The output format is based on sequential fortran binary IO to speed up
% LETKF IO.
% This version takes into account only reflectivity and wind data. Other
% radars may include polarimetric variables.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
rec_length=4;

ref(isnan(ref))=radar.undef;
wind(isnan(wind))=radar.undef;
attenuation(isnan(attenuation))=radar.undef;
qcflag(isnan(qcflag))=radar.undef;

if( isempty(ref) )
    ref=radar.undef*ones(radar.na,radar.nr,radar.ne);
end
if( isempty(wind) )
    wind=radar.undef*ones(radar.na,radar.nr,radar.ne);
end
if( isempty(qcflag) )
    qcflag=radar.undef*ones(radar.na,radar.nr,radar.ne);
end
if( isempty(attenuation) )
    attenuation=radar.undef*ones(radar.na,radar.nr,radar.ne);
end

fid=fopen(filename,'w',endian);

%Write radar characteristics header.
radar.beam_wid_r=100.0;  %Radial resolution in meters.
radar.lambda=4;
radar.nvars=2;  %Currently wind and reflectivity qcflag and attenuation

fwrite(fid,int32(6*rec_length),'int32');%Fortran sequential header
fwrite(fid,radar.year,'single');
fwrite(fid,radar.month,'single');
fwrite(fid,radar.day,'single');
fwrite(fid,radar.hour,'single');
fwrite(fid,radar.minute,'single');
fwrite(fid,radar.second,'single');  
fwrite(fid,int32(6*rec_length),'int32');%Fortran sequential header

fwrite(fid,int32(8*rec_length),'int32');%Fortran sequential header
fwrite(fid,radar.lon,'single');
fwrite(fid,radar.lat,'single');
fwrite(fid,radar.altitude,'single');
fwrite(fid,radar.beam_wid_h,'single');
fwrite(fid,radar.beam_wid_v,'single');
fwrite(fid,radar.beam_wid_r,'single');  
fwrite(fid,radar.lambda,'single');     %Wave length in cm
fwrite(fid,radar.undef,'single');
fwrite(fid,int32(8*rec_length),'int32');%Fortran sequential header

%Write size data header
fwrite(fid,int32(4*rec_length),'int32');%Fortran sequential header
fwrite(fid,int32(radar.na),'int32');
fwrite(fid,int32(radar.nr),'int32');
fwrite(fid,int32(radar.ne),'int32');
fwrite(fid,int32(radar.nvars),'int32');
fwrite(fid,int32(4*rec_length),'int32');%Fortran sequential header

fwrite(fid,int32(radar.na*rec_length),'int32');%Fortran sequential header
fwrite(fid,radar.azimuth,'single');
fwrite(fid,int32(radar.na*rec_length),'int32');%Fortran sequential header

fwrite(fid,int32(radar.nr*rec_length),'int32');%Fortran sequential header
fwrite(fid,radar.range,'single');
fwrite(fid,int32(radar.nr*rec_length),'int32');%Fortran sequential header

fwrite(fid,int32(radar.ne*rec_length),'int32');%Fortran sequential header
fwrite(fid,radar.elevation,'single');
fwrite(fid,int32(radar.ne*rec_length),'int32');%Fortran sequential header

fwrite(fid,int32(1*rec_length),'int32');%Fortran sequential header
fwrite(fid,radar.attenuation_factor,'single');
fwrite(fid,int32(1*rec_length),'int32');%Fortran sequential header

%Write main block of reflectivity data

for ii=1:radar.ne
   fwrite(fid,int32(radar.na*radar.nr*rec_length),'int32');%Fortran sequential header
   fwrite(fid,ref(:,:,ii),'single');
   fwrite(fid,int32(radar.na*radar.nr*rec_length),'int32');%Fortran sequential header
end

%Write main block of wind data

for ii=1:radar.ne
   fwrite(fid,int32(radar.na*radar.nr*rec_length),'int32');%Fortran sequential header
   fwrite(fid,wind(:,:,ii),'single');
   fwrite(fid,int32(radar.na*radar.nr*rec_length),'int32');%Fortran sequential header
end

for ii=1:radar.ne
   fwrite(fid,int32(radar.na*radar.nr*rec_length),'int32');%Fortran sequential header
   fwrite(fid,attenuation(:,:,ii),'single');
   fwrite(fid,int32(radar.na*radar.nr*rec_length),'int32');%Fortran sequential header
end

for ii=1:radar.ne
   fwrite(fid,int32(radar.na*radar.nr*rec_length),'int32');%Fortran sequential header
   fwrite(fid,qcflag(:,:,ii),'single');
   fwrite(fid,int32(radar.na*radar.nr*rec_length),'int32');%Fortran sequential header
end



fclose(fid);

time=toc;
display(['Data was written in ' num2str(time) 'seconds'])

end



