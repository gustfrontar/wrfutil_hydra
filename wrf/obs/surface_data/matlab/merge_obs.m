clear all

close all

inidate='20091117000000';
enddate='20091119000000';

latmin=-36.0;
latmax=-27.1;
lonmin=360-66.1;
lonmax=360-55.1;

endian='b';

obsfreq=3600; %Frequencia de las observaciones en segundos

OBSDATA1='/home/jruiz/datos/DATA/OBS/PREPBUFRSA/';

OBSDATA2='/home/jruiz/datos/DATA/OBS/SURFACESA/';

OBSMERGE='/home/jruiz/datos/DATA/OBS/PREPBUFR_AND_SURFACE/';

cdate=datenum(inidate,'yyyymmddHHMMSS');

edate=datenum(enddate,'yyyymmddHHMMSS');

while( cdate <= edate )

  nobs=0;
  
  display( datestr(cdate,'yyyymmddHHMMSS')  )

  file1=[ OBSDATA1 '/obs_' datestr(cdate,'yyyymmddHHMMSS') '.dat' ];

  file2=[ OBSDATA2 '/obs_' datestr(cdate,'yyyymmddHHMMSS') '.dat' ];

  filemerge=[ OBSMERGE '/obs_' datestr(cdate,'yyyymmddHHMMSS') '.dat' ];

  [obs1]=read_obs(file1,endian);

  [obs2]=read_obs(file2,endian);

  [obsmerge]=[obs1;obs2];

  location_index=( obsmerge(:,2) >= lonmin & obsmerge(:,2) <= lonmax & obsmerge(:,3) >= latmin & obsmerge(:,3) <= latmax );

  obsmerge=obsmerge( location_index , : );

  write_obs(filemerge,obsmerge,endian);

  cdate=cdate+obsfreq/86400;

  cdate=round(cdate*86400)/86400;

  nobs=size(obsmerge,1);
  
  display([num2str(nobs) ' observations written'])

end
