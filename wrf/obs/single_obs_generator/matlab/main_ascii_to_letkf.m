clear all
close all

inidate='20091117180000';
enddate='20091117190000';

obsfreq=300; %Frequencia de las observaciones en segundos

%datapath='/home/jruiz/datos/DATA/OBS/SURFACE_DATA_RAW/';

outpath='/home/jruiz/share/DATA/OBS/SINGLEOBSEXP/';

undef=-9999;

endian='b';

US_OBS=82819;US_OBS_ERR=1.0;
VS_OBS=82820;VS_OBS_ERR=1.0;
TS_OBS=83073;TS_OBS_ERR=1.0;
QS_OBS=83330;QS_OBS_ERR=1.0e-3;
HR_OBS=83331;HR_OBS_ERR=10.0;
PS_OBS=14593;PS_OBS_ERR=200.0;

latmin=-36.0;
latmax=-27.1;
lonmin=360-66.1;
lonmax=360-55.1;

SURFACE_DATA_TYPE=1;

cdate=datenum(inidate,'yyyymmddHHMMSS');
edate=datenum(enddate,'yyyymmddHHMMSS');


while( cdate <= edate )

  outfile=[ outpath '/obs_' datestr(cdate,'yyyymmddHHMMSS') '.dat'  ];

  display( datestr(cdate,'yyyymmddHHMMSS')  )

  nfile=fopen(outfile,'w',endian);

  %Escribo una obs de temperatura.
  wk(1)=QS_OBS            ; %ID
  wk(6)=QS_OBS_ERR        ; %Error
  wk(2)=-60.3             ; %Lon
  wk(3)=-31.7             ; %Lat
  wk(4)=30                ; %Station height
  wk(5)=1.0e-3            ; %Observation
  wk(7)=-9                ; %Code for simulated observations 
  fwrite( nfile , 7*4, 'int32' );
  fwrite( nfile , wk , 'float32' );
  fwrite( nfile , 7*4, 'int32' );

  fclose(nfile);


  cdate=cdate+obsfreq/86400;

  cdate=round(cdate*86400)/86400;

end







