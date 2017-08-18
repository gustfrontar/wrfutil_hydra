clear all
close all

inidate='20091117000000';
enddate='20091119000000';

obsfreq=3600; %Frequencia de las observaciones en segundos

datapath='/home/jruiz/datos/DATA/OBS/SURFACE_DATA_RAW/';

outpath='/home/jruiz/datos/DATA/OBS/SURFACESA/';

undef=-9999;

endian='b';

US_OBS=82819;US_OBS_ERR=2.0;
VS_OBS=82820;VS_OBS_ERR=2.0;
TS_OBS=83073;TS_OBS_ERR=2.0;
QS_OBS=83330;QS_OBS_ERR=2.0e-3;
HR_OBS=83331;HR_OBS_ERR=10.0;
PS_OBS=14593;PS_OBS_ERR=200.0;

latmin=-36.0;
latmax=-27.1;
lonmin=-66.1;
lonmax=-55.1;

SURFACE_DATA_TYPE=1;

[station]=get_station_list_fun([datapath '/station_list']);

files = dir( [ datapath '/*' inidate(1:4) '*' ] );

data=[];

%Leo todos los archivos

for file = files'

    inputdata = load( [ datapath file.name ] )';

    inputdata( inputdata == undef ) = NaN;

    station_code=str2num(file.name(1:5));

    inputdate=datenum( [inputdata(1,:);inputdata(2,:);inputdata(3,:);inputdata(4,:);0*inputdata(4,:);0*inputdata(4,:)]');
    
    %Separamos las diferentes variables.

    temp=inputdata(5,:)/10.0;
    dewt=inputdata(6,:)/10.0;
    pres=inputdata(7,:)/10.0;
    wdir=inputdata(8,:);
    wspd=inputdata(9,:);

    %Transformamos las variables a T, Hr, U10 y V10;
    ep=get_es(dewt+273.16);
    esp=get_es(temp+273.16);
    rh=100*ep./esp;

    [ uwnd vwnd ]=VelDirToUV(wdir,wspd);

    clear tmpdata
   
    tmpdata(1,:)=inputdate;
    
    station_index=find( station.code == station_code );


    if ( ~isempty( station_index ) )
 
      tmpdata(1,:)=round(inputdate*86400)/86400;

      tmpdata(2,:)=station.lon(station_index);
      tmpdata(3,:)=station.lat(station_index);
      tmpdata(4,:)=station.z(station_index);

      tmpdata(5,:)=temp;
      tmpdata(6,:)=rh;
      tmpdata(7,:)=pres;
      tmpdata(8,:)=uwnd;
      tmpdata(9,:)=vwnd;

      data= [ data  tmpdata ] ;

    else

      display(['Station code not found :' num2str(station_code) ])

    end

end


cdate=datenum(inidate,'yyyymmddHHMMSS');
edate=datenum(enddate,'yyyymmddHHMMSS');


while( cdate <= edate )

  nobs=0;

  %Busco los elementos de datoa que corresponden a este tiempo.
  cdata = data( : , ( data(1,:) == cdate) & data(3,:) <= latmax & data(3,:) >= latmin & data(2,:) >= lonmin & data(2,:) <= lonmax );

  outfile=[ outpath '/obs_' datestr(cdate,'yyyymmddHHMMSS') '.dat'  ];

  display( datestr(cdate,'yyyymmddHHMMSS')  )

  nfile=fopen(outfile,'w',endian);

  for it=1:size(cdata,2)
   if( ~isnan( cdata(4,it) )  )
     %Escribo las obs de temperatura
     ivar=5;
     if( ~isnan( cdata(ivar,it) ) )
      wk(1)=TS_OBS            ; %ID
      wk(6)=TS_OBS_ERR        ; %Error
      wk(2)=cdata(2,it)       ; %Lon
      wk(3)=cdata(3,it)       ; %Lat
      wk(4)=cdata(4,it)       ; %Station height
      wk(5)=cdata(ivar,it)    ; %Observation
      wk(7)=SURFACE_DATA_TYPE ;
      fwrite( nfile , 7*4, 'int32' );
      fwrite( nfile , wk , 'float32' );
      fwrite( nfile , 7*4, 'int32' );
      nobs = nobs + 1;
     end
     %Escribo las obs de rh
     ivar=6;
     if( ~isnan( cdata(ivar,it) ) )
      wk(1)=HR_OBS            ; %ID
      wk(6)=HR_OBS_ERR        ; %Error
      wk(2)=cdata(2,it)       ; %Lon
      wk(3)=cdata(3,it)       ; %Lat
      wk(4)=cdata(4,it)       ; %Station height
      wk(5)=cdata(ivar,it)    ; %Observation
      wk(7)=SURFACE_DATA_TYPE ;
      fwrite( nfile , 7*4, 'int32' );
      fwrite( nfile , wk , 'float32' );
      fwrite( nfile , 7*4, 'int32' );
      nobs = nobs + 1;
     end
     %Escribo las obs de U
     ivar=8;
     if( ~isnan( cdata(ivar,it) ) )
      wk(1)=US_OBS            ; %ID
      wk(6)=US_OBS_ERR        ; %Error
      wk(2)=cdata(2,it)       ; %Lon
      wk(3)=cdata(3,it)       ; %Lat
      wk(4)=cdata(4,it)       ; %Station height
      wk(5)=cdata(ivar,it)    ; %Observation
      wk(7)=SURFACE_DATA_TYPE ;
      fwrite( nfile , 7*4, 'int32' );
      fwrite( nfile , wk , 'float32' );
      fwrite( nfile , 7*4, 'int32' );
      nobs = nobs + 1;
     end
     %Escribo las obs de V
     ivar=9;
     if( ~isnan( cdata(ivar,it) ) )
      wk(1)=VS_OBS            ; %ID
      wk(6)=VS_OBS_ERR        ; %Error
      wk(2)=cdata(2,it)       ; %Lon
      wk(3)=cdata(3,it)       ; %Lat
      wk(4)=cdata(4,it)    ; %Station height
      wk(5)=cdata(ivar,it)    ; %Observation
      wk(7)=SURFACE_DATA_TYPE ;
      fwrite( nfile , 7*4, 'int32' );
      fwrite( nfile , wk , 'float32' );
      fwrite( nfile , 7*4, 'int32' );
      nobs = nobs + 1;
     end
     %Escribo las obs de PS
     ivar=7;
     if( ~isnan( cdata(ivar,it) ) & ~isnan( cdata(4,it) ) ) 
      wk(1)=PS_OBS            ; %ID
      wk(6)=PS_OBS_ERR        ; %Error
      wk(2)=cdata(2,it)       ; %Lon
      wk(3)=cdata(3,it)       ; %Lat
      wk(4)=cdata(4,it)       ; %Station height
      wk(5)=cdata(ivar,it)    ; %Observation
      wk(7)=SURFACE_DATA_TYPE ;
      fwrite( nfile , 7*4, 'int32' );
      fwrite( nfile , wk , 'float32' );
      fwrite( nfile , 7*4, 'int32' );
      nobs = nobs + 1;
     end

   end

  end

  fclose(nfile);

  display(['A TOTAL NUMBER OF ' num2str(nobs) ' HAS BEEN WRITTEN TO THE OBSERVATION FILE'])

  cdate=cdate+obsfreq/86400;

  cdate=round(cdate*86400)/86400;

end







