function [ stations ]=station_list_fun( file )

%Esta funcion lee la lista de estaciones que asocia codigo de estacion con lat lon.

stations.n = 0; %Number of stations.

cont=true;

fid = fopen( file );

tline=' ';

while ischar( tline )

    tline = fgetl(fid) ; 

    if ( ischar( tline ) )

       stations.n = stations.n + 1 ;

       stations.code(stations.n)=str2num( tline(1:5) );

       tmplatint=str2num(tline(38:39));
       tmplatdec=str2num(tline(40:41));
       tmplath=tline(42);

       tmplonint=str2num(tline(44:46));
       tmplondec=str2num(tline(47:48));
       tmplonh=tline(49); 

       altura= tline(50:54) ;
       if ( ~isempty( str2num(altura) ) )
          stations.z(stations.n)=str2num( altura );
       else
          stations.z(stations.n)=NaN;
       end

       tmplat=tmplatint+tmplatdec/60;

       tmplon=tmplonint+tmplondec/60;

       if( strcmp( tmplath , 'S' ) )
         tmplat=-tmplat;
       end

       if( strcmp( tmplonh , 'W' ) )
         tmplon=-tmplon;
       end

       stations.lat(stations.n)=tmplat;
       stations.lon(stations.n)=tmplon;

    end

end


