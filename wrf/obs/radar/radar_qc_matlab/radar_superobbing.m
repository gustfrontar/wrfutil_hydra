function [grid_ref , grid_count_ref , grid_dv , grid_count_dv , grid_az_ref , grid_el_ref , grid_ra_ref ]=radar_superobbing( radar , ref , dv , cart , fileout )   

grid_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
grid_dv =zeros(cart.nlon,cart.nlat,cart.nlev);
grid_count_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
grid_count_dv =zeros(cart.nlon,cart.nlat,cart.nlev);

grid_az_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
grid_el_ref=grid_az_ref;
grid_ra_ref=grid_az_ref;

grid_az_dv=grid_az_ref;
grid_el_dv=grid_az_ref;
grid_ra_dv=grid_az_ref;

%Convert from DBZ to power.
ref( ~isnan(ref) )=10.0d0 .^ ( ref( ~isnan(ref) )/10.0d0 ) ;

%AVERAGE DATA AND INCLUDE OBSERVATIONA ERROR.
%We will compute the i,j,k for each radar grid point and box average the data.


    for ia=1:radar.na
     for ir=1:radar.nr
      for ie=1:radar.ne


        %Get i,j,k very simple approahc since we are assuming a regular lat/lon/z grid.
        %[reali , realj , realk ]=lll2ijk(radar.longitude(ia,ir,ie),radar.latitude(ia,ir,ie),radar.Z(ia,ir,ie));
        reali=(radar.longitude(ia,ir,ie)-cart.lon(1,1))/cart.dlon + 1.0d0 ;
        realj=(radar.latitude(ia,ir,ie) -cart.lat(1,1))/cart.dlat + 1.0d0 ;  
        realk=(radar.Z(ia,ir,ie)        -cart.z(1,1,1))/cart.dz   + 1.0d0 ;
        
        
        %Skip data outside the model domain.
        if( reali >= 1 & reali <= cart.nlon & realj >= 1 & realj <= cart.nlat & realk >= 1 & realk <= cart.nlev )

         i=round(reali);
         j=round(realj);
         k=round(realk);

         if ( ~isnan(ref(ia,ir,ie)) ) 

            grid_ref(i,j,k)=grid_ref(i,j,k)+ref(ia,ir,ie) ;
            grid_count_ref(i,j,k)=grid_count_ref(i,j,k)+1.0d0 ;
            
            grid_az_ref(i,j,k)=grid_az_ref(i,j,k)+radar.azimuth(ia);
            grid_el_ref(i,j,k)=grid_el_ref(i,j,k)+radar.elevation(ie);
            grid_ra_ref(i,j,k)=grid_ra_ref(i,j,k)+radar.range(ir);
         end
         if ( ~isnan(dv(ia,ir,ie)) )
            grid_dv(i,j,k)=grid_dv(i,j,k)+dv(ia,ir,ie);
            grid_count_dv(i,j,k)=grid_count_dv(i,j,k)+1.0d0 ;
            grid_az_dv(i,j,k)=grid_az_dv(i,j,k)+radar.azimuth(ia) ;
            grid_el_dv(i,j,k)=grid_el_dv(i,j,k)+radar.elevation(ie) ;
            grid_ra_dv(i,j,k)=grid_ra_dv(i,j,k)+radar.range(ir) ;

         end

        end
     
      end

     end
    end

    
%Compute the average of the observations within each grid box.    
    grid_ref( grid_count_ref > 0 ) = grid_ref( grid_count_ref > 0 ) ./ grid_count_ref( grid_count_ref > 0 );
    grid_az_ref( grid_count_ref > 0 ) = grid_az_ref( grid_count_ref > 0 ) ./ grid_count_ref( grid_count_ref > 0 );
    grid_el_ref( grid_count_ref > 0 ) = grid_el_ref( grid_count_ref > 0 ) ./ grid_count_ref( grid_count_ref > 0 );
    grid_ra_ref( grid_count_ref > 0 ) = grid_ra_ref( grid_count_ref > 0 ) ./ grid_count_ref( grid_count_ref > 0 );

    grid_dv( grid_count_dv > 0 ) = grid_dv( grid_count_dv > 0 ) ./ grid_count_dv( grid_count_dv > 0 );
    grid_az_dv( grid_count_dv > 0 ) = grid_az_dv( grid_count_dv > 0 ) ./ grid_count_dv( grid_count_dv > 0 );
    grid_el_dv( grid_count_dv > 0 ) = grid_el_dv( grid_count_dv > 0 ) ./ grid_count_dv( grid_count_dv > 0 );
    grid_ra_dv( grid_count_dv > 0 ) = grid_ra_dv( grid_count_dv > 0 ) ./ grid_count_dv( grid_count_dv > 0 );

    


 %WRITE DATA IN LETKF FORMAT.
 %WRITE FILE HEADER.
 
 nfile= fopen(fileout,'w','b');
 
 
   %OPEN(UNIT=99,FILE='radarobs.dat',STATUS='unknown',FORM='unformatted')
   %Introduce a small header with the radar possition and two values that might be useful.
   fwrite(nfile,4, 'int32'); 
   fwrite( nfile , radar.lon , 'float32');
   fwrite(nfile,4, 'int32'); 
   fwrite(nfile,4, 'int32'); 
   fwrite( nfile , radar.lat , 'float32');
   fwrite(nfile,4, 'int32'); 
   fwrite(nfile,4, 'int32'); 
   fwrite( nfile , radar.altitude , 'float32' );
   fwrite(nfile,4, 'int32'); 
   
    nobs=0;
    for ii=1:cart.nlon
     for jj=1:cart.nlat
      for kk=1:cart.nlev
         if( grid_count_ref(ii,jj,kk) > 0 )
          
           wk(1)=radar.id_ref_obs         ;
           wk(6)=radar.error_ref          ;
           wk(2)=grid_az_ref(ii,jj,kk)    ;
           wk(3)=grid_el_ref(ii,jj,kk)    ;
           wk(4)=grid_ra_ref(ii,jj,kk)    ;
           wk(5)=grid_ref(ii,jj,kk)       ;
           wk(7)=radar.radar_type ;
           
           fwrite(nfile,7*4, 'int32'); 
           fwrite( nfile , wk ,'float32');
           fwrite(nfile,7*4, 'int32'); 
           nobs = nobs + 1;
          end
          if( grid_count_dv(ii,jj,kk) > 0 )
           wk(1)=radar.id_dv_obs           ;
           wk(6)=radar.error_dv            ;
           wk(2)=grid_az_dv(ii,jj,kk)      ;
           wk(3)=grid_el_dv(ii,jj,kk)      ;
           wk(4)=grid_ra_dv(ii,jj,kk)      ;
           wk(5)=grid_dv(ii,jj,kk)         ;
           wk(7)=radar.radar_type          ;
           
           fwrite(nfile,7*4, 'int32'); 
           fwrite( nfile , wk ,'float32')  ;
           fwrite(nfile,7*4, 'int32'); 
           
           nobs=nobs +1;
          end
      end
     end
    end



display(['A TOTAL NUMBER OF ' num2str(nobs) ' HAS BEEN WRITTEN TO THE OBSERVATION FILE'])


end
