function [ gridrad ]=radar_superobbing( radar , ref , dv , cart , gridrad )   

%grid_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
%grid_dv =zeros(cart.nlon,cart.nlat,cart.nlev);
%grid_count_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
%grid_count_dv =zeros(cart.nlon,cart.nlat,cart.nlev);

%grid_az_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
%grid_el_ref=grid_az_ref;
%grid_ra_ref=grid_az_ref;

%grid_az_dv=grid_az_ref;
%grid_el_dv=grid_az_ref;
%grid_ra_dv=grid_az_ref;

%Convert from DBZ to power.
ref( ~isnan(ref) )=10.0d0 .^ ( ref( ~isnan(ref) )/10.0d0 ) ;


%Convert input averages to total sums.   
    gridrad.grid_ref = gridrad.grid_ref .* gridrad.grid_count_ref;
    gridrad.grid_az_ref = gridrad.grid_az_ref .* gridrad.grid_count_ref;
    gridrad.grid_el_ref = gridrad.grid_el_ref .* gridrad.grid_count_ref;
    gridrad.grid_ra_ref = gridrad.grid_ra_ref .* gridrad.grid_count_ref;
    gridrad.grid_lat_ref = gridrad.grid_lat_ref .* gridrad.grid_count_ref;
    gridrad.grid_lon_ref = gridrad.grid_lon_ref .* gridrad.grid_count_ref;
    gridrad.grid_z_ref = gridrad.grid_z_ref .* gridrad.grid_count_ref;

    gridrad.grid_dv = gridrad.grid_dv .* gridrad.grid_count_dv;
    gridrad.grid_az_dv = gridrad.grid_az_dv .* gridrad.grid_count_dv;
    gridrad.grid_el_dv = gridrad.grid_el_dv .* gridrad.grid_count_dv;
    gridrad.grid_ra_dv = gridrad.grid_ra_dv .* gridrad.grid_count_dv;
    gridrad.grid_lat_dv = gridrad.grid_lat_dv .* gridrad.grid_count_dv;
    gridrad.grid_lon_dv = gridrad.grid_lon_dv .* gridrad.grid_count_dv;
    gridrad.grid_z_dv = gridrad.grid_z_dv .* gridrad.grid_count_dv;


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

            gridrad.grid_ref(i,j,k)=gridrad.grid_ref(i,j,k)+ref(ia,ir,ie) ;
            gridrad.grid_count_ref(i,j,k)=gridrad.grid_count_ref(i,j,k)+1.0d0 ;
            %Azimuth should be carefully averaged, before adding a new value check that the distance
            %between the new value and the average is less than 180
            if( gridrad.grid_count_ref(i,j,k) == 1 )  
               gridrad.grid_az_ref(i,j,k)=gridrad.grid_az_ref(i,j,k)+radar.azimuth(ia);
            else
               tmpvar=gridrad.grid_az_ref(i,j,k)/(gridrad.grid_count_ref(i,j,k)-1) ; %Compute the current average.
               diff=tmpvar - radar.azimuth(ia) ; 
               if( diff > 180 )
                 gridrad.grid_az_ref(i,j,k)=gridrad.grid_az_ref(i,j,k)+radar.azimuth(ia)+360;
               elseif( diff < -180 )
                 gridrad.grid_az_ref(i,j,k)=gridrad.grid_az_ref(i,j,k)+radar.azimuth(ia)-360;
               else
                 gridrad.grid_az_ref(i,j,k)=gridrad.grid_az_ref(i,j,k)+radar.azimuth(ia);
               end
            end
            gridrad.grid_el_ref(i,j,k)=gridrad.grid_el_ref(i,j,k)+radar.elevation(ie);
            gridrad.grid_ra_ref(i,j,k)=gridrad.grid_ra_ref(i,j,k)+radar.range(ir);
            gridrad.grid_lat_ref(i,j,k)=gridrad.grid_lat_ref(i,j,k)+radar.latitude(ia,ir,ie);
            gridrad.grid_lon_ref(i,j,k)=gridrad.grid_lon_ref(i,j,k)+radar.longitude(ia,ir,ie);
            gridrad.grid_z_ref(i,j,k)=gridrad.grid_z_ref(i,j,k)+radar.Z(ia,ir,ie);
         end
         if ( ~isnan(dv(ia,ir,ie)) )
            gridrad.grid_dv(i,j,k)=gridrad.grid_dv(i,j,k)+dv(ia,ir,ie);
            gridrad.grid_count_dv(i,j,k)=gridrad.grid_count_dv(i,j,k)+1.0d0 ;
            if( gridrad.grid_count_dv(i,j,k) == 1 )
               gridrad.grid_az_dv(i,j,k)=gridrad.grid_az_dv(i,j,k)+radar.azimuth(ia);
            else
               tmpvar=gridrad.grid_az_dv(i,j,k)/(gridrad.grid_count_dv(i,j,k)-1) ; %Compute the current average.
               diff=tmpvar - radar.azimuth(ia) ;
               if( diff > 180 )
                 gridrad.grid_az_dv(i,j,k)=gridrad.grid_az_dv(i,j,k)+radar.azimuth(ia)+360;
               elseif( diff < -180 )
                 gridrad.grid_az_dv(i,j,k)=gridrad.grid_az_dv(i,j,k)+radar.azimuth(ia)-360;
               else
                 gridrad.grid_az_dv(i,j,k)=gridrad.grid_az_dv(i,j,k)+radar.azimuth(ia);
               end
            end
            gridrad.grid_el_dv(i,j,k)=gridrad.grid_el_dv(i,j,k)+radar.elevation(ie) ;
            gridrad.grid_ra_dv(i,j,k)=gridrad.grid_ra_dv(i,j,k)+radar.range(ir) ;
            gridrad.grid_lat_dv(i,j,k)=gridrad.grid_lat_dv(i,j,k)+radar.latitude(ia,ir,ie);
            gridrad.grid_lon_dv(i,j,k)=gridrad.grid_lon_dv(i,j,k)+radar.longitude(ia,ir,ie);
            gridrad.grid_z_dv(i,j,k)=gridrad.grid_z_dv(i,j,k)+radar.Z(ia,ir,ie);

         end

        end
     
      end

     end
    end

    
%Compute the average of the observations within each grid box.    
    gridrad.grid_ref( gridrad.grid_count_ref > 0 ) = gridrad.grid_ref( gridrad.grid_count_ref > 0 ) ./ gridrad.grid_count_ref( gridrad.grid_count_ref > 0 );

    gridrad.grid_az_ref( gridrad.grid_count_ref > 0 ) = gridrad.grid_az_ref( gridrad.grid_count_ref > 0 ) ./ gridrad.grid_count_ref( gridrad.grid_count_ref > 0 );
    gridrad.grid_el_ref( gridrad.grid_count_ref > 0 ) = gridrad.grid_el_ref( gridrad.grid_count_ref > 0 ) ./ gridrad.grid_count_ref( gridrad.grid_count_ref > 0 );
    gridrad.grid_ra_ref( gridrad.grid_count_ref > 0 ) = gridrad.grid_ra_ref( gridrad.grid_count_ref > 0 ) ./ gridrad.grid_count_ref( gridrad.grid_count_ref > 0 );
    gridrad.grid_lat_ref( gridrad.grid_count_ref > 0 ) = gridrad.grid_lat_ref( gridrad.grid_count_ref > 0 ) ./ gridrad.grid_count_ref( gridrad.grid_count_ref > 0 );
    gridrad.grid_lon_ref( gridrad.grid_count_ref > 0 ) = gridrad.grid_lon_ref( gridrad.grid_count_ref > 0 ) ./ gridrad.grid_count_ref( gridrad.grid_count_ref > 0 );
    gridrad.grid_z_ref( gridrad.grid_count_ref > 0 ) = gridrad.grid_z_ref( gridrad.grid_count_ref > 0 ) ./ gridrad.grid_count_ref( gridrad.grid_count_ref > 0 );

    gridrad.grid_dv( gridrad.grid_count_dv > 0 ) = gridrad.grid_dv( gridrad.grid_count_dv > 0 ) ./ gridrad.grid_count_dv( gridrad.grid_count_dv > 0 );
    gridrad.grid_az_dv( gridrad.grid_count_dv > 0 ) = gridrad.grid_az_dv( gridrad.grid_count_dv > 0 ) ./ gridrad.grid_count_dv( gridrad.grid_count_dv > 0 );
    gridrad.grid_el_dv( gridrad.grid_count_dv > 0 ) = gridrad.grid_el_dv( gridrad.grid_count_dv > 0 ) ./ gridrad.grid_count_dv( gridrad.grid_count_dv > 0 );
    gridrad.grid_ra_dv( gridrad.grid_count_dv > 0 ) = gridrad.grid_ra_dv( gridrad.grid_count_dv > 0 ) ./ gridrad.grid_count_dv( gridrad.grid_count_dv > 0 );
    gridrad.grid_lat_dv( gridrad.grid_count_dv > 0 ) = gridrad.grid_lat_dv( gridrad.grid_count_dv > 0 ) ./ gridrad.grid_count_dv( gridrad.grid_count_dv > 0 );
    gridrad.grid_lon_dv( gridrad.grid_count_dv > 0 ) = gridrad.grid_lon_dv( gridrad.grid_count_dv > 0 ) ./ gridrad.grid_count_dv( gridrad.grid_count_dv > 0 );
    gridrad.grid_z_dv( gridrad.grid_count_dv > 0 ) = gridrad.grid_z_dv( gridrad.grid_count_dv > 0 ) ./ gridrad.grid_count_dv( gridrad.grid_count_dv > 0 );

    %Check that averaged azimuth falls within the 0-360 range.
    for i=1:size(gridrad.grid_ref,1)
      for j=1:size(gridrad.grid_ref,2)
        for k=1:size(gridrad.grid_ref,3)
            if( gridrad.grid_az_ref(i,j,k) > 360 )
              gridrad.grid_az_ref(i,j,k)=gridrad.grid_az_ref(i,j,k)-360;
            end
            if( gridrad.grid_az_ref(i,j,k) < 0 )
              gridrad.grid_az_ref(i,j,k)=gridrad.grid_az_ref(i,j,k)+360;
            end

            if( gridrad.grid_az_dv(i,j,k) > 360 )
              gridrad.grid_az_dv(i,j,k)=gridrad.grid_az_dv(i,j,k)-360;
            end
            if( gridrad.grid_az_dv(i,j,k) < 0 )
              gridrad.grid_az_dv(i,j,k)=gridrad.grid_az_dv(i,j,k)+360;
            end
        end
      end
    end




end
