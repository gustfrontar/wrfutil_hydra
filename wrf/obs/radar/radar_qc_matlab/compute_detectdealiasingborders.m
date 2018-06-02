function [ dvc ] = compute_detectdealiasingborders( radar , dv , dvda , nx , ny , nz , nh_thresh )

%dv  doppler velocity without unfolding.
%dvda doppler velocity with unfolding.
%nx,ny,nz size of the neighborhood to be explored in the border between
%dealiased and undealiased pixels.
%nh_thresh threshold which will be used to detect strong wind jumps in this
%regions.

%dvc output corrected wind fild.


diff = dv-dvda;

dvc=dvda;
       
%nh_size=2;    %Size of the neighborhood.
%nh_thresh=10; %Umbral por encima del cual si la diferencia es mayor se descarta el pixel.
       
%%Loop over elevations
for ie=1:radar.ne
  %%Loop over azimuths
  for ia=1:radar.na
    %%Loop over ranges
    for ir=1:radar.nr
        
        if ( diff(ia,ir,ie) ~= 0 && ~isnan(diff(ia,ir,ie)) )
               
           %Look for neighbors.
           for iia=-nx:nx
              for iir=-ny:ny
                 for iie=-nz:nz
                    mya=iia+ia;
                    myr=iir+ir;
                    mye=iie+ie;
                    if( mya < 1 );mya=1;end
                    if( mya > radar.na);mya=radar.na;end
                    if( myr < 1 );myr=1;end
                    if( myr > radar.nr);myr=radar.nr;end
                    if( mye < 1 );mye=1;end
                    if( mye > radar.ne);mye=radar.ne;end
                     
                    %If a neighbor has not been dealiased and its doppler
                    %velocity difference with respect to the current pixel
                    %is larger than the threshold, then delete that pixel.
                    if( diff(mya,myr,mye) == 0 && abs( dvda(mya,myr,mye) - dvda(ia,ir,ie) ) > nh_thresh )
                        dvc(mya,myr,mye)=NaN;
                    end
                    
                    %If this pixel is in the border (it has a non dealiased
                    %neighbor) then it is also elminated.
                    if( abs(iia) <= 1 && abs(iir) <= 1 && abs(iie) <= 1 && diff(mya,myr,mye) == 0)
                        dvc(ia,ir,ie) = NaN;
                    end
                     
                 end
              end
           end
                

        end %End if 
     end  %End loop over ranges 
   end  %End loop over azimuths
end  %End loop over elevations


end %End function

