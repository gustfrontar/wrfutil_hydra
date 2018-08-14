function [ref2]=compute_detectmissing(radar,ref,ref_thresh,minref,nmissing_max)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% data, a 3D array (usually azimuth,r,elev) containing a radar field
% (reflectivity or wind). In the case of wind take into account that
% missing values are filled with 0. That migth not be the best choice for
% computing wind texture.
% nx,ny,nz are the sizes of the local box over which texture will be
% computed. The box can be 1D, 2D or 3D. Size of the box in each direction
% is computed as: 2nx+1, 2ny+1 and 2nz+1 respectively.

ref2=ref;

for ie=1:radar.ne
  for ia=1:radar.na
     ismissing=false;
     nmissing=1;
     for ir=2:radar.nr 
         if( abs( ref(ia,ir-1,ie) - ref(ia,ir,ie) ) > ref_thresh & ref(ia,ir,ie) <= minref ) 
             ismissing=true;
             nmissing=1;
         end
         if( ismissing & ref(ia,ir,ie) > minref )
             ismissing=false;
         end
         if( nmissing > nmissing_max )
             ismissing=false;
         end
         if( ismissing )
            ref2(ia,ir,ie)=NaN;
            nmissing=nmissing+1;
         end
     end

  end
end


time=toc;

display(['Detect Missing computed in ' num2str(time) ' seconds']);






