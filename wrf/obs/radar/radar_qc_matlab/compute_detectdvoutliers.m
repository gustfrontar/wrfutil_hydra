function [ref2]=compute_detectdvoutliers(radar,dv,dv_thresh,nx,ny,nz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% data, a 3D array (usually azimuth,r,elev) containing a radar field
% (reflectivity or wind). In the case of wind take into account that
% missing values are filled with 0. That migth not be the best choice for
% computing wind texture.
% nx,ny,nz are the sizes of the local box over which texture will be
% computed. The box can be 1D, 2D or 3D. Size of the box in each direction
% is computed as: 2nx+1, 2ny+1 and 2nz+1 respectively.

tic
% Juan Ruiz 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic


ref2=ref;

for ie=1:radar.ne
  for ia=1:radar.na
     last_valid_val=NaN;
     for ir=2:radar.nr
         if( ~isnan( dv(ia,ir-1,ie) ) )
            last_valid_val = dv(ia,ir-1,ie) ;
         end
         %First step, detect sharp transitions between DV in adyacent
         %poinst along a ray.
         if( abs( last_valid_val - dv(ia,ir,ie) ) > dv_thresh & ~isnan( last_valid_val) ) 
             ismissing=true;
             nmissing=1;
             
             %If there is a jump, which pixel is more coherent with the
             %sourroundings?
             
             
             
             
         end
         
     end
       


  end
end


time=toc;

display(['Detect Missing computed in ' num2str(time) ' seconds']);






