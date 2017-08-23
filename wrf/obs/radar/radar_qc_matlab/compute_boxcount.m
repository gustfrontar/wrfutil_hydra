function [count]=compute_boxcount(data,nx,ny,nz,threshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% data, a 3D array (usually azimuth,r,elev) containing a radar field
% (reflectivity or wind). In the case of wind take into account that
% missing values are filled with 0. That migth not be the best choice for
% computing wind texture.
% nx,ny,nz are the sizes of the local box over which median, mean and std
% will be computed. 
% Threshold is the reflectivity threshold. Data above this threshold will
% be counted. Tipically 5 dbz for reflectivity count (Lakshmaman et al
% 2006).
% This parameter is similar to what is known as echo size. It can also be
% applied to echo top instead of reflectivity. 
% To avoid dependence with box size the count relative to the box size is
% provided (percentage of points with variable values above the desired
% threshold, of course NaN values are never above any threshold).
% Juan Ruiz 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

[na nr ne]=size(data);

%data(isnan(data))=0;

nx=squeeze(nx);
ny=squeeze(ny);
nz=squeeze(nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATISTICS OVER THE SQUARE LOCAL BOX (CAN BE A 1D, 2D or 3D BOX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%To reduce memory size the computation will be done level by level.
box_size=(2*nx+1)*(2*ny+1)*(2*nz+1);

count=NaN(size(data));

for kk=1:ne
   tmp_field=NaN(na,nr,box_size);
   tmp_index=1;
   for j=-ny:ny
       for i=-nx:nx
         if(j <= 0)
           jmax=nr;jmin=1-j;
         else
           jmax=nr-j;jmin=1;
         end
          for k=-nz:nz
             if( kk+k >= 1 && kk+k <= ne)
              tmp_field(:,jmin:jmax,tmp_index)= ...
              circshift(data(:,jmin+j:jmax+j,kk+k),[i 0]);
             end
            tmp_index=tmp_index+1;
          end
       end
   end

 count(:,:,kk)=sum(tmp_field > threshold,3)/box_size;

end

count(isnan(count))=0;

time=toc;

display(['Count computed in ' num2str(time) ' seconds']);






