function [mean]=compute_boxmean(data,nx,ny,nz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% In this version of boxmean no cyclic conditions are assumed.

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

mean=NaN(size(data));

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
         if(i <= 0)
            imax=na;imin=1-i;
         else
            imax=na-i;imin=1;
         end
          for k=-nz:nz
             if( kk+k >= 1 && kk+k <= ne)
              tmp_field(imin:imax,jmin:jmax,tmp_index)= ...
              data(imin+i:imax+i,jmin+j:jmax+j,kk+k);
             end
            tmp_index=tmp_index+1;
          end
       end
   end

 mean(:,:,kk)=nanmean(tmp_field,3);

end

time=toc;

display(['Mean no cyclic computed in ' num2str(time) ' seconds']);






