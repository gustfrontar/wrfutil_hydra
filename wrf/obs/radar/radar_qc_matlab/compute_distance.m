function [dist]=compute_distance(data,data2,nx,ny,nz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data is the array whose consistency we want to evaluate.
%data2 is a masked array where some pixels has been removed (those which are
%inconsistent). Where data2 is NaN this point is not going to be considered
%for dist computation.
%The idea is to perform this filter in an iteratively way where a first
%selection of the data is performed, but points that are good can be
%reclasified as that by a second pass of the filter.


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

dist=NaN(size(data));

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
             if( kk+k >= 1 && kk+k <= ne )
              tmp_field(:,jmin:jmax,tmp_index)= ...
              circshift(data2(:,jmin+j:jmax+j,kk+k),[i 0])-data(:,jmin:jmax,kk);
              tmp_index=tmp_index+1;
             end
          end
       end
     end
     
   
%    for j=-ny:ny
%        for i=-nx:nx
%          if(j <= 0)
%            jmax=nr;jmin=1-j;
%          else
%            jmax=nr-j;jmin=1;
%          end
%          if(i <= 0)
%             imax=na;imin=1-i;
%          else
%             imax=na-i;imin=1;
%          end
%           for k=-nz:nz
%              if( kk+k >= 1 && kk+k <= ne)
%              
%               tmp_field(imin:imax,jmin:jmax,tmp_index)= ...
%               data2(imin+i:imax+i,jmin+j:jmax+j,kk+k)-data(imin:imax,jmin:jmax,kk);
%              end
%             tmp_index=tmp_index+1;
%           end
%        end
%    end

 dist(:,:,kk)=sqrt(nanmean((tmp_field).^2,3));

end






time=toc;

display(['Distance computed in ' num2str(time) ' seconds']);






