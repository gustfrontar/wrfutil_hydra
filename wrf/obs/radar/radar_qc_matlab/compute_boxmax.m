function [maxfield]=compute_boxmax(data,nx,ny,nz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% data, a 3D array
% nx,ny,nz are the sizes of the local box over which max will be computed
% Some of these parameters are used as part of the REC algorithm described
% in:
% THE RADAR ECHO CLASSIFIER: A FUZZY LOGIC ALGORITHM FOR THE WSR-88D
%Cathy Kessinger, Scott Ellis, and Joseph Van Andel
%National Center for Atmospheric Research
%Boulder, CO 80305
%In this paper a value of nx=2 and ny=20 (for Z) and 10 (for Wind) are
%used (nz is set to 0).

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

maxfield=NaN(size(data));

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

 maxfield(:,:,kk)=max(tmp_field,[],3);

end

time=toc;

display(['Mean computed in ' num2str(time) ' seconds']);






