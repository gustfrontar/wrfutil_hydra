function [meanv, stdv]=compute_boxmean(data,nx,ny,nz,nx2,ny2,nz2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% data, a 3D array (usually azimuth,r,elev) containing a radar field
% (reflectivity or wind). In the case of wind take into account that
% missing values are filled with 0. That migth not be the best choice for
% computing wind texture.
% nx,ny,nz are the outer sizes of the local box over which median, mean and
% std will be computed.
% nx2,ny2, nz2 are the outer sizes of the local box that won't be used for
% computation.
% Some of these parameters are used as part of the REC algorithm described
% in:
% THE RADAR ECHO CLASSIFIER: A FUZZY LOGIC ALGORITHM FOR THE WSR-88D
%Cathy Kessinger, Scott Ellis, and Joseph Van Andel
%National Center for Atmospheric Research
%Boulder, CO 80305
%In this paper a value of nx=2 and ny=20 (for Z) and 10 (for Wind) are
%used (nz is set to 0).
%This version omits the value of the central point.

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
box_size=(2*nx+1)*(2*ny+1)*(2*nz+1)-(2*nx2+1)*(2*ny2+1)*(2*nz2+1);

meanv=NaN(size(data));
stdv=NaN(size(data));

if box_size>0

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
             if( kk+k >= 1 && kk+k <= ne & ~(abs(j)<=ny2 & abs(i)<=nx2 & abs(k)<=nz2 ))
              tmp_field(:,jmin:jmax,tmp_index)= ...
              circshift(data(:,jmin+j:jmax+j,kk+k),[i 0]);
              tmp_index=tmp_index+1;
             end
          end
       end
     end
   
  
%     for ia=1:na
%      for ir=1:nr
%        meanv(ia,ir,kk)=nanmean(squeeze(tmp_field(ia,ir,:)));
%        stdv(ia,ir,kk)=nanstd(squeeze(tmp_field(ia,ir,:)));
%      end
%    end

      meanv(:,:,kk)=nanmean(tmp_field,3);
      stdv(:,:,kk)=nanstd(tmp_field,[],3);

   end

end

time=toc;

display(['Mean computed in ' num2str(time) ' seconds']);






