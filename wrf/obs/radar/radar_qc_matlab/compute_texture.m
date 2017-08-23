function [texture]=compute_texture(radar,data,nx,ny,nz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% data, a 3D array (usually azimuth,r,elev) containing a radar field
% (reflectivity or wind). In the case of wind take into account that
% missing values are filled with 0. That migth not be the best choice for
% computing wind texture.
% nx,ny,nz are the sizes of the local box over which texture will be
% computed. The box can be 1D, 2D or 3D. Size of the box in each direction
% is computed as: 2nx+1, 2ny+1 and 2nz+1 respectively.
% The box is always centered at each grid point.
% The output is an array the same size of radar.data with the texture
% computed at each location.
% The texture parameter is defined in 
%
% THE RADAR ECHO CLASSIFIER: A FUZZY LOGIC ALGORITHM FOR THE WSR-88D
%Cathy Kessinger, Scott Ellis, and Joseph Van Andel
%National Center for Atmospheric Research
%Boulder, CO 80305
%In this paper a value of nx=2 and ny=20 (for Z) and 10 (for Wind) are
%used (nz is set to 0).

texture_dim=1;  %Dimension of texture computation (independent of box dimension)
                %When texture_dim = 1 then classical texture is computed.


% Juan Ruiz 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

%data(isnan(data))=radar.replacerefmissing;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF RADIAL BEAM DIFFERENCE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_data=NaN(size(data));

tmp_data(:,1:end-1,:)=diff(data,1,2);

tmp_data(:,end,:)=tmp_data(:,end-1,:);

if (texture_dim == 2 )
    
     tmp_data(:,1:end-1,:)=tmp_data(1:end,:,:) + diff(data,1,1);

     tmp_data(end,:,:)=tmp_data(end-1,:,:)+data(end,:,:)+data(1,:,:);
    
    if ( texture_dim == 3 )
    
     tmp_data(:,:,1:end-1)=tmp_data(:,:,1:end-1) + diff(data,1,3);

     tmp_data(:,:,1:end-1)=tmp_data(:,:,1:end-1);   
           
    end
end

tmp_data=(tmp_data./texture_dim).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUM DIFERENCE OVER THE SQUARE LOCAL BOX (CAN BE A 1D, 2D or 3D BOX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%To reduce memory size the computation will be done level by level.
box_size=(2*nx+1)*(2*ny+1)*(2*nz+1);

texture=NaN(size(data));
for kk=1:radar.ne
   tmp_field=NaN(radar.na,radar.nr,box_size);
   tmp_index=1;
   for j=-ny:ny
       for i=-nx:nx
         if(j <= 0)
           jmax=radar.nr;jmin=1-j;
         else
           jmax=radar.nr-j;jmin=1;
         end
          for k=-nz:nz
             if( kk+k >= 1 && kk+k <= radar.ne)
              tmp_field(:,jmin:jmax,tmp_index)= ...
              circshift(tmp_data(:,jmin+j:jmax+j,kk+k),[i 0]);
             end
            tmp_index=tmp_index+1;
          end
       end
   end

   
 texture(:,:,kk)=nanmean(tmp_field,3);


end

time=toc;

display(['Texture computed in ' num2str(time) ' seconds']);






