function [ ]=radar_superobbing( radar , gridrad , fileout , debug_output )   

 fileout_letkf=[fileout '.dat'];
 fileout_debug=[fileout '.grd'];

 nx=size(gridrad.grid_ref,1);
 ny=size(gridrad.grid_ref,2);
 nz=size(gridrad.grid_ref,3);

 
 nfile= fopen(fileout_letkf,'w','b');
 
   fwrite(nfile,4, 'int32'); 
   fwrite( nfile , radar.lon , 'float32');
   fwrite(nfile,4, 'int32'); 
   fwrite(nfile,4, 'int32'); 
   fwrite( nfile , radar.lat , 'float32');
   fwrite(nfile,4, 'int32'); 
   fwrite(nfile,4, 'int32'); 
   fwrite( nfile , radar.altitude , 'float32' );
   fwrite(nfile,4, 'int32'); 
   
    nobs=0;
    for ii=1:nx
     for jj=1:ny
      for kk=1:nz
         if( gridrad.grid_count_ref(ii,jj,kk) > 0 )
          
           wk(1)=radar.id_ref_obs         ;
           wk(6)=radar.error_ref          ;
           wk(2)=gridrad.grid_az_ref(ii,jj,kk)    ;
           wk(3)=gridrad.grid_el_ref(ii,jj,kk)    ;
           wk(4)=gridrad.grid_ra_ref(ii,jj,kk)    ;
           wk(5)=gridrad.grid_ref(ii,jj,kk)       ;
           wk(7)=radar.radar_type ;
           
           fwrite(nfile,7*4, 'int32'); 
           fwrite( nfile , wk ,'float32');
           fwrite(nfile,7*4, 'int32'); 
           nobs = nobs + 1;
          end
          if( gridrad.grid_count_dv(ii,jj,kk) > 0 )
           wk(1)=radar.id_dv_obs           ;
           wk(6)=radar.error_dv            ;
           wk(2)=gridrad.grid_az_dv(ii,jj,kk)      ;
           wk(3)=gridrad.grid_el_dv(ii,jj,kk)      ;
           wk(4)=gridrad.grid_ra_dv(ii,jj,kk)      ;
           wk(5)=gridrad.grid_dv(ii,jj,kk)         ;
           wk(7)=radar.radar_type          ;
           
           fwrite(nfile,7*4, 'int32'); 
           fwrite( nfile , wk ,'float32')  ;
           fwrite(nfile,7*4, 'int32'); 
           
           nobs=nobs +1;
          end
      end
     end
    end



display(['A TOTAL NUMBER OF ' num2str(nobs) ' HAS BEEN WRITTEN TO THE OBSERVATION FILE'])

if ( debug_output )

      
  nfile= fopen(fileout_debug,'w','b');

  %Ref
  for iz=1:nz
      fwrite(nfile,nx*ny, 'int32');
      fwrite(nfile,gridrad.grid_ref(:,:,iz)','float32')  ;
      fwrite(nfile,nx*ny, 'int32');
  end
  for iz=1:nz
      fwrite(nfile,nx*ny, 'int32');
      fwrite(nfile,gridrad.grid_count_ref(:,:,iz)','float32')  ;
      fwrite(nfile,nx*ny, 'int32');
  end
  for iz=1:nz
      fwrite(nfile,nx*ny, 'int32');
      fwrite(nfile,gridrad.grid_dv(:,:,iz)','float32')  ;
      fwrite(nfile,nx*ny, 'int32');
  end
  for iz=1:nz
      fwrite(nfile,nx*ny, 'int32');
      fwrite(nfile,gridrad.grid_count_dv(:,:,iz)','float32')  ;
      fwrite(nfile,nx*ny, 'int32');
  end
  for iz=1:nz
      fwrite(nfile,nx*ny, 'int32');
      fwrite(nfile,gridrad.grid_az_ref(:,:,iz)','float32')  ;
      fwrite(nfile,nx*ny, 'int32');
  end
  for iz=1:nz
      fwrite(nfile,nx*ny, 'int32');
      fwrite(nfile,gridrad.grid_ra_ref(:,:,iz)','float32')  ;
      fwrite(nfile,nx*ny, 'int32');
  end
  for iz=1:nz
      fwrite(nfile,nx*ny, 'int32');
      fwrite(nfile,gridrad.grid_el_ref(:,:,iz)','float32')  ;
      fwrite(nfile,nx*ny, 'int32');
  end
  for iz=1:nz
      fwrite(nfile,nx*ny, 'int32');
      fwrite(nfile,gridrad.grid_az_dv(:,:,iz)','float32')  ;
      fwrite(nfile,nx*ny, 'int32');
  end
  for iz=1:nz
      fwrite(nfile,nx*ny, 'int32');
      fwrite(nfile,gridrad.grid_ra_dv(:,:,iz)','float32')  ;
      fwrite(nfile,nx*ny, 'int32');
  end
  for iz=1:nz
      fwrite(nfile,nx*ny, 'int32');
      fwrite(nfile,gridrad.grid_el_dv(:,:,iz)','float32')  ;
      fwrite(nfile,nx*ny, 'int32');
  end




end


end
