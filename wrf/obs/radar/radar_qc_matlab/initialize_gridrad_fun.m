
function gridrad = initialize_gridrad_fun( cart )

       gridrad.grid_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_count_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_dv=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_count_dv=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_az_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_el_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_ra_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_lat_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_lon_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_z_ref=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_az_dv=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_el_dv=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_ra_dv=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_lat_dv=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_lon_dv=zeros(cart.nlon,cart.nlat,cart.nlev);
       gridrad.grid_z_dv=zeros(cart.nlon,cart.nlat,cart.nlev);

end
