clear all
close all

filename='/home/jruiz/WRF/WPS/geo_em.d01.nc';

topo = read_netcdf_var(filename,'HGT_M',false);
lon = read_netcdf_var(filename,'XLONG_M',false);
lat = read_netcdf_var(filename,'XLAT_M',false);

save ('topo_wrf.mat','lon','lat','topo');