
function var = read_netcdf_var(filename,variable,unpack_flag)
%Unpack flag, if true then search for offset and scale factors and unpack
%data.

ncid = netcdf.open(filename,'nowrite');


varid = netcdf.inqVarID(ncid,variable);


var =  netcdf.getVar(ncid,varid) ;


 if( unpack_flag )
 var=double(var);
 fill_value = netcdf.getAtt(ncid,varid,'_FillValue');
 var ( var==fill_value )=NaN;
 add_offset = netcdf.getAtt(ncid,varid,'add_offset');
 scale_factor = netcdf.getAtt(ncid,varid,'scale_factor');
 var = var * scale_factor + add_offset ;
 end

netcdf.close(ncid)


end
