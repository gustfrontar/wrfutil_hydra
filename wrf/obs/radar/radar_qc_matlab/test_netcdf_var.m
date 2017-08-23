function [var_test]=test_netcdf_var(file,var)

%Test if a particular variable is included in a netcdf file.
%var_test = false means that the variable is not present, var_test=true
%means that the variable is present.

ncid = netcdf.open(file,'nowrite');

var_test=true;

try
    ID = netcdf.inqVarID(ncid,var);
catch %exception
    %if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
    warning(['Variable ' var ' could not be read'])
    var_test=false;

    %    str = 'bad';
    %end
end
netcdf.close(ncid)



end
