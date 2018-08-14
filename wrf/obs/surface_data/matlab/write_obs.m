
function []=read_obs( file , obs , endian )

nfile=fopen(file,'w',endian);


nobs=size(obs,1);


for io=1:nobs

 fwrite( nfile , 7*4 , 'int32' );
 fwrite( nfile , obs(io,:) , 'float32' );
 fwrite( nfile , 7*4 , 'int32' );

end

end
