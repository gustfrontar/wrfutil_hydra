
function [obs,nobs]=read_obs( file , endian )

nfile=fopen(file,'r',endian);

obs=[];

nobs=0;

cont=true;

while( cont )

 fread( nfile , 1  , 'int32' );
 tmp=fread( nfile , 7 , 'float32' );
 fread( nfile , 1  , 'int32' );

 if( ~isempty( tmp ) )
   nobs=nobs+1;
   obs(nobs,:)=tmp;
 else
   cont=false;
 end

end

end
