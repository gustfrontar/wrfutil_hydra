function [order_var , order_azimuth , levels , order_time , azimuth_exact ]= order_variable(var,azimuth,elevation,time,azimuth_res)

%order_var es la variable ordenada con los azimuths entre 0 y 360 (si hay rayos repetidos se promedian).
%order_azimuth es el azimuth "aproximado" utilizando 0 como azimuth inicial y avanzando en intervalos regulares e iguales a la resolucion
%del azimuth en grados.
%levels son los angulos de elevacion.
%azimuth_exact es un array que contiene para cada nivel el valor exacto del azimuth que corresponde a cada rayo. 

levels=unique(elevation);
order_azimuth=0:azimuth_res:360-azimuth_res;  %Asumo que el rango de los azimuths se define en forma regular.
naz=length(order_azimuth);
nel=length(levels);

order_var=NaN(naz,size(var,1),nel);
order_time=NaN(naz,nel);

azimuth_exact=NaN(naz,nel);

for ilev=1:nel

  levmask= elevation == levels(ilev) ;
  %Busco los azimuth que corresponden a este nivel.
  azlev=azimuth( levmask );
  timelev=time( levmask );
  %Busco los valores de la variable que corresponden a este nivel
  varlev=var(:,levmask);

  %Para el primer azimuth (que es el que contiene el salto en angulo de 0 a 360 y por eso es un caso especial).
  az_index=( azlev <= azimuth_res/2.0 | azlev >= 360 - azimuth_res/2.0 );
  if( sum(az_index) > 0 )
    order_var(1,:,ilev) = nanmean( varlev(:,az_index) , 2);
    order_time(1,ilev) = nanmean( timelev( az_index ) );

    azimuth_exact(1,ilev) = nanmean( azlev( az_index ) );
  end

  %Para los que vienen despues.
  for iaz=2:naz
    %Busco todos los rayos que para este nivel estan alrededor del valor de cada uno de los azimuths presentes en la variable.
    %order azimuth.
    az_index=( azlev <= order_azimuth(iaz) + azimuth_res/2.0 & azlev >= order_azimuth(iaz) - azimuth_res/2.0 );
    %Si hay mas de un rayo que cae dentro del mismo rango de azimuth los promedio y si no hay ninguno (gap), entonces queda el NaN con
    %el que fue originalmente definida la variable.
    if( sum(az_index) > 0 )
      order_var(iaz,:,ilev) = nanmean( varlev(:,az_index) , 2 ); 
      order_time(iaz,ilev) = nanmean( timelev( az_index ) );

      azimuth_exact(iaz,ilev) = nanmean( azlev( az_index ) );
      
    end
   
  end 
    
end













end
