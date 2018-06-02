%Traduce la SLP a la presion de superficie medida por la estacion segun la cuenta sugerida por la WMO.
%Esta cuenta debe usarse para estaciones con alturas sobre el nivel del mar menores a 700 metros!!.

function sp = SlpToSP( slp , t , td , z )

%slp presion a nivel del mar (Pa)
%t temperatura en K
%q humedad espcifica
%z altura de la estacion

R=287; %Constante de los gases
a=0.0065; %Lapse rate K/m (atmosfera standard)
ch=0.12 ; %Coeficiente
g=9.81;   %Aceleracion de la gravedad

[e]=get_es(td) ; %Presion de vapor
e=e/100;         %Paso la presion de vapor a hPa.

sp = slp .* exp( -(  g*z/R  ) ./ ( t + a * z / 2 + e * ch )  );

end

