
function [ cart ] = define_cartesian_grid( radar , dx , dz , maxz ) 

%-----------------------------------------------------------------------
% Define grid
%-----------------------------------------------------------------------
   %Translate DX into an appropiate DLON and DLAT.
   %Hopefully nobody will put a radar at the pole.
   cdeg2rad=pi/180;
   crad2deg=180/pi;
   re=6371e3;
   
   DLON=crad2deg*dx/(cos(radar.lat*cdeg2rad)*re);
   DLAT=crad2deg*dx/re;

   %Compute possible value for NLON in order to cover the maximum radar range.
   MAXRANGE=max( radar.range ) ; 
   MAXRANGE=2.0*MAXRANGE;
   NLON=ceil( MAXRANGE / dx );
   NLAT=NLON;
   
   %maxz=squeeze( max(max(max(radar.Z))) );

   NLEV=ceil(maxz/dz);

   %Force grid dimensions to be odd
   if( mod( NLON , 2 )== 0);NLON=NLON+1;end
   if( mod( NLAT , 2 )== 0);NLAT=NLAT+1;end

   cart.lon=NaN(NLON,NLAT);
   cart.lat=NaN(NLON,NLAT);
   cart.z=NaN(NLON,NLAT,NLEV);
   %DEFINE LAT AND LON
   for i=1:NLON
     for j=1:NLAT
       cart.lon(i,j)=radar.lon + DLON*( -1.0-(NLON-1.0)/2.0 + i );
       cart.lat(i,j)=radar.lat + DLAT*( -1.0-(NLAT-1.0)/2.0 + j );
     end
   end

   %DEFINE Z
   for k=1:NLEV
      cart.z(:,:,k)=dz*(k-1);
   end

   cart.nlon=NLON;
   cart.nlat=NLAT;
   cart.nlev=NLEV;
   
   cart.dlon=DLON;
   cart.dlat=DLAT;
   
   cart.dz=dz;
   cart.dx=dx;
   
end
