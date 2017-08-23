
function [lon lat]=ll_arc_distance(lon0,lat0, arc_dist, az)

% NAME:
%        LL_ARC_DISTANCE
%
% PURPOSE:
%         This function returns the longitude and latitude [lon, lat] of
%        a point a given arc distance (-pi <= Arc_Dist <= pi), and azimuth (Az),
%        from a specified location Lon_lat0.
%
% CATEGORY:
%        Mapping, geography.
%
% CALLING SEQUENCE:
%        Result = LL_ARC_DISTANCE(Lon_lat0, Arc_Dist, Az)
%
% INPUTS:
%            Lon_lat0: A 2-element vector containing the longitude and latitude
%                  of the starting point. Values are assumed to be in radians
%                  unless the keyword DEGREES is set.
%            Arc_Dist: Distance in Km.    
%            Az:          The azimuth from the center. The value is assumed to be in
%                  degrees.
%
%
% OUTPUTS:
%        This function returns lat and lon in degrees.
%
% PROCEDURE:
%        Formula from Map Projections - a working manual.  USGS paper
%        1395.  Equations (5-5) and (5-6).
%
%
R=12742/2;     %Earth radius


if (arc_dist == 0)
    lon=lon0;
    lat=lat0;

else
    
    arc_dist=arc_dist/R;
    

 cdist = cos(arc_dist);               %
 sdist = sin(arc_dist);

  %az = az*pi/180;
  
  sinll1 = sind(lat0);
  cosll1 = cosd(lat0);

  lat=asind(sinll1*cdist + cosll1*sdist.*cosd(az));
  lon=lon0 + (180/pi)*atan2(sdist.*sind(az),cosll1.*cdist - sinll1.*sdist.*cosd(az));
  
end