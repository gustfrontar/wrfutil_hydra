
function [U,V]=VelDirToUV(Vel,Dir)

U=-Vel.*sind(Dir);

V=-Vel.*cosd(Dir);

end
