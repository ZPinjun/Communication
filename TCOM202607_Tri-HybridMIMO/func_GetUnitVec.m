function u = func_GetUnitVec(theta,phi)
% This function converts spatial angles to a 3D unit directional vector
% The unit of theta and phi are degree

theta = reshape(theta,1,[]);
phi = reshape(phi,1,[]);
u = [sind(theta).*cosd(phi); sind(theta).*sind(phi); cosd(theta)];

end

