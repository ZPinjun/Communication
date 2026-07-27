function AMG = func_GetFixedAMG(RP_fixed,az,el)
% This function returns the antenna magnitude gain given a direction and a
% radiation pattern
% RP_fixed: [ phi(deg), theta(deg), power gain ], az: phi in [-pi,pi]; el: theta in [0,pi]

% Find the closest angle based on the input
az_snapped = interp1(-180:0.5:180, -180:0.5:180, az, 'nearest');
el_snapped = interp1(0:0.5:180, 0:0.5:180, el, 'nearest');
match_index = RP_fixed(:,1:2) == [az_snapped,el_snapped];
AMG = RP_fixed(match_index(:,1)&match_index(:,2),3);

end

