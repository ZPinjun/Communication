% References:
% [1] R. Wang, P. Zheng et al., “Electromagnetically reconfigurable fluid
% antenna system for wireless communications: Design, modeling, algorithm,
% fabrication, and experiment,” 2025. [Online]. Available:
% https://arxiv.org/abs/2502.19643 

function G = func_GetFixedPattern(az,el)
% This function returns the antenna magnitude gain given a direction using
% a fixed symmetric radiation pattern in [1]

% Load the radiation patterns obatained through HFSS full-wave simulation in [1]
load('PRs_data.mat') % [ phi(deg), theta(deg), power gain(dB) ], az: phi in [-pi,pi]; el: theta in [0,pi]

% Select a radiation pattern and convert to magnitude (or field) pattern
ind = 19;
RP_magnitude = dataMatrix(:,[1,2,ind+2]);
RP_magnitude(:,3) = 10.^(dataMatrix(:,ind+2)/20);

% Find the closest angle based on the input
az_snapped = interp1(-180:0.5:180, -180:0.5:180, az, 'nearest');
el_snapped = interp1(0:0.5:180, 0:0.5:180, el, 'nearest');
match_index = RP_magnitude(:,1:2) == [az_snapped,el_snapped];
G = RP_magnitude(match_index(:,1)&match_index(:,2),3);

end

