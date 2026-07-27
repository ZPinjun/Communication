% References:
% [1] R. Wang, P. Zheng et al., “Electromagnetically reconfigurable fluid
% antenna system for wireless communications: Design, modeling, algorithm,
% fabrication, and experiment,” 2025. [Online]. Available:
% https://arxiv.org/abs/2502.19643   

function vec_basis = func_GetPatternBasis(model,az,el)

if strcmp(model,'model I')
    % Load the radiation patterns obatained through HFSS full-wave simulation in [1]
    load('PRs_data.mat') % [ phi(deg), theta(deg), power gain(dB) ], az: phi in [-pi,pi]; el: theta in [0,pi]

    % Convert to magnitude (or field) pattern
    RP_magnitude = dataMatrix;
    RP_magnitude(:,3:end) = nan;
    for ind = 1:64
        RP_magnitude(:,ind+2) = 10.^(dataMatrix(:,ind+2)/20);
    end

    % Find the closest angle based on the input
    az_snapped = interp1(-180:0.5:180, -180:0.5:180, az, 'nearest');
    el_snapped = interp1(0:0.5:180, 0:0.5:180, el, 'nearest');
    match_index = RP_magnitude(:,1:2) == [az_snapped,el_snapped];
    vec_basis = RP_magnitude(match_index(:,1)&match_index(:,2),3:end).';

elseif strcmp(model,'model II')

end

end

