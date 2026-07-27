function c = func_GetSHcoef(T)

c = nan(T,1);
Lmax = 0;
while 1
    if Lmax^2+2*Lmax+1 < T
        Lmax = Lmax + 1;
    elseif Lmax^2+2*Lmax+1 == T
        break
    else
        error('Invalid T !!!')
    end
end

%% Generate spherical coordinates
N = 64;                                             % sampling number
theta_range = linspace(0, pi, N);
phi_range = linspace(-pi, pi, 2*N);
[theta, phi] = meshgrid(theta_range, phi_range);    % theta in [0,pi], phi in [0,2pi]


%% Pick a radiation pattern
load('data_patterns.mat')
% Convert to magnitude (or field) pattern 
RP_magnitude = dataMatrix;
RP_magnitude(:,3:end) = nan;
for ind = 1:64
    RP_magnitude(:,ind+2) = 10.^(dataMatrix(:,ind+2)/20);
end
ind = 19;
RP = RP_magnitude(:,[1:2,ind+2]);
% Get radiation matrix
E = nan(length(phi_range),length(theta_range));
for ind_theta = 1:length(theta_range)
    theta_i = theta_range(ind_theta)*180/pi;
    theta_snapped = interp1(0:0.5:180, 0:0.5:180, theta_i, 'nearest');
    for ind_phi = 1:length(phi_range)
        phi_i = phi_range(ind_phi)*180/pi;
        phi_snapped = interp1(-180:0.5:180, -180:0.5:180, phi_i, 'nearest');
        match_index = RP(:,1:2) == [phi_snapped,theta_snapped];
        E(ind_phi,ind_theta) = RP(match_index(:,1)&match_index(:,2),3);
    end
end


%% Perform spherical harmonic decomposition
alm = zeros(Lmax+1, 2*Lmax+1);  % Coefficent matrix
for l = 0:Lmax
    for m = -l:l
        % Compute harmonics coefficents
        Y_lm = harmonicY(l, m, theta, phi,'type','real','phase',false);
        % Compute projection coefficient a_lm (inner product)
        c(l^2 + l + m + 1) = sum(sum(E .* conj(Y_lm) .* sin(theta))) * (pi/N) * (2*pi/(2*N));
    end
end

end

