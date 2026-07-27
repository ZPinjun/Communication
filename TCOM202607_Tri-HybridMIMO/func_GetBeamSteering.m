function H_out = func_GetBeamSteering(phi_all,theta_all,type,Fanten,paras)

N = size(Fanten,2);

% Get positions of all antennas at BS's local coordinate system, pn (Nx3)
file_path_antenna_local = fullfile('data_RayTracing', 'pos_BS_antenna_local.xlsx');
p_n = readmatrix(file_path_antenna_local);
k0 = 2*pi/paras.lambda;
% Get unit directional vectors (3xNumAoDs)
U = [sind(theta_all).*cosd(phi_all); sind(theta_all).*sind(phi_all); cosd(theta_all)];
% Get radiation pattern over all directions and antennas
G = nan(length(phi_all),N);
if strcmp(type,'Model I')

    load('data_patterns.mat')
    S = paras.S;
    Fsel = Fanten;
    % Convert to magnitude (or field) pattern
    RP_magnitude = dataMatrix;
    RP_magnitude(:,3:end) = nan;
    for ind = 1:64
        RP_magnitude(:,ind+2) = 10.^(dataMatrix(:,ind+2)/20);
    end
    for ind_n = 1:N
        b_vec = Fsel((1:S)+(ind_n-1)*S,ind_n);
        ind = find(b_vec==1);
        RP = RP_magnitude(:,[1:2,ind+2]);
        % Get radiation matrix
        Gn_vec = zeros(length(phi_all),1);
        az_snapped = interp1(-180:0.5:180, -180:0.5:180, phi_all, 'nearest');
        el_snapped = interp1(0:0.5:180, 0:0.5:180, theta_all, 'nearest');
        for ind_az = 1:length(phi_all)
            az_i_snapped = az_snapped(ind_az);
            el_i_snapped = el_snapped(ind_az);
            match_index = RP(:,1:2) == [az_i_snapped,el_i_snapped];
            Gn_vec(ind_az) = RP(match_index(:,1)&match_index(:,2),3);
        end
        G(:,ind_n) = Gn_vec;
    end

elseif strcmp(type,'Model II')
    T = paras.T;
    LT = paras.Rmax;
    Fcof = Fanten;
    for ind_n = 1:N
        c_vec = Fcof((1:T)+(ind_n-1)*T,ind_n);
        Gn_vec = zeros(length(phi_all),1);
        ind = 1;
        for l = 0:LT
            for m = -l:l
                Y_lm = harmonicY(l, m, deg2rad(theta_all).', deg2rad(phi_all).','type','real','phase',false);
                Gn_vec = Gn_vec + c_vec(ind) * Y_lm;
                ind = ind + 1;
            end
        end
        G(:,ind_n) = Gn_vec;
    end
else
    error('Incorrect input!')
end
% Form the beam-steering matrix
H_out = G.*exp(1j*k0*U.'*p_n.');

% figure
% plot(phi_all,G(:,61).'); hold on

end
