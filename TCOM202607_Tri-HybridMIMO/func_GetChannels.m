% References:
% [1] J. Hoydis, S. Cammerer et al., "Sionna: An open-source library for
% next-generation physical layer research," 2022. [Online]. Available:
% https://nvlabs.github.io/sionna/ 
% [2] R. Wang, P. Zheng et al., “Electromagnetically reconfigurable fluid
% antenna system for wireless communications: Design, modeling, algorithm,
% fabrication, and experiment,” 2025. [Online]. Available:
% https://arxiv.org/abs/2502.19643 

function [H_cv,H_sel_hardware,H_sel_fictitious,H_cof,paras] = func_GetChannels(N)
% This function generates wireless channels and effective channels based on
% the ray-tracing simulation results using Sionna [1]. Input N is the
% total antenna number at the BS (N <= 100)

%% Set system parameters 
Mk = [4,4,4].';             % Number of antennas at each UE
Lk = [11,6,4].';            % Number of paths of each channel
S = 64;                     % Total number of available radiation patterns
Rmax = 6;              
T = Rmax^2 + 2*Rmax + 1;    % Spherical harmonics truncation length
fc = 30e9;                  % Signal frequency (30 GHz)


%% Load the radiation patterns obatained through HFSS full-wave simulation in [2]
load('data_patterns.mat') % [ phi(deg), theta(deg), power gain(dB) ], az: phi in [-pi,pi]; el: theta in [0,pi]
% - All available patterns
RP_all = dataMatrix;
RP_all(:,3:end) = 10.^(RP_all(:,3:end)/20);
% - Select a fixed pattern
ind = 19;
RP_fixed = dataMatrix(:,[1,2,ind+2]);
RP_fixed(:,3) = 10.^(dataMatrix(:,ind+2)/20);


%% Generate the fictitious pencil beams for benchmark
parameters.sigmae2 = 0.8;
parameters.theta_bar = linspace(90,180,8);
parameters.phi_bar = linspace(-90,90,8);
f = @(theta) exp(-(2*theta.^2)/parameters.sigmae2).*sin(theta);
result = integral(f, 0, pi);
parameters.Gmax = sqrt(2/result);

% el = -90:1:90;
% az = 0*ones(size(el));
% theta_bar = 0;
% phi_bar = 0;
% u_bar = func_GetUnitVec(theta_bar,phi_bar);
% u_all = func_GetUnitVec(el,az);
% y = parameters.Gmax*exp(-acos(u_all.'*u_bar).^2/parameters.sigmae2);
% 
% figure
% plot(el,y.')


%% --- UE 1
ind_k = 1;
% Read ray-tracing simulation results
file_path = fullfile('data_RayTracing', 'UE1_a.xlsx'); 
Chan_gain = readmatrix(file_path);
file_path = fullfile('data_RayTracing', 'UE1_tau.xlsx');
Delay = PerSheetRead(file_path);
file_path = fullfile('data_RayTracing', 'UE1_AoD_el.xlsx');
AoD_el = PerSheetRead(file_path);
file_path = fullfile('data_RayTracing', 'UE1_AoD_az.xlsx');
AoD_az = PerSheetRead(file_path);
file_path = fullfile('data_RayTracing', 'UE1_AoA_el.xlsx');
AoA_el = PerSheetRead(file_path);
file_path = fullfile('data_RayTracing', 'UE1_AoA_az.xlsx');
AoA_az = PerSheetRead(file_path);

% process AoD/AoA data
for ind_l = 1:Lk(ind_k)
    % el
    AoD_el{ind_l}(:,3) = rad2deg(AoD_el{ind_l}(:,3));
    AoA_el{ind_l}(:,3) = rad2deg(AoA_el{ind_l}(:,3));
    % az
    AoD_az{ind_l}(:,3) = rad2deg(AoD_az{ind_l}(:,3))-90;
    AoD_az{ind_l}(:,3) = wrapTo180(AoD_az{ind_l}(:,3));
    AoA_az{ind_l}(:,3) = rad2deg(AoA_az{ind_l}(:,3))-90;
    AoA_az{ind_l}(:,3) = wrapTo180(AoA_az{ind_l}(:,3));
end

% Generate wireless channel 
H1_cv = zeros(Mk(ind_k),N);
H1_sel_hardware = zeros(Mk(ind_k),N*S);
H1_sel_fictitious = zeros(Mk(ind_k),N*S);
H1_cof = zeros(Mk(ind_k),N*T);
for ind_l = 1:Lk(ind_k)
    C = ones(Mk(ind_k),N)*Chan_gain(ind_l,end);
    A = nan(Mk(ind_k),N);
    GUE = nan(Mk(ind_k),N);
    GBS_fixed = nan(Mk(ind_k),N);
    GBS_sel_hardware = nan(Mk(ind_k),N*S);
    GBS_sel_fictitious = nan(Mk(ind_k),N*S);
    GBS_cof = nan(Mk(ind_k),N*T);
    for ind_n = 1:N
        for ind_m = 1:Mk(ind_k)
            match_index = Delay{ind_l}(:,1:2) == [ind_n-1,ind_m-1];
            % Propagation delay (matrix A)
            delay_mn = Delay{ind_l}(match_index(:,1)&match_index(:,2),end).';
            A(ind_m,ind_n) = (1/sqrt(N*Mk(ind_k)))*exp(-1j*delay_mn*fc*2*pi);
            % Antenna magnitude gain at UE (matrix GUE)
            AoA_az_mn = AoA_az{ind_l}(match_index(:,1)&match_index(:,2),end).';
            AoA_el_mn = AoA_el{ind_l}(match_index(:,1)&match_index(:,2),end).';
            GUE(ind_m,ind_n) = func_GetFixedAMG(RP_fixed,AoA_az_mn,AoA_el_mn);
            % Antenna magnitude gain at BS (matrix GBS_fixed, GBS_sel, and GBS_cof)
            AoD_az_mn = AoD_az{ind_l}(match_index(:,1)&match_index(:,2),end).';
            AoD_el_mn = AoD_el{ind_l}(match_index(:,1)&match_index(:,2),end).';
            % - Fixed pattern
            GBS_fixed(ind_m,ind_n) = func_GetFixedAMG(RP_fixed,AoD_az_mn,AoD_el_mn);
            % - Model I
            GBS_sel_hardware(ind_m,(ind_n-1)*S+(1:S)) = func_GetAllAMG(RP_all,AoD_az_mn,AoD_el_mn).';
            GBS_sel_fictitious(ind_m,(ind_n-1)*S+(1:S)) = func_GetAllFMG(parameters,AoD_az_mn,AoD_el_mn).';
            % - Model II
            gamma = nan(1,T);
            for r = 0:Rmax
                for m = -r:r
                    ind_xi = r^2 + r + m + 1;
                    gamma(ind_xi) = harmonicY(r,m,deg2rad(AoD_el_mn),deg2rad(AoD_az_mn),'type','real','phase',false);
                end
            end
            GBS_cof(ind_m,(ind_n-1)*T+(1:T)) = gamma;
        end
    end
    % Compute conventional channel
    H1_cv = H1_cv + sqrt(N*Mk(ind_k)/Lk(ind_k))*C.*A.*GUE.*GBS_fixed;
    H1_sel_hardware = H1_sel_hardware + sqrt(N*Mk(ind_k)/Lk(ind_k))*kron((C.*A.*GUE),ones(1,S)).*GBS_sel_hardware;
    H1_sel_fictitious = H1_sel_fictitious + sqrt(N*Mk(ind_k)/Lk(ind_k))*kron((C.*A.*GUE),ones(1,S)).*GBS_sel_fictitious;
    H1_cof = H1_cof + sqrt(N*Mk(ind_k)/Lk(ind_k))*kron((C.*A.*GUE),ones(1,T)).*GBS_cof;
end


%% --- UE 2
ind_k = 2;
% Read ray-tracing simulation results
file_path = fullfile('data_RayTracing', 'UE2_a.xlsx'); 
Chan_gain = readmatrix(file_path);
file_path = fullfile('data_RayTracing', 'UE2_tau.xlsx');
Delay = PerSheetRead(file_path);
file_path = fullfile('data_RayTracing', 'UE2_AoD_el.xlsx');
AoD_el = PerSheetRead(file_path);
file_path = fullfile('data_RayTracing', 'UE2_AoD_az.xlsx');
AoD_az = PerSheetRead(file_path);
file_path = fullfile('data_RayTracing', 'UE2_AoA_el.xlsx');
AoA_el = PerSheetRead(file_path);
file_path = fullfile('data_RayTracing', 'UE2_AoA_az.xlsx');
AoA_az = PerSheetRead(file_path);

% process AoD/AoA data
for ind_l = 1:Lk(ind_k)
    % el
    AoD_el{ind_l}(:,3) = rad2deg(AoD_el{ind_l}(:,3));
    AoA_el{ind_l}(:,3) = rad2deg(AoA_el{ind_l}(:,3));
    % az
    AoD_az{ind_l}(:,3) = rad2deg(AoD_az{ind_l}(:,3))-90;
    AoD_az{ind_l}(:,3) = wrapTo180(AoD_az{ind_l}(:,3));
    AoA_az{ind_l}(:,3) = rad2deg(AoA_az{ind_l}(:,3))-90;
    AoA_az{ind_l}(:,3) = wrapTo180(AoA_az{ind_l}(:,3));
end

% Generate wireless channel 
H2_cv = zeros(Mk(ind_k),N);
H2_sel_hardware = zeros(Mk(ind_k),N*S);
H2_sel_fictitious = zeros(Mk(ind_k),N*S);
H2_cof = zeros(Mk(ind_k),N*T);
for ind_l = 1:Lk(ind_k)
    C = ones(Mk(ind_k),N)*Chan_gain(ind_l,end);
    A = nan(Mk(ind_k),N);
    GUE = nan(Mk(ind_k),N);
    GBS_fixed = nan(Mk(ind_k),N);
    GBS_sel_hardware = nan(Mk(ind_k),N*S);
    GBS_cof = nan(Mk(ind_k),N*T);
    for ind_n = 1:N
        for ind_m = 1:Mk(ind_k)
            match_index = Delay{ind_l}(:,1:2) == [ind_n-1,ind_m-1];
            % Propagation delay (matrix A)
            delay_mn = Delay{ind_l}(match_index(:,1)&match_index(:,2),end).';
            A(ind_m,ind_n) = (1/sqrt(N*Mk(ind_k)))*exp(-1j*delay_mn*fc*2*pi);
            % Antenna magnitude gain at UE (matrix GUE)
            AoA_az_mn = AoA_az{ind_l}(match_index(:,1)&match_index(:,2),end).';
            AoA_el_mn = AoA_el{ind_l}(match_index(:,1)&match_index(:,2),end).';
            GUE(ind_m,ind_n) = func_GetFixedAMG(RP_fixed,AoA_az_mn,AoA_el_mn);
            % Antenna magnitude gain at BS (matrix GBS_fixed, GBS_sel, and GBS_cof)
            AoD_az_mn = AoD_az{ind_l}(match_index(:,1)&match_index(:,2),end).';
            AoD_el_mn = AoD_el{ind_l}(match_index(:,1)&match_index(:,2),end).';
            % - Fixed pattern
            GBS_fixed(ind_m,ind_n) = func_GetFixedAMG(RP_fixed,AoD_az_mn,AoD_el_mn);
            % - Model I
            GBS_sel_hardware(ind_m,(ind_n-1)*S+(1:S)) = func_GetAllAMG(RP_all,AoD_az_mn,AoD_el_mn).';
            GBS_sel_fictitious(ind_m,(ind_n-1)*S+(1:S)) = func_GetAllFMG(parameters,AoD_az_mn,AoD_el_mn).';
            % - Model II
            gamma = nan(1,T);
            for r = 0:Rmax
                for m = -r:r
                    ind_xi = r^2 + r + m + 1;
                    gamma(ind_xi) = harmonicY(r,m,deg2rad(AoD_el_mn),deg2rad(AoD_az_mn),'type','real','phase',false);
                end
            end
            GBS_cof(ind_m,(ind_n-1)*T+(1:T)) = gamma;
        end
    end
    % Compute conventional channel
    H2_cv = H2_cv + sqrt(N*Mk(ind_k)/Lk(ind_k))*C.*A.*GUE.*GBS_fixed;
    H2_sel_hardware = H2_sel_hardware + sqrt(N*Mk(ind_k)/Lk(ind_k))*kron((C.*A.*GUE),ones(1,S)).*GBS_sel_hardware;
    H2_sel_fictitious = H2_sel_fictitious + sqrt(N*Mk(ind_k)/Lk(ind_k))*kron((C.*A.*GUE),ones(1,S)).*GBS_sel_fictitious;
    H2_cof = H2_cof + sqrt(N*Mk(ind_k)/Lk(ind_k))*kron((C.*A.*GUE),ones(1,T)).*GBS_cof;
end


%% --- UE 3
ind_k = 3;
% Read ray-tracing simulation results
file_path = fullfile('data_RayTracing', 'UE3_a.xlsx'); 
Chan_gain = readmatrix(file_path);
file_path = fullfile('data_RayTracing', 'UE3_tau.xlsx');
Delay = PerSheetRead(file_path);
file_path = fullfile('data_RayTracing', 'UE3_AoD_el.xlsx');
AoD_el = PerSheetRead(file_path);
file_path = fullfile('data_RayTracing', 'UE3_AoD_az.xlsx');
AoD_az = PerSheetRead(file_path);
file_path = fullfile('data_RayTracing', 'UE3_AoA_el.xlsx');
AoA_el = PerSheetRead(file_path);
file_path = fullfile('data_RayTracing', 'UE3_AoA_az.xlsx');
AoA_az = PerSheetRead(file_path);

% process AoD/AoA data
for ind_l = 1:Lk(ind_k)
    % el
    AoD_el{ind_l}(:,3) = rad2deg(AoD_el{ind_l}(:,3));
    AoA_el{ind_l}(:,3) = rad2deg(AoA_el{ind_l}(:,3));
    % az
    AoD_az{ind_l}(:,3) = rad2deg(AoD_az{ind_l}(:,3))-90;
    AoD_az{ind_l}(:,3) = wrapTo180(AoD_az{ind_l}(:,3));
    AoA_az{ind_l}(:,3) = rad2deg(AoA_az{ind_l}(:,3))-90;
    AoA_az{ind_l}(:,3) = wrapTo180(AoA_az{ind_l}(:,3));
end

% Generate wireless channel 
H3_cv = zeros(Mk(ind_k),N);
H3_sel_hardware = zeros(Mk(ind_k),N*S);
H3_sel_fictitious = zeros(Mk(ind_k),N*S);
H3_cof = zeros(Mk(ind_k),N*T);
for ind_l = 1:Lk(ind_k)
    C = ones(Mk(ind_k),N)*Chan_gain(ind_l,end);
    A = nan(Mk(ind_k),N);
    GUE = nan(Mk(ind_k),N);
    GBS_fixed = nan(Mk(ind_k),N);
    GBS_sel_hardware = nan(Mk(ind_k),N*S);
    GBS_sel_fictitious = nan(Mk(ind_k),N*S);
    GBS_cof = nan(Mk(ind_k),N*T);
    for ind_n = 1:N
        for ind_m = 1:Mk(ind_k)
            match_index = Delay{ind_l}(:,1:2) == [ind_n-1,ind_m-1];
            % Propagation delay (matrix A)
            delay_mn = Delay{ind_l}(match_index(:,1)&match_index(:,2),end).';
            A(ind_m,ind_n) = (1/sqrt(N*Mk(ind_k)))*exp(-1j*delay_mn*fc*2*pi);
            % Antenna magnitude gain at UE (matrix GUE)
            AoA_az_mn = AoA_az{ind_l}(match_index(:,1)&match_index(:,2),end).';
            AoA_el_mn = AoA_el{ind_l}(match_index(:,1)&match_index(:,2),end).';
            GUE(ind_m,ind_n) = func_GetFixedAMG(RP_fixed,AoA_az_mn,AoA_el_mn);
            % Antenna magnitude gain at BS (matrix GBS_fixed, GBS_sel, and GBS_cof)
            AoD_az_mn = AoD_az{ind_l}(match_index(:,1)&match_index(:,2),end).';
            AoD_el_mn = AoD_el{ind_l}(match_index(:,1)&match_index(:,2),end).';
            % - Fixed pattern
            GBS_fixed(ind_m,ind_n) = func_GetFixedAMG(RP_fixed,AoD_az_mn,AoD_el_mn);
            % - Model I
            GBS_sel_hardware(ind_m,(ind_n-1)*S+(1:S)) = func_GetAllAMG(RP_all,AoD_az_mn,AoD_el_mn).';
            GBS_sel_fictitious(ind_m,(ind_n-1)*S+(1:S)) = func_GetAllFMG(parameters,AoD_az_mn,AoD_el_mn).';
            % - Model II
            gamma = nan(1,T);
            for r = 0:Rmax
                for m = -r:r
                    ind_xi = r^2 + r + m + 1;
                    gamma(ind_xi) = harmonicY(r,m,deg2rad(AoD_el_mn),deg2rad(AoD_az_mn),'type','real','phase',false);
                end
            end
            GBS_cof(ind_m,(ind_n-1)*T+(1:T)) = gamma;
        end
    end
    % Compute conventional channel
    H3_cv = H3_cv + sqrt(N*Mk(ind_k)/Lk(ind_k))*C.*A.*GUE.*GBS_fixed;
    H3_sel_hardware = H3_sel_hardware + sqrt(N*Mk(ind_k)/Lk(ind_k))*kron((C.*A.*GUE),ones(1,S)).*GBS_sel_hardware;
    H3_sel_fictitious = H3_sel_fictitious + sqrt(N*Mk(ind_k)/Lk(ind_k))*kron((C.*A.*GUE),ones(1,S)).*GBS_sel_fictitious;
    H3_cof = H3_cof + sqrt(N*Mk(ind_k)/Lk(ind_k))*kron((C.*A.*GUE),ones(1,T)).*GBS_cof;
end


%% Summarize necessary system parameters
% System parameters
paras.K = 3;
paras.S = S;
paras.T = T;
paras.Rmax = Rmax;
% Signal parameters
paras.fc = 30e9;                            % Center Frequency [Hz]
paras.c0 = 299792458;                       % Speed of Light [m/s]
paras.lambda = paras.c0/paras.fc;           % Wavelength
paras.Dk = 2*ones(paras.K,1);               % Number of data streams at each UE
paras.sigma2 = 10^(-90/10);                 % Noise power at each UE       
paras.Pn = 10^(0/10);                       % Power constraint per antenna
paras.NRF = sum(paras.Dk)+3;                % Number of RFCs
paras.beta = ones(paras.K,1)/paras.K;       % Weight assigned to each user's rate


%% Save the generated channels
H_cv = {H1_cv,H2_cv,H3_cv};         % Fixed-pattern antennas
H_sel_hardware = {H1_sel_hardware,H2_sel_hardware,H3_sel_hardware};             % Pattern-reconfigurable antennas (Model I)
H_sel_fictitious = {H1_sel_fictitious,H2_sel_fictitious,H3_sel_fictitious};     % Pattern-reconfigurable antennas (Model I)
H_cof = {H1_cof,H2_cof,H3_cof};     % Pattern-reconfigurable antennas (Model II)

end

%% Functions
function data = PerSheetRead(filename)
% Obtain sheet names
sheets = sheetnames(filename);
% Create a cell to store all sheet contents
data = cell(length(sheets), 1);
% Iterate sheets
for i = 1:length(sheets)
    data{i} = readmatrix(filename, 'Sheet', sheets{i});
end
end

function angle_out = wrapTo180(angle_in)
% wrapTo180  Wraps any angle to the range (-180, 180], mapping -180 to 180.
%
%   angle_out = wrapTo180(angle_in) takes an input angle_in (in degrees)
%   and returns angle_out wrapped to the range (-180, 180], including 180.
%
%   Example:
%       wrapTo180(190)              % returns -170
%       wrapTo180(-190)             % returns 170
%       wrapTo180(-180)             % returns 180
%       wrapTo180([360 -270 45])    % returns [0 90 45]
    angle_out = mod(angle_in + 180, 360) - 180;
    angle_out(angle_out == -180) = 180;
end

