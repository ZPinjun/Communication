%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Pinjun Zheng                          
% Last Modified: July, 2026
%
% If you use this code or any (modified) part of it in any publication,
% please cite the paper: 
%
% P. Zheng, Y. Zhang, T. Y. Al-Naffouri, M. J. Hossain and A. Chaaban,
% "Tri-Hybrid Multi-User Precoding Using Pattern-Reconfigurable Antennas:
% Fundamental Models and Practical Algorithms," in IEEE Transactions on
% Communications, doi: 10.1109/TCOMM.2026.3714350.  
%
% Contact email: pinjun.zheng@ubc.ca; pinjun.zheng@alumni.kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates Fig. 7 in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all

seed = 1;

%% Obtain the optimized radiation pattern
load('result_sim00.mat')


%% Obtain multipaths AoDs
ind_anten = 90;     % Choose an antenna
num_users = 3;
Npath = [11,6,4];
% Read ray-tracing data
az_path1 = fullfile('data_RayTracing', 'UE1_AoD_az.xlsx');
el_path1 = fullfile('data_RayTracing', 'UE1_AoD_el.xlsx');
a_path1 = fullfile('data_RayTracing', 'UE1_a.xlsx');
az_path2 = fullfile('data_RayTracing', 'UE2_AoD_az.xlsx');
el_path2 = fullfile('data_RayTracing', 'UE2_AoD_el.xlsx');
a_path2 = fullfile('data_RayTracing', 'UE2_a.xlsx');
az_path3 = fullfile('data_RayTracing', 'UE3_AoD_az.xlsx');
el_path3 = fullfile('data_RayTracing', 'UE3_AoD_el.xlsx');
a_path3 = fullfile('data_RayTracing', 'UE3_a.xlsx');
az_path = {az_path1,az_path2,az_path3};
el_path = {el_path1,el_path2,el_path3};
a_path = {a_path1,a_path2,a_path3};
az_all = cell(3,1);
el_all = cell(3,1);
loss_dB_all = cell(3,1);
vec_minPL = nan(3,1);
for user = 1:num_users
    % AoD and path loss data
    az = nan(1, Npath(user));         % AoD azimuth
    el = nan(1, Npath(user));         % AoD elevation
    loss_dB = nan(1, Npath(user));    % path loss dB
    
    AoD_az = func_PerSheetRead(az_path{user});
    AoD_el = func_PerSheetRead(el_path{user});
    PL_a = func_PerSheetRead(a_path{user});
    for ind_l = 1:Npath(user)
        az(ind_l) = rad2deg(AoD_az{ind_l}(ind_anten,3))-90;
        el(ind_l) = rad2deg(AoD_el{ind_l}(ind_anten,3));
        loss_dB(ind_l) = 10*log10(abs(PL_a{1}(ind_l,2)));
    end

    az_all{user} = az;
    el_all{user} = el;
    loss_dB_all{user} = loss_dB;
    vec_minPL(user) = min(loss_dB);
end


%% Calculate beam-steering vectors
phi_all = linspace(-90,90,100);
theta_test = linspace(94,121,10);
FD1_CV = FRF_CV*FBB_CV(:,1:2);
FD2_CV = FRF_CV*FBB_CV(:,sum(paras.Dk(1))+(1:2));
FD3_CV = FRF_CV*FBB_CV(:,sum(paras.Dk(1:2))+(1:2));
FD1_MI = FRF_MI_hw*FBB_MI_hw(:,1:2);
FD2_MI = FRF_MI_hw*FBB_MI_hw(:,sum(paras.Dk(1))+(1:2));
FD3_MI = FRF_MI_hw*FBB_MI_hw(:,sum(paras.Dk(1:2))+(1:2));
FD1_MII = FRF_MII*FBB_MII(:,1:2);
FD2_MII = FRF_MII*FBB_MII(:,sum(paras.Dk(1))+(1:2));
FD3_MII = FRF_MII*FBB_MII(:,sum(paras.Dk(1:2))+(1:2));
E1_CV = nan(length(theta_test),length(phi_all));
E2_CV = nan(length(theta_test),length(phi_all));
E3_CV = nan(length(theta_test),length(phi_all));
E1_MI = nan(length(theta_test),length(phi_all));
E2_MI = nan(length(theta_test),length(phi_all));
E3_MI = nan(length(theta_test),length(phi_all));
E1_MII = nan(length(theta_test),length(phi_all));
E2_MII = nan(length(theta_test),length(phi_all));
E3_MII = nan(length(theta_test),length(phi_all));
sel_vec = zeros(paras.S,1);
sel_vec(19) = 1;
FCV = kron(eye(size(FRF_CV,1)),sel_vec);
for ind_theta = 1:length(theta_test)

    disp(['Progress: ', num2str(100*ind_theta/length(theta_test)) '%']);

    theta_all = theta_test(ind_theta)*ones(size(phi_all));
    H_out_CV = func_GetBeamSteering(phi_all,theta_all,'Model I',FCV,paras);
    H_out_MI = func_GetBeamSteering(phi_all,theta_all,'Model I',Fsel_hw,paras);
    H_out_MII = func_GetBeamSteering(phi_all,theta_all,'Model II',Fcof,paras);
    E1_CV(ind_theta,:) = sqrt(sum(abs((H_out_CV*FD1_CV).').^2));
    E2_CV(ind_theta,:) = sqrt(sum(abs((H_out_CV*FD2_CV).').^2));
    E3_CV(ind_theta,:) = sqrt(sum(abs((H_out_CV*FD3_CV).').^2));
    E1_MI(ind_theta,:) = sqrt(sum(abs((H_out_MI*FD1_MI).').^2));
    E2_MI(ind_theta,:) = sqrt(sum(abs((H_out_MI*FD2_MI).').^2));
    E3_MI(ind_theta,:) = sqrt(sum(abs((H_out_MI*FD3_MI).').^2));
    E1_MII(ind_theta,:) = sqrt(sum(abs((H_out_MII*FD1_MII).').^2));
    E2_MII(ind_theta,:) = sqrt(sum(abs((H_out_MII*FD2_MII).').^2));
    E3_MII(ind_theta,:) = sqrt(sum(abs((H_out_MII*FD3_MII).').^2));
end

save result_sim03B.mat

figure
subplot(3,1,1)
user_colors = [
    65 105 225;     % #4169E1 - Royal Blue
    255 153 0;      % #FF9900 - Orange
    34 139 34       % #228B22 - Forest Green
] / 255;
y = 10*log10(max(abs(E1_CV)));
y = y - max(y);
h1 = plot(phi_all,y,':','Color',user_colors(1,:), 'LineWidth', 2); hold on
y = 10*log10(max(abs(E2_CV)));
y = y - max(y);
h1 = [h1, plot(phi_all,y,':','Color',user_colors(2,:), 'LineWidth', 2)]; hold on
y = 10*log10(max(abs(E3_CV)));
y = y - max(y);
h1 = [h1, plot(phi_all,y,':','Color',user_colors(3,:), 'LineWidth', 2)]; hold on
grid on; box on
% Plot multipath
h2 = []; 
minPL = min(vec_minPL);
for ind_k = 1:num_users
    az_user = az_all{ind_k};
    PL_a = func_PerSheetRead(a_path{ind_k});
    for ind_l = 1:length(az_user)
        x = az_user(ind_l); 
        scale = 10*log10(abs(PL_a{1}(ind_l,2))) - minPL + 1; 
        y_low = -17; 
        y_high = y_low + scale;           
        h2 = [h2, plot([x x], [y_low y_high], 'Color', user_colors(ind_k, :), 'LineWidth', 1)];
    end
end
xlim([-90,90])
ylim([-17,5])
xlabel('$\phi$ ($^\circ$)','Interpreter','latex');
ylabel('Normalized Magnitude Gain (dB)','Interpreter','latex');
legend([h1(1), h1(2), h1(3) h2(1), h2(12), h2(18)], ...
    {'User 1 beampattern', 'User 2 beampattern', 'User 3 beampattern', ...
    'User 1 multipath', 'User 2 multipath', 'User 3 multipath'});
title('Array Beampattern Optimized via WMMSE Using Fixed Antennas')


subplot(3,1,2)
user_colors = [
    65 105 225;     % #4169E1 - Royal Blue
    255 153 0;      % #FF9900 - Orange
    34 139 34       % #228B22 - Forest Green
] / 255;
y = 10*log10(max(abs(E1_MI)));
y = y - max(y);
plot(phi_all,y,':','Color',user_colors(1,:), 'LineWidth', 2); hold on
y = 10*log10(max(abs(E2_MI)));
y = y - max(y);
plot(phi_all,y,':','Color',user_colors(2,:), 'LineWidth', 2); hold on
y = 10*log10(max(abs(E3_MI)));
y = y - max(y);
plot(phi_all,y,':','Color',user_colors(3,:), 'LineWidth', 2); hold on
grid on; box on
% Plot multipath
h2 = []; 
minPL = min(vec_minPL);
for ind_k = 1:num_users
    az_user = az_all{ind_k};
    PL_a = func_PerSheetRead(a_path{ind_k});
    for ind_l = 1:length(az_user)
        x = az_user(ind_l); 
        scale = 10*log10(abs(PL_a{1}(ind_l,2))) - minPL + 1; 
        y_low = -17; 
        y_high = y_low + scale;           
        h2 = [h2, plot([x x], [y_low y_high], 'Color', user_colors(ind_k, :), 'LineWidth', 1)];
    end
end
xlim([-90,90])
ylim([-17,5])
xlabel('$\phi$ ($^\circ$)','Interpreter','latex');
ylabel('Normalized Magnitude Gain (dB)','Interpreter','latex');
title('Array Beampattern Optimized via Model I Using Reconfigurable Antennas')


subplot(3,1,3)
user_colors = [
    65 105 225;     % #4169E1 - Royal Blue
    255 153 0;      % #FF9900 - Orange
    34 139 34       % #228B22 - Forest Green
] / 255;
y = 10*log10(max(abs(E1_MII)));
y = y - max(y);
plot(phi_all,y,':','Color',user_colors(1,:), 'LineWidth', 2); hold on
y = 10*log10(max(abs(E2_MII)));
y = y - max(y);
plot(phi_all,y,':','Color',user_colors(2,:), 'LineWidth', 2); hold on
y = 10*log10(max(abs(E3_MII)));
y = y - max(y);
plot(phi_all,y,':','Color',user_colors(3,:), 'LineWidth', 2); hold on
grid on; box on
% Plot multipath
h2 = []; 
minPL = min(vec_minPL);
for ind_k = 1:num_users
    az_user = az_all{ind_k};
    PL_a = func_PerSheetRead(a_path{ind_k});
    for ind_l = 1:length(az_user)
        x = az_user(ind_l); 
        scale = 10*log10(abs(PL_a{1}(ind_l,2))) - minPL + 1; 
        y_low = -17; 
        y_high = y_low + scale;           
        h2 = [h2, plot([x x], [y_low y_high], 'Color', user_colors(ind_k, :), 'LineWidth', 1)];
    end
end
xlim([-90,90])
ylim([-17,5])
xlabel('$\phi$ ($^\circ$)','Interpreter','latex');
ylabel('Normalized Magnitude Gain (dB)','Interpreter','latex');
title('Array Beampattern Optimized via Model II Using Reconfigurable Antennas')



