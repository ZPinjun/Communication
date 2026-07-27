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
% This script generates Fig. 6 in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all

seed = 1;

%% Load the optimized patterns
load('result_sim00.mat')
ind_anten = 35;     % Choose an antenna


%% Visualize multipaths
user_colors = [
    65 105 225;     % #4169E1 - Royal Blue
    255 153 0;      % #FF9900 - Orange
    34 139 34       % #228B22 - Forest Green
] / 255;
num_users = 3;
Npath = [11,6,4];

figure
centerX = 450;
centerY = 450;
width = 1000;
height = 400;
set(gcf,'position',[centerX, centerY,width, height])
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'tight');
view_az = 45;
view_el = 45;

nexttile
hold on;
axis equal;
grid on;
box on;
xlabel('$x$','Interpreter','latex', 'FontSize', 13);
ylabel('$y$','Interpreter','latex', 'FontSize', 13);
zlabel('$z$','Interpreter','latex', 'FontSize', 13);
title('\textbf{Multipath Departure Directions}','\textbf{at A Selected BS Antenna}','Interpreter','latex', 'FontSize', 13);
view(view_az, view_el)
origin = [0, 0, 0];
legend_handles = gobjects(num_users, 1);
legend_labels = {'User 1', 'User 2', 'User 3'};
for user = 1:num_users
    legend_handles(user) = plot3(NaN, NaN, NaN, '-', ...
        'Color', user_colors(user, :), 'LineWidth', 2);
end
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

% Plot arrows
minPL = min(vec_minPL);
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

    for i = 1:Npath(user)
        az_rad = deg2rad(az(i));
        el_rad = deg2rad(el(i));
        x = sin(el_rad) * cos(az_rad);
        y = sin(el_rad) * sin(az_rad);
        z = cos(el_rad);
        dir = [x, y, z];
        scale = loss_dB(i) - minPL + 1; 
        vec = dir * scale;
        quiver3(origin(1), origin(2), origin(3), ...
                vec(1), vec(2), vec(3), 0, ...
                'Color', user_colors(user, :), ...
                'LineWidth', 1.5, ...
                'MaxHeadSize', 0.3);
    end
end
legend(legend_handles, legend_labels, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 13);
xticklabels([]);
yticklabels([]);
zticklabels([]);


%% Visualize the formed radiation pattern based on Model I
nexttile
PlotModelI(Fsel_hw,paras,ind_anten,view_az,view_el)


%% Visualize the formed radiation pattern based on Model II
nexttile
PlotModelII(Fcof,paras,ind_anten,view_az,view_el)


%% Add subfig labels
yPos = 0.05;
annotation('textbox', [0.17, yPos, 0.05, 0.05], 'String', '(a)', ...
    'EdgeColor', 'none', 'FontSize', 15, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
annotation('textbox', [0.475, yPos, 0.05, 0.05], 'String', '(b)', ...
    'EdgeColor', 'none', 'FontSize', 15, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
annotation('textbox', [0.8, yPos, 0.05, 0.05], 'String', '(c)', ...
    'EdgeColor', 'none', 'FontSize', 15, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');


%% Functions
function PlotModelI(Fsel,paras,ind_anten,view_az,view_el)

S = paras.S;
N = 64;     % sampling number
theta_range = linspace(0, pi, N);
phi_range = linspace(-pi, pi, 2*N);
[theta, phi] = meshgrid(theta_range, phi_range);    % theta in [0,pi], phi in [0,2pi]

% Pick a radiation pattern
load('data_patterns.mat')
% Convert to magnitude (or field) pattern 
RP_magnitude = dataMatrix;
RP_magnitude(:,3:end) = nan;
for ind = 1:64
    RP_magnitude(:,ind+2) = 10.^(dataMatrix(:,ind+2)/20);
end
b_vec = Fsel((1:S)+(ind_anten-1)*S,ind_anten);
ind = find(b_vec==1);
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

[X, Y, Z] = sph2cart(phi, pi/2-theta,abs(E));         % Convert to Cartesian coordinate
surf(X, Y, Z, E, 'EdgeColor', 'none');
title({'\textbf{Selected Optimal Pattern}','\textbf{of the Antenna via Model I}'},'Interpreter','latex','FontSize',13);
view(view_az, view_el)
shading interp;     % Smooth shading to remove blocky colors
colormap(gca, 'jet');    
clim([0 3]);    
% Add lighting to improve 3D perception
light('Position', [1 1 1], 'Style', 'infinite'); % Light source
lighting phong;     % Use Phong lighting for smooth shading
material dull;      % Adjust material properties for a natural look
% Labels and title
xlabel('$x$','Interpreter','latex', 'FontSize', 13);
ylabel('$y$','Interpreter','latex', 'FontSize', 13);
zlabel({'$\quad$','$\quad$','$z$'},'Interpreter','latex', 'FontSize', 13);
yticks([-1,0,1]);
axis equal; grid on; box on;
ax = gca;  
set(ax, 'TickLabelInterpreter', 'latex');
set(ax, 'FontSize', 13); 
end


function PlotModelII(Fcof,paras,ind_anten,view_az,view_el)
T = paras.T;
LT = paras.Rmax;
N = 64;     % sampling number
theta_range = linspace(0, pi, N);
phi_range = linspace(-pi, pi, 2*N);
[theta, phi] = meshgrid(theta_range, phi_range);    % theta in [0,pi], phi in [0,2pi]

c_vec = Fcof((1:T)+(ind_anten-1)*T,ind_anten);
F_reconstructed = zeros(length(phi_range),length(theta_range));
ind = 1;
for l = 0:LT
    for m = -l:l
        Y_lm = harmonicY(l, m, theta, phi,'type','real','phase',false);
        F_reconstructed = F_reconstructed + c_vec(ind) * Y_lm;
        ind = ind + 1;
    end
end
% Visualize
[X, Y, Z] = sph2cart(phi, pi/2-theta,abs(F_reconstructed)); % Convert to Cartesian coordinate
surf(X, Y, Z, F_reconstructed, 'EdgeColor', 'none');
title({'\textbf{Arbitrarily Optimized Pattern}','\textbf{of the Antenna via Model II}'}','Interpreter','latex','FontSize',13);
view(view_az, view_el)
shading interp;    % Smooth shading to remove blocky colors
colormap(gca, 'jet');
clim([0 3]);
% Add lighting to improve 3D perception
light('Position', [1 1 1], 'Style', 'infinite'); % Light source
lighting phong;    % Use Phong lighting for smooth shading
material dull;     % Adjust material properties for a natural look
% Labels and title
xlabel('$x$','Interpreter','latex', 'FontSize', 13);
ylabel('$y$','Interpreter','latex', 'FontSize', 13);
zlabel({'$\quad$','$\quad$','$z$'},'Interpreter','latex', 'FontSize', 13);
axis equal; grid on; box on;
ax = gca;  
set(ax, 'TickLabelInterpreter', 'latex');
set(ax, 'FontSize', 13); 
end
