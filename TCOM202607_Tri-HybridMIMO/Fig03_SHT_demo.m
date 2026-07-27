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
% This script generates Fig. 3 in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all

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
ind = 36;
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


%% Plot the original 3D radiation pattern
figure
centerX = 450;
centerY = 450;
width = 1500;
height = 400;
set(gcf,'position',[centerX, centerY,width, height])
tiledlayout(1, 5, 'Padding', 'compact', 'TileSpacing', 'compact');
%
view_az = 30;
view_el = 45;
[X, Y, Z] = sph2cart(phi, pi/2-theta,abs(E));         % Convert to Cartesian coordinate
nexttile;
surf(X, Y, Z, E, 'EdgeColor', 'none');
title('\textbf{Original Pattern}','Interpreter','latex','FontSize',13);
view(view_az, view_el)
shading interp;     % Smooth shading to remove blocky colors
colormap(gca, 'jet');    
clim([0 2.5]);    
% Add lighting to improve 3D perception
light('Position', [1 1 1], 'Style', 'infinite'); % Light source
lighting phong;     % Use Phong lighting for smooth shading
material dull;      % Adjust material properties for a natural look
% Labels and title
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
yticks([-1,0,1]);
axis equal;         % Keep aspect ratio equal
grid on;         
box on;
ax = gca;  
set(ax, 'TickLabelInterpreter', 'latex');
set(ax, 'FontSize', 13); 


%% Perform spherical harmonic decomposition
Lmax = 10;                      % Highest degree
alm = zeros(Lmax+1, 2*Lmax+1);  % Coefficent matrix
for l = 0:Lmax
    for m = -l:l
        % Compute harmonics coefficents
        Y_lm = harmonicY(l, m, theta, phi,'type','real','phase',false);
        % Compute projection coefficient a_lm (inner product)
        alm(l+1, m+Lmax+1) = sum(sum(E .* conj(Y_lm) .* sin(theta))) * (pi/N) * (2*pi/(2*N));
    end
end


%% Visualize harmonic coefficents
nexttile;
b = bar3(alm); 
for k = 1:length(b)
    ZData = get(b(k), 'ZData'); 
    set(b(k), 'CData', ZData, 'FaceColor', 'interp'); 
end
title('\textbf{Harmonic Coefficients}','Interpreter','latex','FontSize',13);
xlabel('Order $q$','Interpreter','latex');
ylabel('Degree $u$','Interpreter','latex');
colormap(gca,flipud(hot)); clim([0 2.5]); %colorbar;
xticks([1, 6, 11, 16, 21]);
xticklabels({'-10', '-5', '0', '5', '10'});
yticks([1, 6, 11]);
yticklabels({'0', '5', '10'});
ax = gca;  
set(ax, 'TickLabelInterpreter', 'latex');
set(ax, 'FontSize', 13); 
box on


%% Reconstruct radiation pattern
%
LT = 2;
Alm = alm(1:LT+1,:);
C1 = norm(Alm(1,:),2)^2;
C2 = 4*pi - C1;
Alm(2:end,:) = sqrt(C2)*Alm(2:end,:)/norm(Alm(2:end,:),'fro');
F_reconstructed = zeros(size(E));
for l = 0:LT
    for m = -l:l
        Y_lm = harmonicY(l, m, theta, phi,'type','real','phase',false);
        F_reconstructed = F_reconstructed + Alm(l+1, m+Lmax+1) * Y_lm;
    end
end
% Visualize
view_az = 30;
view_el = 45;
[X, Y, Z] = sph2cart(phi, pi/2-theta,abs(F_reconstructed)); % Convert to Cartesian coordinate
nexttile;
surf(X, Y, Z, F_reconstructed, 'EdgeColor', 'none');
title({'\textbf{Truncated Reconstruction}','\textbf{($U=2$)}'},'Interpreter','latex','FontSize',13);
view(view_az, view_el)
shading interp;    % Smooth shading to remove blocky colors
colormap(gca, 'jet');     
clim([0 2.5]);    
% Add lighting to improve 3D perception
light('Position', [1 1 1], 'Style', 'infinite'); % Light source
lighting phong;    % Use Phong lighting for smooth shading
material dull;     % Adjust material properties for a natural look
% Labels and title
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
axis equal; grid on; box on;
ax = gca;  
set(ax, 'TickLabelInterpreter', 'latex');
set(ax, 'FontSize', 13); 
%
LT = 4;
Alm = alm(1:LT+1,:);
C1 = norm(Alm(1,:),2)^2;
C2 = 4*pi - C1;
Alm(2:end,:) = sqrt(C2)*Alm(2:end,:)/norm(Alm(2:end,:),'fro');
F_reconstructed = zeros(size(E));
for l = 0:LT
    for m = -l:l
        Y_lm = harmonicY(l, m, theta, phi,'type','real','phase',false);
        F_reconstructed = F_reconstructed + Alm(l+1, m+Lmax+1) * Y_lm;
    end
end
% Visualize
view_az = 30;
view_el = 45;
[X, Y, Z] = sph2cart(phi, pi/2-theta,abs(F_reconstructed)); % Convert to Cartesian coordinate
nexttile;
surf(X, Y, Z, F_reconstructed, 'EdgeColor', 'none');
title({'\textbf{Truncated Reconstruction}','\textbf{($U=4$)}'},'Interpreter','latex','FontSize',13);
view(view_az, view_el)
shading interp;    % Smooth shading to remove blocky colors
colormap(gca, 'jet');
clim([0 2.5]);    
% Add lighting to improve 3D perception
light('Position', [1 1 1], 'Style', 'infinite'); % Light source
lighting phong;    % Use Phong lighting for smooth shading
material dull;     % Adjust material properties for a natural look
% Labels and title
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
axis equal; grid on; box on;
ax = gca;  
set(ax, 'TickLabelInterpreter', 'latex');
set(ax, 'FontSize', 13); 
%
LT = 8;
Alm = alm(1:LT+1,:);
C1 = norm(Alm(1,:),2)^2;
C2 = 4*pi - C1;
Alm(2:end,:) = sqrt(C2)*Alm(2:end,:)/norm(Alm(2:end,:),'fro');
F_reconstructed = zeros(size(E));
for l = 0:LT
    for m = -l:l
        Y_lm = harmonicY(l, m, theta, phi,'type','real','phase',false);
        F_reconstructed = F_reconstructed + Alm(l+1, m+Lmax+1) * Y_lm;
    end
end
% Visualize
view_az = 30;
view_el = 45;
[X, Y, Z] = sph2cart(phi, pi/2-theta,abs(F_reconstructed)); % Convert to Cartesian coordinate
nexttile;
surf(X, Y, Z, F_reconstructed, 'EdgeColor', 'none');
title({'\textbf{Truncated Reconstruction}','\textbf{($U=8$)}'},'Interpreter','latex','FontSize',13);
view(view_az, view_el)
shading interp;    % Smooth shading to remove blocky colors
colormap(gca, 'jet');
clim([0 2.5]);    
% Add lighting to improve 3D perception
light('Position', [1 1 1], 'Style', 'infinite'); % Light source
lighting phong;    % Use Phong lighting for smooth shading
material dull;     % Adjust material properties for a natural look
% Labels and title
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
axis equal; grid on; box on;
ax = gca;  
set(ax, 'TickLabelInterpreter', 'latex');
set(ax, 'FontSize', 13); 


%% Add subfig labels
yPos = 0.12;
annotation('textbox', [0.07, yPos, 0.05, 0.05], 'String', '(a)', ...
    'EdgeColor', 'none', 'FontSize', 17, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
annotation('textbox', [0.28, yPos, 0.05, 0.05], 'String', '(b)', ...
    'EdgeColor', 'none', 'FontSize', 17, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
annotation('textbox', [0.475, yPos, 0.05, 0.05], 'String', '(c)', ...
    'EdgeColor', 'none', 'FontSize', 17, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
annotation('textbox', [0.67, yPos, 0.05, 0.05], 'String', '(d)', ...
    'EdgeColor', 'none', 'FontSize', 17, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
annotation('textbox', [0.875, yPos, 0.05, 0.05], 'String', '(e)', ...
    'EdgeColor', 'none', 'FontSize', 17, 'HorizontalAlignment', 'center', 'Interpreter', 'latex');