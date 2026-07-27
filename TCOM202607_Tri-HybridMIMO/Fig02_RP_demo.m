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
% This script generates Fig. 2 in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all

%% Read data
file_path = fullfile('data_measurement', 'measurement.xlsx'); 
data = readmatrix(file_path);


%% Process data 
% Convert degree to rad
azimuth = data(:,1);
phi = deg2rad(azimuth);
% Convert antenna gain to maganitude gain in scale
data_sim = data(:,2:4);
data_mea = data(:,6:8);
data_sim = 10.^(data_sim/20);
data_mea = 10.^(data_mea/20);

pHFSS_1 = data_sim(:,3); 
pHFSS_2 = data_sim(:,2); 
pHFSS_3 = data_sim(:,1); 
pMEA_1 = data_mea(:,1); 
pMEA_2 = data_mea(:,2); 
pMEA_3 = data_mea(:,3);


%% Plot radiation patterns
figure;
% font: Times New Roman, size 14
p = pHFSS_1;
polarplot(phi, p,'--', 'Color', '#4169E1' ,'LineWidth', 1.5); hold on
p = pHFSS_2;
polarplot(phi, p,'--', 'Color', '#FF9900' , 'LineWidth', 1.5); hold on
p = pHFSS_3;
polarplot(phi, p,'--', 'Color', '#228B22' , 'LineWidth', 1.5); hold on
p = pMEA_1;
polarplot(phi, p,'-', 'Color', '#4169E1' , 'LineWidth', 1.5); hold on
p = pMEA_2;
polarplot(phi, p,'-', 'Color', '#FF9900' , 'LineWidth', 1.5); hold on
p = pMEA_3;
polarplot(phi, p,'-', 'Color', '#228B22' , 'LineWidth', 1.5); hold on
ax = gca; 
ax.ThetaZeroLocation = 'top';
ax.RTickLabel = [];
grid on;
legend({'HFSS (state 1)', 'HFSS (state 2)', 'HFSS (state 3)', ...
    'Measured (state 1)', 'Measured (state 2)', 'Measured (state 3)'}, ...
    'NumColumns', 2,'Interpreter','latex');