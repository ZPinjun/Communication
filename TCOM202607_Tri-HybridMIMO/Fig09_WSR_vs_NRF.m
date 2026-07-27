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
% This script generates Fig. 9 in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all

seed = 1;

% Generate channels
load('data_channels.mat')
paras.Pn = 10^(10/10);
% Optimization using pattern fixed antennas (Classical WMMSE optimizer)
[FD_CV,~] = func_OptCV(H_CV,paras,seed);
% Optimization using pattern reconfigurable antennas (Model I)
[Fsel,FD_MI,~] = func_OptSel(FD_CV,Heff_MI_hw,paras,seed);
% Optimization using pattern reconfigurable antennas (Model II)
[Fcof1,FD_MII1,~] = func_OptCof(FD_CV,Heff_MII,paras,0.7,seed);
[Fcof2,FD_MII2,~] = func_OptCof(FD_CV,Heff_MII,paras,0.9,seed);
NRF_test = 0:1:5;
WSR_CV_all = nan(length(NRF_test),1);
WSR_MI_all = nan(length(NRF_test),1);
WSR_MII1_all = nan(length(NRF_test),1);
WSR_MII2_all = nan(length(NRF_test),1);

save data_sim05.mat


%% We load the generated data directly
load('data_sim05.mat')
h = waitbar(0.1, 'Iterating over various number of RFCs at BS...');
for ind_NRF = 1:length(NRF_test)
    paras.NRF = sum(paras.Dk) + NRF_test(ind_NRF);
    [FRF_CV,FBB_CV] = func_HybridDecomp(FD_CV,paras,seed);
    [FRF_MI,FBB_MI] = func_HybridDecomp(FD_MI,paras,seed);
    [FRF_MII1,FBB_MII1] = func_HybridDecomp(FD_MII1,paras,seed);
    [FRF_MII2,FBB_MII2] = func_HybridDecomp(FD_MII2,paras,seed);
    WSR_CV = func_CalWSR(eye(size(H_CV{1},2)),FRF_CV,FBB_CV,H_CV,paras);
    WSR_MI = func_CalWSR(Fsel,FRF_MI,FBB_MI,Heff_MI_hw,paras);
    WSR_MII1 = func_CalWSR(Fcof1,FRF_MII1,FBB_MII1,Heff_MII,paras);
    WSR_MII2 = func_CalWSR(Fcof2,FRF_MII2,FBB_MII2,Heff_MII,paras);
    
    % Collect results
    WSR_CV_all(ind_NRF) = WSR_CV;
    WSR_MI_all(ind_NRF) = WSR_MI;
    WSR_MII1_all(ind_NRF) = WSR_MII1;
    WSR_MII2_all(ind_NRF) = WSR_MII2;

    waitbar(ind_NRF/length(NRF_test), h, sprintf('Progress：%d%%', round(100*ind_NRF/length(NRF_test))));  
end
close(h);

save result_sim05.mat

colors = {[65, 105, 225]/255, ...   % Royal Blue
          [255, 127, 14]/255, ...   % Burnt Orange
          [34, 139, 34]/255, ...    % Forest Green
          [153, 50, 204]/255};      % Dark Orchid
figure
plot(NRF_test,WSR_MII1_all,'-^','Color',colors{3},'LineWidth',1.5); hold on
plot(NRF_test,WSR_MII2_all,'-s','Color',colors{3},'LineWidth',1.5); hold on
plot(NRF_test,WSR_MI_all,'-o','Color',colors{2},'LineWidth',1.5); hold on
plot(NRF_test,WSR_CV_all,'-d','Color',colors{1},'LineWidth',1.5); hold on
grid on; box on
xlabel('Number of BS RFCs');
ylabel('Weighted Sum-Rate (bps/Hz)');
legend({'Tri-hybrid (Model II, $\rho=0.7$)', 'Tri-hybrid (Model II, $\rho=0.9$)', ...
    'Tri-hybrid (Model I, hardware)', 'Hybrid (WMMSE)'}, 'NumColumns', 1, 'Location', 'best','Interpreter','latex');

