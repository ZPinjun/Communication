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
% This script generates Fig. 8 in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all

seed = 1;

Pn_test = linspace(-20,10,6).';
WSR_CV_FD_all = nan(length(Pn_test),1);
WSR_CV_all = nan(length(Pn_test),1);
WSR_ZF_FD_all = nan(length(Pn_test),1);
WSR_ZF_all = nan(length(Pn_test),1);
WSR_MI_FD_hw_all = nan(length(Pn_test),1);
WSR_MI_hw_all = nan(length(Pn_test),1);
WSR_MI_FD_ft_all = nan(length(Pn_test),1);
WSR_MI_ft_all = nan(length(Pn_test),1);
WSR_MII_FD_all = nan(length(Pn_test),1);
WSR_MII_all = nan(length(Pn_test),1);

h = waitbar(0, 'Iterating over various per-antenna power budgets...');  
for ind_Pn = 1:length(Pn_test)

    % Generate channels
    load('data_channels.mat')
    paras.Pn = 10^(Pn_test(ind_Pn)/10);

    % Optimization using pattern fixed antennas (Classical WMMSE optimizer)
    [FD_CV,MSE_seq_CV] = func_OptCV(H_CV,paras,seed);
    [FRF_CV,FBB_CV] = func_HybridDecomp(FD_CV,paras,seed);
    WSR_CV_FD = func_CalWSR(eye(size(H_CV{1},2)),eye(size(FD_CV,1)),FD_CV,H_CV,paras);
    WSR_CV = func_CalWSR(eye(size(H_CV{1},2)),FRF_CV,FBB_CV,H_CV,paras);

    % Optimization using pattern fixed antennas (ZF benchmark)
    P_total = norm(FD_CV,'fro')^2;
    [FD_ZF, Ak] = func_SumRateBD(H_CV, P_total, paras);
    [FRF_ZF,FBB_ZF] = func_HybridDecomp(FD_ZF,paras,seed);
    WSR_ZF_FD = func_CalWSR(eye(size(H_CV{1},2)),eye(size(FD_CV,1)),FD_ZF,H_CV,paras);
    WSR_ZF = func_CalWSR(eye(size(H_CV{1},2)),FRF_ZF,FBB_ZF,H_CV,paras);

    % Optimization using pattern reconfigurable antennas (Model I, hardware)
    [Fsel_hw,FD_MI_hw,MSE_seq_MI_hw] = func_OptSel(FD_CV,Heff_MI_hw,paras,seed);
    [FRF_MI_hw,FBB_MI_hw] = func_HybridDecomp(FD_MI_hw,paras,seed);
    WSR_MI_FD_hw = func_CalWSR(Fsel_hw,eye(size(FD_MI_hw,1)),FD_MI_hw,Heff_MI_hw,paras);
    WSR_MI_hw = func_CalWSR(Fsel_hw,FRF_MI_hw,FBB_MI_hw,Heff_MI_hw,paras);

    % Optimization using pattern reconfigurable antennas (Model I, fictitious)
    [Fsel_ft,FD_MI_ft,MSE_seq_MI_ft] = func_OptSel(FD_CV,Heff_MI_ft,paras,seed);
    [FRF_MI_ft,FBB_MI_ft] = func_HybridDecomp(FD_MI_ft,paras,seed);
    WSR_MI_FD_ft = func_CalWSR(Fsel_ft,eye(size(FD_MI_ft,1)),FD_MI_ft,Heff_MI_ft,paras);
    WSR_MI_ft = func_CalWSR(Fsel_ft,FRF_MI_ft,FBB_MI_ft,Heff_MI_ft,paras);

    % Optimization using pattern reconfigurable antennas (Model II, arbitrary)
    [Fcof,FD_MII,MSE_seq_MII] = func_OptCof(FD_CV,Heff_MII,paras,0.7,seed);
    [FRF_MII,FBB_MII] = func_HybridDecomp(FD_MII,paras,seed);
    WSR_MII_FD = func_CalWSR(Fcof,eye(size(FD_MII,1)),FD_MII,Heff_MII,paras);
    WSR_MII = func_CalWSR(Fcof,FRF_MII,FBB_MII,Heff_MII,paras);

    % Collect results
    WSR_CV_FD_all(ind_Pn) = WSR_CV_FD;
    WSR_CV_all(ind_Pn) = WSR_CV;
    WSR_ZF_FD_all(ind_Pn) = WSR_ZF_FD;
    WSR_ZF_all(ind_Pn) = WSR_ZF;
    WSR_MI_FD_hw_all(ind_Pn) = WSR_MI_FD_hw;
    WSR_MI_hw_all(ind_Pn) = WSR_MI_hw;
    WSR_MI_FD_ft_all(ind_Pn) = WSR_MI_FD_ft;
    WSR_MI_ft_all(ind_Pn) = WSR_MI_ft;
    WSR_MII_FD_all(ind_Pn) = WSR_MII_FD;
    WSR_MII_all(ind_Pn) = WSR_MII;
    
    waitbar(ind_Pn/length(Pn_test), h, sprintf('Progress：%d%%', round(100*ind_Pn/length(Pn_test))));  
end
close(h);

save result_sim04.mat

colors = {[65, 105, 225]/255, ...   % Royal Blue
          [255, 127, 14]/255, ...   % Burnt Orange
          [34, 139, 34]/255, ...    % Forest Green
          [153, 50, 204]/255, ...   % Dark Orchid
          [218, 165, 32]/255};      % Goldenrod
figure
plot(Pn_test,WSR_MII_FD_all,'--^','Color',colors{3},'LineWidth',1.5); hold on
plot(Pn_test,WSR_MI_FD_ft_all,'--s','Color',colors{5},'LineWidth',1.5); hold on
plot(Pn_test,WSR_MI_FD_hw_all,'--o','Color',colors{2},'LineWidth',1.5); hold on
plot(Pn_test,WSR_CV_FD_all,'--d','Color',colors{1},'LineWidth',1.5); hold on
plot(Pn_test,WSR_ZF_FD_all,'-->','Color',colors{4},'LineWidth',1.5); hold on

plot(Pn_test,WSR_MII_all,'-^','Color',colors{3},'LineWidth',1.5); hold on
plot(Pn_test,WSR_MI_ft_all,'-s','Color',colors{5},'LineWidth',1.5); hold on
plot(Pn_test,WSR_MI_hw_all,'-o','Color',colors{2},'LineWidth',1.5); hold on
plot(Pn_test,WSR_CV_all,'-d','Color',colors{1},'LineWidth',1.5); hold on
plot(Pn_test,WSR_ZF_all,'->','Color',colors{4},'LineWidth',1.5); hold on
grid on; box on
xlabel('Per-Antenna Power Budget (dBm)');
ylabel('Weighted Sum-Rate (bps/Hz)');
ylim([0,21])
legend({'$\mathbf{F}_\mathrm{D}^\mathrm{opt},\mathbf{F}_\mathrm{cof}^\mathrm{opt}$ (Model II)', ...
    '$\mathbf{F}_\mathrm{D}^\mathrm{opt},\mathbf{F}_\mathrm{sel}^\mathrm{opt}$ (Model I, fiticious)', ...
    '$\mathbf{F}_\mathrm{D}^\mathrm{opt},\mathbf{F}_\mathrm{sel}^\mathrm{opt}$ (Model I, hardware)', ...
    '$\mathbf{F}_\mathrm{D}^\mathrm{opt}$ (WMMSE) + fixed antennas', ...
    '$\mathbf{F}_\mathrm{D}^\mathrm{opt}$ (ZF) + fixed antennas', ...
    '$\mathbf{F}_\mathrm{BB}^\mathrm{opt},\mathbf{F}_\mathrm{RF}^\mathrm{opt},\mathbf{F}_\mathrm{cof}^\mathrm{opt}$ (Model II)', ...
    '$\mathbf{F}_\mathrm{BB}^\mathrm{opt},\mathbf{F}_\mathrm{RF}^\mathrm{opt},\mathbf{F}_\mathrm{sel}^\mathrm{opt}$ (Model I, fiticious)', ...
    '$\mathbf{F}_\mathrm{BB}^\mathrm{opt},\mathbf{F}_\mathrm{RF}^\mathrm{opt},\mathbf{F}_\mathrm{sel}^\mathrm{opt}$ (Model I, hardware)', ...
    '$\mathbf{F}_\mathrm{BB}^\mathrm{opt},\mathbf{F}_\mathrm{RF}^\mathrm{opt}$ (WMMSE) + fixed antennas', ...
    '$\mathbf{F}_\mathrm{BB}^\mathrm{opt},\mathbf{F}_\mathrm{RF}^\mathrm{opt}$ (ZF) + fixed antennas'}, ...
       'NumColumns', 2, 'Location', 'best','Interpreter','latex');