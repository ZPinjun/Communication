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
% This script shows a demo of the core algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all

seed = 1;

%% Generate channels
% The following function generates wireless channels based on ray-tracing results
% [H_CV,Heff_MI_hw,Heff_MI_ft,Heff_MII,paras] = func_GetChannels(100);
% save data_channels.mat

% We load the generated channels directly 
load('data_channels.mat')

%% Optimization using pattern fixed antennas (Classic WMMSE optimizer)
[FD_CV,MSE_seq_CV] = func_OptCV(H_CV,paras,seed);
[FRF_CV,FBB_CV] = func_HybridDecomp(FD_CV,paras,seed);
[WSR_CV_FD,rates_CV_FD] = func_CalWSR(eye(size(H_CV{1},2)),eye(size(FD_CV,1)),FD_CV,H_CV,paras);
[WSR_CV,rates_CV] = func_CalWSR(eye(size(H_CV{1},2)),FRF_CV,FBB_CV,H_CV,paras);

%% Optimization using pattern fixed antennas (ZF benchmark)
P_total = norm(FD_CV,'fro')^2;
[FD_ZF, Ak] = func_SumRateBD(H_CV, P_total, paras);
[FRF_ZF,FBB_ZF] = func_HybridDecomp(FD_ZF,paras,seed);
[WSR_ZF_FD,rates_ZF_FD] = func_CalWSR(eye(size(H_CV{1},2)),eye(size(FD_CV,1)),FD_ZF,H_CV,paras);
[WSR_ZF,rates_ZF] = func_CalWSR(eye(size(H_CV{1},2)),FRF_ZF,FBB_ZF,H_CV,paras);

%% Optimization using pattern reconfigurable antennas (Model I, hardware) 
[Fsel_hw,FD_MI_hw,MSE_seq_MI_hw] = func_OptSel(FD_CV,Heff_MI_hw,paras,seed);
[FRF_MI_hw,FBB_MI_hw] = func_HybridDecomp(FD_MI_hw,paras,seed);
[WSR_MI_FD_hw,rates_MI_FD_hw] = func_CalWSR(Fsel_hw,eye(size(FD_MI_hw,1)),FD_MI_hw,Heff_MI_hw,paras);
[WSR_MI_hw,rates_MI_hw] = func_CalWSR(Fsel_hw,FRF_MI_hw,FBB_MI_hw,Heff_MI_hw,paras);

%% Optimization using pattern reconfigurable antennas (Model I, fictitious) 
[Fsel_ft,FD_MI_ft,MSE_seq_MI_ft] = func_OptSel(FD_CV,Heff_MI_ft,paras,seed);
[FRF_MI_ft,FBB_MI_ft] = func_HybridDecomp(FD_MI_ft,paras,seed);
[WSR_MI_FD_ft,rates_MI_FD_ft] = func_CalWSR(Fsel_ft,eye(size(FD_MI_ft,1)),FD_MI_ft,Heff_MI_ft,paras);
[WSR_MI_ft,rates_MI_ft] = func_CalWSR(Fsel_ft,FRF_MI_ft,FBB_MI_ft,Heff_MI_ft,paras);

%% Optimization using pattern reconfigurable antennas (Model II, arbitrary) 
[Fcof,FD_MII,MSE_seq_MII] = func_OptCof(FD_CV,Heff_MII,paras,0.7,seed);
[FRF_MII,FBB_MII] = func_HybridDecomp(FD_MII,paras,seed);
[WSR_MII_FD,rates_MII_FD] = func_CalWSR(Fcof,eye(size(FD_MII,1)),FD_MII,Heff_MII,paras);
[WSR_MII,rates_MII] = func_CalWSR(Fcof,FRF_MII,FBB_MII,Heff_MII,paras);

save result_sim00.mat


%% Display weighted sum-rate evaluations
disp('--- Fixed-pattern antennas optimized using classical WMMSE ---')
disp(['Fully digital WSR: ',num2str(WSR_CV_FD),' bps/Hz; ','Hybrid WSR: ',num2str(WSR_CV), ' bps/Hz.']);
disp('--- Fixed-pattern antennas optimized using benchmark ZF ---')
disp(['Fully digital WSR: ',num2str(WSR_ZF_FD),' bps/Hz; ','Hybrid WSR: ',num2str(WSR_ZF), ' bps/Hz.']);
disp('--- Pattern reconfigurable antennas optimized based on Model I (real hardware patterns) ---')
disp(['Fully digital WSR: ',num2str(WSR_MI_FD_hw),' bps/Hz; ','Hybrid WSR: ',num2str(WSR_MI_hw), ' bps/Hz.']);
disp('--- Pattern reconfigurable antennas optimized based on Model I (fictitious patterns) ---')
disp(['Fully digital WSR: ',num2str(WSR_MI_FD_ft),' bps/Hz; ','Hybrid WSR: ',num2str(WSR_MI_ft), ' bps/Hz.']);
disp('--- Pattern reconfigurable antennas optimized based on Model II (arbitrary optimization) ---')
disp(['Fully digital WSR: ',num2str(WSR_MII_FD),' bps/Hz; ','Hybrid WSR: ',num2str(WSR_MII), ' bps/Hz.']);



