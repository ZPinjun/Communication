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
% This script generates Fig. 5 in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all

%% User settings
seed = 1;
rho_model_ii = 0.7;
include_model_i_fictitious = false;

% Practical defaults so the script finishes in a reasonable time.
max_iter_wmmse = 150;
max_iter_model_i = 150;
max_iter_model_ii = 150;

%% Load channels and common initialization
load('data_channels.mat', 'H_CV', 'Heff_MI_hw', 'Heff_MI_ft', 'Heff_MII', 'paras');
FD_init_common = generate_common_fd_init(size(H_CV{1}, 2), sum(paras.Dk), paras.Pn, seed);

%% Run convergence traces
[FD_CV_trace, WSR_seq_CV, MSE_seq_CV] = trace_fixed_wmmse( ...
    FD_init_common, H_CV, paras, seed, max_iter_wmmse);

[Fsel_hw, FD_MI_hw, WSR_seq_MI_hw, MSE_seq_MI_hw] = trace_model_i( ...
    FD_init_common, Heff_MI_hw, paras, seed, max_iter_model_i);

[Fcof, FD_MII, WSR_seq_MII, MSE_seq_MII] = trace_model_ii( ...
    FD_init_common, Heff_MII, paras, rho_model_ii, seed, max_iter_model_ii);

if include_model_i_fictitious
    [Fsel_ft, FD_MI_ft, WSR_seq_MI_ft, MSE_seq_MI_ft] = trace_model_i( ...
        FD_init_common, Heff_MI_ft, paras, seed, max_iter_model_i);
else
    Fsel_ft = [];
    FD_MI_ft = [];
    WSR_seq_MI_ft = [];
    MSE_seq_MI_ft = [];
end

%% Save results
save('result_sim08.mat', ...
    'seed', 'rho_model_ii', 'include_model_i_fictitious', ...
    'max_iter_wmmse', 'max_iter_model_i', 'max_iter_model_ii', ...
    'FD_init_common', ...
    'FD_CV_trace', 'WSR_seq_CV', 'MSE_seq_CV', ...
    'Fsel_hw', 'FD_MI_hw', 'WSR_seq_MI_hw', 'MSE_seq_MI_hw', ...
    'Fcof', 'FD_MII', 'WSR_seq_MII', 'MSE_seq_MII', ...
    'Fsel_ft', 'FD_MI_ft', 'WSR_seq_MI_ft', 'MSE_seq_MI_ft');

%% Plot WSR vs. iteration
colors = [...
    65, 105, 225; ...   % Royal Blue
    255, 127, 14; ...   % Burnt Orange
    34, 139, 34; ...    % Forest Green
    153, 50, 204] / 255; % Dark Orchid

figure
plot(0:(length(WSR_seq_CV) - 1), WSR_seq_CV, '-d', ...
    'Color', colors(1, :), 'LineWidth', 1.8, 'MarkerSize', 5); hold on
plot(0:(length(WSR_seq_MI_hw) - 1), WSR_seq_MI_hw, '-o', ...
    'Color', colors(2, :), 'LineWidth', 1.8, 'MarkerSize', 5); hold on
plot(0:(length(WSR_seq_MII) - 1), WSR_seq_MII, '-s', ...
    'Color', colors(3, :), 'LineWidth', 1.8, 'MarkerSize', 5); hold on

legend_entries = {'Fixed antennas (WMMSE)', 'Model I (hardware)', 'Model II'};
if include_model_i_fictitious
    plot(0:(length(WSR_seq_MI_ft) - 1), WSR_seq_MI_ft, '-d', ...
        'Color', colors(4, :), 'LineWidth', 1.8, 'MarkerSize', 5); hold on
    legend_entries{end + 1} = 'Model I (fictitious)'; 
end

grid on; box on
xlabel('Outer Iteration');
ylabel('Weighted Sum-Rate (bps/Hz)');
title('Convergence of the alternating optimization algorithms');
legend(legend_entries, 'Location', 'best', 'Interpreter', 'latex');
set(gcf, 'PaperPositionMode', 'auto');


%% Local functions
function FD_init = generate_common_fd_init(N, D, Pn, seed)
rng(seed)

FD_init = randn(N, D) + 1j * randn(N, D);
for ind_n = 1:N
    FD_init(ind_n, :) = Pn * FD_init(ind_n, :) / norm(FD_init(ind_n, :), 2);
end
end


function [FD_opt, WSR_seq, MSE_seq] = trace_fixed_wmmse(FD_init, H, paras, seed, max_iter)
rng(seed)

K = paras.K;
N = size(H{1}, 2);
M = nan(1, K);
for ind_k = 1:K
    M(ind_k) = size(H{ind_k}, 1);
end
D_all = paras.Dk;
D = sum(D_all);
sigma2 = paras.sigma2;
Pn = paras.Pn;
beta = paras.beta;

FD_opt = FD_init;

W_opt = cell(K, 1);
for ind_k = 1:K
    W_opt{ind_k} = eye(D_all(ind_k));
end

MSE_seq = nan(max_iter, 1);
WSR_seq = nan(max_iter + 1, 1);
WSR_seq(1) = func_CalWSR(eye(N), eye(N), FD_opt, H, paras);

for ind_iter = 1:max_iter
    FD_part = cell(K, 1);
    for ind_k = 1:K
        FD_part{ind_k} = FD_opt(:, (sum(D_all(1:ind_k)) - D_all(ind_k) + 1):sum(D_all(1:ind_k)));
    end

    U_opt = cell(K, 1);
    for ind_k = 1:K
        Hk = H{ind_k};
        Mk = M(ind_k);
        Matrix1 = zeros(Mk, Mk);
        for ind_j = 1:K
            Matrix1 = Matrix1 + Hk * (FD_part{ind_j} * FD_part{ind_j}') * Hk';
        end
        U_opt{ind_k} = (Matrix1 + sigma2 * eye(Mk))^(-1) * Hk * FD_part{ind_k};
        W_opt{ind_k} = (eye(D_all(ind_k)) - U_opt{ind_k}' * Hk * FD_part{ind_k})^(-1);
    end

    A = cal_A_fixed_local(U_opt, W_opt, beta, K);
    B = cal_B_fixed_local(U_opt, W_opt, beta, K);
    Hconc = [];
    for ind_k = 1:K
        Hconc = [Hconc; H{ind_k}]; 
    end

    for ind_n = 1:N
        hn = Hconc(:, ind_n);
        an = hn' * A * hn;
        C = zeros(D, sum(M));
        for ind_l = 1:N
            if ind_l ~= ind_n
                hl = Hconc(:, ind_l);
                C = C + FD_opt(ind_l, :)' * hl';
            end
        end
        bn = -B * hn + C * A * hn;
        pn = -bn * min([1 / an, sqrt(Pn) / norm(bn, 2)]);
        FD_opt(ind_n, :) = pn';
    end

    MSE_seq(ind_iter) = cal_MSE_fixed_local(U_opt, W_opt, FD_opt, H, D_all, sigma2, beta);
    WSR_seq(ind_iter + 1) = func_CalWSR(eye(N), eye(N), FD_opt, H, paras);

    if ind_iter > 1 && abs(MSE_seq(ind_iter) - MSE_seq(ind_iter - 1)) < 2.5e-3
        MSE_seq = MSE_seq(1:ind_iter);
        WSR_seq = WSR_seq(1:(ind_iter + 1));
        break
    end
end

if any(isnan(MSE_seq))
    last_valid = find(~isnan(MSE_seq), 1, 'last');
    if isempty(last_valid)
        last_valid = 0;
    end
    MSE_seq = MSE_seq(1:last_valid);
    WSR_seq = WSR_seq(1:(last_valid + 1));
end
end


function [Fsel_opt, FD_opt, WSR_seq, MSE_seq] = trace_model_i(FD_init, Heff, paras, seed, max_iter)
rng(seed)

S = paras.S;
K = paras.K;
N = size(Heff{1}, 2) / S;
M = nan(1, K);
for ind_k = 1:K
    M(ind_k) = size(Heff{ind_k}, 1);
end
D_all = paras.Dk;
D = sum(D_all);
sigma2 = paras.sigma2;
Pn = paras.Pn;
beta = paras.beta;

W = cell(K, 1);
for ind_k = 1:K
    W{ind_k} = eye(D_all(ind_k));
end

B = zeros(S, N);
B(min(19, S), :) = 1;

W_opt = W;
B_opt = B;
FD_opt = FD_init;
MSE_seq = nan(max_iter, 1);
WSR_seq = nan(max_iter + 1, 1);
WSR_seq(1) = func_CalWSR(func_B2Fsel(B_opt), eye(N), FD_opt, Heff, paras);

for ind_iter = 1:max_iter
    H = Heff;
    for ind_k = 1:K
        H{ind_k} = H{ind_k} * func_B2Fsel(B_opt);
    end

    FD_part = cell(K, 1);
    for ind_k = 1:K
        FD_part{ind_k} = FD_opt(:, (sum(D_all(1:ind_k)) - D_all(ind_k) + 1):sum(D_all(1:ind_k)));
    end

    U_opt = cell(K, 1);
    for ind_k = 1:K
        Hk = H{ind_k};
        Mk = M(ind_k);
        Matrix1 = zeros(Mk, Mk);
        for ind_i = 1:K
            FDi = FD_part{ind_i};
            Matrix1 = Matrix1 + Hk * (FDi * FDi') * Hk';
        end
        U_opt{ind_k} = (Matrix1 + sigma2 * eye(Mk))^(-1) * Hk * FD_part{ind_k};
        W_opt{ind_k} = (eye(D_all(ind_k)) - U_opt{ind_k}' * Hk * FD_part{ind_k})^(-1);
    end

    BB = zeros(N * S, N * S);
    DD = [];
    for ind_k = 1:K
        BB = BB + beta(ind_k) * Heff{ind_k}' * U_opt{ind_k} * W_opt{ind_k} * U_opt{ind_k}' * Heff{ind_k};
        DD = [DD; beta(ind_k) * W_opt{ind_k} * U_opt{ind_k}' * Heff{ind_k}]; %#ok<AGROW>
    end

    for ind_n = 1:N
        Bnn = BB((1:S) + (ind_n - 1) * S, (1:S) + (ind_n - 1) * S);
        Qn = get_Qn_local(ind_n, FD_opt, B_opt, S, D, N, BB);
        Dn = DD(:, (1:S) + (ind_n - 1) * S);

        cost_opt = inf;
        fn_opt = [];
        bn_opt = [];
        for ind_s = 1:S
            bn = zeros(S, 1);
            bn(ind_s) = 1;
            v1 = bn' * Bnn * bn;
            v2 = (Qn - Dn) * bn;
            fn = -v2 * min(1 / v1, sqrt(Pn) / norm(v2, 2));
            cost = norm(fn, 2)^2 * v1 + 2 * real(fn' * v2);
            if cost < cost_opt
                cost_opt = cost;
                fn_opt = fn;
                bn_opt = bn;
            end
        end
        FD_opt(ind_n, :) = fn_opt';
        B_opt(:, ind_n) = bn_opt;
    end

    MSE_seq(ind_iter) = cal_MSE_local(U_opt, W_opt, FD_opt, B_opt, Heff, D_all, sigma2, beta);
    WSR_seq(ind_iter + 1) = func_CalWSR(func_B2Fsel(B_opt), eye(N), FD_opt, Heff, paras);

    if ind_iter > 1 && abs(MSE_seq(ind_iter) - MSE_seq(ind_iter - 1)) < 2.8e-3
        MSE_seq = MSE_seq(1:ind_iter);
        WSR_seq = WSR_seq(1:(ind_iter + 1));
        break
    end
end

if any(isnan(MSE_seq))
    last_valid = find(~isnan(MSE_seq), 1, 'last');
    if isempty(last_valid)
        last_valid = 0;
    end
    MSE_seq = MSE_seq(1:last_valid);
    WSR_seq = WSR_seq(1:(last_valid + 1));
end

Fsel_opt = func_B2Fsel(B_opt);
end


function [Fsel_opt, FD_opt, WSR_seq, MSE_seq] = trace_model_ii(FD_init, Heff, paras, rho, seed, max_iter)
rng(seed)

T = paras.T;
K = paras.K;
N = size(Heff{1}, 2) / T;
M = nan(1, K);
for ind_k = 1:K
    M(ind_k) = size(Heff{ind_k}, 1);
end
D_all = paras.Dk;
D = sum(D_all);
sigma2 = paras.sigma2;
Pn = paras.Pn;
beta = paras.beta;

W = cell(K, 1);
for ind_k = 1:K
    W{ind_k} = eye(D_all(ind_k));
end

c_vec = func_GetSHcoef(T);
c_vec = [2 * sqrt(rho * pi); (2 * sqrt((1 - rho) * pi)) * c_vec(2:end) / norm(c_vec(2:end), 2)];
C_opt = c_vec * ones(1, N);
FD_opt = FD_init;
W_opt = W;

MSE_seq = nan(max_iter, 1);
WSR_seq = nan(max_iter + 1, 1);
WSR_seq(1) = func_CalWSR(func_B2Fsel(C_opt), eye(N), FD_opt, Heff, paras);

for ind_iter = 1:max_iter
    H = Heff;
    for ind_k = 1:K
        H{ind_k} = H{ind_k} * func_B2Fsel(C_opt);
    end

    FD_part = cell(K, 1);
    for ind_k = 1:K
        FD_part{ind_k} = FD_opt(:, (sum(D_all(1:ind_k)) - D_all(ind_k) + 1):sum(D_all(1:ind_k)));
    end

    U_opt = cell(K, 1);
    for ind_k = 1:K
        Hk = H{ind_k};
        Mk = M(ind_k);
        Matrix1 = zeros(Mk, Mk);
        for ind_i = 1:K
            FDi = FD_part{ind_i};
            Matrix1 = Matrix1 + Hk * (FDi * FDi') * Hk';
        end
        U_opt{ind_k} = (Matrix1 + sigma2 * eye(Mk))^(-1) * Hk * FD_part{ind_k};
        W_opt{ind_k} = (eye(D_all(ind_k)) - U_opt{ind_k}' * Hk * FD_part{ind_k})^(-1);
    end

    BB = zeros(N * T, N * T);
    DD = [];
    for ind_k = 1:K
        BB = BB + beta(ind_k) * Heff{ind_k}' * U_opt{ind_k} * W_opt{ind_k} * U_opt{ind_k}' * Heff{ind_k};
        DD = [DD; beta(ind_k) * W_opt{ind_k} * U_opt{ind_k}' * Heff{ind_k}]; %#ok<AGROW>
    end

    for ind_n = 1:N
        Bnn = BB((1:T) + (ind_n - 1) * T, (1:T) + (ind_n - 1) * T);
        Qn = get_Qn_local(ind_n, FD_opt, C_opt, T, D, N, BB);
        Dn = DD(:, (1:T) + (ind_n - 1) * T);

        cn = C_opt(:, ind_n);
        v1 = cn' * Bnn * cn;
        v2 = (Qn - Dn) * cn;
        fn = -v2 * min(1 / v1, sqrt(Pn) / norm(v2, 2));
        cn = func_RieOptMII(cn, Bnn, Qn, Dn, fn);

        FD_opt(ind_n, :) = fn';
        C_opt(:, ind_n) = cn;
    end

    MSE_seq(ind_iter) = cal_MSE_local(U_opt, W_opt, FD_opt, C_opt, Heff, D_all, sigma2, beta);
    WSR_seq(ind_iter + 1) = func_CalWSR(func_B2Fsel(C_opt), eye(N), FD_opt, Heff, paras);

    if ind_iter > 1 && abs(MSE_seq(ind_iter) - MSE_seq(ind_iter - 1)) < 2.8e-3
        MSE_seq = MSE_seq(1:ind_iter);
        WSR_seq = WSR_seq(1:(ind_iter + 1));
        break
    end
end

if any(isnan(MSE_seq))
    last_valid = find(~isnan(MSE_seq), 1, 'last');
    if isempty(last_valid)
        last_valid = 0;
    end
    MSE_seq = MSE_seq(1:last_valid);
    WSR_seq = WSR_seq(1:(last_valid + 1));
end

Fsel_opt = func_B2Fsel(C_opt);
end


function Qn = get_Qn_local(n, FD, B, S, D, N, BB)
Qn = zeros(D, S);
for ind_q = 1:N
    if ind_q ~= n
        Qn = Qn + FD(ind_q, :)' * B(:, ind_q)' * BB((1:S) + (ind_q - 1) * S, (1:S) + (n - 1) * S);
    end
end
end


function A = cal_A_fixed_local(U, W, beta, K)
A = [];
for ind_k = 1:K
    A = blkdiag(A, beta(ind_k) * U{ind_k} * W{ind_k} * U{ind_k}');
end
end


function B = cal_B_fixed_local(U, W, beta, K)
B = [];
for ind_k = 1:K
    B = blkdiag(B, beta(ind_k) * W{ind_k} * U{ind_k}');
end
end


function MSE = cal_MSE_fixed_local(U, W, FD, H, D_all, sigma2, beta)
K = length(D_all);
MSE = 0;
for ind_k = 1:K
    Ek = cal_Ek_fixed_local(ind_k, U, FD, H, D_all, sigma2);
    Wk = W{ind_k};
    MSE = MSE + beta(ind_k) * (trace(Wk * Ek) - log2(det(Wk)));
end
MSE = real(MSE);
end


function Ek = cal_Ek_fixed_local(k, U, FD, H, D_all, sigma2)
K = length(D_all);
Mk = size(H{k}, 1);
FD_part = cell(K, 1);
for ind_k = 1:K
    FD_part{ind_k} = FD(:, (sum(D_all(1:ind_k)) - D_all(ind_k) + 1):sum(D_all(1:ind_k)));
end

Hk = H{k};
Matrix1 = eye(D_all(k)) - U{k}' * Hk * FD_part{k};
matrix2 = sigma2 * eye(Mk);
for ind_i = 1:K
    if ind_i ~= k
        matrix2 = matrix2 + (Hk * FD_part{ind_i}) * (Hk * FD_part{ind_i})';
    end
end
Ek = Matrix1 * Matrix1' + U{k}' * matrix2 * U{k};
end


function MSE = cal_MSE_local(U, W, FD, B, Heff, D_all, sigma2, beta)
K = length(D_all);
MSE = 0;
for ind_k = 1:K
    Ek = cal_Ek_local(ind_k, U, FD, B, Heff, D_all, sigma2);
    Wk = W{ind_k};
    MSE = MSE + beta(ind_k) * (trace(Wk * Ek) - log2(det(Wk)));
end
MSE = real(MSE);
end


function Ek = cal_Ek_local(k, U, FD, B, Heff, D_all, sigma2)
K = length(D_all);
Mk = size(Heff{k}, 1);
FD_part = cell(K, 1);
for ind_k = 1:K
    FD_part{ind_k} = FD(:, (sum(D_all(1:ind_k)) - D_all(ind_k) + 1):sum(D_all(1:ind_k)));
end

Fsel = func_B2Fsel(B);
Hk = Heff{k} * Fsel;
Matrix1 = eye(D_all(k)) - U{k}' * Hk * FD_part{k};
matrix2 = sigma2 * eye(Mk);
for ind_i = 1:K
    if ind_i ~= k
        matrix2 = matrix2 + (Hk * FD_part{ind_i}) * (Hk * FD_part{ind_i})';
    end
end
Ek = Matrix1 * Matrix1' + U{k}' * matrix2 * U{k};
end
