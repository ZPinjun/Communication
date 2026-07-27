% Reference: 
% [1] J. Choi, S. Han and J. Joung, "Low-Complexity Multiuser MIMO Precoder
% Design Under Per-Antenna Power Constraints," in IEEE Transactions on
% Vehicular Technology, vol. 67, no. 9, pp. 9011-9015, Sept. 2018, doi:
% 10.1109/TVT.2018.2849021.   
% [2] Q. H. Spencer, A. L. Swindlehurst and M. Haardt, "Zero-forcing
% methods for downlink spatial multiplexing in multiuser MIMO channels," in
% IEEE Transactions on Signal Processing, vol. 52, no. 2, pp. 461-471, Feb.
% 2004, doi: 10.1109/TSP.2003.821107.   

function [FD_ZF, Ak] = func_SumRateBD(H_all, P_total, paras)
% This function implements the weighted sum-capacity block diagonalization
% in [1],[2]
%
% Inputs:
%   H_all           - Cell array of user channels: H_all{k} is [Nr x Nt]
%   P_total         - Total transmit power
%   noisePower      - Noise variance at each user
%   weights         - User weight vector [w1, w2, ..., wK]
%
% Outputs:
%   FD_ZF           - ZF precoder
%   Ak              - Ak defined in (4) in [1]
%   C_weighted_sum  - Total weighted sum capacity

noisePower = paras.sigma2;
weights = paras.beta;
Dk = paras.Dk;

K = length(H_all);
V_all = cell(K,1);
P_all = cell(K,1);

eig_list = [];      % Will store all subchannel gains
weight_list = [];   % Associated user weights
map_list = [];      % Mapping from global subchannel index to (user, stream)

null_spaces = cell(K,1);
sing_vals = cell(K,1);

% === Step 1: Block Diagonalization + effective channel decomposition ===
for k = 1:K
    % Construct interference channel (excluding user k)
    H_bar = [];
    for j = [1:k-1, k+1:K]
        H_bar = [H_bar; H_all{j}];
    end
    [~, ~, Vh] = svd(H_bar);
    null_space = Vh(:, size(H_bar,1)+1:end);
    null_spaces{k} = null_space;

    % Effective channel
    H_eff = H_all{k} * null_space;
    [~, S, ~] = svd(H_eff);
    d_k = Dk(k);
    sing_vals{k} = diag(S(1:d_k, 1:d_k)).^2;

    % Store mapping info
    eig_list = [eig_list; sing_vals{k}];
    weight_list = [weight_list; weights(k) * ones(d_k,1)];
    map_list = [map_list; [k * ones(d_k,1), (1:d_k)']];
end

% === Step 2: Global Weighted Waterfilling ===
p_all = weightedWaterfilling(eig_list, weight_list, P_total, noisePower);

% === Step 3: Map power back to users and build precoders ===
Ak = cell(K,1);
for k = 1:K
    d_k = length(sing_vals{k});
    p_k = zeros(d_k,1);

    for i = 1:d_k
        idx = find(map_list(:,1) == k & map_list(:,2) == i);
        p_k(i) = p_all(idx);
    end

    P_all{k} = diag(p_k);
    H_eff = H_all{k} * null_space;
    [~, ~, V_svd] = svd(H_eff);
    Ak{k} = V_svd(:, 1:d_k)*sqrt(P_all{k});
    V_all{k} = null_spaces{k} * Ak{k};
end

FD_ZF = [];
for ind = 1:length(V_all)
    FD_ZF = [FD_ZF,V_all{ind}];
end

end

function p = weightedWaterfilling(g, w, P, sigma2)
% Weighted waterfilling across subchannels with weights w and gains g

N = length(g);
mu_low = 0;
mu_high = 1e6;
delta = 1e-6;
max_iter = 1000;

for iter = 1:max_iter
    mu = (mu_low + mu_high) / 2;
    p = max(0, (w ./ (mu*log(2)) - sigma2 ./ g));
    if sum(p) > P
        mu_low = mu;
    else
        mu_high = mu;
    end
    if abs(sum(p) - P) < delta
        break;
    end
end

end