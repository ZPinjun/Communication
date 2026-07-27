function [FD_opt,MSE_seq] = func_OptCV(H,paras,seed)

rng(seed)
    
%% Initialization
% Get parameters
K = paras.K;            % Total number of users    
N = size(H{1},2);       % Number of antennas at BS
M = nan(1,K);
for ind_k = 1:K         % Number of antennas at each UE
    M(ind_k) = size(H{ind_k},1);
end
D_all = paras.Dk;       % Number of data stream per user
D = sum(D_all);         % Total data stream at BS
sigma2 = paras.sigma2;  % Noise power
Pn = paras.Pn;          % Per-antenna power constraint
beta = paras.beta;      % Sum-rate weights
% Fully digital precoder
FD = randn(N,D) + 1j*randn(N,D);
for ind_n = 1:N
    FD(ind_n,:) = Pn*FD(ind_n,:)/norm(FD(ind_n,:),2);
end
% Weight matrix
W = cell(K,1);
for ind_k = 1:K
    W{ind_k} = eye(D_all(ind_k));
end


%% WMMSE optimization
MaxIter = 400;
epsilon = 2.5e-3;         % Convergence threshold
W_opt = W;
FD_opt = FD;
MSE_seq = nan(MaxIter,1);
for ind_iter = 1:MaxIter

    % Partition FD
    P = cell(K,1);
    for ind_k = 1:K
        P{ind_k} = FD_opt(:,(sum(D_all(1:ind_k))-D_all(ind_k)+1):sum(D_all(1:ind_k)));
    end

    % Update U & W
    U_opt = cell(K,1);
    for ind_k = 1:K
        Hk = H{ind_k};
        % Uk
        Mk = M(ind_k);
        Matrix1 = zeros(Mk,Mk);
        for ind_j = 1:K
            Pj = P{ind_j};
            Matrix1 = Matrix1 + Hk*(Pj*Pj')*Hk';
        end
        U_opt{ind_k} = (Matrix1 + sigma2*eye(Mk))^(-1)*Hk*P{ind_k};
        % Wk
        W_opt{ind_k} = (eye(D_all(ind_k)) - U_opt{ind_k}'*Hk*P{ind_k})^(-1);
    end

    % Per-antenna precoders update
    A = cal_A(U_opt,W_opt,beta,K);
    B = cal_B(U_opt,W_opt,beta,K);
    Hconc = [];
    for ind_k = 1:K
        Hconc = [Hconc;H{ind_k}];
    end
    for ind_n = 1:N
        % hn
        hn = Hconc(:,ind_n);
        % an
        an = hn'*A*hn;
        % C
        C = zeros(D,sum(M));
        for ind_l = 1:N
            if ind_l ~= ind_n
                hl = Hconc(:,ind_l);
                C = C + FD_opt(ind_l,:)'*hl';
            end
        end
        % bn
        bn = -B*hn + C*A*hn;
        % update
        pn = -bn*min([1/an,sqrt(Pn)/norm(bn,2)]);
        FD_opt(ind_n,:) = pn';
    end

    % Evaluate objective function
    MSE_seq(ind_iter) = cal_MSE(U_opt,W_opt,FD_opt,H,D_all,sigma2,beta);
    if ind_iter > 1 && abs(MSE_seq(ind_iter) - MSE_seq(ind_iter-1)) < epsilon
        break
    end
end

end


%% Functions
function A = cal_A(U,W,beta,K)
A = [];
for ind_k = 1:K
    A = blkdiag(A,beta(ind_k)*U{ind_k}*W{ind_k}*U{ind_k}');
end
end

function B = cal_B(U,W,beta,K)
B = [];
for ind_k = 1:K
    B = blkdiag(B,beta(ind_k)*W{ind_k}*U{ind_k}');
end
end

function MSE = cal_MSE(U,W,FD,H,D_all,sigma2,beta)
K = length(D_all);
MSE = 0;
for ind_k = 1:K
    Ek = cal_Ek(ind_k,U,FD,H,D_all,sigma2);
    Wk = W{ind_k};
    MSE = MSE + beta(ind_k)*(trace(Wk*Ek)-log2(det(Wk)));
end
MSE = real(MSE);
end

function Ek = cal_Ek(k,U,FD,H,D_all,sigma2)
K = length(D_all);
Mk = size(H{k},1);
% Partition FD
FD_part = cell(K,1);
for ind_k = 1:K
    FD_part{ind_k} = FD(:,(sum(D_all(1:ind_k))-D_all(ind_k)+1):sum(D_all(1:ind_k)));
end
% Term 1
Hk = H{k};
Matrix1 = eye(D_all(k)) - U{k}'*Hk*FD_part{k};
% Term 2
matrix2 = sigma2*eye(Mk);
for ind_i = 1:K
    if ind_i ~= k
        matrix2 = matrix2 + (Hk*FD_part{ind_i})*(Hk*FD_part{ind_i})';
    end
end
%
Ek = Matrix1*Matrix1' + U{k}'*matrix2*U{k};
end
