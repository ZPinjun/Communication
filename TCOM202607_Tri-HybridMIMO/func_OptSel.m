function [Fsel_opt,FD_opt,MSE_seq] = func_OptSel(FD_init,Heff,paras,seed)

rng(seed)
    
%% Initialization
% Get parameters
S = paras.S;            % Total number of antenna states
K = paras.K;            % Total number of users    
N = size(Heff{1},2)/S;  % Number of antennas at BS
M = nan(1,K);
for ind_k = 1:K         % Number of antennas at each UE
    M(ind_k) = size(Heff{ind_k},1);
end
D_all = paras.Dk;       % Number of data stream per user
D = sum(D_all);         % Total data stream at BS
sigma2 = paras.sigma2;  % Noise power
Pn = paras.Pn;          % Per-antenna power constraint
beta = paras.beta;      % Sum-rate weights
% Weight matrix
W = cell(K,1);
for ind_k = 1:K
    W{ind_k} = eye(D_all(ind_k));
end
% Selection matrix
B = zeros(S,N);
B(19,:) = 1;


%% WMMSE optimization
MaxIter = 400;
epsilon = 2.8e-3;         % Convergence threshold
% Partition Heff
Heff_part = Heff;
for ind_k = 1:K
    Heff_part{ind_k} = reshape(Heff_part{ind_k},[M(ind_k),S,N]);
end
W_opt = W;
B_opt = B;
FD_opt = FD_init;
MSE_seq = nan(MaxIter,1);
for ind_iter = 1:MaxIter

    % Process Heff
    H = Heff;
    for ind_k = 1:K
        H{ind_k} = H{ind_k}*func_B2Fsel(B_opt);
    end

    % Partition FD
    FD_part = cell(K,1);
    for ind_k = 1:K
        FD_part{ind_k} = FD_opt(:,(sum(D_all(1:ind_k))-D_all(ind_k)+1):sum(D_all(1:ind_k)));
    end

    % Update U & W
    U_opt = cell(K,1);
    for ind_k = 1:K
        Hk = H{ind_k};
        % Uk
        Mk = M(ind_k);
        Matrix1 = zeros(Mk,Mk);
        for ind_i = 1:K
            FDi = FD_part{ind_i};
            Matrix1 = Matrix1 + Hk*(FDi*FDi')*Hk';
        end
        U_opt{ind_k} = (Matrix1 + sigma2*eye(Mk))^(-1)*Hk*FD_part{ind_k};
        % Wk
        W_opt{ind_k} = (eye(D_all(ind_k)) - U_opt{ind_k}'*Hk*FD_part{ind_k})^(-1);
    end

    % Per-antenna precoders update
    BB = zeros(N*S,N*S);
    DD = [];
    for ind_k = 1:K
        BB = BB + beta(ind_k)*Heff{ind_k}'*U_opt{ind_k}*W_opt{ind_k}*U_opt{ind_k}'*Heff{ind_k};
        DD = [DD; beta(ind_k)*W_opt{ind_k}*U_opt{ind_k}'*Heff{ind_k}];
    end
    for ind_n = 1:N
        Bnn = BB((1:S)+(ind_n-1)*S,(1:S)+(ind_n-1)*S);
        Qn = get_Qn(ind_n,FD_opt,B_opt,S,D,N,BB);
        Dn = DD(:,(1:S)+(ind_n-1)*S);
        % Search the optimal state
        cost_opt = inf;
        fn_opt = [];
        bn_opt = [];
        for ind_s = 1:S
            bn = zeros(S,1);
            bn(ind_s) = 1;
            v1 = bn'*Bnn*bn;
            v2 = (Qn-Dn)*bn;
            fn = -v2*min(1/v1,sqrt(Pn)/norm(v2,2));
            cost = norm(fn,2)^2*v1 + 2*real(fn'*v2);
            if cost < cost_opt
                cost_opt = cost;
                fn_opt = fn;
                bn_opt = bn;
            end
        end
        FD_opt(ind_n,:) = fn_opt';
        B_opt(:,ind_n) = bn_opt;
    end

    % Evaluate objective function
    MSE_seq(ind_iter) = cal_MSE(U_opt,W_opt,FD_opt,B_opt,Heff,D_all,sigma2,beta);
    disp(['Optimizing based on Model I ... Iter = ', num2str(ind_iter) '; MSE =',num2str(MSE_seq(ind_iter))]);
    if ind_iter > 1 && abs(MSE_seq(ind_iter) - MSE_seq(ind_iter-1)) < epsilon
        break
    end
end


%% Assign output
Fsel_opt = func_B2Fsel(B_opt);

end


%% Functions
function Qn = get_Qn(n,FD,B,S,D,N,BB)
Qn = zeros(D,S);
for ind_q = 1:N
    if ind_q ~= n
        Qn = Qn + FD(ind_q,:)'*B(:,ind_q)'*BB((1:S)+(ind_q-1)*S,(1:S)+(n-1)*S);
    end
end
end

function MSE = cal_MSE(U,W,FD,B,Heff,D_all,sigma2,beta)
K = length(D_all);
MSE = 0;
for ind_k = 1:K
    Ek = cal_Ek(ind_k,U,FD,B,Heff,D_all,sigma2);
    Wk = W{ind_k};
    MSE = MSE + beta(ind_k)*(trace(Wk*Ek)-log2(det(Wk)));
end
MSE = real(MSE);
end

function Ek = cal_Ek(k,U,FD,B,Heff,D_all,sigma2)
K = length(D_all);
Mk = size(Heff{k},1);
% Partition FD
FD_part = cell(K,1);
for ind_k = 1:K
    FD_part{ind_k} = FD(:,(sum(D_all(1:ind_k))-D_all(ind_k)+1):sum(D_all(1:ind_k)));
end
% Term 1
Fsel = func_B2Fsel(B);
Hk = Heff{k}*Fsel;
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
