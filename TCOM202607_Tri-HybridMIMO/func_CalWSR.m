function [WSR,SRs] = func_CalWSR(Fsel,FRF,FBB,Heff,paras)
% This script calculates the weighted sum-rate given a set of precoders
% {Fsel,FRF,FBB} 

%% Get parameters
K = paras.K;            % Number of users
beta = paras.beta;      % Sum-rate weights
D_all = paras.Dk;       % Number of data stream per user
sigma2 = paras.sigma2;  % Noise power

%% Partition FRF*FBB = FD
FD = FRF*FBB;
FD_part = {};
for ind_k = 1:K
    FD_part{ind_k} = FD(:,(sum(D_all(1:ind_k))-D_all(ind_k)+1):sum(D_all(1:ind_k)));
end

%% Compute weighted sum-rate
SRs = nan(K,1);
for ind_k = 1:K
    Hk = Heff{ind_k}*Fsel;
    Mk = size(Hk,1);
    Matrix1 = sigma2*eye(Mk);
    for ind_i = 1:K
        if ind_i ~= ind_k
            Matrix1 = Matrix1 + (Hk*FD_part{ind_i})*(Hk*FD_part{ind_i})';
        end
    end
    Matrix2 = eye(Mk) + (Hk*FD_part{ind_k})*(Hk*FD_part{ind_k})'*Matrix1^(-1);
    SRs(ind_k) = beta(ind_k)*log2(det(Matrix2));
end

SRs = real(SRs);
WSR = sum(SRs);
WSR = real(WSR);

end

