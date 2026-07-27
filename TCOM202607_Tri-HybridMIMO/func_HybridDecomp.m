% [1] H. Huang et al., "Hybrid Precoder Design for Angle-of-Departure
% Estimation with Limited-Resolution Phase Shifters," in IEEE Transactions
% on Communications, doi: 10.1109/TCOMM.2024.3487517.

function [FRF_opt,FBB_opt] = func_HybridDecomp(FD,paras,seed)

rng(seed);
NRF = paras.NRF;
Pn = paras.Pn;
N = size(FD,1);
D = sum(paras.Dk);

%% Generate a random initialization FRF_init
FRF_init = (1/sqrt(N))*exp(1j*2*pi*randn(N,NRF));
FRF_init(:,1:D) = (1/sqrt(N))*exp(1j*angle(FD));

%% Decompose FD into FRF and FBB (sum-power constraint) according to Algo. 1 in [1]
Pmax = norm(FD,'fro')^2;
FRFdcom = FRF_init;
Imax = 2000;
Kmax = 200;
found = false;
RelativeError_old = inf;
for i = 1:Imax
    FRF_SeudoInv = (FRFdcom'*FRFdcom)^(-1)*FRFdcom';
    FBBdcom = sqrt(Pmax)*FRF_SeudoInv*FD/norm(FRFdcom*FRF_SeudoInv*FD,'fro');
    rho = max([sqrt(2)*norm(FBBdcom*FBBdcom','fro'),norm(FBBdcom,'fro')^2]);
    FRF_tilde = FRFdcom;
    U = zeros(size(FRFdcom));
    for k = 1:Kmax
        FRFdcom = (1/sqrt(N))*exp(1j*angle(FRF_tilde+U));
        FRF_tilde = (FD*FBBdcom'+rho*(FRFdcom-U))*(FBBdcom*FBBdcom'+rho*eye(NRF))^(-1);
        U = U + FRF_tilde - FRFdcom;
    end
    FRFdcom = (1/sqrt(N))*exp(1j*angle(FRF_tilde));
    RelativeError = norm(FD-FRFdcom*FBBdcom,'fro')^2/norm(FD,'fro')^2;
    disp(['Fully digital to hybrid decomposition ... Relative Error =',num2str(RelativeError)]);
    if abs(RelativeError-RelativeError_old) < 1e-6 || (RelativeError-RelativeError_old) > 1e-4 
        break;
    end
    RelativeError_old = RelativeError;
end

FRF_SeudoInv = (FRFdcom'*FRFdcom)^(-1)*FRFdcom';
FBBdcom = sqrt(Pmax)*FRF_SeudoInv*FD/norm(FRFdcom*FRF_SeudoInv*FD,'fro');

%% Normalization to stay within per-antenna power constraint
v = nan(N,1);
Matrix1 = FRFdcom*FBBdcom;
Matrix2 = Matrix1*Matrix1';
for ind_n = 1:N
    v(ind_n) = sqrt(Pn/Matrix2(ind_n,ind_n));
end
v = [v; 1];
FBB_opt = FBBdcom*min(v);
FRF_opt = FRFdcom;

end

