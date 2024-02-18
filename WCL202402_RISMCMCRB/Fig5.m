clc; clear; close all

rng(0)

%% simulation parameters
dtest = [0.02, 0.1, 0.5];           % different inter-distances of RIS elements to be tested
Ptest = -10:10:80;                  % different signal power (SNR) to be tested
G = 256;                            % times of transmissions

%% preparatory calculation
% mutual impedance zST & zRS (we fix zST & zRS for the LB/RMSE test to evaluate the impact of mutual coupling only)
sp = gen_defaultSetup();
sp = gen_updateSetup(sp);
NRIS = sp.RIS_dim(1)*sp.RIS_dim(2);
zST = zeros(NRIS,1);    % zST: mutual impedance between Tx and RIS
zRS = zeros(NRIS,1);    % zRS: mutual impedance between RIS and Rx
for i = 1:NRIS
    p_q = sp.RIS_G(:,i);
    zST(i,1) = func_MutuImp_antenna(sp.p_T, p_q, sp);
    zRS(i,1) = func_MutuImp_antenna(sp.p_R, p_q, sp);
end
% RIS mutual impedance for different inter-distances
ZSS_all = zeros(NRIS,NRIS,length(dtest));
z_self = func_MutuImp_antenna(zeros(3,1), zeros(3,1), sp);
ZSS_self = diag(z_self*ones(NRIS,1));
for i = 1:length(dtest)
    sp = gen_defaultSetup();
    sp.RIS_spacing = dtest(i);  
    sp = gen_updateSetup(sp);
    ZSS_mutual = zeros(NRIS, NRIS);
    for r = 1:NRIS
        p_p = sp.RIS_G(:,r);
        for c = r+1:NRIS
            p_q = sp.RIS_G(:,c);
            z_qp = func_MutuImp_antenna(p_p, p_q, sp);
            ZSS_mutual(r,c) = z_qp;
        end
    end
    ZSS_mutual = ZSS_mutual + ZSS_mutual.';     % by reciprocity
    ZSS_all(:,:,i) = ZSS_self + ZSS_mutual;                
end
% RIS tunable loads
RSmn = 0.1 + 10*rand(NRIS,G);              % randomly generated between [0.1, 10.1] (Ohm)
LSmn = 1e-9 * (0.1 + 10*rand(NRIS,G));     % randomly generated between [0.1, 10.1] (nH)
zRIS_all = RSmn + 1j*2*pi*sp.f*LSmn;
% quantities for MCRB calculation
x_bar = [real(zST); imag(zST)];
B_true = zeros(NRIS,G,length(dtest));
B_tilde = zeros(NRIS,G);
for g = 1:G
    % for dtest(1)
    B_true(:,g,1) = ( zRS.'*(ZSS_all(:,:,1) + diag(zRIS_all(:,g)))^(-1) ).';
    % for dtest(2)
    B_true(:,g,2) = ( zRS.'*(ZSS_all(:,:,2) + diag(zRIS_all(:,g)))^(-1) ).';
    % for dtest(3)
    B_true(:,g,3) = ( zRS.'*(ZSS_all(:,:,3) + diag(zRIS_all(:,g)))^(-1) ).';
    % for mismatched case
    B_tilde(:,g) = ( zRS.'*(ZSS_all(:,:,1).*diag(ones(NRIS,1)) + diag(zRIS_all(:,g)))^(-1) ).'; 
end
D_true = zeros(2*G,2*NRIS,3);
for i = 1:3
    B = B_true(:,:,i);
    D_true(:,:,i) = [real(B).', -imag(B).'; imag(B).', real(B).'];
end
D_tilde = [real(B_tilde).', -imag(B_tilde).'; imag(B_tilde).', real(B_tilde).'];


%% compute lower bounds
sigma2 = sp.sigma^2; 
LB_nomismatch = zeros(3,length(Ptest));
LB_mismatched = zeros(3,length(Ptest));
Bias = zeros(3,1);
for i = 1:length(Ptest)

    PT = db2pow(Ptest(i));
    SNR = PT/sigma2;
    LBM_nomismatch = zeros(2*NRIS, 2*NRIS,length(dtest));
    LBM_mismatched = zeros(2*NRIS, 2*NRIS,length(dtest));
    
    % --------- case 1: mismatch-free (classical CRLB)
    % for dtest(1)
    D = D_true(:,:,1);
    LBM_nomismatch(:,:,1) = (D.'*D)^(-1)/(2*SNR);
    % for dtest(2)
    D = D_true(:,:,2);
    LBM_nomismatch(:,:,2) = (D.'*D)^(-1)/(2*SNR);
    % for dtest(3)
    D = D_true(:,:,3);
    LBM_nomismatch(:,:,3) = (D.'*D)^(-1)/(2*SNR);
    % --------- case 2: mismatched model (derived LB)
    % for dtest(1)
    D = D_true(:,:,1);
    x0 = (D_tilde.'*D_tilde)^(-1)*D_tilde.'*D*x_bar;
    Bias(1) = (x_bar-x0).'*(x_bar-x0);
    LBM_mismatched(:,:,1) = (D_tilde.'*D_tilde)^(-1)/(2*SNR) + (x_bar-x0)*(x_bar-x0).';
    % for dtest(2)
    D = D_true(:,:,2);
    x0 = (D_tilde.'*D_tilde)^(-1)*D_tilde.'*D*x_bar;
    Bias(2) = (x_bar-x0).'*(x_bar-x0);
    LBM_mismatched(:,:,2) = (D_tilde.'*D_tilde)^(-1)/(2*SNR) + (x_bar-x0)*(x_bar-x0).';
    % for dtest(3)
    D = D_true(:,:,3);
    x0 = (D_tilde.'*D_tilde)^(-1)*D_tilde.'*D*x_bar;
    Bias(3) = (x_bar-x0).'*(x_bar-x0);
    LBM_mismatched(:,:,3) = (D_tilde.'*D_tilde)^(-1)/(2*SNR) + (x_bar-x0)*(x_bar-x0).';

    LB_nomismatch(1,i) = sqrt(trace(LBM_nomismatch(:,:,1)));
    LB_nomismatch(2,i) = sqrt(trace(LBM_nomismatch(:,:,2)));
    LB_nomismatch(3,i) = sqrt(trace(LBM_nomismatch(:,:,3)));
    LB_mismatched(1,i) = sqrt(trace(LBM_mismatched(:,:,1)));
    LB_mismatched(2,i) = sqrt(trace(LBM_mismatched(:,:,2)));
    LB_mismatched(3,i) = sqrt(trace(LBM_mismatched(:,:,3)));
end

% plot derived lower bounds
figure(1)
semilogy(Ptest,LB_mismatched(1,:),'-sb'); hold on
semilogy(Ptest,LB_mismatched(2,:),'-og'); hold on
semilogy(Ptest,LB_mismatched(3,:),'-dr'); hold on
semilogy(Ptest,LB_nomismatch(1,:),'--b'); hold on
semilogy(Ptest,LB_nomismatch(2,:),'--g'); hold on
semilogy(Ptest,LB_nomismatch(3,:),'-.r'); hold on
xlabel('Transmit power [dBm]','Interpreter','latex')
ylabel('RMSE','Interpreter','latex')
grid on


%% test RMSE of ML estimator
wb = waitbar(0,'Test RMSE using maximum likelihood estimator ...');
data_RMSE = zeros(3,length(Ptest));
repeatNum = 1000;
for i = 1:length(Ptest)
    waitbar(i/length(Ptest),wb);
    PT = db2pow(Ptest(i));
    Noise = randn(2*G,repeatNum) * sqrt(0.5*sigma2);
    y1_pure = sqrt(PT)*D_true(:,:,1)*x_bar;
    y2_pure = sqrt(PT)*D_true(:,:,2)*x_bar;
    y3_pure = sqrt(PT)*D_true(:,:,3)*x_bar;
    data_error = zeros(repeatNum,3);
    for j = 1:repeatNum
        y1 = y1_pure + Noise(:,j);
        y2 = y2_pure + Noise(:,j);
        y3 = y3_pure + Noise(:,j);
        x_hat = (D_tilde.'*D_tilde)^(-1)*D_tilde.'*[y1,y2,y3]./sqrt(PT);
        data_error(j,1) = norm(x_hat(:,1) - x_bar,2)^2;
        data_error(j,2) = norm(x_hat(:,2) - x_bar,2)^2;
        data_error(j,3) = norm(x_hat(:,3) - x_bar,2)^2;
    end
    data_RMSE(:,i) = sqrt(mean(data_error)).';
end
close(wb);

semilogy(Ptest,data_RMSE(1,:),'-+b'); hold on
semilogy(Ptest,data_RMSE(2,:),'-xg'); hold on
semilogy(Ptest,data_RMSE(3,:),'-*r'); hold on
semilogy(Ptest,sqrt(Bias(1))*ones(size(Ptest)),'-.b'); hold on
semilogy(Ptest,sqrt(Bias(2))*ones(size(Ptest)),'-.g'); hold on
semilogy(Ptest,sqrt(Bias(3))*ones(size(Ptest)),'-.r'); hold on
legend('LB, $d=0.02\lambda$','LB, $d=0.1\lambda$','LB, $d=0.5\lambda$',...
    'CRLB, $d=0.02\lambda$','CRLB, $d=0.1\lambda$','CRLB, $d=0.5\lambda$',...
    'RMSE, $d=0.02\lambda$','RMSE, $d=0.1\lambda$','RMSE, $d=0.5\lambda$','Interpreter','latex');
set(gca,'GridLineStyle','--')




    