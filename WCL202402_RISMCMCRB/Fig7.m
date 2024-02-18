clc; clear; close all

rng(0);

dtest = 1./[500,200,100,50,20,10,5,2,1,0.4];        % different inter-distances of RIS elements to be tested
sizetest = [4,8,12];                                % different RIS size to be tested
G = 256;                                            % number of transmissions
PT = db2pow(40);                                    % fix the transmission power
data = zeros(length(sizetest),length(dtest));

wb = waitbar(0,'Evaluate mismatch-free CRLB under different inter-distances and RIS size ...');
for s = 1:length(sizetest)
    sp = gen_defaultSetup();
    sp.RIS_dim = [sizetest(s), sizetest(s)];
    sp = gen_updateSetup(sp);
    % RIS tunable loads
    NRIS = sp.RIS_dim(1)*sp.RIS_dim(2);
    RSmn = 0.1 + 10*rand(NRIS,G);              % randomly generated between [0.1, 10.1] (Ohm)
    LSmn = 1e-9 * (0.1 + 10*rand(NRIS,G));     % randomly generated between [0.1, 10.1] (nH)
    zRIS = RSmn + 1j*2*pi*sp.f*LSmn;
    for d = 1:length(dtest)
        waitbar((d+(s-1)*length(dtest))/(length(sizetest)*length(dtest)),wb);

        sp = gen_defaultSetup();
        sp.RIS_dim = [sizetest(s), sizetest(s)];
        sp.RIS_spacing = dtest(d);
        sp = gen_updateSetup(sp);
        % mutual impedance zST & zRS
        zST = zeros(NRIS,1);    % zST: mutual impedance between Tx and RIS
        zRS = zeros(NRIS,1);    % zRS: mutual impedance between RIS and Rx
        for i = 1:NRIS
            p_q = sp.RIS_G(:,i);
            zST(i,1) = func_MutuImp_antenna(sp.p_T, p_q, sp);
            zRS(i,1) = func_MutuImp_antenna(sp.p_R, p_q, sp);
        end
        % RIS mutual impedances
        z_self = func_MutuImp_antenna(zeros(3,1), zeros(3,1), sp);
        ZSS_self = diag(z_self*ones(NRIS,1));
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
            ZSS = ZSS_self + ZSS_mutual;
        % quantities for MCRB calculation
        x_bar = [real(zST); imag(zST)];
        B_true = zeros(NRIS,G);
        for g = 1:G
            B_true(:,g) = ( zRS.'*(ZSS + diag(zRIS(:,g)))^(-1) ).';
        end
        D_true = [real(B_true).', -imag(B_true).'; imag(B_true).', real(B_true).'];

        sigma2 = sp.sigma^2;
        SNR = PT/sigma2;
        data(s,d) = sqrt( trace((D_true.'*D_true)^(-1)/(2*SNR)) );
    end
end
close(wb);

figure(1)
loglog(dtest,data(1,:),'-sb'); hold on
loglog(dtest,data(2,:),'-go'); hold on
loglog(dtest,data(3,:),'-rv'); hold on
xlabel('Antenna distance [$\lambda$]','Interpreter','latex')
ylabel('$\mathrm{Bias}$','Interpreter','latex')
legend('RIS size $4\times 4$','RIS size $8\times 8$','RIS size $12\times 12$','Interpreter','latex')
grid on

