function paras = func_GetParas(H_cv,H_sel)
% This function collects and defines all the necessary system parameters

% System parameters
paras.K = length(H_cv);
paras.S = size(H_sel{1},2)/size(H_cv{1},2);

% Signal parameters
paras.fc = 30e9;                            % Center Frequency [Hz]
paras.c0 = 299792458;                       % Speed of Light [m/s]
paras.lambda = paras.c0/paras.fc;           % Wavelength
paras.Dk = 2*ones(paras.K,1);               % Number of data streams at each UE
paras.sigma2 = 10^(-90/10);                 % Noise power at each UE       
paras.Pn = 10^(0/10);                       % Power constraint per antenna
paras.NRF = sum(paras.Dk)+3;                % Number of RFCs
paras.beta = ones(paras.K,1)/paras.K;       % Weight assigned to each user's rate



end

