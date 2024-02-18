% This function updates the system parameters 
function sp = gen_updateSetup(sp)

sp.lambda = sp.c/sp.f;          % signal wavelength [m]
sp.k = 2*pi/sp.lambda;          % wavenumber
sp.h = sp.lambda*sp.h_sca;      % antenna half-lenth 
sp.a = sp.lambda*sp.a_sca;      % antenna radius

%% generate coordinates of RIS elements:
% arrays are deployed on the X-O-Y plane, Z-axis is the normal direction
RIS_spacing = sp.RIS_spacing;
RIS_dim = sp.RIS_dim;
lambda = sp.lambda;
xrange = (  (0:RIS_dim(1)-1) - 0.5*(RIS_dim(1)-1)  ) * (RIS_spacing*lambda);
yrange = (  (0:RIS_dim(2)-1) - 0.5*(RIS_dim(2)-1)  ) * (RIS_spacing*lambda);
pRIS_L = [kron(xrange, ones(1,RIS_dim(1)));
        kron(ones(1,RIS_dim(1)), yrange);
        zeros(1,RIS_dim(1)*RIS_dim(2))];
pRIS_G = sp.p_RIS + pRIS_L;
sp.RIS_G = pRIS_G;

%% compute thermal noise power
sp.BW = 1;
N0 = sp.Kb*sp.T*1000;                               % thermal noise PSD [mW/Hz]
Pn = N0*sp.BW;                                      % thermal noise power [mW]
sigma_in = sqrt(Pn);                                % input noise standard deviation
sp.sigma = sqrt(10^(sp.NoiseFigure/10))*sigma_in;   % thermal noise standard deviation at the receiver

end

