% This function generates the default system parameters
function sp = gen_defaultSetup()

% system parameters
sp.p_T = [5; -5; 3];        % position of transmitter
sp.p_R = [5; 5; 1];         % position of reveiver
sp.p_RIS = [0; 0; 0];       % position of RIS

% signal parameters
sp.f = 28e9;                % signal frequency [Hz]

% antenna parameters
sp.h_sca = 1/64;            % antenna half-length = wavelenth*sp.h_sca
sp.a_sca = 1/500;           % antenna radius = wavelenth*sp.a_sca
sp.RIS_spacing = 0.2;       % RIS element spacing = lambda*sp.RIS_spacing
sp.RIS_dim = [4, 4];        % RIS antenna dimension

% environment parameters
sp.T = 298.15;              % System temperature (Kelvin), 298.15 Kelvin = 25 celsius
sp.NoiseFigure = 10;         % noise figure in [dB]

% constants
sp.c = 3e8;                 % speed of light [m/s]
sp.mu = 4*pi*10^(-7);       % the magnetic permeability of free space [N/A^2]
sp.epsilon = 8.854e-12;     % the electric permittivity of free space [F/m]
sp.eta = 377;               % the intrinsic impedance of free space [Ohm]
sp.Kb = 1.3806e-23;         % Boltzmann constant


end

