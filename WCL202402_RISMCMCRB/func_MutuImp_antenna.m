% This function computes the mutual impedance of two antennas according to Eq. (2) & (3) of [2]
function z_qp = func_MutuImp_antenna(pp, pq, sp)
% ----- Input:
% pp: center position of the antenna p
% pq: center position of the antenna q
% sp: other parameters and constants
% ----- Output:
% z_qp: mutual/self impedance between antennas p & q

% get parameters
eta = sp.eta;       % the intrinsic impedance of free space
k0 = sp.k;          % wavenumber
hp = sp.h;          % half-length of antenna p
hq = sp.h;          % half-length of antenna q
aq = sp.a;          % antenna radius

% some preparatory calculations
if norm(pp-pq,2) == 0       % if p = q
    rho1 = aq;
    rho2 = 0;
else                        % if p != q
    rho1 = sqrt( (pp(1) - pq(1))^2 + (pp(2) - pq(2))^2 );   
    rho2 = pp(3) - pq(3);
end

% compute integral Eq. (2) of [2]
R  = @(xi,z) sqrt(rho1^2 + (z - xi + rho2).^2);
f1 = @(xi,z) exp(-1j*k0*R(xi,z)).*sin(k0*(hp-abs(xi))).*sin(k0*(hq-abs(z))) ./ ( R(xi,z)*sin(k0*hp)*sin(k0*hq) );
f2 = @(xi,z) k0^2 - 1j*k0./R(xi,z) - (k0^2*(z-xi+rho2).^2+1)./R(xi,z).^2 + 3*1j*k0*(z-xi+rho2).^2./R(xi,z).^3 + 3*(z-xi+rho2).^2./R(xi,z).^4;
f  = @(xi,z) f1(xi,z) .* f2(xi,z);
z_qp = 1j*eta/(4*pi*k0) * integral2(f, -hp, hp, -hq, hq);

end

% [1] Gradoni, Gabriele, and Marco Di Renzo. "End-to-end mutual coupling aware communication model for reconfigurable intelligent surfaces: An electromagnetic-compliant approach based on mutual impedances." IEEE Wireless Communications Letters 10.5 (2021): 938-942.
% [2] Di Renzo, Marco, Vincenzo Galdi, and Giuseppe Castaldi. "Modeling the Mutual Coupling of Reconfigurable Metasurfaces." 17th European Conference on Antennas and Propagation (EuCAP). IEEE, 2023.
