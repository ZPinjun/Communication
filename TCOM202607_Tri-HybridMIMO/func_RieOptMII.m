function cn_opt = func_RieOptMII(cn_init,Bnn,Qn,Dn,fn)

% Decompose to fix the zero-order component
rho = cn_init(1)^2/norm(cn_init,2)^2;
X_init = cn_init(2:end)/sqrt(4*(1-rho)*pi);
if abs(norm(X_init,2) - 1) < 1e-3
    X_init = X_init/norm(X_init,2);
else
    error('Incorrect initialization!!!')
end
b2 = Bnn(1,2:end).';
b3 = Bnn(2:end,1);
Bnn_tilde = Bnn(2:end,2:end);
v1 = norm(fn,2)^2*4*pi*sqrt(rho*(1-rho))*(b2+b3);
s0 = norm(fn,2)^2*4*pi*(1-rho);
v0 = 2*real(fn'*(Qn-Dn));
v2 = 2*sqrt((1-rho)*pi)*v0(2:end).';
% % Check derivation
% cost1 = norm(fn,2)^2*cn_init.'*Bnn*cn_init + 2*real(fn'*(Qn-Dn))*cn_init;
% cn = X_init;
% cost2 = v1.'*cn + 4*rho*pi*Bnn(1,1) + s0*cn.'*Bnn_tilde*cn + v2.'*cn + v0(1)*2*sqrt(rho*pi);



%% =================== use manopt toolbox
% Create the problem structure
Manifold = spherefactory(length(X_init));
problem.M = Manifold;
% Define the problem cost function and its Euclidean derivatives
problem.cost = @(X) cost_function(X,Bnn_tilde,v1,v2,s0);
problem.egrad = @(X) egrad_function(X,Bnn_tilde,v1,v2,s0);
% Solve the problem
warning('off', 'manopt:getHessian:approx')
warning('off', 'manopt:getGradient:approx')
opt.minstepsize = 1e-20;
opt.maxiter = 10;
% opt.tolgradnorm = 1e-10;
opt.verbosity = 0;
[X_opt,~] = trustregions(problem, X_init, opt);
cn_opt = [sqrt(4*rho*pi);sqrt(4*(1-rho)*pi)*X_opt];
cn_opt = real(cn_opt);

end

function f_cost = cost_function(cn,Bnn_tilde,v1,v2,s0)
f_cost = (v1+v2).'*cn + s0*cn.'*Bnn_tilde*cn;
end

function egrad = egrad_function(cn,Bnn_tilde,v1,v2,s0)
egrad = v1 + v2 + s0*(Bnn_tilde+Bnn_tilde.')*cn;
end

