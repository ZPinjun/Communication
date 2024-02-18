clc; clear; close all

sp = gen_defaultSetup();
sp = gen_updateSetup(sp);

% ----- mutual coupling illustration: |z_qp| vs. d
dtest = sp.lambda./[500,200,100,50,20,10,9,8,7,6,5,4,3,2,1,0.4];
data = zeros(1,length(dtest));
for i = 1:length(dtest)
    d = dtest(i);
    p_p = [0; 0; 0];
    p_q = p_p + [d; 0; 0];
    z_qp = func_MutuImp_antenna(p_p, p_q, sp);
    data(1,i) = abs(z_qp);
end

figure(1)
semilogy(dtest/sp.lambda,data);
xlabel('Distance between two radiators [$\lambda$]','Interpreter','latex')
ylabel('$|z_{qp}|\ (\Omega)$','Interpreter','latex')
grid on



    