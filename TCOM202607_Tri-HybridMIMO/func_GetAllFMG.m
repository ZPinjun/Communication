function vec_basis = func_GetAllFMG(parameters,AoD_az,AoD_el)
% This function returns all available antenna magnitude gains given a
% direction and a radiation pattern 

Gmax = parameters.Gmax;
sigmae2 = parameters.sigmae2;
theta_bar = parameters.theta_bar;
phi_bar = parameters.phi_bar;

vec_basis = nan(length(theta_bar)*length(phi_bar),1);
i = 1;
for ind_theta = 1:length(theta_bar)
    th_bar = theta_bar(ind_theta);
    for ind_phi = 1:length(phi_bar)
        ph_bar = phi_bar(ind_phi);
        u_bar = func_GetUnitVec(th_bar,ph_bar);
        u = func_GetUnitVec(AoD_el,AoD_az);
        vec_basis(i) = parameters.Gmax*exp(-acos(u.'*u_bar).^2/parameters.sigmae2);
        i = i + 1;
    end
end

end