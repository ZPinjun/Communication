function Fsel = func_B2Fsel(B)
Numc = size(B,2);
Fsel = [];
for ind_c = 1:Numc
    Fsel = blkdiag(Fsel,B(:,ind_c));
end
end

