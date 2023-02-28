function [WP] = FMKKM_update_Wp(Hstar,HP)

k = size(HP,2);
numker = size(HP,3);
WP = zeros(k,k,numker);
for p = 1 : numker
    Tp = HP(:,:,p)'*Hstar;
    [Up,~,Vp] = svd(Tp,'econ');
    WP(:,:,p) = Up*Vp';
end
end


