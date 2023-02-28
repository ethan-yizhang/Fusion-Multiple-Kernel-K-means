function [Beta] = FMKKM_update_beta(Hstar,HP,WP)
%OPM 此处显示有关此函数的摘要
%   此处显示详细说明
numker = size(HP,3);
A = zeros(numker,1);
for  p=1:numker
    A(p) = trace(Hstar'*(HP(:,:,p)*WP(:,:,p)));
end
Beta = A./norm(A);
Beta((Beta<eps))=0;
Beta = Beta./norm(Beta);
end

