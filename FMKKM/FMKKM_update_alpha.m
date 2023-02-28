function [Alpha] = FMKKM_update_alpha(KH,HP)
%OPM 此处显示有关此函数的摘要
%   此处显示详细说明
num = size(KH,1);
numker = size(KH,3);
A = zeros(numker,1);
for  p=1:numker
    U0 = eye(num)-HP(:,:,p)*HP(:,:,p)';
    A(p) = trace(KH(:,:,p)*U0);
end
% Gamma = A./norm(A); %max Sigma r A
% Gamma((Gamma<eps))=0;
% Gamma = Gamma./norm(Gamma);
Alpha = (1./A)/sum(1./A);%min Sigma r^2 A
Alpha((Alpha<eps))=0;
Alpha = Alpha./sum(Alpha);
end


