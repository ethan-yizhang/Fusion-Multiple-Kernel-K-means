function [Gamma] = FMKKM_update_gamma(KH,Hstar)
%OPM 此处显示有关此函数的摘要
%   此处显示详细说明
num = size(KH,1);
numker = size(KH,3);
U0 = eye(num)-Hstar*Hstar';
A = zeros(numker,1);
for  p=1:numker
    A(p) = trace(KH(:,:,p)*U0);
end
% Gamma = A./norm(A); %max Sigma r A
% Gamma((Gamma<eps))=0;
% Gamma = Gamma./norm(Gamma);
Gamma = (1./A)/sum(1./A);%min Sigma r^2 A
Gamma((Gamma<eps))=0;
Gamma = Gamma./sum(Gamma);
end


