function KH3 = algorithm_laplacian_imputation(KH,S,numclass,alpha0)

num = size(KH,1);
numker = size(KH,3);
gamma0 = ones(numker,1)/numker;
KH3 = zeros(num,num,numker);
for p =1:numker
    indx = setdiff(1:num,S{p}.indx');
    KAp = KH(indx,indx,p);
    KH3(indx,indx,p) = (KAp+KAp')/2;
end
clear KH
Kmatrix  = mycombFun(KH3,gamma0.^2);
H = mykernelkmeans(Kmatrix,numclass);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kx = eye(num) - H*H';
%% 
mis_indx = S{p}.indx;
obs_indx = setdiff(1:num,mis_indx);
for p =1:numker
    KH3(:,:,p) = absentKernelImputation(Kx,KH3(obs_indx,obs_indx,p),mis_indx,alpha0);
end