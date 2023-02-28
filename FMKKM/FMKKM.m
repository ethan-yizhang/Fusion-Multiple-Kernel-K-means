function [Hstar,Alpha,Beta,Gamma,WP,obj,obj1,obj2,obj3] = FMKKM(KH,lambda1,lambda2,numclass,option)

% lambda1 = 1;
% lambda2 = 1;
num = size(KH,1);
numker = size(KH,3);

if ~isfield(option,'goldensearch_deltmax')
    option.goldensearch_deltmax=5e-2;
end
if ~isfield(option,'goldensearchmax')
    option.goldensearchmax=1e-8;
end
if ~isfield(option,'firstbasevariable')
    option.firstbasevariable='first';
end


%% Initializing
Alpha = ones(numker,1)/numker;
Beta = ones(numker,1)/sqrt(numker);
Gamma = ones(numker,1)/numker;
[HP,WP] = myInitialiHp(KH,numclass);
Kmatrix = sumKbeta(KH,Gamma.^2);
[Hstar,~]= mykernelkmeans(Kmatrix,numclass);

iter = 1;
flag = 1;

[obj(iter),obj1(iter),obj2(iter),obj3(iter)] = FMKKM_Objectve(KH,Hstar,HP,WP,Alpha,Beta,Gamma,lambda1,lambda2);

while flag
    iter = iter + 1;
    [WP] = FMKKM_update_Wp(Hstar,HP);

    A = zeros(num,num);
    C = zeros(num,numclass);
    for p = 1:numker
        A = A + Gamma(p).^2 * KH(:,:,p);
        C = C + lambda2 * Beta(p) * HP(:,:,p) * WP(:,:,p);
    end 
    [Hstar, ~]= FMKKM_Opt_orth_st(Hstar, @FMKKM_Opt_fun, option, A,C);
    clear A C


    for p = 1:numker
        A = lambda1 * Alpha(p).^2 * KH(:,:,p);
        B = lambda2 * Hstar' * Beta(p);
        C = WP(:,:,p);
        [HP(:,:,p), ~]= FMKKM_Opt_orth_st(Hstar, @FMKKM_Opt_fun2, option, A,B,C);
    end
    clear A B C
    
    [Gamma] = FMKKM_update_gamma(KH,Hstar);
    [Alpha] = FMKKM_update_alpha(KH,HP);
    [Beta] = FMKKM_update_beta(Hstar,HP,WP);
    
    [obj(iter),obj1(iter),obj2(iter),obj3(iter)] = FMKKM_Objectve(KH,Hstar,HP,WP,Alpha,Beta,Gamma,lambda1,lambda2);
    
    
    if iter >= 2 && abs((obj(iter)-obj(iter-1))/obj(iter))<1e-5
        flag = 0;
    end
end
end

