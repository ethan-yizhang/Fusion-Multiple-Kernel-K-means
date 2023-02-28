clear
clc
warning off;

addpath(genpath('./'));



dataName = 'flower17';
load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
numclass = length(unique(Y));
Y(Y<1)=numclass;
numker = size(KH,3);
num = size(KH,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KH = kcenter(KH);
KH = knorm(KH);
options.seuildiffsigma=1e-5;        % stopping criterion for weight variation
%------------------------------------------------------
% Setting some numerical parameters
%------------------------------------------------------
options.goldensearch_deltmax=1e-1; % initial precision of golden section search
options.numericalprecision=1e-16;   % numerical precision weights below this value are set to zero
%------------------------------------------------------
% some algorithms paramaters
%------------------------------------------------------
options.firstbasevariable='first'; % tie breaking method for choosing the base variable in the reduced gradient method
options.nbitermax=500;             % maximal number of iteration
options.seuil=0;                   % forcing to zero weights lower than this
options.seuilitermax=10;           % value, for iterations lower than this one
options.miniter=0;                 % minimal number of iterations
options.threshold = 1e-4;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qnorm = 2;



%% %%---MKKM--Coordinate --%%%
tic;
[H_normalized2,Sigma2,obj2] = mkkmeans_train(KH,numclass);
[res_mean(:,2),res_std(:,2)] = myNMIACCV2(H_normalized2,Y,numclass);
timecost(2) = toc;


%% FMKKM
tic;
rhoset4 = 2.^[1:1:9];
lambdaset4 = 2.^[3:1:10];

res_allmean4 = zeros(4,length(rhoset4),length(lambdaset4));
res_allstd4 = zeros(4,length(rhoset4),length(lambdaset4));
Sigma1_all4 = zeros(numker,length(rhoset4),length(lambdaset4));
Sigma2_all4 = zeros(numker,length(rhoset4),length(lambdaset4));
Sigma3_all4 = zeros(numker,length(rhoset4),length(lambdaset4));
for ir =1:length(rhoset4)
    for il = 1:length(lambdaset4)
        
        [H_normalized4,Alpha4,Beta4,Gamma4,WP4,obj4] = FMKKM(KH,rhoset4(ir),lambdaset4(il),numclass,options);
        [res_mean4,res_std4]= myNMIACCV2(H_normalized4,Y,numclass);
        res_allmean4(:,ir,il)=res_mean4;
        res_allstd4(:,ir,il)=res_std4;
        Sigma1_all4(:,ir,il)=Alpha4;
        Sigma2_all4(:,ir,il)=Beta4;
        Sigma3_all4(:,ir,il)=Gamma4;        
        fprintf('lambda1: %.2f lambda2: %.2f   Proposed: %.4f  \n',rhoset4(ir),lambdaset4(il),res_mean4(1));
        
    end
end

[~,max_idx]=max(res_allmean4(1,:),[],'all','linear');
res_mean(:,4) = res_allmean4(:,max_idx);
res_std(:,4) = res_allstd4(:,max_idx);
Sigma4 = Sigma3_all4(:,max_idx);
timecost(4) = toc;




