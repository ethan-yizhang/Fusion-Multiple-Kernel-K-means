clear
clc
warning off;

path = 'D:\work2015\';
addpath(genpath('./'));



dataNameSet =  {'3Sources_3view','bbcsport_2view','Cal-5views','Caltech101-7','football_3view','football_9view',...
    'LandUse-21','MF-first3view','MSRA-6view','olympics_3view','olympics_9view','politicsie_3view','politicsie_9view',...
    'YALE','washington','wpbc','wdbc','bbcsport2view','bbcsport', 'AR10P','cornell',...
    'politicsuk_3view','politicsuk_9view','rugby_3view','rugby_9view','scene15_mtv','still-2-mtv','willow-mtv',...
    'texas','wisconsin','heart','ionosphere','sonar','ORL','pima','liver',...
    'citeseer','cora','Caltech101-20','4Area-Tao','flower17','proteinFold','SensITVehicle_1500sample_2view_3cluster','UCI_DIGIT','caltech101_nTrain5_48',...
    'caltech101_nTrain10_48','caltech101_nTrain15_48',...
    'plant','psortPos','psortNeg','nonpl','segment_2310',...
    'CCV','caltech101_nTrain20_48','caltech101_nTrain25_48','caltech101_nTrain30_48','flower102',...
    'ALOI-100', 'caltech101', 'caltech101_mit', 'Handwritten_numerals', 'MNIST_5000' ,'Reuters','YoutubeFace_sel_fea_10','NUS-WIDE-OBJECT_fea','SUNRGBD_fea','AwA_fea'};

resultpath = [path,'work2022\FMKKM\result2023\'];
mkdir(resultpath);
% logname = ['./V2demo.txt'];
% diary(logname);
% diary on;

for id = [1]% :length(dataNameSet)
    
    dataName = dataNameSet{id}
%     load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
    if id<65%length(dataNameset)-1
        load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
    else
        load([path,'datasets\',dataName,'_Kmatrix'],'K','Y');
        for i = 1:size(K,1)
            KH(:,:,i) = K{i};
        end
        clear K
    end
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
    options.numericalprecision=1e-16;   % numerical precision weights below this value
    % are set to zero
    %------------------------------------------------------
    % some algorithms paramaters
    %------------------------------------------------------
    options.firstbasevariable='first'; % tie breaking method for choosing the base
    % variable in the reduced gradient method
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
    rhoset4 = 2.^[-5:1:10];
    lambdaset4 = 2.^[0:1:10];
    rhoset4 = 2.^[1:1:9];
    lambdaset4 = 2.^[3:1:10];
    rhoset4 = 2.^[1:1:2];
    lambdaset4 = 2.^[3:1:4];
    
    res_allmean4 = zeros(4,length(rhoset4),length(lambdaset4));
    res_allstd4 = zeros(4,length(rhoset4),length(lambdaset4));
    Sigma1_all4 = zeros(numker,length(rhoset4),length(lambdaset4));
    Sigma2_all4 = zeros(numker,length(rhoset4),length(lambdaset4));
    Sigma3_all4 = zeros(numker,length(rhoset4),length(lambdaset4));
    for ir =1:length(rhoset4)
        for il = 1:length(lambdaset4)[H_normalized4,Alpha4,Beta4,Gamma4,WP4,obj4] = FMKKM(KH,rhoset4(ir),lambdaset4(il),numclass,options);
    
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
    


    save([resultpath,dataName,'_FMKKM_clustering.mat'],'res_mean','res_std','timecost',...
        'Sigma4','res_allmean4','res_allstd4','Sigma1_all4','Sigma2_all4','Sigma3_all4','max_idx');
    
    clear KH Y WP3 H_normalized4 Sigma4 H_normalized5 WP5 Sigma5
end