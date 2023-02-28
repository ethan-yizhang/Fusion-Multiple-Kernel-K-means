function ImputedKM = algorithm_knn_imputation(K,S,k)
% KM is the Missing kernels
% S is the index of missing samples
% k is the k of k-nn
% ImputedKM is the imputed kernels
% ImputedKM is of the same size of KM
KM = algorithm_zero_imputation(K,S);
ImputedKM = KM;
numker = size(KM,3);
for p = 1 : numker
    % ��Ե�i��kernel���в�ȫ
    Indx = setdiff(1:numker,p);
    % ��i��kernel��ȱʧ�����
    MisingIndex = S{p}.indx;
    % ����knn
    TempK = sum(KM(:,:,Indx),3)/(numker-1);
    IDX = knnsearch(TempK,TempK(MisingIndex ,:) ,'k',k+1);
    % ����ȫ
    for j = 1 : length(MisingIndex)
        heheTemp = mean(TempK(IDX(j,2:k+1) , :) , 1);
        ImputedKM(MisingIndex(j),:,p) = heheTemp;
        ImputedKM(:,MisingIndex(j),p) = heheTemp;
    end
    ImputedKM(:,:,p) = (ImputedKM(:,:,p)+ImputedKM(:,:,p)')/2;
end