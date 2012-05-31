function [Sigma,X,M,Clusters]=CreateInitializationClusters(X,M,p)
S=length(X);
n=size(X{1},1);
Sigma=cell(S,1);
SigmaVec=zeros(S,n^2);

for k=1:S
    X{k}=X{k}-repmat(mean(X{k},2),1,p);
    M{k}=M{k}-repmat(mean(M{k},2),1,p);
    Sigma{k}=X{k}*X{k}';
    SigmaVec(k,:)=Sigma{k}(:)';
end
if S>1
    d = pdist(SigmaVec);
    z = linkage(d);
    c = cluster(z,'maxclust',round(S/4));
    Nclust=length(unique(c));
    MaxSize=3;
    kClust=1;
    Clusters=zeros(length(c),MaxSize);
    for k=1:Nclust
        ClustIndx=find(c==k);
        if length(ClustIndx)>MaxSize
            Nparts=ceil(length(ClustIndx)/MaxSize);
           for i=1:(Nparts-1)
               PartIndx=ClustIndx((i-1)*MaxSize+1:i*MaxSize);
                Clusters(kClust,1:MaxSize)=PartIndx;
               kClust=kClust+1;
           end
           PartIndx=ClustIndx((Nparts-1)*MaxSize+1:end);
           Clusters(kClust,1:length(PartIndx))=PartIndx;
           kClust=kClust+1;
        else
            Clusters(kClust,1:length(ClustIndx))=ClustIndx;
            kClust=kClust+1;
        end
    end
    NClustNew=kClust-1;
Clusters=Clusters(1:NClustNew,:);
else
    MaxSize=1;
    Clusters=1;
end