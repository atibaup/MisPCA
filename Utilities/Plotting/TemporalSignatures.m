function [Signature,ClusterMembers,Xr,AlignedData,T]=TemporalSignatures(Clusters,F,X,Masks,dopt_MisPCA,Lambdas,Genes,GeneNames,GenesAffy,pValues,Options,TimeVec,ClusterLoadingsOrResponses,FixedDelaysOnOff,ClustersToShow)
n=size(Masks{1},1);
p=length(Genes);
k=size(F,2);
Nclusters=length(unique(Clusters));
Nclust=Nclusters;
S=length(X);
Loadings=zeros(k,p);
Data=zeros(n,S*p);
DataR=zeros(n-2*Options.Padding,S*p);
AlignedData=zeros(n,p);
ClustersAux=[];
nPlots=1;
TimeVecAux=TimeVec;
TimeVec=cell(S,1);
%figure;
Xr=cell(S,1);
Last=0;
FigOnOff=1;
%figure
for i=1:S
    ni=length(find(Masks{i}(:,1)==1));
    TimeIndx{i}=Options.Padding+1:(ni-Options.Padding);
    Indxtime=(Masks{i}(:,1)==1);
    
    TimeVec{i}=TimeVecAux(Indxtime);
    TimeVec{i}=TimeVec{i}(TimeIndx{i});
    Fadj=zeros(n,k);
    SubjDelays=[];
    for j=1:k
       SubjDelays=[SubjDelays,dopt_MisPCA{j}(i)]; 
    end
    Faligned=zeros(n,k);
    for j=1:k
         Fadj(:,j)= circshift(F(:,j),-dopt_MisPCA{j}(i));
         Faligned(:,j)= circshift(F(:,j),-(dopt_MisPCA{j}(1)));%(-dopt_MisPCA{j}(i)+dopt_MisPCA{1}(i))
    end

    Fadj=Fadj*diag(Lambdas.^(1/2));
    Faligned=Faligned*diag(Lambdas.^(1/2));
    FadjOld=Fadj;
    Fadj=Fadj(Indxtime,:);
    Faligned=Faligned(Indxtime,:);
    if i==1
         Fref=FadjOld;
    end    
    PFX=(Fadj'*Fadj)\(Fadj'*X{i}(Indxtime,Genes));
    
    %Loadings=[Loadings,PFX];
    Loadings=Loadings+1/S*PFX;

    if FixedDelaysOnOff
        Xr{i}=Fref(Indxtime,:)*PFX;
    else        
        Xr{i}=Faligned*PFX;
    end
        
    Xr{i}=1/max(max(abs(Xr{i})))*Xr{i};
    
    if FigOnOff
        figure;
        subplot(3,1,1)
        plot(Faligned)
        axis tight;
        subplot(3,1,2)
        pshow=300;
        plot(1/max(max(abs(Xr{i}(:,1:pshow))))*Xr{i}(:,1:pshow))
        axis tight;
        subplot(3,1,3)
        plot(1/max(max(abs(X{i}(Indxtime,1:pshow))))*X{i}(Indxtime,1:pshow)); hold on;
        plot(1/max(max(abs(Fadj)))*Fadj,'LineWidth',3);
        axis tight;
    end
    
    Data(Indxtime,p*(i-1)+1:p*i)=[X{i}(Indxtime,Genes)];
    DataR(Indxtime,p*(i-1)+1:p*i)=[Xr{i}];
    
    AlignedData(Indxtime,:)=AlignedData(Indxtime,:)+1/S*Xr{i};
    
    %AlignedData(Last+1:(Last+length(find(Indxtime))),:)=[Xr{i}];
    %Last=Last+length(find(Indxtime));
    
    nPlots=nPlots+2;
    ClustersAux=[ClustersAux;Clusters];
end

%compute clusters from aligned subspace
if ClusterLoadingsOrResponses==2
    linkage='complete';
    distance='correlation';
    Nclust=Nclusters;  
    T2 = clusterdata(AlignedData','maxclust',Nclust,'distance',distance,'linkage',linkage); % 'median' seuclidean correlation
    [SortedT,IndxSorted]=sort(T2);    
    figure
    subplot(2,1,1)
    imagesc(AlignedData(:,IndxSorted))
    subplot(2,1,2)
    imagesc(SortedT')
    
    ClustersAux2=[];
    for i=1:S
        ClustersAux2=[ClustersAux2;T2];
    end
    ClustersAux=ClustersAux2;
    T=T2;
else
    T=Clusters;
end


fid1=fopen(strcat('Results/Clustering/BackgroundGenes.txt'),'w');
for i=1:length(GenesAffy)
    fprintf(fid1,'%s \n',strrep(GenesAffy{i},'_at',''));
end
fclose(fid1);

ClusterMembers=cell(Nclust,1);
fid=zeros(Nclust,1);
for i=1:Nclust
    fid(i)=fopen(strcat('Results/Clustering/Cluster_',num2str(i),'.txt'),'w');
    
    GeneIndx=find(T==i);
    [orderedpVals,orderedIndx]=sort(pValues(GeneIndx),'ascend');
    ClusterMembers{i}.Ids=GeneIndx(orderedIndx);
    ClusterMembers{i}.Genes=GeneNames(GeneIndx(orderedIndx));
    ClusterMembers{i}.GenesAffy=GenesAffy(GeneIndx(orderedIndx));
    ClusterMembers{i}.pVals=orderedpVals;
    for l=1:length( ClusterMembers{i}.GenesAffy)
        fprintf(fid(i),'%s \n',strrep(ClusterMembers{i}.GenesAffy{l},'_at',''));
    end
    fclose(fid(i));
end


figure('Name','Clusters on Aligned subspace')
Signature=cell(length(ClustersToShow),1);
TimeVecAux=Options.Padding+1:(Options.n-Options.Padding);
for i=1:length(ClustersToShow)
    
    Signature{i}.Data=DataR(:,ClustersAux==ClustersToShow(i));
    Signature{i}.DataRaw=Data(:,ClustersAux==ClustersToShow(i));
    subplot(length(ClustersToShow),1,i)
    plot(Signature{i}.Data(TimeVecAux,:),'.r','LineWidth',1,'MarkerSize',1); hold on;
    
    % stats computation (there're missing vals so we need to take this
    % into account)
    NZ_Points=Signature{i}.Data(:,abs(Signature{i}.Data(k,:))>0)';
    AverageSign=mean(NZ_Points,1);
    NZ_PointsRaw=Signature{i}.DataRaw(:,abs(Signature{i}.Data(k,:))>0)';
    AverageSignRaw=mean(NZ_PointsRaw,1);

    StdDev=std(NZ_Points,0,1);
    StdDevRaw=std(NZ_PointsRaw,0,1);
    UCi=AverageSign(:)+StdDev(:);
    LCi=AverageSign(:)-StdDev(:);
    Signature{i}.AverageSign=AverageSign;
    Signature{i}.UCi=UCi;
    Signature{i}.LCi=LCi;
    Signature{i}.AverageSignRaw=AverageSignRaw;
    Signature{i}.UCiRaw=AverageSignRaw(:)+StdDevRaw(:);
    Signature{i}.LCiRaw=AverageSignRaw(:)-StdDevRaw(:);

    plot(AverageSign(TimeVecAux),'-k','LineWidth',3); hold on;
    plot(UCi(TimeVecAux),'--^b','LineWidth',1); hold on;
    plot(LCi(TimeVecAux),'--vb','LineWidth',1); hold on;
    minY=1.3*min(min(Signature{i}.LCi));
    maxY=1.3*max(max(Signature{i}.UCi));
    axis([1, (Options.n-2*Options.Padding),minY,maxY]);
    ClusterTitle='';
    for j=1:min(10,length(ClusterMembers{i}.Genes))
        if j==1
            ClusterTitle=[ClusterMembers{i}.Genes{j}];
        else
            ClusterTitle=[ClusterTitle,'-',ClusterMembers{i}.Genes{j}];
            if mod(j,4)==0
                ClusterTitle=[ClusterTitle,10];
            end
        end
    end
    if length(ClusterMembers{i}.Genes)>10
        ClusterTitle=[ClusterTitle,'- And ',num2str(length(ClusterMembers{i}.Genes)-10),' other'];
    end
    ClusterTitle='';
    title(strcat('Cl-',num2str(i),' (',ClusterTitle,')'),'FontSize',10,'FontWeight','bold')
end
% SignatureO=cell(Nclust,1);
% figure('Name','Clusters on Observed space')
% for i=1:Nclust
%     SignatureO{i}.Data=Data(:,ClustersAux==i);
%     subplot(Nclust,1,i)
%     AverageSign=zeros(n,1);
%     VarSignature=zeros(n,1);
%     StdSignature=zeros(n,1);
%     for k=1:n
%         NZ_Points=Signature{i}.Data(k,abs(SignatureO{i}.Data(k,:))>0);
%         AverageSign(k)=mean(NZ_Points);
%         VarSignature(k)=var(NZ_Points);
%         StdSignature(k)=sqrt(VarSignature(k));
%     end
%     plot(TimeVecAux,SignatureO{i}.Data,'.r','LineWidth',1,'MarkerSize',1); hold on;
%     plot(TimeVecAux,AverageSign,'-r','LineWidth',3); hold on;
%     plot(TimeVecAux,AverageSign+StdSignature,'.-b','LineWidth',3); hold on;
%     plot(TimeVecAux,AverageSign-StdSignature,'.-g','LineWidth',3); hold on;
%     axis([min(TimeVecAux), max(TimeVecAux), min(AverageSign-StdSignature)*1.1,max(AverageSign+StdSignature)*1.1]);
% end