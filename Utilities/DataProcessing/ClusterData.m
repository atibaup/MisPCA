function [T,Centroid,GeneClusters,GeneRowIds,h,Results]=ClusterData(Nclust,GeneNames,X,linkage,distance,TimeVec,figOnOff)

S=length(X);
h=[];
p=size(X{1},2);
n=size(X{1},1);
Data=[];
for i=1:S
    Data=[Data, 1/max(max(X{i}))*X{i}'];
    X{i}=1/max(max(X{i}))*X{i};
end

%% Clustering

T = clusterdata(Data,'maxclust',Nclust,'distance',distance,'linkage',linkage); % 'median' seuclidean correlation
clims=[0,1];
[Ts SrtIndx]=sort(T);
if figOnOff
h(1)=figure('Name','Loading Clusters');
ax=subplot(1,2,1);
imagesc(Ts);
ax=subplot(1,2,2);
imagesc(Data(SrtIndx,:))
end

%% Waveforms for each cluster
Centroid=cell(1,Nclust);
CentroidL=cell(1,Nclust);
CentroidU=cell(1,Nclust);
CentroidOr=cell(1,Nclust);
ax=cell(1,Nclust);
if figOnOff
h(2)=figure('Name','Cluster Signatures');
end
Maxy=0;
Miny=+inf;
for i=1:Nclust
    Indx=find(T==i);
    Centroid{i}=zeros(n,1);
    CentroidL{i}=zeros(n,1);
    CentroidU{i}=zeros(n,1);
    CentroidOr{i}=zeros(n,1);
    Aux=[];
    for j=1:S
%         if i==3 
%             AdjustedDelays=dhat{j}; %mod(dhat{j}-min(dhat{j})+1,16); %ceil(mean(dhat{j}))
%             NAdjX=1/max(max((RecoveredX(MeanA_o,MeanA_no,F,M_no,AdjustedDelays))))*RecoveredX(MeanA_o,MeanA_no,F,M_no,AdjustedDelays);
%            	AdjustedDelays=dhat{j}-min(dhat{j})+1; %mod(dhat{j}-min(dhat{j})+1,16); %ceil(mean(dhat{j}))
%             AdjX=1/max(max((RecoveredX(MeanA_o,MeanA_no,F,M_no,AdjustedDelays))))*RecoveredX(MeanA_o,MeanA_no,F,M_no,AdjustedDelays);
%             figure(h6)
%             subplot(2,1,1)
%             plot(NAdjX); hold on;
%             title('Recovered')
%             subplot(2,1,2)
%             plot(AdjX); hold on;
%             title('Time-Aligned Recovered')
%         end
        Aux=[Aux 1/(max(mean(X{j}(:,Indx),2)))*mean(X{j}(:,Indx),2)];
%         figure(h6)
%         subplot(Nclust,1,i)
%         plot(AdjX); hold on;
%         CentroidOr{i}=CentroidOr{i}+1/(max(mean(X{j}(:,Indx),2)))*mean(X{j}(:,Indx),2);
    end
    Centroid{i}=mean(Aux,2);
    VarR{i}=std(Aux,0,2);
    CentroidU{i}=Centroid{i}+VarR{i};
    CentroidL{i}=Centroid{i}-VarR{i};
    CentroidOr{i}=1/(S)*CentroidOr{i};
    Maxy=max(Maxy,max(CentroidU{i}));
    Miny=min(Miny,min(min(CentroidL{i})));
    str2{i}(1)={''};
    if figOnOff
    ax2{i}=subplot(2,ceil(Nclust/2),i);
    plot(TimeVec,Centroid{i},'-r','LineWidth',2); hold on;
    %plot(CentroidOr{i},'-.b','LineWidth',.5); hold on;
    xlabel('Time','FontSize',12,'FontWeight','bold')
    ylabel('Expression level','FontSize',12,'FontWeight','bold')
    plot(TimeVec,CentroidU{i},'o-g'); hold on;
    hplot{i}=plot(TimeVec,CentroidL{i},'.-b'); hold on;
    end
    kl=1;
    for l=1:min(length(Indx),8)
        if length(strfind(GeneNames{Indx(l)},'AFFX'))==0
            str2{i}(l)={GeneNames{Indx(l)}};
        end       
    end
    if length(Indx)>5
       str2{i}(l+1)={strcat('And +',num2str(length(Indx)-5),' more.')};
    end
    l=l;
    if figOnOff
    title(strcat('Cluster:',num2str(i),' (R) \sigma=',num2str(mean(VarR{i}),1)),'FontSize',12,'FontWeight','bold');
    x = get(hplot{i},'XData'); 		% Get the plotted data
    y = get(hplot{i},'YData');
    imin = find(min(y) == y);		% Find the index of the min and max
    imax = find(max(x) == x);
    text(max(x),min(y),str2{i},'HorizontalAlignment','right','FontSize',8,'FontWeight','bold');
    end
end

Results.Centroid=Centroid;
Results.CentroidL=CentroidL;
Results.CentroidU=CentroidU;
Results.VarR=VarR;
Miny=Miny;
Maxy=Maxy;
if figOnOff
for i=1:Nclust
    subplot(ax2{i});
    axis([min(TimeVec) max(TimeVec) .95*Miny 1.05*Maxy]);

end
suptitle('Mean Recovered Expression Signatures for each cluster');
end
%% Output gene clusters
for i=1:Nclust
    GeneClusters{i}=GeneNames(T==i);
    GeneRowIds{i}=find(T==i);
end
    

