function h=Plot_MisPCAClusters(TMisPCA,TPCA,ResultsAlignedData,ResultsMisAlignedData,A_misPCA,A_PCA,TimeVec,OrderClustersMisPCA,OrderClustersPCA,ClusterMembersMisPCA,ClusterMembersPCA,Options)
h=figure;
nclust=length(OrderClustersMisPCA);
Col1=cool(nclust);
Col2=hot(nclust);
colors=[Col1(1:floor(nclust/2),:);  Col2(1:ceil(nclust/2),:);];
%colors=hsv(nclust);
MarkerChoice=['^^^vvv'];

Indxplot1=[];
Indxplot2=[];
for i=1:nclust
    Indxplot1=[Indxplot1, 2+(i-1)*4];
    Indxplot2=[Indxplot2, 3+(i-1)*4];
end

%% Compute 2-D MDS embeddings for 
D_misPCA=pdist(A_misPCA);
D_PCA=pdist(A_PCA);
A_MisPCA_mds = MDSCALE(D_misPCA,2);
A_PCA_mds = MDSCALE(D_PCA,2);
A_MisPCA_mds =A_MisPCA_mds*diag(1./max(abs(A_MisPCA_mds)));
A_PCA_mds =A_PCA_mds*diag(1./max(abs(A_PCA_mds)));

minA1=1.1*min(min(A_MisPCA_mds(:,1)),min(A_PCA_mds(:,1)));
maxA1=1.1*max(max(A_MisPCA_mds(:,1)),max(A_PCA_mds(:,1)));
minA2=1.1*min(min(A_MisPCA_mds(:,2)),min(A_PCA_mds(:,2)));
maxA2=1.1*max(max(A_MisPCA_mds(:,2)),max(A_PCA_mds(:,2)));
%% Plot centroids and MDS coordinates
ObsWndw=Options.Padding+1:(Options.n-Options.Padding);
Maxy=zeros(nclust,1);
Miny=inf*ones(nclust,1);
Maxy2=zeros(nclust,1);
Miny2=inf*ones(nclust,1);
ax1=cell(nclust,1);
ax2=cell(nclust,1);
for i=1:nclust
    figure(h)
    ax1{i}=subplot(nclust,4,1+(i-1)*4);
    
    %% Plot centroids os MisPCA-based clusters
    Mean1=ResultsAlignedData{OrderClustersMisPCA(i)}.AverageSign(ObsWndw);
    U1=ResultsAlignedData{OrderClustersMisPCA(i)}.UCi(ObsWndw);
    L1=ResultsAlignedData{OrderClustersMisPCA(i)}.LCi(ObsWndw);
    plot(TimeVec(ObsWndw),Mean1,'-','LineWidth',2,'Color',colors(i,:),'Marker',MarkerChoice(i),'MarkerSize',7); hold on;
    plot(TimeVec(ObsWndw),U1,'^-','Color',colors(i,:)); hold on;
    plot(TimeVec(ObsWndw),L1,'v-','Color',colors(i,:)); hold on;
    VarclustMisPCA(i)=mean(var(ResultsAlignedData{OrderClustersMisPCA(i)}.Data,1,2));
    text(max(TimeVec(ObsWndw))-50,max(Mean1),['{\sigma}^2_',num2str(i),'=',num2str(VarclustMisPCA(i),3)],'FontSize',10)
    if i==nclust
    xlabel('Time (hours)','FontSize',14,'FontWeight','bold')
    end
    axis([min(TimeVec(ObsWndw)) max(TimeVec(ObsWndw)) .5 1.1]);
        
    %% Plot coordinates on 2-D MDS space for MisPCA and PCA coordinates
    
    IndxClustMisPCA=(TMisPCA==OrderClustersMisPCA(i));
    IndxClustPCA=(TPCA==OrderClustersPCA(i));
    figure(h);
    subplot(nclust,4,Indxplot1);
    plot(A_MisPCA_mds(IndxClustMisPCA,1),A_MisPCA_mds(IndxClustMisPCA,2),'.','Color',colors(i,:),'Marker',MarkerChoice(i),'MarkerSize',7); hold on;
    axis([minA1,maxA1,minA2,maxA2]);
    xlabel('1st MDS PC','FontSize',14,'FontWeight','bold');
    ylabel('2nd MDS PC','FontSize',14,'FontWeight','bold');
    
    grid on;
    subplot(nclust,4,Indxplot2)
    plot(A_PCA_mds(IndxClustPCA,1),A_PCA_mds(IndxClustPCA,2),'.','Color',colors(i,:),'Marker',MarkerChoice(i),'MarkerSize',7); hold on;
    axis([minA1,maxA1,minA2,maxA2]);
    xlabel('1st MDS PC','FontSize',14,'FontWeight','bold');
    ylabel('2nd MDS PC','FontSize',14,'FontWeight','bold');
    grid on;
    
    %% Plot Centroids of PCA-based clustering
    ax2{i}=subplot(nclust,4,4+(i-1)*4);
    Mean2=ResultsMisAlignedData{OrderClustersPCA(i)}.AverageSignRaw(ObsWndw);
    U2=ResultsMisAlignedData{OrderClustersPCA(i)}.UCiRaw(ObsWndw);
    L2=ResultsMisAlignedData{OrderClustersPCA(i)}.LCiRaw(ObsWndw);

    plot(TimeVec(ObsWndw),Mean2,'-','LineWidth',2,'Color',colors(i,:),'Marker',MarkerChoice(i),'MarkerSize',7); hold on;
    plot(TimeVec(ObsWndw),U2,'^-','Color',colors(i,:)); hold on;
    plot(TimeVec(ObsWndw),L2,'v-','Color',colors(i,:)); hold on;
    VarclustPCA(i)=mean(var(ResultsMisAlignedData{OrderClustersPCA(i)}.Data,1,2));
    text(max(TimeVec(ObsWndw))-50,max(Mean1),['{\sigma}^2_',num2str(i),'=',num2str(VarclustPCA(i),3)],'FontSize',10)

    if i==nclust
    xlabel('Time (hours)','FontSize',14,'FontWeight','bold')
    end
    axis([min(TimeVec(ObsWndw)) max(TimeVec(ObsWndw)) .5 1.1]);

    Maxy(i)=max(Maxy(i),max(max(U1),max(U2)));
    Miny(i)=min(Miny(i),min(min(L1),min(L2)));
    Maxy2(i)=max(Maxy(i),max(U2));
    Miny2(i)=min(Miny(i),min(L2));
end


%% Adjust axis
figure(h)
for i=1:nclust
    subplot(ax1{i});
    axis([min(TimeVec) max(TimeVec) 1.05*Miny(i) 1.05*Maxy(i)]);
    %% Annotate MisPCA centroids
    Ngenes=length(ClusterMembersMisPCA{OrderClustersMisPCA(i)}.Genes);
    text(min(TimeVec)-100,1.1*Maxy(i),[num2str(Ngenes),' genes:'],'FontSize',10);
    for j=1:min(5,length(ClusterMembersMisPCA{OrderClustersMisPCA(i)}.Genes))
         txtstr=[ClusterMembersMisPCA{OrderClustersMisPCA(i)}.Genes{j},'(',num2str(ClusterMembersMisPCA{OrderClustersMisPCA(i)}.pVals(j),'%0.1g'),')'];
         text(min(TimeVec)-99,1.05*Maxy(i)-j/6*(1.05*Maxy(i)-1.05*Miny(i)),txtstr,'FontSize',10);
    end
    text(min(TimeVec)-99,1.05*Maxy(i)-6/6*(1.05*Maxy(i)-1.05*Miny(i)),'...','FontSize',10);
    subplot(ax2{i});
    axis([min(TimeVec) max(TimeVec) .95*Miny2(i) 1.05*Maxy2(i)]);
    %% Annotate PCA centroids
    Ngenes=length(ClusterMembersPCA{OrderClustersPCA(i)}.Genes);
    text(max(TimeVec)+10,1.1*Maxy(i),[num2str(Ngenes),' genes:'],'FontSize',10);
    for j=1:min(5,length(ClusterMembersPCA{OrderClustersPCA(i)}.Genes))
         txtstr=[ClusterMembersPCA{OrderClustersPCA(i)}.Genes{j},'(',num2str(ClusterMembersPCA{OrderClustersPCA(i)}.pVals(j),'%0.1g'),')'];
         text(max(TimeVec)+11,1.05*Maxy(i)-j/6*(1.05*Maxy(i)-1.05*Miny(i)),txtstr,'FontSize',10);
    end
    text(max(TimeVec)+11,1.05*Maxy(i)-6/6*(1.05*Maxy(i)-1.05*Miny(i)),'...','FontSize',10);
end

subplot(ax1{nclust});
text(mean(TimeVec), 1.05*Miny(i)-10,['Average Centroid Variance=',num2str(mean(VarclustMisPCA),1)],'FontSize',10);

subplot(ax2{nclust});
text(mean(TimeVec), 1.05*Miny(i)-10,['Average Centroid Variance=',num2str(mean(VarclustPCA),1)],'FontSize',10);

CentroidVariances=['Average Centroid Variance: MisPCA=',num2str(mean(VarclustMisPCA),1),'/ PCA=',num2str(mean(VarclustPCA),3)]


