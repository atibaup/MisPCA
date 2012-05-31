function []=Plot_MisPCAvsPCA_Clusters(TMisPCA,TRaw,ResultsAlignedData,ResultsData,A_o,A_PCA,TimeVec,OrderClustersOPFA,OrderClustersRaw,Options)
h=figure;
nclust=length(OrderClustersOPFA);
Col1=winter(4);
Col2=autumn(4);
colors=[Col1(1:ceil(nclust/2),:);  Col2(1:(nclust/2),:);];
MarkerChoice=['xso*+v^'];

Indxplot1=[];
Indxplot2=[];
for i=1:nclust
    Indxplot1=[Indxplot1, 2+(i-1)*4];
    Indxplot2=[Indxplot2, 3+(i-1)*4];
end
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
    %%
    Mean1=ResultsAlignedData{OrderClustersOPFA(i)}.AverageSign(ObsWndw);
    U1=ResultsAlignedData{OrderClustersOPFA(i)}.UCi(ObsWndw);
    L1=ResultsAlignedData{OrderClustersOPFA(i)}.LCi(ObsWndw);
    plot(TimeVec(ObsWndw),Mean1,'-','LineWidth',2,'Color',colors(i,:),'Marker',MarkerChoice(i),'MarkerSize',7); hold on;
    plot(TimeVec(ObsWndw),U1,'^-','Color',colors(i,:)); hold on;
    plot(TimeVec(ObsWndw),L1,'v-','Color',colors(i,:)); hold on;
    if i==nclust
    xlabel('Time (hours)','FontSize',14,'FontWeight','bold')
    end
    axis([min(TimeVec(ObsWndw)) max(TimeVec(ObsWndw)) .5 1.1]);
    %title(strcat('Cluster:',num2str(i),' '),'FontSize',12,'FontWeight','bold');
    
    %%
    IndxClustOPFA=(TMisPCA==OrderClustersOPFA(i));
    figure(h)
    subplot(nclust,4,Indxplot1);
    %p1=plot3(y_OPFA(IndxClustOPFA,1),y_OPFA(IndxClustOPFA,2),y_OPFA(IndxCl
    %ustOPFA,3),'.','Color',colors(i,:),'MarkerSize',7); hold on;

    for s=1:Options.S
    %plot(A_o{s}(1,IndxClustOPFA),A_o{s}(2,IndxClustOPFA),'.','Color',colors(i,:),'Marker',MarkerChoice(i),'MarkerSize',7); hold on;  
        plot3(A_o{s}(1,IndxClustOPFA),A_o{s}(2,IndxClustOPFA),A_o{s}(3,IndxClustOPFA),'.','Color',colors(i,:),'Marker',MarkerChoice(i),'MarkerSize',7); hold on;
    end

    %plot(GeneTraj(1),GeneTraj(2),'.','Color',colors(i,:),'Marker',MarkerChoice(i),'MarkerSize',7); hold on;
    xlabel('1st Factor','FontSize',14,'FontWeight','bold')
    ylabel('2nd Factor','FontSize',14,'FontWeight','bold')
    grid on
    subplot(nclust,4,Indxplot2)
    IndxClust=(TRaw==OrderClustersRaw(i));
    %p2=plot3(y(IndxClust,1),y(IndxClust,2),y(IndxClust,3),'.','Color',colors(i,:),'MarkerSize',7); hold on;
    %p2=plot(y(IndxClust,1),y(IndxClust,2),'.','Color',colors(i,:),'MarkerSize',7); hold on;
    for s=1:Options.S
        %p1=plot(y_OPFA(IndxClustOPFA,1),y_OPFA(IndxClustOPFA,2),'.','Color',colors(i,:),'MarkerSize',7); hold on;
        %p1=plot3(A_PCA{s}(Factors(1),IndxClust),A_PCA{s}(Factors(2),IndxClust),A_PCA{s}(Factors(3),IndxClust),'.','Color',colors(i,:),'MarkerSize',7); hold on;
        %plot(A_PCA{s}(1,IndxClust),A_PCA{s}(2,IndxClust),'.','Color',colors(i,:),'Marker',MarkerChoice(i),'MarkerSize',7); hold on;
        plot3(A_PCA{s}(1,IndxClust),A_PCA{s}(2,IndxClust),A_PCA{s}(3,IndxClust),'.','Color',colors(i,:),'Marker',MarkerChoice(i),'MarkerSize',7); hold on;
    end
    xlabel('1st PC','FontSize',14,'FontWeight','bold')
    ylabel('1st PC','FontSize',14,'FontWeight','bold')
    grid on
    %%
    ax2{i}=subplot(nclust,4,4+(i-1)*4);
    Mean2=ResultsData{OrderClustersRaw(i)}.AverageSign(ObsWndw);
    U2=ResultsData{OrderClustersRaw(i)}.UCi(ObsWndw);
    L2=ResultsData{OrderClustersRaw(i)}.LCi(ObsWndw);

    plot(TimeVec(ObsWndw),Mean2,'-','LineWidth',2,'Color',colors(i,:),'Marker',MarkerChoice(i),'MarkerSize',7); hold on;
    plot(TimeVec(ObsWndw),U2,'^-','Color',colors(i,:)); hold on;
    plot(TimeVec(ObsWndw),L2,'v-','Color',colors(i,:)); hold on;
    if i==nclust
    xlabel('Time (hours)','FontSize',14,'FontWeight','bold')
    end
    axis([min(TimeVec(ObsWndw)) max(TimeVec(ObsWndw)) .5 1.1]);
    %title(strcat('Cluster:',num2str(i),' '),'FontSize',12,'FontWeight','bold');
    Maxy(i)=max(Maxy(i),max(U1));
    Miny(i)=min(Miny(i),min(L1));
    Maxy2(i)=max(Maxy(i),max(U2));
    Miny2(i)=min(Miny(i),min(L2));
%%       
    %%

end
figure(h)
for i=1:nclust
    subplot(ax1{i});
    axis([min(TimeVec) max(TimeVec) .95*Miny(i) 1.05*Maxy(i)]);
    subplot(ax2{i});
    axis([min(TimeVec) max(TimeVec) .95*Miny2(i) 1.05*Maxy2(i)]);
end

