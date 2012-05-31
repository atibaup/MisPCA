function [ImageLoadings,GenesCluster,Err]=ClusterGeneResponses(T,dopt_MisPCA,F_MisPCA,ChosenFactors,Load,Genes,GenesPv,X,Masks,Options,Lambdas,TimeVec,FigOnOff)
Nclust=length(unique(T));
ImageLoadings=[];
ClusterLegend=[];
GenesCluster=cell(Nclust,1);
for i=1:Nclust
    GenesCluster{i}.Ids=find(T==i);
    GenesCluster{i}.Names=Genes(GenesCluster{i}.Ids);
    GenesShown=GenesCluster{i}.Ids(1:min(100,end));
    if FigOnOff
        [Err,h]=OriginalVSReconstructed({dopt_MisPCA{ChosenFactors}},F_MisPCA(:,ChosenFactors),X,X,Masks,...
            Options.Padding,GenesShown,1:Options.S,1,Options.SxSubjectIDs,strcat('Cluster #',num2str(i),'- ',num2str(length(GenesCluster{i}.Ids)),' members'),...
            Genes(GenesShown),ones(Options.S,1),Lambdas(ChosenFactors),TimeVec,GenesPv(GenesShown));
    else
        Err=-1;
    end
    ImageLoadings=[ImageLoadings Load(:,(T==i))];
    ClusterLegend=[ClusterLegend i*ones(1,length(find(T==i)))];
end
MinVal=min(min(ImageLoadings));
MaxVal=max(max(ImageLoadings));

ImageLoadings=[ImageLoadings; (MaxVal-MinVal)*1/max(ClusterLegend-mean(ClusterLegend))*(ClusterLegend-mean(ClusterLegend))];

figure
imagesc(ImageLoadings);