function [IntGenes]=InterestingGenesFromClusters(IntGenesNames,GeneClusters,nclust,Genes)

NIntGenes=length(IntGenesNames);
k=1;
for i=1:nclust
    for j=1:min(10,length(GeneClusters{i}))
        IntGenesNames{NIntGenes+k}=GeneClusters{i}{j};
        k=k+1;
    end
end
[IntGenes]=FindGenesRows(IntGenesNames,Genes);
IntGenes=IntGenes(IntGenes>0);