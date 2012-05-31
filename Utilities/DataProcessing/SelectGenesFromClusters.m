function [Genes]=SelectGenesFromClusters(T,ChosenClusters,Genes)
GeneList=[];
T=T(:);
for i=1:length(ChosenClusters)
    GeneList=[GeneList; find(T==ChosenClusters(i))];
end
Genes=Genes(GeneList);