function [RowIds,Genes]=FindGenesRows(GenesOfInterest,rownames)

RowIds=[];
Nfound=1;
Genes=[];
for i=1:length(GenesOfInterest)
    for j=1:length(rownames)
        if strcmp(GenesOfInterest{i},rownames{j})
            RowIds(Nfound)=j;
            Genes{Nfound}=GenesOfInterest{i};
            Nfound=Nfound+1;
            
        end
    end
end
