function [A_misPCA,A_PCA,X,Genes,GenesAffy,SortedPv]=ConcatLoadings(Results,PCAResults,Options,X,Genes,GenesAffy,SortedPv,newp,AverageOnOff)

A_misPCA=[];
A_PCA=[];
AverageA_misPCA=zeros(size(Results.Aaligned{2},1),Options.p);
AverageA_PCA=zeros(size(PCAResults.Aaligned{2},1),Options.p);
for i=1:Options.S
    
    % normalize
    Aligned1=[];
    Aligned2=[];
    
    for j=1:Options.p
        Aligned1=[Aligned1,1/max(abs(Results.Aaligned{i}(:,j)))*Results.Aaligned{i}(:,j)];
        Aligned2=[Aligned2,1/max(abs(PCAResults.Aaligned{i}(:,j)))*PCAResults.Aaligned{i}(:,j)];
    end
    AverageA_misPCA=AverageA_misPCA+1/Options.S*Aligned1;
    AverageA_PCA=AverageA_PCA+1/Options.S*Aligned2;
    
    
    A_misPCA=[A_misPCA,Aligned1'];
    A_PCA=[A_PCA,Aligned2'];
end

NormAvloadings=diag(AverageA_misPCA'*AverageA_misPCA);

[vals,indxGenes]=sort(NormAvloadings,'descend');

figure; plot(vals); hold on;
stem(newp,max(vals),'-ok');


indxGenes=indxGenes(1:newp);

for i=1:Options.S    
    X{i}=X{i}(:,indxGenes);    
end
if AverageOnOff
 A_misPCA=AverageA_misPCA(:,indxGenes)';
A_PCA=AverageA_PCA(:,indxGenes)';   
else
A_misPCA=A_misPCA(indxGenes,:);
A_PCA=A_PCA(indxGenes,:);
end
Genes=Genes(indxGenes);
GenesAffy=GenesAffy(indxGenes);
SortedPv=SortedPv(indxGenes);