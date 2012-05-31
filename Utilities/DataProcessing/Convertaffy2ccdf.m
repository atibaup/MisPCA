function [GeneListccdf,GeneListaffynew]=Convertaffy2ccdf(GeneListaffy,DataDir)
load(strcat(DataDir,'GeneNames'));
Ngenes=length(GeneListaffy);
k=1;
for i=1:Ngenes
    Aux=find(strcmp(names_RNAaffy,GeneListaffy(i))==1);
    if length(Aux)==1
        GeneListccdf(k)=names_RNAccdf(Aux);
        GeneListaffynew(k)=GeneListaffy(i);
        k=k+1;
    elseif length(Aux)>1
        Aux=Aux
        GeneListccdf(k)=names_RNAccdf(Aux(1));
        GeneListaffynew(k)=GeneListaffy(i);
        k=k+1;
        disp(strcat('Gene found twice: ',GeneListaffy(i),' rows:',num2str(Aux)))
    elseif length(Aux)==0
        disp(strcat('Gene not found: ',GeneListaffy(i)))
    end
end