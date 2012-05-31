function [GeneListaffy,GeneListccdfnew]=Convertccdf2affy(GeneListccdf,DataDir)
load(strcat(DataDir,'GeneNames'));
Ngenes=length(GeneListccdf);
k=1;
for i=1:Ngenes
    Aux=find(strcmp(names_RNAccdf,GeneListccdf(i))==1);
    if length(Aux)==1
        GeneListaffy(k)=names_RNAaffy(Aux);
        GeneListccdfnew(k)=GeneListccdf(i);
        k=k+1;
    elseif length(Aux)>1
        Aux=Aux
        GeneListaffy(k)=names_RNAaffy(Aux(1));
        GeneListccdfnew(k)=GeneListccdf(i);
        k=k+1;
        disp(strcat('Gene found twice: ',GeneListccdf(i),' rows:',num2str(Aux)))
    elseif length(Aux)==0
        disp(strcat('Gene not found: ',GeneListccdf(i)))
    end
end