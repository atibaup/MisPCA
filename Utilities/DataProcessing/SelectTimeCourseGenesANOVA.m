function [Genes,GenesAffy,RowIds,X,Masks,SortedPv]=SelectTimeCourseGenesANOVA(Genes,GenesAffy,RowIds,Times,X,Masks,Options)

S=length(X);
nvars=size(X{1},2);
n=size(Masks{1},1); 
p=zeros(nvars,1);
pTh=Options.p;
F=zeros(length(Genes),1);

fprintf('\n ANOVA gene selection ...')
for i=1:nvars
    if mod(i,500)==0
        fprintf('\n \t  Processing gene %d/%d ...',i,nvars)
    end
    Xconcat=zeros(S*size(Times,2),size(Times,1));
    IndxPadding=zeros(n,1);
    IndxPadding((Options.Padding+1):(n-Options.Padding))=1;
    for k=1:size(Times,1)
        aux=[];
        for s=1:S
            Signal=zeros(n,1);
            Indx=((Masks{s}(:,1)==1)&IndxPadding);
            Signal(Indx,:)=X{s}(Indx,i);
            aux=[aux;Signal(Times(k,:))];
        end
        Xconcat(:,k)=aux;
    end    
    [p(i),table]=anova1(Xconcat,[],'off');
    F(i)=table{2,5};
end

[SortedPv,I] = sort(p,'ascend');
GenesOld=Genes;
GenesAffyOld=GenesAffy;
RowIdsOld=RowIds;
I=I(1:pTh);
clear Genes GenesAffy RowIds
for i=1:pTh
    Genes{i}=GenesOld{I(i)};
    GenesAffy{i}=GenesAffyOld{I(i)};
    RowIds(i)=RowIdsOld(I(i));
end
SortedPv=SortedPv(1:pTh);

for j=1:S
    X{j}=X{j}(:,I);
    X{j}=1/max(max(abs(X{j})))*X{j};
    Masks{j}=Masks{j}(:,I);
end
