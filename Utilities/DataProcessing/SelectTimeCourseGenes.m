function [Genes,GenesAffy,RowIds,X,Indx,SortedPv]=SelectTimeCourseGenes(Genes,GenesAffy,RowIds,sample,X,y,pTh)

S=length(X);
p=length(Genes);
if p <=5000
    pv=pv_behrens_fisher2(sample.blood.RNA.data(RowIds,:),[2:5],'minT',...
        sample.blood.RNA.imtrx,sample.blood.RNA.excelmat,y,sample.blood.RNA.alltimes,0);
    
    [SortedPv, Indx]=sort(pv,'ascend');
    figure
    plot(SortedPv); hold on;
    %plot(find(SortedPv=pTh,max(SortedPv),'xr');
    title('p-values')
    if pTh <1
        Indx=Indx(SortedPv<pTh);
    else
        Indx=Indx(1:pTh);
    end
    RowIds=RowIds(Indx);
    Genes=Genes(Indx);
    GenesAffy=GenesAffy(Indx);
    for i=1:S
        X{i}=X{i}(:,Indx);
    end
else
    disp('Can t do behren s for more than 1000 genes')
    Indx=1:p;
    SortedPv=(1:p)/p;
end