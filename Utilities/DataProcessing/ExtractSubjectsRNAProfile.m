function [X,Mask,AllTimes,Samples,y,RowIds,Genes,GenesAffy]=ExtractSubjectsRNAProfile(sample,subjectIDs,Virus,Genes,FigureOnOff,Padding,SubsMeanOnOff,InterpOnOff)
% Same function as ExtractSubjectsRNAProfile but extracts
% a mask of non-missing vlaues to deal with incomplete data.
AllTimes=unique(sample.excelmat);
p=length(Genes);
if p<length(sample.blood.RNA.rownames)-62
    [RowIds,Genes]=FindGenesRows(Genes,sample.blood.RNA.rownames);
else
    Genes=Genes;
    RowIds=1:p;
end
GenesAffy=sample.rownames_affy(RowIds);
for i=1:length(subjectIDs)
    [SamplesAux]=ExtractSubjectRNASamples(sample,subjectIDs(i),Virus);
    AllTimes=union(AllTimes,SamplesAux.SampleTimes);
end
AllTimes=AllTimes;
Ntimes=length(AllTimes);
k=1;
for i=1:length(subjectIDs)
    [SamplesAux]=ExtractSubjectRNASamples(sample,subjectIDs(i),Virus);
    ColIds=SamplesAux.Samples;
    X{k}=sample.blood.RNA.data(RowIds,ColIds)';
    
    [CommonTimes,ObsIndx{k},IB]=intersect(AllTimes,SamplesAux.SampleTimes);
    X{k}=X{k}(IB,:);
    n=size(X{k},1);
    if n>0
        if SubsMeanOnOff
            X{k}=X{k}-repmat(mean(X{k},1),n,1);
        end
        %X{k}=1/max(max(abs(X{k})))*X{k};
        Samples{k}=SamplesAux;
        if FigureOnOff==1
            figure
            plot(X{k});
        end
        y(k)=SamplesAux.Label(length(SamplesAux.Label));
        k=k+1;
    end
end

for i=1:(k-1)
    Aux=zeros(Ntimes,p);
    Aux(ObsIndx{i},:)=X{i};
    X{i}=[repmat(Aux(1,:),Padding,1); Aux; repmat(Aux(end,:),Padding,1)];
    
    Mask{i}=zeros(Ntimes,p);
    Mask{i}(ObsIndx{i},:)=ones(length(ObsIndx{i}),p);
    Mask{i}=[ones(Padding,p);Mask{i};ones(Padding,p)];   
    
    n=length(Mask{i});
    if InterpOnOff
        x=find(Mask{i}(:,1)==1);
        if length(x)<n
            xi=find(Mask{i}(:,1)==0);
            Xinterp{i}=zeros(n,p);
             for g=1:p
                 Y=X{i}(:,g);
                 Xinterp{i}(x,g)=Y;
                 Xinterp{i}(xi,g)=interp1(x,Y,xi);
             end
            Mask{i}=ones(n,p);
            X{i}=Xinterp{i};
        end
    end
end