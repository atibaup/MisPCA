function [Load,IntGenes,AverageLoad,A_o]=findLoadingsOnAlignedSubspace(X,Masks,F_MisPCA,Lambdas,dopt_MisPCA,Options,ChosenFactors,Genes,Threshold,NormalizeOnOff,FigOnOff)
f=length(ChosenFactors);
p=length(Genes);
S=length(X);
Load=0;
Load=[];
AverageLoad=zeros(f,p);
A_o=cell(1,Options.S);
for i=1:Options.S
    Fadj=zeros(Options.n,f);
    for j=1:f
         [T]=CreateTranslationMatrix(dopt_MisPCA{j}(i),Options.n);
         Fadj(:,j)= T*F_MisPCA(:,ChosenFactors(j));
    end
    Fadj=Fadj(Masks{i}==1,:)*diag(sqrt(Lambdas(ChosenFactors)));
    ProjOnF=(Fadj'*Fadj)\Fadj';
    %Load=Load+ 1/Options.S*diag(Lambdas(ChosenFactors).^(-1))*ProjOnF*X{i};
    if ~NormalizeOnOff
        Aux=ProjOnF*X{i};
        Load=[Load; Aux];
        AverageLoad=AverageLoad+1/S*Aux;
    else
        Aux=ProjOnF*X{i};
        A_o{i}=Aux;
        %Aux=1/sum(Aux)*Aux;
        Load=[Load; Aux];
        AverageLoad=AverageLoad+1/S*Aux;
    end
    AverageLoad=AverageLoad+1/S*Aux;
end
if ~NormalizeOnOff
    for j=1:p
        AverageLoad(:,j)=1/max(max(abs(AverageLoad(:,j))))*AverageLoad(:,j);
        Load(:,j)=1/max(max(abs(Load(:,j))))*Load(:,j);
    end
end

Norms=diag(Load'*Load);
IntGenes.Ids=find(Norms>=Threshold*max(Norms));
IntGenes.Names=Genes(IntGenes.Ids);
if FigOnOff
    figure;subplot(2,1,1);
    imagesc(Load);
    subplot(2,1,2);
    plot(Norms); hold on;
    plot(1:length(Genes),Threshold*max(Norms),'-r');
    title(strcat(num2str(length(IntGenes.Ids)),' genes selected.'))
end