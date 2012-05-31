function [Xhat,A,FactorsDistance,MError]=FMisPCA_Performance(X,D,Masks,F,FMisPC,d,dhat,Lambdas,FigOnOff)

S=length(X);
n=size(Masks{1},1);
p=size(Masks{1},2);

fo=size(F,2);
f=size(FMisPC,2);
A=cell(S,1);
Err=zeros(S,1);
FactorsDistance=zeros(S,1);
Xhat=cell(S,1);
FMisPCadj=cell(S,1);
Fadj=cell(S,1);
for i=1:S
    FMisPCadj{i}=zeros(n,f);
    for j=1:f
        FMisPCadj{i}(:,j)=circshift(FMisPC(:,j),-dhat{i}(j));
    end

    FMisPCadj{i}=FMisPCadj{i}*diag(sqrt(Lambdas));
    
    Fadj{i}=zeros(n,fo);
    for j=1:fo
        Fadj{i}(:,j)=circshift(F(:,j),d{i}(j));
    end

    Indx=(Masks{i}(:,1)==1);

    Fp=FMisPCadj{i}(Indx,:);
    A{i}=(Fp'*Fp)\(Fp'*X{i}(Indx,:));
    
    Xhat{i}=FMisPCadj{i}*A{i};   
    Err(i)=norm(Xhat{i}-D{i},'fro')/norm(D{i},'fro');
    FactorsDistance(i)=FactorsSubspaceDistance(FMisPCadj{i},Fadj{i});

end
MError=mean(Err);
FactorsDistance=mean(FactorsDistance);


if FigOnOff
    h=figure;
    h2=figure;
    kPlot=1;
    for i=1:S
        figure(h)
        subplot(S,3,kPlot)
        plot(D{i});
        title('Original')
        
        subplot(S,3,kPlot+1)
        plot(Xhat{i});
        title('Reconstructed')
        
        subplot(S,3,kPlot+2)
        plot(Xhat{i}-D{i});
        title('Residual')
        
        kPlot=kPlot+3;
        figure(h2);
        subplot(2,1,1);
        plot(FMisPCadj{i}); hold on;
        title('Recovered')
        subplot(2,1,2);
        plot(Fadj{i}); hold on;
        title('Original')
    end
    
end