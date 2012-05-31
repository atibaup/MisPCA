function [Xhat,A,Err,MError]=RecoveredX(X,Masks,F,dhat,Lambdas)

S=length(X);
n=size(Masks{1},1);
p=size(Masks{1},2);

f=size(F,2);
A=cell(S,1);
Err=zeros(S,1);

Xhat=cell(S,1);
for i=1:S
    Fadj=zeros(n,f);

    for j=1:f
        Fadj(:,j)=circshift(F(:,j),-dhat{i}(j));
    end
    Fadj=Fadj*diag(sqrt(Lambdas));

    Indx=(Masks{i}(:,1)==1);

    Fp=Fadj(Indx,:);
    A{i}=(Fp'*Fp)\(Fp'*X{i}(Indx,:));
    
    Xhat{i}=Fadj*diag(sqrt(Lambdas))*A{i};   
    Err(i)=norm(Xhat{i}(Indx,:)-X{i}(Indx,:),'fro')^2;
end
MError=mean(Err);