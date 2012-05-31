function [Res,Sigma,SigmaAv,Fadj]=ComputeResidual(dopt,F,X)
S=size(X,2);
n=size(X{1},1);
p=size(X{1},2);
k=size(F,2);
%figure
Res=cell(S,1);
for i=1:S
    Fadj=zeros(n,k);
    for j=1:k
         [T]=CreateTranslationMatrix(dopt{j}(i),n);
         Fadj(:,j)= T*F(:,j);
    end
    warning off;
    P=eye(n)-Fadj*((Fadj'*Fadj)\Fadj');
    Res{i}=P*X{i};
end
Sigma=cell(S,1);
SigmaAv=zeros(n);
for i=1:S
    Sigma{i}=cov(Res{i}');%(Res{i}-repmat(mean(Res{i},2),1,p))*(Res{i}-repmat(mean(Res{i},2),1,p))';
    Sigma{i}=1/2*(Sigma{i}+Sigma{i}');
    SigmaAv=SigmaAv+1/S*Sigma{i};
end