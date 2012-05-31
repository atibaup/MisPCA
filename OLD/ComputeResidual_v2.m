function [Sigma,SigmaAv,FadjAv]=ComputeResidual_v2(dopt,F,SigmaOr,lambdamax)
S=length(SigmaOr);
n=size(SigmaOr{1},1);
k=size(F,2);
Fadj=zeros(n,k);
SigmaAv=zeros(n);
Sigma=cell(S,1);
FadjAv=zeros(n,k);
for i=1:S
    for j=1:k
         [T]=CreateTranslationMatrix(dopt{j}(i),n);
         Fadj(:,j)= T*F(:,j);
    end
    Sigma{i}= SigmaOr{i}-Fadj*diag(lambdamax)*Fadj';
    Sigma{i}=1/2*(Sigma{i}+Sigma{i}');
    SigmaAv=SigmaAv+1/S*Sigma{i};
    FadjAv=FadjAv+1/S*Fadj;
end