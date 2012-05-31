function [d,Fhat,Lambdas,sigma_h,sigma_n]=OraclePCA(Sigma,F,dtrue)
S=length(Sigma);
n=size(Sigma{1},1);
SigmaAv=zeros(n);
d=cell(S,1);
for i=1:S
    Aux=circshift(Sigma{i}',-dtrue{i}(1))';
    Aux=circshift(Aux,-dtrue{i}(1));
    SigmaAv=SigmaAv+1/S*Aux;
    d{i}=-dtrue{i};
end


[Fhat,Lambdas]=eigs(SigmaAv,F);
Lambdas=diag(Lambdas);
sigma_n=1/(n-1)*(trace(SigmaAv)-sum(Lambdas));
sigma_h=sum(Lambdas)-sigma_n;

