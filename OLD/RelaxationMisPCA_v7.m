function [dhat,Fhat,sigma_h,sigma_n,lambdamax]=RelaxationMisPCA_v7(X,Masks,hini,Set,Niter,tol,NrndmIni,FigOnOff,beta,GridSize)
S=length(X);
n=size(X{1},1);
p=size(X{1},2);
Sigma=cell(S,1);

for i=1:S
    if p==1
         Sigma{i}=1/p*(X{i})*(X{i})';
    else
        Sigma{i}=cov(X{i}');
    end
end
[dhat,Fhat,sigma_h,sigma_n,lambdamax]=RelaxationMisPCA_v7_fromSigma(Sigma,Masks,hini,Set,Niter,tol,NrndmIni,FigOnOff,beta,GridSize);