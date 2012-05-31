function [dhat,Fhat,sigma_h,sigma_n,lambdamax]=RelaxationMisPCA_v8_fromSigma(Sigma,Masks,hini,Set,Niter,tol,NrndmIni,FigOnOff,beta,GridSize)
S=length(Sigma);
n=size(Masks{1},1);


SigmaAvIni=zeros(n);

for i=1:S
    SigmaAvIni(Masks{i}==1,Masks{i}==1)=SigmaAvIni(Masks{i}==1,Masks{i}==1)+1/S* Sigma{i};
end

cost=zeros(Niter,1);
eigmax=zeros(Niter,1);

%% Obtain PCA estimate to initialize



sigma_n=1/(n-1)*(trace(SigmaAv)-lambdamax);
sigma_h=lambdamax-sigma_n;
Fhat=Fhat;
lambdamax=lambdamax;