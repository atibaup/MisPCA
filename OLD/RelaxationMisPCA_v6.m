function [dhat,Fhat,sigma_h,sigma_n,lambdamax]=RelaxationMisPCA_v6(X,FigOnOff)
S=length(X);
n=size(X{1},1);
p=size(X{1},2);
Sigma=cell(S,1);
SigmaAv=zeros(n);

opts.disp=0;
for i=1:S
    if p==1
         Sigma{i}=1/p*(X{i})*(X{i})';
    else
        Sigma{i}=cov(X{i}');
    end
end
z=cell(S,n)
for i=1:S
  for j=1:n
      T=CreateTranslationMatrix(j-1,n);
      X{i,j}=T'*Sigma{i}*T;
  end
end
cvx_begin
    variable H(n,n) symmetric;
    expression c(n*S,1);

    maximize(max(c));
    subject to
        H== semidefinite(n);
        trace(H)==1;    
cvx_end

C=-cvx_optval;
[hest,d]=eigs(H,1,'lm',opts);

%% Obtain delays from initial estimate

SigmaAligned=zeros(n);
dhat=zeros(S,1);
for i=1:S
    d=zeros(n,1);
    for k=1:n
        T=CreateTranslationMatrix(k-1,n);
        d(k)=hest'*T'*Sigma{i}*T*hest;
    end
    [maxVal indx]=max(d);
    dhat(i)=indx-1;
    [T]=CreateTranslationMatrix(dhat(i),n);
    SigmaAligned=SigmaAligned+1/S*T'*Sigma{i}*T;
end
%
[Fhat,lambdamax]=eigs(SigmaAligned,1,'lm',opts);%
%Fhat=hest;

sigma_n=1/(n-1)*(trace(SigmaAv)-lambdamax);
sigma_h=lambdamax-sigma_n;
Fhat=Fhat;
lambdamax=lambdamax;