function [dhat,Fhat,sigma_h,sigma_n,lambdamax]=RelaxationMisPCA_v4(X,FigOnOff)
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

cvx_begin
    variable Z(n,n) symmetric;
    variable a(n*S,1);
    expression c(n*S,1);
        for i=1:S
            for j=1:n
                T=CreateTranslationMatrix(j-1,n);
                c((i-1)*n+j)=trace(Z*T'*Sigma{i}*T);
            end
        end
    minimize(a'*c);
    subject to
        for i=1:S
            norm(a((i-1)*n+1:i*n),1) <= 1;
            a((i-1)*n+1:i*n)>=zeros(n,1);
        end
        Z== semidefinite(n);
        trace(Z)==1;    

cvx_end

for i=1:S
    Shownorma=  norm(a((i-1)*n+1:i*n),1)
end

C=-cvx_optval;

for i=1:S
    for j=1:n
        [T]=CreateTranslationMatrix(j,n);
        SigmaAv=SigmaAv+1/S*a((i-1)*n+j)* T*Sigma{i}*T;
    end
end
[hest,lambdamax]=eigs(SigmaAv,1,'lm',opts);%

if FigOnOff
    figure
    subplot(3,1,1)
    plot(a); hold on;
    %    stem(a,'r')
    subplot(3,1,2)
    imagesc(SigmaAv);
    subplot(3,1,3)
    %    [v,d]=eigs(H,1,'lm');
    plot(hest)
end


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