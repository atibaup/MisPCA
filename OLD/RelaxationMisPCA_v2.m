function [dhat,Fhat,sigma_h,sigma_n,lambdamax]=RelaxationMisPCA_v2(X,hini,Niter,tol,NrndmIni,FigOnOff)
S=length(X);
n=size(X{1},1);
p=size(X{1},2);

Sigma=cell(S,1);
SigmaAv=zeros(n);
MaxEig=zeros(S,1);
opts.disp=0;
for i=1:S
    if p==1
         Sigma{i}=1/p*(X{i})*(X{i})';
    else
        Sigma{i}=cov(X{i}');
    end
    SigmaAv=SigmaAv+1/S* Sigma{i};
    MaxEig(i)=eigs(Sigma{i},1,'lm',opts);
end

cost=zeros(Niter,1);
eigmax=zeros(Niter,1);

%% Obtain PCA estimate to initialize
if length(hini)>1
    h=hini;
else
    [h, eigmax(1)]=eigs(SigmaAv,1,'lm',opts);
end

costVec=zeros(NrndmIni,1);
eigmaxVec=zeros(NrndmIni,1);
lambdamaxVec=zeros(NrndmIni,1);
hVec=cell(NrndmIni,1);
for r=1:NrndmIni
    h=h+1/n*randn(n,1);
    %h=randn(n,1);
    h=h/norm(h);
    eigmax(1)=h'*SigmaAv*h;
    cost(1)=0;
    %% Solve relaxed problem
    convergence=0;
    H=h*h';
    k=2;   
    while ~convergence && k<= Niter
        [c]=constructC(Sigma,H);
        [a]=Estimatea(c,n,S);
        [H,cost(k),SigmaAv]=EstimateHnucNorm_v2(a,Sigma);
        [hest, eigmax(k)]=eigs(SigmaAv,1,'lm',opts);       
        if FigOnOff
            figure
            subplot(5,1,1)
            plot(c); hold on;
            stem(max(c)/max(a)*a,'r')
            subplot(5,1,2)
            imagesc(SigmaAv);
            subplot(5,1,3)
            %[h,d]=eigs(H,1,'lm');
            %H=h*h';
            plot(hest);
            %eigmax(k)=h'*SigmaAv*h;
            subplot(5,1,4)
            imagesc(H)
            subplot(5,1,5)
            [D]=eig(H);
            plot(D);
        end        
        if k>2
            if abs(cost(k)-cost(k-1))<=tol
                convergence=1;
            end
        end
        k=k+1;
    end
    
    costVec(r)=cost(k-1);
    eigmaxVec(r)=eigmax(k-1);
    lambdamaxVec(r)=eigmax(k-1);
    hVec{r}=hest;
    
    if FigOnOff
        figure
        plot(cost(1:k-1)); hold on;
        plot(eigmax(1:k-1),'-or');
        legend('Cost','Eigmax')
    end
end

if FigOnOff
    figure
    plot(costVec); hold on;
    plot(eigmaxVec,'-or');
    title('Costs different initialization')
    legend('Cost','Eigmax')
end

[maxVal,maxIndx]=max(eigmaxVec);
maxIndx=maxIndx(1);
hest=hVec{maxIndx};

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