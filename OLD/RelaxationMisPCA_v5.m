function [dhat,Fhat,sigma_h,sigma_n,lambdamax]=RelaxationMisPCA_v5(X,hini,Niter,tol,NrndmIni,GridSize,dmax,FigOnOff,beta)
S=length(X);
n=size(X{1},1);
p=size(X{1},2);

Grid=0:GridSize:dmax;
nGrid=length(Grid);

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
Set=[zeros(S,1),(n-1)*ones(S,1)];
%% Obtain PCA estimate to initialize
if length(hini)>1
    h=hini;
else
    [h, eigmax(1)]=eigs(SigmaAv,1,'lm',opts);
end
if nargin < 8
    beta=0;
end
costVec=zeros(NrndmIni,1);
eigmaxVec=zeros(NrndmIni,1);
lambdamaxVec=zeros(NrndmIni,1);
hVec=cell(NrndmIni,1);
hVec2=cell(NrndmIni,1);
aEst=cell(NrndmIni,1);
delta=0e-1;
disp('RelaxMisPCA-start...')
for r=1:NrndmIni
    if FigOnOff; fprintf('... %d/%d ... ',r,NrndmIni); end;
    %h=h+1/n*randn(n,1);
    h=randn(n,1);
    h=h/norm(h);
    eigmax(1)=h'*SigmaAv*h;
    cost(1)=0;
    %% Solve relaxed problem
    convergence=0;
    H=h*h';
    k=2;  
    while ~convergence && k<= Niter

        [c]=constructC(Sigma,H,Set,GridSize,dmax);
        [a]=EstimateaWuncertainty(c,nGrid,S,delta);
        [H,cost(k),SigmaAv]=EstimateHnucNorm_v2(a,Sigma,beta,GridSize,dmax);
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
    [hVec2{r},d]=eigs(H,1,'lm');
    hVec{r}=hest;
    aEst{r}=a;
    if FigOnOff
        figure
        plot(cost(1:k-1)); hold on;
        plot(eigmax(1:k-1),'-or');
        legend('Cost','Eigmax')
    end
end

disp('...RelaxMisPCA-END')
if FigOnOff
    figure
    plot(costVec); hold on;
    plot(eigmaxVec,'-or');
    title('Costs different initialization')
    legend('Cost','Eigmax')
end

[lambdamax,maxIndx]=max(eigmaxVec);
maxIndx=maxIndx(1);
lambdamax=lambdamax(1);
Fhat=hVec{maxIndx};
a=aEst{maxIndx};

%% Obtain delays from initial estimate
% 
% SigmaAligned=zeros(n);
dhat=zeros(S,1);
for i=1:S
    d=zeros(n,1);
    for k=1:nGrid
        T=CreateTranslationMatrix(Grid(k),n);
        d(k)=Fhat'*T'*Sigma{i}*T*Fhat;
    end
    [~, indx]=max(d);
    dhat(i)=Grid(indx);
% %     [T]=CreateTranslationMatrix(dhat(i),n);
% %     SigmaAligned=SigmaAligned+1/S*T'*Sigma{i}*T;
%         for j=1:nGrid
%             if  a((i-1)*nGrid+j)>0
%                 dhat(i)=Grid(j);
%                 %[T]=CreateTranslationMatrix(dhat(i),n);
%                 %SigmaAligned=SigmaAligned+1/S*a((i-1)*n+j)*T'*Sigma{i}*T;
%             end
%         end
end
%
%[Fhat,lambdamax]=eigs(SigmaAligned,1,'lm',opts);
%Fhat=hest;

sigma_n=1/(n-1)*(trace(SigmaAv)-lambdamax);
sigma_h=lambdamax-sigma_n;
Fhat=Fhat;
lambdamax=lambdamax;