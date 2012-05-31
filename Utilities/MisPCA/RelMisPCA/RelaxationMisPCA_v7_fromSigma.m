function [dhat,Fhat,sigma_h,sigma_n,lambdamax]=RelaxationMisPCA_v7_fromSigma(Sigma,Masks,hini,Set,Niter,tol,NrndmIni,FigOnOff,beta,GridSize)
S=length(Sigma);
n=size(Masks{1},1);
SigmaAvIni=sparse(n,n);

for i=1:S
    Indx=(Masks{i}(:,1)==1);
    SigmaAvIni(Indx,Indx)=SigmaAvIni(Indx,Indx)+1/S* Sigma{i}(Indx,Indx);
end

%% Obtain PCA estimate to initialize

if nargin < 8
    beta=0;
end
costVec=zeros(NrndmIni,1);

hVec2=cell(NrndmIni,1);
aEst=cell(NrndmIni,1);
SigmaAvVec=cell(NrndmIni,1);

%disp('RelaxMisPCA-start...')
% Parallelizable for (deactivated for OCTAVE compatibility) 
 for r=1:NrndmIni
    %fprintf('... %d/%d ... ',r,NrndmIni)

    if length(hini)<=1
        h=randn(n,1);
    end
    h=h/norm(h);
    eigmax=zeros(Niter,1);
    cost=zeros(Niter,1);
    cost(1)=h'*SigmaAvIni*h;
    %% Solve relaxed problem
    convergence=0;
    H=h*h';
    k=2;  
    %tic
    while ~convergence && k<= Niter
        %disp('Iter')
        %tic
        [c,a,SigmaAvTmp,Weights]=constructC_v3(Sigma,Masks,h,Set,GridSize);
        %[a]=EstimateaWuncertainty(c,n,S,delta);
        %Time1=toc
        %tic
        [H,cost(k),SigmaAvTmp,h]=EstimateHnucNorm_v5(SigmaAvTmp,Weights,beta);
%        [hest, eigmax(k)]=eigs(SigmaAv,1,'lm',opts);       
        %Time2=toc
        if FigOnOff
            figure
            subplot(5,1,1)
            plot(1/max(c)*c); hold on;
            stem(1/max(a)*a,'r')
            subplot(5,1,2)
            imagesc(SigmaAvTmp);
            subplot(5,1,3)
            hest=eigs(H,1,'lm');
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
    %TimeLoop=toc
    costVec(r)=cost(k-1);

    hVec2{r}=h;
    eigmaxVec(r)=hVec2{r}'*SigmaAvTmp*hVec2{r};
    %hVec{r}=hest;
    aEst{r}=a;
    SigmaAvVec{r}=SigmaAvTmp;
    if FigOnOff
        figure
        plot(cost(1:k-1)); hold on;
        %plot(eigmax(1:k-1),'-or');
        legend('Objective');%,'Eigmax')
    end
end

%disp('...RelaxMisPCA-END')
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
Fhat=hVec2{maxIndx};
a=aEst{maxIndx};
SigmaAv=SigmaAvVec{maxIndx};
%% Obtain delays from initial estimate
% 
% SigmaAligned=zeros(n);
dhat=zeros(S,1);
for i=1:S
%     d=zeros(n,1);
%     for k=1:n
%         T=CreateTranslationMatrix(k-1,n);
%         d(k)=hest'*T'*Sigma{i}*T*hest;
%     end
%     [maxVal indx]=max(d);
%     dhat(i)=indx-1;
%     [T]=CreateTranslationMatrix(dhat(i),n);
%     SigmaAligned=SigmaAligned+1/S*T'*Sigma{i}*T;
        for j=1:GridSize:(Set(i,2)-Set(i,1))
            if  a((i-1)*n+j)>0
                dhat(i)=Set(i,1)+j-1; %Set(i,1)+
                %[T]=CreateTranslationMatrix(dhat(i),n);
                %SigmaAligned=SigmaAligned+1/S*a((i-1)*n+j)*T'*Sigma{i}*T;
            end
        end
end
%
%[Fhat,lambdamax]=eigs(SigmaAligned,1,'lm',opts);
%Fhat=hest;
sigma_n=1/(n-1)*(trace(SigmaAv)-lambdamax);
sigma_h=lambdamax-sigma_n;
