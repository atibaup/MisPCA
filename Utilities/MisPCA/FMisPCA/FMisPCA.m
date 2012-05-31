function [H,dhat,Lambdas,SigmaAv,sigma_h,sigma_n]=FMisPCA(Sigma,Masks,Options)
S=length(Sigma);
p=size(Sigma{1},1);
F=Options.k;
beta=Options.beta;
dmax=Options.dmax;

Niter=10;
tol=1e-3;
NrndmIni=Options.Nrndm;
Set=[-round(dmax/2)*ones(S,1),round(dmax/2)*ones(S,1)];
S=length(Sigma);

costVec=zeros(NrndmIni,1);
Hest=cell(NrndmIni,1);
aEst=cell(NrndmIni,1);
SigmaAvVec=cell(NrndmIni,1);
Lambdas=cell(NrndmIni,1);

e = ones(p,1);
W = spdiags([-e 2*e -e], -1:1, p, p);
        
%disp('RelaxMisPCA-start...')
% Parallelizable for (deactivated for OCTAVE compatibility) 
 for r=1:NrndmIni
    %fprintf('... %d/%d ... ',r,NrndmIni)
    if isempty(Options.Hini)
        H=randn(p,F);
        [H,D]=eigs(H*H',F);
    else
        H=Options.Hini;
        D=[];
    end

    HHt=H*H';    
    [dummy,a,SigmaAvTmp,Weights]=constructC_FMisPCA(Sigma,Masks,H,Set,Options.GridSize);
    B=diag(sqrt(Weights))*SigmaAvTmp*diag(sqrt(Weights))-beta*(W'*W);
    cost=zeros(Niter,1);
    cost(1)=trace(B*HHt);
    
    %% A-MisPCA
    convergence=0;
    k=2;  
    %tic
    while ~convergence && k<= Niter
        %disp('Iter')
        %tic
        [c,a,SigmaAvTmp,Weights]=constructC_FMisPCA(Sigma,Masks,H,Set,Options.GridSize);
        %Time1=toc
        
        %tic
        B=diag(sqrt(Weights))*SigmaAvTmp*diag(sqrt(Weights))-beta*(W'*W);
        [H,D]=eigs(B,F);
        HHt=H*H';
        cost(k)=sum(diag(D));
        %Time2=toc
        
        if Options.FiguresOnOff
            figure
            subplot(5,1,1)
            plot(1/max(c)*c); hold on;
            stem(1/max(a)*a,'r')
            subplot(5,1,2)
            imagesc(SigmaAvTmp);
            subplot(5,1,3)
            plot(H);
            subplot(5,1,4)
            imagesc(HHt)
            subplot(5,1,5)
            plot(diag(D));
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
    aEst{r}=a;
    Hest{r}=H;
    Lambdas{r}=diag(D);
    SigmaAvVec{r}=SigmaAvTmp;
    if  Options.FiguresOnOff
        figure
        plot(cost(1:k-1)); hold on;
        legend('Objective');%,'Eigmax')
    end
end
%fprintf('\n ');
%disp('...RelaxMisPCA-END')
if  Options.FiguresOnOff
    figure
    plot(costVec); hold on;
    title('Costs different initialization')
    legend('Cost')
end

[dummy,maxIndx]=max(costVec);
maxIndx=maxIndx(1);
Lambdas=Lambdas{maxIndx};
SigmaAv=SigmaAvVec{maxIndx};
a=aEst{maxIndx};
H=Hest{maxIndx};

%% Obtain delays from initial estimate

d=zeros(S,1);
dhat=cell(S,1); % cell format, compatible with OPFA
for i=1:S
        for j=1:Options.GridSize:(Set(i,2)-Set(i,1))
            if  a((i-1)*p+j)>0
                d(i)=Set(i,1)+j-1;
                dhat{i}=d(i)*ones(F,1);
            end
        end
end

sigma_n=1/(p-1)*(trace(SigmaAv)-sum(Lambdas));
sigma_h=sum(Lambdas)-sigma_n;
