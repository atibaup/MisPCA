function [H,dhat,Lambdas,SigmaAv,sigma_h,sigma_n]=SeqFMisPCA(Sigma,Options)

S=length(Sigma);
p=size(Sigma{1},1);
F=Options.k;
beta=Options.beta;
dmax=Options.dmax;
Set=[-round(dmax/2)*ones(S,1),round(dmax/2)*ones(S,1)];

e = ones(p,1);
W = spdiags([-e 2*e -e], -1:1, p, p);
WtW=-beta*(W'*W);

dhat=cell(S,1);
%% 1st iteration

SigmaAv=Sigma{1};   
B=SigmaAv-WtW;
[H,Lambdas]=eigs(B,F);
dhat{1}=zeros(F,1); 

%% 2nd and more iterations
if S>1
    Left=2:S;
    while(~isempty(Left))
        Next=Left(1);
        CostNext=zeros(round((Set(Next,2)-Set(Next,1))/Options.GridSize),1);
        for j=1:Options.GridSize:(Set(Next,2)-Set(Next,1))
            d=-(Set(Next,1)+j-1);
            Hc=circshift(H,d);
            CostNext(j)=trace(Hc'*Sigma{Next}*Hc);
        end
        [dummy,IndxMax]=max(CostNext);
        IndxMax=IndxMax(1);
        d=-(Set(Next,1)+IndxMax-1);
        dhat{Next}=-d*ones(F,1);
        Sigmac=circshift(Sigma{Next}',-d)';
        Sigmac=circshift(Sigmac,-d);
        SigmaAv=Next/S*SigmaAv+ 1/S*Sigmac;
        B=SigmaAv-WtW;
        [H,Lambdas]=eigs(B,F);
        Left=setdiff(Left,Next);
    end
end
        
Lambdas=diag(Lambdas);
sigma_n=1/(p-1)*(trace(SigmaAv)-sum(Lambdas));
sigma_h=sum(Lambdas)-sigma_n;

