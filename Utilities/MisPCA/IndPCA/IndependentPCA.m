function [dhat,Fhat,sigma_h,sigma_n,lambdamax]=IndependentPCA(X,GridSize,dmax)
S=length(X);
n=size(X{1},1);
p=size(X{1},2);

Grid=round(-round(dmax/2):GridSize:round(dmax/2));
nGrid=length(Grid);

Sigma=cell(S,1);
SigmaAv=zeros(n);
Xs=[];
for i=1:S
    if p==1
        Xs=[Xs X{i}];
        Sigma{i}=X{i}*X{i}';
    else
        Sigma{i}=cov(X{i}');
        SigmaAv=SigmaAv+1/S*Sigma{i};
    end
end
if p==1
    SigmaAv=cov(Xs');
end

opts.disp=0;
V=cell(S,1);
eigInd=zeros(S,1);
if p>1
    for i=1:S
        [V{i},eigInd(i)]=eigs(Sigma{i},1,'lm',opts);
    end
else
    for i=1:S
        V{i}=X{i}/norm(X{i});
        eigInd(i)=norm(X{i});
    end
end

dhat=zeros(S,1);
F=cell(1,n);
d=zeros(1,n);
Rest=1/S*V{1};
if S>1
    Left=2:S;
    while(~isempty(Left))
        Next=Left(1);
        dplus=zeros(nGrid,1);
        dminus=zeros(nGrid,1);
        for k=1:nGrid
            Vshift=circshift(V{Next},-Grid(k));
            dplus(k)=norm(Rest-Vshift,2);
            dminus(k)=norm(Rest+Vshift,2);
        end
        [minvalplus, indxpl]=min(dplus);
        [minvalminus, indxminus]=min(dminus);
        if minvalplus <= minvalminus
            sign=1;
            indx=indxpl;
        else
             sign=-1;
            indx=indxminus;
        end
        dhat(Next)=-Grid(indx);
        Vshift=circshift(V{Next},-Grid(indx));
        Rest=Rest+ sign*1/S*Vshift;
        Left=setdiff(Left,Next);
    end
end
Fhat=1/norm(Rest)*Rest;
maxVal=mean(eigInd);
sigma_n=1/(n-1)*(trace(SigmaAv)-maxVal);
sigma_h=maxVal-sigma_n;
Fhat=Fhat;
lambdamax=maxVal;
