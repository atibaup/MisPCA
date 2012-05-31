function [H,dhat,Lambdas]=BF_FMisPCA(Sigma,Options)

S=length(Sigma);
p=size(Sigma{1},1);
F=Options.k;
beta=Options.beta;
dmax=Options.dmax;
GridSize=Options.GridSize;

Set=[(-round(dmax/2)+1)*ones(S,1),round(dmax/2)*ones(S,1)];

e = ones(p,1);
W = spdiags([-e 2*e -e], -1:1, p, p);
WtW=-beta*(W'*W);

[Grid]=GridSize*CreateNdimGrid(round(Set/GridSize));
Ngrid=length(Grid);

dhat=cell(S,1);
opts.disp=0;
opts.issym=1;
opts.tol=1e-2;

if S>1
    C=zeros(Ngrid,1);
    for i=1:Ngrid
        if mod(i,round(Ngrid/3))==0;    disp(['Iter ',num2str(i),'/',num2str(Ngrid)]); end;
        SigmaAv=zeros(p);
        for k=1:S
            Aux=circshift(Sigma{k},[-Grid(i,k),-Grid(i,k)]);
            SigmaAv=SigmaAv+1/S*Aux;
        end
        Lambdas=eigs(SigmaAv-WtW,F,'lm',opts);
        C(i)=sum(Lambdas);
    end
    [dummy, maxIndx]=max(C);
    maxIndx=maxIndx(1);
    d=Grid(maxIndx,:);
    SigmaAv=zeros(p);
    for k=1:S
        Aux=circshift(Sigma{k}',-d(k))';
        Aux=circshift(Aux,-d(k));
        SigmaAv=SigmaAv+1/S*Aux;
        dhat{k}=-d(k)*ones(F,1);
    end
    [H Lambdas]=eigs(SigmaAv-WtW,F,'lm',opts);
else
    [H Lambdas]=eigs(Sigma{1}-WtW,F,'lm',opts);
end

Lambdas=diag(Lambdas);
