function [c,a,SigmaAv,Weights]=constructC_FMisPCA(Sigma,Masks,H,Set,GridSize)
n=size(Sigma{1},1);
S=length(Sigma);
%tic
Niter=0;
a=zeros(n*S,1);
c=zeros(n*S,1);

SigmaAv=zeros(n);
Weights=zeros(n,1);
HHt=H*H';
for i=1:S
    %ShowGrid= Set(i,1)+1:GridSize:Set(i,2)+1
    maxSoFar=0;  
    Aux=circshift(Sigma{i}',Set(i,1))';
    Aux=circshift(Aux,Set(i,1));
    for j= 1:GridSize:(Set(i,2)-Set(i,1))
        
        d=-(Set(i,1)+j-1);
        SigmaA=circshift(Sigma{i}',-d)';
        SigmaA=circshift(SigmaA,-d);
        Indx=(Masks{i}(:,1)==1);
        c((i-1)*n+j)=trace(SigmaA(Indx,Indx)*HHt(Indx,Indx));
        
        if c((i-1)*n+j)>=maxSoFar
           a((i-1)*n+1:i*n)=zeros(n,1);
           a((i-1)*n+j)=1;           
           Aux=SigmaA;
           maxSoFar=c((i-1)*n+j);
        end
        Niter=Niter+1;
    end
    SigmaAv=SigmaAv+1/S*Aux;
    Indx=(Masks{i}(:,1)==1);
    Weights=Weights+double(Indx);
end
Weights=1/S*Weights;
%Niter=Niter
%Timeloop=toc