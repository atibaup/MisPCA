function [H,eta,SigmaAv]=EstimateHnucNorm_v4(a,Sigma,Masks,beta)
n=size(Masks{1},1);
S=length(Sigma);
SigmaAv=zeros(n);
for i=1:S
    for j=1:n
        if a((i-1)*n+j)>0
            T=CreateTranslationMatrix(j-1,n);
            Aux=zeros(n);
            Aux(Masks{i}==1,Masks{i}==1)=Sigma{i};
            Aux=T'*Aux*T;
            %Aux=circshift((circshift(Aux',-(j-1)))',-(j-1));
            SigmaAv=SigmaAv+1/S*a((i-1)*n+j)*Aux;
        end
    end
end
delta=1/2;
W=eye(n)+delta*(diag(-ones(n-1,1),+1)+diag(-ones(n-1,1),-1));
    
B=SigmaAv+beta*W;

eta=eigs(B,1,'lm');

W=eta*speye(n)-B;
%tic
%V = null(W);
[V]=ComputeNullSpaceRankr(W,6);
%TimeNullspace=toc
H=V*diag(1/size(V,2)*ones(size(V,2),1))*V';
