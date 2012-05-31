function [C]=TraceObjective_v2(H,Sigma,Masks,a)
warning off;
n=size(Masks{1},1);
S=length(Sigma);
C=0;
cvx_begin
expression C;
for i=1:S
    for j=1:n
        if a((i-1)*n+j)>0
            T=CreateTranslationMatrix(j-1,n);
%             Aux=zeros(n);
%             Aux(Masks{i}==1,Masks{i}==1)=Sigma{i};
%             Aux=T'*Aux*T;
%             TrsMask=T*Masks{i};
%             TrsMask=(TrsMask==1);
%             Aux=Aux(TrsMask,TrsMask);
%           C=C-1/S*a((i-1)*n+j)*trace(H(TrsMask,TrsMask)*Aux);
            TransH=T*H*T';
            Aux=Sigma{i};
            C=C-1/S*a((i-1)*n+j)*trace(TransH(Masks{i}==1,Masks{i}==1)*Aux);
        end
    end
end
cvx_end