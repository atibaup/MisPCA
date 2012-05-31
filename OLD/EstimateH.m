function [H,C,SigmaAv]=EstimateH(a,Sigma,beta)
n=size(Sigma{1},1);
S=length(Sigma);
C=0;
H=sdpvar(n,n);
SigmaAv=zeros(n);
for i=1:S
    for j=1:n
        if a((i-1)*n+j)>0
            T=CreateTranslationMatrix(j-1,n);
            C=C-1/S*trace(H*T'*Sigma{i}*T);
            SigmaAv=SigmaAv+1/S*T'*Sigma{i}*T;
        end
    end
end
opts = sdpsettings('verbose',0);
F=set(H>0)+set(trace(H'*H) <=1);
sol=solvesdp(F,C,opts); %+beta*trace(H'*H)
%sol=sol

C=-double(C);
H=double(H);
