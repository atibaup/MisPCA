function [X,SigmaAv,H,dini,Scell,D]=GenerateFMisPCAData(p,n,S,F,SNR,r,dmax,Spacing,sigma,H)

if isempty(sigma)
   sigma=ones(F,1); 
end

if isempty(H)
    H=zeros(p,F);
    H(:,1)=[ones(r,1);zeros(p-r,1);]; % upregulation
    for i=2:F
        H(:,i)=circshift(smoothSig(H(:,1)',6,'triang'),round(rand*dmax));
    end
    [H,dummy]=eigs(H*H',F);
end

d=Spacing*(unidrnd(round(dmax/Spacing)+1,S,1)-1);

X=cell(S,1);
Scell=cell(S,1);
D=cell(S,1);
dini=cell(S,1);
SigmaAv=zeros(p);
for i=1:S
    D{i}=sqrt(10^(SNR/10))*circshift(H*diag(sigma),d(i))*randn(F,n);
    X{i}=D{i} + randn(p,n);
    Scell{i}=X{i}*X{i}';
    SigmaAv=SigmaAv+1/S*Scell{i};
    dini{i}=d(i)*ones(F,1);
end
