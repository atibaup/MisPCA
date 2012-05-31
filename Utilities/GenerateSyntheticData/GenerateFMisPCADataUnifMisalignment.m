function [X,S,H,d,Scell]=GenerateFMisPCADataUnifMisalignment(p,n,F,SNR,r,dmax,Spacing,sigma,H)

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


NperDelay=round(n/(round(dmax/Spacing)+1))

n=NperDelay*(round(dmax/Spacing)+1);
d=[];
for i=1:round(dmax/Spacing)+1
  d=[d;(i-1)*Spacing*ones(NperDelay,1)];
end


A=randn(F,n);
X=zeros(p,n);
Scell=cell(n,1);
for i=1:n
    X(:,i)=sqrt(10^(SNR/10))*circshift(H*diag(sqrt(sigma)),d(i))*A(:,i) + randn(p,1);
    Scell{i}=X(:,i)*X(:,i)';
end
S=cov(X');
