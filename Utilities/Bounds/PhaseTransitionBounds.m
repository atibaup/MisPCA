function [B,PT,PT2,PTUB,PTLB,lambda,lambda2,lambda3L,lambda3U]=PhaseTransitionBounds(H,d,s_d,sigma)
n=length(d);
F=size(H,2);
p=size(H,1);

dmax=max(d);
%dsorted=sort(unique([0;d]),'ascend');
dsorted=sort(unique(d),'ascend');
dPoints=length(dsorted);

%compute r_h
r_h=zeros(F,F*dPoints);
for i=1:dPoints
    r_h(:,(i-1)*F+1:i*F)=circshift(H,-dsorted(1))'*circshift(H,-dsorted(i));
end
% generate block_toeplitz matrix

R_H=sparse(F*dPoints,F*dPoints);
for i=1:dPoints
    for j=i:dPoints
        Indx=(abs(i-j))*F+1:((abs(i-j)+1)*F);
        Indxi=(j-1)*F+1:j*F;
        Indxj=(i-1)*F+1:i*F;
        R_H(Indxi,Indxj)=r_h(:,Indx);
        R_H(Indxj,Indxi)=r_h(:,Indx)';
    end
end

% compute: lambda1(diag(s_d \kron sigma)^(1/2)*R_H*diag(s_d \kron sigma)^(1/2))

D=diag(kron(s_d(:),sigma(:)).^(1/2));
B=D*R_H*D;
lambda=eigs(B,F);

% Bound using Toeplitz Bounds if F=1: lambda1(R_H)
if F==1
    [lambdaT,lambda1]=ToeplitzBounds(H,dsorted,F);
    lambdaLB=lambda1(:,1);
    lambdaUB=lambda1(:,2);
    lambdaT=lambdaT;
end
% Phase Transition Point, knowledge of h and d:

PT=sqrt(p/n)./lambda;

% Phase Transition Point, knowledge of h, uniform delays:

D2=diag(kron(ones(dPoints,1),sigma(:)).^(1/2));
B2=D2*R_H*D2;
lambda2=1/dPoints*eigs(B2,F);
PT2=sqrt(p/n)./lambda2;

% Phase Transition Point, knowledge of r_h, uniform delays:
if F==1
    lambda3L=1/dPoints*lambdaLB;
    lambda3U=1/dPoints*lambdaUB;
    PTUB=sqrt(p/n)./lambda3L;
    PTLB=sqrt(p/n)./lambda3U;
else
    PTUB=NaN;
    PTLB=NaN;
end

% figure;
% plot(PT,'-ko'); hold on;
% plot(PT2,'-rs'); hold on;
% plot(PTUB,'-m^'); hold on;
% plot(PTLB,'-mv'); hold on;
% legend('PT','PT uniform dets','PT Toeplitz bounds');

% figure; 
% subplot(3,1,1)
% imagesc(r_h)
% subplot(3,1,2)
% imagesc(R_H)
% subplot(3,1,3)
% imagesc(B)
