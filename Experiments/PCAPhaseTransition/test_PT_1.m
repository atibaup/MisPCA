
F=3;
p=100;
n=p;

%d=Spacing*round(rand(n,1)*dmax/Spacing); % uniform
d=Spacing*round((1/2+1/10*randn(n,1))*dmax/Spacing); % gaussian
dsorted=sort(unique(d),'ascend');

SpacingSig=round(p/3);

H1=zeros(p,F);
H1(:,1)=[ones(r,1);zeros(p-r,1)];
H1(:,2)=[zeros(SpacingSig,1);triang(r);zeros(p-SpacingSig-r,1)];
H1(:,3)=[zeros(2*SpacingSig,1);ones(2*r,1);zeros(p-2*SpacingSig-2*r,1)];
H1=H1*diag(1./norm(H1,2));


if dmax>0
    s_d=1/n*hist(d,dsorted);
else
    s_d=1;
end


H_=zeros(p,p);
for i=1:length(s_d)
    H_=H_+s_d(i)*circshift(H1,dsorted(i))*diag(sigmabar)*circshift(H1,dsorted(i))';
end
Sigma=H_+eye(p);
[V,D]=eigs(Sigma,F);

FigOnOff=1;
SNR=-10:2:30;
[PredictedLambda1,PredictedCorr1,B]=AsymptoticPrediction(H1,d,s_d,sigmabar,SNR,0);


Corr1=zeros(length(SNR),F,Nrndm);
Lambda1=zeros(length(SNR),F,Nrndm);
for k=1:length(SNR)
    k=k
    for r=1:Nrndm     
        A=randn(F,n);
        X=zeros(p,n);
        for i=1:n
            X(:,i)=sqrt(10^(SNR(k)/10))*circshift(H1*diag(sigmabar.^(1/2)),d(i))*A(:,i) + randn(p,1);
        end
        S=cov(X');
        
        [VS,DS]=eigs(S,F);
        
        for f=1:F
            Corr1(k,f,r)=max([VS(:,f)'*V(:,f),VS(:,f)'*(-V(:,f))])^2;
            Lambda1(k,f,r)=DS(f,f);
        end
    end
end


