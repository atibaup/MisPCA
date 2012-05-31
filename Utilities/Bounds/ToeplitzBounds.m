function [lambda,lambda1,lambda2,rh_F,rh_F_12]=ToeplitzBounds(h,dsorted,F)

p=length(h);
dPoints=length(dsorted);
rh=zeros(1,dPoints);
dsorted=dsorted;
for i=1:dPoints
    rh(i)=h'*circshift(h,-dsorted(i));
end
rh=rh;
Rh=toeplitz(rh);
figure; imagesc(Rh)
W=dftmtx(dPoints);

%% Using autocorrelation function
rh_F=real(W*rh(:));

n=1:dPoints;
d=exp(1i*(2*pi*(n-1))/(2*dPoints));
rh_F_12=real(W*diag(d)*rh(:));

omega=sort(rh_F,'descend');

deltabar=min(rh_F_12);
deltahat=max(rh_F_12);

lambda=eigs(Rh,F);
lambda=sort(lambda,'descend');
lambda1(:,1)=(omega(1:F)+deltabar(1)-1);
lambda1(:,2)=(omega(1:F)+deltahat(1)-1);

%% Using ferreira's result directly

c=[rh,0,fliplr(rh(2:end))];
W2=dftmtx(2*(dPoints));

c_F=real(W2*c');

mu2k1=c_F(1:2:end);
mu2k=c_F(2:2:end);

omega=sort(mu2k1,'descend');
omega=omega(1:F);
deltabar=min(mu2k);
deltahat=max(mu2k);
lambda2(:,1)=1/2*(omega+deltabar);
lambda2(:,2)=1/2*(omega+deltahat);
