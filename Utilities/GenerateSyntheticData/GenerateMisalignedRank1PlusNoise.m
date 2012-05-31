function [S,d,s_d,V,D]=GenerateMisalignedRank1PlusNoise(h,n,dmax,Spacing)


p=length(h);
F=1;

d=Spacing*round(rand(n,1)*dmax/Spacing);
dsorted=sort(unique(d),'ascend');
dPoints=length(dsorted);

if dmax>0
    s_d=1/n*hist(d,unique([0;dsorted]));
else
    s_d=1;
end

A=randn(F,n);
X=zeros(p,n);
for i=1:n
    X(:,i)=circshift(h,d(i))*A(:,i) + randn(p,1);
end
S=cov(X');

H_=zeros(p,p);
for i=1:dPoints
    H_=H_+s_d(i)*circshift(h,dsorted(i))*circshift(h,dsorted(i))';
end

Sigma=H_+eye(p);

[V,D]=eigs(Sigma,F);