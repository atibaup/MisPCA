function [T]=CreateTranslationMatrix(d,n)

v=sparse(n,1);
d=mod(d,n);
v(d+1)=1;
T=circulant(v);
