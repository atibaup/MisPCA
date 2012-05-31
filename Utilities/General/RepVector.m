function [V]=RepVector(v,r)
Nv=length(v);
V=zeros(Nv*r,1);
for i=1:Nv
    V((i-1)*r+1:i*r)=v(i)*ones(r,1);
end