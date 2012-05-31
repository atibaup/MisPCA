function [a]=Estimatea(c,n,S)
% opts = sdpsettings('verbose',0);
% a=sdpvar(n*S,1);
% F=set(sum(abs(a(1:n))) <= 1);
% for i=1:S
%    F=F+set(sum(abs(a((i-1)*n+1:i*n))) <= 1);
% end
% %F=sum(abs(a)) <= lambda;
% sol=solvesdp(F,-a'*c,opts);
% sol=sol
a=zeros(n*S,1);
for i=1:S
    v=zeros(n,1);
    [mval,mindx]=max(c((i-1)*n+1:i*n));
    v(mindx)=1;
    a((i-1)*n+1:i*n)=v;
end
