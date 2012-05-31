function [C]=TraceObjectiveJoint(H,Sigma,a)
warning off;
n=size(Sigma{1},1);
S=length(Sigma);
C=0;
cvx_begin
expression C(1);
for i=1:S
    for j=1:n
            T=CreateTranslationMatrix(j-1,n);
            C=C-1/S*a((i-1)*n+j)*trace(H*T'*Sigma{i}*T);
    end
end
cvx_end