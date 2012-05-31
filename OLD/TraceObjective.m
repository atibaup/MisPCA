function [C]=TraceObjective(H,Sigma,a,GridSize,dmax)
warning off;
n=size(Sigma{1},1);
S=length(Sigma);
Grid=0:GridSize:dmax;
nGrid=length(Grid);
C=0;
cvx_begin
expression C;
for i=1:S
    for j=1:nGrid
        if a((i-1)*nGrid+j)>0
            T=CreateTranslationMatrix(Grid(j),n);
            C=C-1/S*a((i-1)*nGrid+j)*trace(H*T'*Sigma{i}*T);
        end
    end
end
cvx_end